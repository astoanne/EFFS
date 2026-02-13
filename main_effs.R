#Begin
#spline.fit.X is the function which used to get the spline fitting of the data.
spline.fit.X<-function( X , nbasis.X , orders , Time.t ){
  Basis.spline<-create.bspline.basis( c(0,1) , nbasis=nbasis.X , norder=orders ) 
  Spline.results<-smooth.basis(Time.t, t(X), Basis.spline)
  X.fd<-Spline.results$fd
  coeff.X<-t(Spline.results$fd$coefs)
  predict.X<-t(eval.fd(Time.t, X.fd))
  errors.X<-rowMeans((predict.X-X)^2) 
  return(list(coeff.X=coeff.X, X.fd=X.fd, predict.X=predict.X, errors.X=errors.X))
}

#############
spline.fit.vector<-function( X , nbasis.X , orders , Time.t ){
  Basis.spline<-create.bspline.basis( c(0,1) , nbasis=nbasis.X , norder=orders ) 
  Spline.results<-smooth.basis(Time.t, X, Basis.spline)
  X.fd<-Spline.results$fd
  coeff.X<-Spline.results$fd$coefs
  predict.X<-eval.fd(Time.t, X.fd)
  errors.X<-mean((predict.X-X)^2) 
  return(list(coeff.X=coeff.X, X.fd=X.fd, predict.X=predict.X, errors.X=errors.X))
}

#####
#L^2 distance between two function formed data or the original data.
Dist.L2<-function(Xa, Time.t, Xb){
  Isfd.a<-is.fd(Xa)
  Isfd.b<-is.fd(Xb)
  delta.X<-Time.t[2]-Time.t[1]
  
  L2.trapzoid<-function(Xa,Xb,Time.t){
    L2<-rep(0,times=length(Xa))
    for (i in 1:(length(Xa)-1)){
      L2[i]<-((Xa[i]-Xb[i])^2+(Xa[i+1]-Xb[i+1])^2)*(Time.t[i+1]-Time.t[i])/2
    }
    L2<-sum(L2)
    return(L2)
  }

  if ((Isfd.a==0)&(Isfd.b==0)){
    na.inda<-which(is.na(Xa)==1)
    na.indb<-which(is.na(Xb)==1)
    na.ind<-unique(c(na.inda,na.indb))
    if (isempty(na.ind)==0){
      Xa<-Xa[-na.ind]
      Xb<-Xb[-na.ind]
      Time.t<-Time.t[-na.ind]
    }
    L2<- L2.trapzoid(Xa,Xb,Time.t)
  } else if ((Isfd.a==0)&(Isfd.b==1)){
    Xb.value<-eval.fd(Time.t, Xb)
    na.ind<-which(is.na(Xa)==1)
    if (isempty(na.ind)==0){
      Xa<-Xa[-na.ind]
      Xb.value<-Xb.value[-na.ind]
      Time.t<-Time.t[-na.ind]
    }
    L2<- L2.trapzoid(Xa,Xb.value,Time.t)
    
  } else if ((Isfd.a==1)&(Isfd.b==0)){
    Xa.value<-eval.fd(Time.t, Xa)
    
    na.ind<-which(is.na(Xb)==1)
    if (isempty(na.ind)==0){
      Xb<-Xb[-na.ind]
      Xa.value<-Xa.value[-na.ind]
      Time.t<-Time.t[-na.ind]
    }
    L2<- L2.trapzoid(Xa.value,Xb,Time.t)
    
  }else{
    Xa.value<-eval.fd(Time.t, Xa)
    Xb.value<-eval.fd(Time.t, Xb)
    L2<- L2.trapzoid(Xa.value,Xb.value,Time.t)
  }
  return(L2)
}  

#####################
#geometric similarity measure
similarity<-function( cen.a , cen.b , Time.t, radius.a , radius.b){
  sim<-exp(-Dist.L2( cen.a, Time.t, cen.b ) - abs( radius.a - radius.b ) )
  return(sim)
}
###############################
main.STR.FFS.visual <- function(add.thr, merge.thr, Time.t, Time.t.len,
                                nbasis.input, nbasis.output,
                                Time.input.in, Time.output.in,
                                NN, data.X.Y.d, orders,
                                ncol.input, ncol.output,
                                M.clu,
                                pruning.strategy = c("none", "usage", "age", "importance", "utility","ema"),
                                prune.threshold = list(
                                  usage = 2,
                                  age = 10,
                                  importance = 0.2,
                                  activation_ratio = 0.4,
                                  utility_sum = 2.0,     # example sum threshold (tunable)
                                  utility_avg = 0.02,     # example average threshold (tunable)
                                  ema = 0.02   # example cutoff
                                )) {
  pruning.strategy <- match.arg(pruning.strategy)
  
  ## --- Initialization --- ##
  PN       <- matrix(list(diag(nbasis.output)), nrow = M.clu, ncol = 1)
  center   <- matrix(0, nrow = M.clu, ncol = nbasis.input)
  cons.para<- matrix(0, nrow = M.clu, ncol = nbasis.output)
  radius   <- numeric(M.clu)
  c.num    <- numeric(M.clu)
  membership <- numeric(M.clu)
  
  # For age-based pruning, track last time each rule was winner
  last_active <- rep(0, M.clu)
  # For importance-based pruning, track cumulative activation
  importance  <- rep(0, M.clu)
  # For utility-based pruning
  birth_t    <- rep(0, M.clu)   # time step when the rule was created
  util_sum   <- rep(0, M.clu)   # sum of normalized firing strengths since birth (LEOA-style)
  util_avg   <- rep(0, M.clu)   # average normalized firing strength since birth
  act_ema <- rep(0, M.clu)   # EMA of normalized firing strength
  ema_beta <- 0.95   
  
  ensure_capacity <- function(n_needed) {
    cur <- nrow(center)
    if (n_needed <= cur) return(invisible(NULL))
    grow_by <- max(1L, ceiling(cur * 0.5))
    center    <<- rbind(center,    matrix(0, nrow = grow_by, ncol = ncol(center)))
    cons.para <<- rbind(cons.para, matrix(0, nrow = grow_by, ncol = ncol(cons.para)))
    radius    <<- c(radius, rep(0, grow_by))
    c.num     <<- c(c.num, rep(0, grow_by))
    last_active <<- c(last_active, rep(0, grow_by))
    importance  <<- c(importance,  rep(0, grow_by))
    PN        <<- c(PN, replicate(grow_by, diag(nbasis.output), simplify = FALSE))
    birth_t  <<- c(birth_t,  rep(0, grow_by))
    util_sum <<- c(util_sum, rep(0, grow_by))
    util_avg <<- c(util_avg, rep(0, grow_by))
    act_ema <<- c(act_ema, rep(0, grow_by))
  }
  # After first data point and initial cluster creation:
  n.clu <- 1
  c.num[1] <- 1
  last_active[1] <- 1
  importance[1]  <- 0          # start from 0 so sums reflect from birth
  birth_t[1]     <- 1
  util_sum[1]    <- 0
  util_avg[1]    <- 0
  act_ema[1] <- 0
  
  gc()

  data.X.Y<-as.data.frame(data.X.Y.d)
  data.X.Y<-datastream_dataframe(data=data.X.Y)
  
  data.real<-data.X.Y$get_points(1)
  data.real<-as.numeric(data.real)
  x.data<-data.real[1:ncol.input]
  y.data<-data.real[(ncol.input+1):(ncol.input+ncol.output)]
  t<-1
  
  Time.input<-Time.input.in#The time inputs for Spline fitting of the inputs data X.
  NA.ind<-which(is.na(x.data)==1)
  if(isempty(NA.ind)==0){
    Time.input<-Time.input[-NA.ind]
    x.data<-x.data[ -NA.ind]
  }
  Time.output<-Time.output.in#The time inputs for Spline fitting of the inputs data X.
  NA.ind<-which(is.na(y.data)==1)
  if(isempty(NA.ind)==0){
    Time.output<-Time.output[-NA.ind]
    y.data<-y.data[ -NA.ind]
  }
  
  x.spline<-spline.fit.vector( x.data , nbasis.input , orders , Time.input )
  y.spline<-spline.fit.vector( y.data , nbasis.output , orders , Time.output )
  
  x.para<-x.spline$X.fd$coefs
  y.para<-y.spline$X.fd$coefs
  
  n.clu<-1#Count the number of clusters
  rule.monitor<-1
  #rm(data.real, x.data)
  
  center[1, ]<-matrix(x.spline$coeff.X ,  nrow = 1) #The cluster center of the first cluster.
  cons.para[1, ]<-matrix( y.spline$coeff.X , nrow = 1)#The consequent parameters of the first cluster.
  cons.fd<-y.spline$X.fd
  #####################################################
  #\{\int B(t)B(t)^{T} dt\}^{-1} only needs to be calculated for one time.
  B.spline.Y<-create.bspline.basis(c(0,1), nbasis = nbasis.output, norder = orders)
  B.t<-eval.basis(Time.t, B.spline.Y)
  B.prime<-matrix(list(), 1, Time.t.len)
  for (i in 1:Time.t.len){
    B.prime[[i]]<-matrix(B.t[i,], nbasis.output, 1)%*%matrix(B.t[i,], 1, nbasis.output)*(1/Time.t.len)
  }
  B.prime<-Reduce('+', B.prime)

  ###############################################
  radius[1]<-0 #The initial value. 
  c.num[1]<-1 #The number of the data within this cluster.
  y.train.func.mse<-Dist.L2(y.spline$X.fd, Time.t, cons.fd)
  predict.y.spline<-eval.fd(Time.output, cons.fd)
  y.train.real.mse<-mean((predict.y.spline-y.data)^2) 
  
  ############################################################################
  ############################################################################
  #The second data comes.
  data.real<-data.X.Y$get_points(1)#The second data.
  data.real<-as.numeric(data.real)
  x.data<-data.real[1:ncol.input]
  y.data<-data.real[(ncol.input+1):(ncol.input+ncol.output)]
  t<-t+1
  
  Time.input<-Time.input.in#The time inputs for Spline fitting of the inputs data X.
  NA.ind<-which(is.na(x.data)==1)
  if(isempty(NA.ind)==0){
    Time.input<-Time.input[-NA.ind]
    x.data<-x.data[ -NA.ind]
  }
  Time.output<-Time.output.in#The time inputs for Spline fitting of the inputs data X.
  NA.ind<-which(is.na(y.data)==1)
  if(isempty(NA.ind)==0){
    Time.output<-Time.output[-NA.ind]
    y.data<-y.data[ -NA.ind]
  }
  
  x.spline<-spline.fit.vector( x.data , nbasis.input , orders , Time.input )
  y.spline<-spline.fit.vector( y.data , nbasis.output , orders , Time.output)
  
  x.para<-cbind(x.para, x.spline$X.fd$coefs)
  y.para<-cbind(y.para, y.spline$X.fd$coefs)
  
  ################################
  #Compute the testing error.
  y.spline.test<-cons.para[1,]
  y.test.fd<-y.spline$X.fd
  y.test.fd$coefs<-as.matrix(y.spline.test)
  
  y.test.func.mse<-Dist.L2(y.spline$X.fd, Time.t, y.test.fd)
  predict.y.spline<-eval.fd(Time.output, y.test.fd)
  y.test.real.mse<-( mean((predict.y.spline-y.data)^2) )/(t-1)
  ###############################
  #Add new cluster, compute the cluster radiuses.
  n.clu<-n.clu+1
  ensure_capacity(n.clu)
  center[n.clu, ]<-matrix(x.spline$coeff.X, nrow = 1 )
  
  cons.para[n.clu, ]<-matrix(y.spline$coeff.X, nrow = 1) 
  
  cen.a<-x.spline$X.fd
  cen.b<-x.spline$X.fd
  cen.a$coefs<-as.matrix(center[1,])
  cen.b$coefs<-as.matrix(center[2,])
  radius[1:n.clu]<-rep(sqrt(Dist.L2(cen.a, Time.t, cen.b))/2 , times=2) 
  c.num[n.clu]<-1
  
  birth_t[n.clu]  <- t
  importance[n.clu] <- 0   # start accumulating from now
  util_sum[n.clu]  <- 0
  util_avg[n.clu]  <- 0
  act_ema[n.clu] <- 0
  
  cc.fd<-x.spline$X.fd
  
  ###############################
  #Compute the membership of the new data regarding to each cluster.
  #Estimate the spline coefficients of the output of the new data (x.data, y.data). 
  #Compute the training results and training error.
  for (i in 1:n.clu){
    #Compute the membership
    cc.fd$coefs<-as.matrix(center[i,])
    den <- 2 * (pmax(radius[i], 1e-8)^2)  # Guard against zero radius
    membership[i] <- exp(-Dist.L2(x.spline$X.fd, Time.t, cc.fd) / den)
  }
  if (all(membership[1:n.clu]==0)){
    membership[1:n.clu]<-rep(1/n.clu, times=n.clu)
  }else{
    membership[1:n.clu]<-membership[1:n.clu]/sum(membership[1:n.clu])
  }
  
  ####
  #Update the consequent parameters using the real output
  #B.t<-eval.basis(Time.output, B.spline.Y)
  
  Time.output.len<-length(Time.output)
  PN[[1]]<-PN[[1]]-(membership[1]*PN[[1]]%*%B.prime)%*%inv(membership[1]*PN[[1]]%*%B.prime+diag(nbasis.output))%*%PN[[1]]
  
  cons.para[1,]<-matrix(cons.para[1,],nrow = 1)+(matrix(y.spline$coeff.X, nrow = 1)-matrix(cons.para[1,],nrow = 1))%*%(PN[[1]]*membership[1])%*%B.prime
  
  y.spline.train<-colsums(do.call(cbind, replicate(nbasis.output, matrix(membership[1:n.clu], nrow = n.clu, ncol=1), simplify = FALSE))*cons.para[1:n.clu, ])
  y.train.fd<-y.spline$X.fd
  y.train.fd$coefs<-as.matrix(y.spline.train)
  
  y.train.func.mse<-( (t-1)*y.train.func.mse + Dist.L2(y.spline$X.fd, Time.t, y.train.fd) ) / t
  predict.y.spline<-eval.fd(Time.output, y.train.fd)
  y.train.real.mse<-((t-1)*y.train.real.mse + mean((predict.y.spline-y.data)^2) ) / t
  
  
  ########################################################################################################
  ########################################################################################################
  #New data comes.
  data.real<-data.X.Y$get_points(1)
  data.real<-as.numeric(data.real)
  x.data<-data.real[1:ncol.input]
  y.data<-data.real[(ncol.input+1):(ncol.input+ncol.output)]
  t<-t+1
  
  Time.input<-Time.input.in#The time inputs for Spline fitting of the inputs data X.
  NA.ind<-which(is.na(x.data)==1)
  if(isempty(NA.ind)==0){
    Time.input<-Time.input[-NA.ind]
    x.data<-x.data[ -NA.ind]
  }
  Time.output<-Time.output.in#The time inputs for Spline fitting of the inputs data X.
  NA.ind<-which(is.na(y.data)==1)
  if(isempty(NA.ind)==0){
    Time.output<-Time.output[-NA.ind]
    y.data<-y.data[ -NA.ind]
  }
  
  x.spline<-spline.fit.vector( x.data , nbasis.input , orders , Time.input )
  y.spline<-spline.fit.vector( y.data , nbasis.output , orders , Time.output )
  
  x.para<-cbind(x.para, x.spline$X.fd$coefs)
  y.para<-cbind(y.para, y.spline$X.fd$coefs)
  
  null.x<-is.null(x.data)
  
  rule.ave<-1
  rule.max<-2
  rule.monitor<-c(rule.monitor, n.clu)
  
  gc()
  gcinfo(verbose = FALSE)
  
  while(null.x==0){
    ###################################################################################################
    #Compute the membership of the new data regarding to each cluster.
    #Compute the testing results and testing error.
    for (i in 1:n.clu){
      #Compute the membership
      cc.fd$coefs<-as.matrix(center[i,])
      den <- 2 * (pmax(radius[i], 1e-8)^2)  # Guard against zero radius
      membership[i] <- exp(-Dist.L2(x.spline$X.fd, Time.t, cc.fd) / den)
      
    }
    activation<-membership[1:n.clu]
    
    if (all(membership[1:n.clu]==0)){
      membership[1:n.clu]<-rep(1/n.clu, times=n.clu)
    }else{
      membership[1:n.clu]<-membership[1:n.clu]/sum(membership[1:n.clu])
    }
    
    max.act <- which.max(membership[1:n.clu])
    last_active[max.act] <- t
    
    y.test<-colsums(do.call(cbind, replicate(nbasis.output, matrix(membership[1:n.clu], nrow = n.clu, ncol=1), simplify = FALSE))*cons.para[1:n.clu, ])
    y.test.fd<-y.spline$X.fd
    y.test.fd$coefs<-as.matrix(y.test)
    
    y.test.func.mse<-( (t-1)*y.test.func.mse + Dist.L2(y.spline$X.fd, Time.t, y.test.fd)) / t
    predict.y.spline<-eval.fd(Time.output, y.test.fd)
    y.test.real.mse<-((t-1)*y.test.real.mse + mean((predict.y.spline-y.data)^2) ) / t
    
    
    if (t<=NN){
      y.test.func.mse.c<-0
      y.test.real.mse.c<-0
    }else if (t>NN){
      
      y.test.func.mse.c<-( (t-1-NN)*y.test.func.mse.c + Dist.L2(y.spline$X.fd, Time.t, y.test.fd)) / (t-NN)
      #predict.y.spline<-eval.fd(Time.output, y.test.fd)
      y.test.real.mse.c<-((t-1-NN)*y.test.real.mse.c + mean((predict.y.spline-y.data)^2) ) / (t-NN)
    }
    ################################################################################################
    #Fuzzy rule adding. The threshold is set according to the \alpha-cut of clusters. \alpha is determined by exp(-(2*sigma)^2 / 2*sigma). 
    #Compute the activation degrees of the input regarding to all the clusters. 
    #A new fuzzy rule should be added when all the activation degrees are lower than the threshold.
    
    add.ind<-0
    if ((max(activation)<add.thr)&&(n.clu<M.clu)){
      
      add.ind<-1
      n.clu<-n.clu+1
      ensure_capacity(n.clu)
      center[n.clu, ]<-matrix(x.spline$coeff.X, nrow = 1 )
      cons.para[n.clu, ]<-matrix(y.spline$coeff.X, nrow = 1) 
      c.num[n.clu]<-1
      birth_t[n.clu]  <- t
      importance[n.clu] <- 0   # start accumulating from now
      util_sum[n.clu]  <- 0
      util_avg[n.clu]  <- 0
      act_ema[n.clu] <- 0
      
      max.act<-which(activation==max(activation))
      if (length(max.act)==1){
        cen.a<-x.spline$X.fd
        cen.a$coefs<-as.matrix(center[max.act, ])
        radius[n.clu]<-sqrt(Dist.L2(cen.a, Time.t, x.spline$X.fd))/2
        
      }else{
        cen.a<-x.spline$X.fd
        radius.poss<-rep( 0 , times=length(max.act) )
        for ( i in 1:length(max.act) ){
          cen.a$coefs<-as.matrix( center[max.act[i], ] )
          radius.poss[i]<-sqrt( Dist.L2( cen.a, Time.t, x.spline$X.fd ) )/2
        }
        radius[n.clu]<-min( radius.poss ) 
        
      }
      #gc()
    }else{
      add.ind<-0
      max.act<-which(activation==max(activation))
      if (length(max.act)==1){
        center.win<-( c.num[max.act]*center[max.act, ]+ x.spline$coeff.X)/( c.num[max.act] + 1 )
        cen.a$coefs<-as.matrix( center.win )
        cen.b$coefs<-as.matrix( center[max.act, ] )
        radius[max.act]<-sqrt( ( 1/(c.num[max.act]+1) ) * ( c.num[max.act]*( radius[max.act]^2 + Dist.L2( cen.a, Time.t, cen.b ) ) 
                                                            + Dist.L2( cen.a, Time.t, x.spline$X.fd ) )   )
        
        center[max.act, ]<-center.win
        
        c.num[max.act]<-c.num[max.act]+1
        
      }else{
        
        for (i in 1:length(max.act) ){
          center.win<-( c.num[ max.act[i] ] * center[ max.act[i] , ] + x.spline$coeff.X )/( c.num[ max.act[i] ] + 1 )
          cen.a$coefs<-as.matrix( center.win )
          cen.b$coefs<-as.matrix( center[max.act[i], ] )
          radius[ max.act[i] ]<-sqrt( ( 1/(c.num[ max.act[i] ]+1) ) * ( c.num[max.act[i]]*( radius[max.act[i]]^2 + Dist.L2( cen.a, Time.t, cen.b ) ) 
                                                                        + Dist.L2( cen.a, Time.t, x.spline$X.fd ) )   )
          center[max.act[i], ]<-center.win
        }
        
        c.num[max.act]<-c.num[max.act]+rep( 1 , times=length(max.act) )
      }
      
      #remove(center.win)
      gc()
    }
    ##################################################################################################
    # If no rule added, consider merging or pruning.
    # Merging.
    
    merged <- FALSE
    if (add.ind==0){
      if (length(max.act)>1){
        max.act<-max.act[1]
      }
      if (n.clu>1){
        cen.a$coefs<-as.matrix( center[max.act, ] )
        sim<-rep( 0 , times=n.clu )
        for (i in 1:n.clu){
          cen.b$coefs<-as.matrix( center[i, ] )
          sim[i]<-similarity( cen.a , cen.b , Time.t , radius[max.act] , radius[i] )
        }
        sim[max.act]=0
        if (max(sim)>merge.thr){
          merged <- TRUE
          max.sim<-which(sim==max(sim))
          max.sim<-max.sim[1]
          #Merge cluster max.sim with max.act.
          cen.b$coefs<-as.matrix( center[max.sim, ] )
          cen.c<-cen.a
          center.new<-( c.num[max.act]*center[max.act, ] + c.num[max.sim]*center[max.sim, ] )/( c.num[max.act]+c.num[max.sim] )
          cen.c$coefs<-as.matrix(center.new)
          radius[max.act]<-sqrt( (1/(c.num[max.act]+c.num[max.sim])) * ( c.num[max.act] * radius[max.act]^2 + c.num[max.sim] * radius[max.sim]^2 
                                                                         + c.num[max.sim] * Dist.L2(cen.b, Time.t, cen.c)  )  ) 
          center[max.act, ]<-center.new
          c.num[max.act]<-c.num[max.act]+c.num[max.sim]
          cons.para[max.act, ]<-( c.num[max.act]*cons.para[max.act, ] + c.num[max.sim]*cons.para[max.sim, ] )/(c.num[max.act]+c.num[max.sim])
          PN[[max.act]]<-(c.num[max.act]*PN[[max.act]]+c.num[max.sim]*PN[[max.sim]])/(c.num[max.act]+c.num[max.sim])
          
          # --- merge tracking stats too ---
          # Utility/importance are additive
          importance[max.act] <- importance[max.act] + importance[max.sim]
          util_sum[max.act]   <- util_sum[max.act]   + util_sum[max.sim]
          
          # Birth time should reflect the *earliest* origin of the merged rule
          birth_t[max.act] <- min(birth_t[max.act], birth_t[max.sim])
          
          # Recompute util_avg with updated life length
          life_len_maxact <- max(1, t - birth_t[max.act] + 1)
          util_avg[max.act] <- util_sum[max.act] / life_len_maxact
          
          # EMA isn’t strictly additive; a reasonable choice is importance-weighted average
          w1 <- importance[max.act]
          w2 <- importance[max.sim]
          w  <- max(1e-12, w1 + w2)   # avoid divide by zero
          act_ema[max.act] <- (w1 * act_ema[max.act] + w2 * act_ema[max.sim]) / w
          # --- end merge tracking stats ---
          center.tem<-center[1:n.clu, ]
          cons.para.tem<-cons.para[1:n.clu, ]
          radius.tem<-radius[1:n.clu]
          c.num.tem<-c.num[1:n.clu]
          PN.tem<-PN[1:n.clu]
          
          n.clu<-n.clu-1
          center[1:n.clu, ]<-center.tem[-max.sim, ]
          cons.para[1:n.clu, ]<-cons.para.tem[-max.sim, ]
          radius[1:n.clu]<-radius.tem[-max.sim]
          c.num[1:n.clu]<-c.num.tem[-max.sim]
          PN[1:n.clu]<-PN.tem[-max.sim]
          gc()
        }
        #rm(sim)
        gc()
      }
    }
    if (merged) {
      # Recompute membership if clusters were merged, so pruning sees updated structure
      for (i in 1:n.clu){
        cc.fd$coefs <- as.matrix(center[i,])
        membership[i] <- exp(-Dist.L2(x.spline$X.fd, Time.t, cc.fd) / (2*(radius[i]^2)))
      }
      if (all(membership[1:n.clu]==0)){
        membership[1:n.clu] <- rep(1/n.clu, times=n.clu)
      } else {
        membership[1:n.clu] <- membership[1:n.clu] / sum(membership[1:n.clu])
      }
    }
    
    # Keep your existing cumulative importance
    importance[1:n.clu] <- importance[1:n.clu] + membership[1:n.clu]
    
    # Utility (LEOA-style sum and average since birth)
    util_sum[1:n.clu] <- util_sum[1:n.clu] + membership[1:n.clu]
    life_len <- pmax(1, t - birth_t[1:n.clu] + 1)
    util_avg[1:n.clu] <- util_sum[1:n.clu] / life_len
    act_ema[1:n.clu] <- ema_beta * act_ema[1:n.clu] + (1 - ema_beta) * membership[1:n.clu]
    
    ####
    #Pruning.
    min_rules_keep <- 2      # at least 2
    warmup_steps   <- 10     # first 10 no prune
    if (pruning.strategy != "none" && n.clu > min_rules_keep && t > warmup_steps) {
      # Identify least-used rule
      min_usage <- min(c.num[1:n.clu])
      min_idx   <- which(c.num[1:n.clu] == min_usage)[1]
      prune_now <- FALSE
      
      if (pruning.strategy == "usage" &&
          min_usage < prune.threshold$usage &&
          membership[min_idx] < prune.threshold$activation_ratio * add.thr) {
        prune_now <- TRUE
        
      } else if (pruning.strategy == "age" &&
                 (t - last_active[min_idx]) > prune.threshold$age &&
                 membership[min_idx] < prune.threshold$activation_ratio * add.thr) {
        prune_now <- TRUE
      } else if (pruning.strategy == "importance" &&
                 importance[min_idx] < prune.threshold$importance) {
        prune_now <- TRUE
      } else if (pruning.strategy == "utility") {
        # Choose one (or both) criteria. Commonly the *average* is used (Angelov 2011).
        # Example: prune if average utility is too low
        if (util_avg[min_idx] < prune.threshold$utility_avg) {
          prune_now <- TRUE
        }
        # Or LEOA-like sum criterion:
        # if (util_sum[min_idx] < prune.threshold$utility_sum) prune_now <- TRUE
      } else if (pruning.strategy == "ema" &&
                 act_ema[min_idx] < prune.threshold$ema) {
        prune_now <- TRUE
      }
      if (prune_now) {
        cat("Pruning rule:", min_idx, "count=", min_usage, "\n")
        center      <- center[-min_idx, , drop = FALSE]
        cons.para   <- cons.para[-min_idx, , drop = FALSE]
        radius      <- radius[-min_idx]
        c.num       <- c.num[-min_idx]
        PN          <- PN[-min_idx]
        last_active <- last_active[-min_idx]
        importance  <- importance[-min_idx]
        birth_t   <- birth_t[-min_idx]
        util_sum  <- util_sum[-min_idx]
        util_avg  <- util_avg[-min_idx]
        act_ema   <- act_ema[-min_idx]
        n.clu       <- n.clu - 1
      } 
    }
    ####
    #Update the consequent parameters and compute the training error.
    
    for (i in 1:n.clu){
      #Compute the membership
      cc.fd$coefs<-as.matrix(center[i,])
      den <- 2 * (pmax(radius[i], 1e-8)^2)  # Guard against zero radius
      membership[i] <- exp(-Dist.L2(x.spline$X.fd, Time.t, cc.fd) / den)
    }
    
    activation<-membership[1:n.clu]
    
    
    if (all(membership[1:n.clu]==0)){
      membership[1:n.clu]<-rep(1/n.clu, times=n.clu)
    }else{
      membership[1:n.clu]<-membership[1:n.clu]/sum(membership[1:n.clu])
    }
    for (k in 1:n.clu){
      PN[[k]]<-PN[[k]]-(membership[k]*PN[[k]]%*%B.prime)%*%inv(membership[k]*PN[[k]]%*%B.prime+diag(nbasis.output))%*%PN[[k]]
      
      cons.para[k,]<-matrix(cons.para[k,],nrow=1)+membership[k]*(matrix(y.spline$coeff.X, nrow = 1)-matrix(cons.para[k,],nrow=1))%*%PN[[k]]%*%B.prime
    }
    ####
    #y.spline.train --- the output estimated after using (x.data , y.data) for training.
    y.spline.train<-colsums(do.call(cbind, replicate(nbasis.output, matrix(membership[1:n.clu], nrow = n.clu, ncol=1), simplify = FALSE))*cons.para[1:n.clu, ])
    y.train.fd<-y.spline$X.fd
    y.train.fd$coefs<-as.matrix(y.spline.train)
    
    y.train.func.mse<-( (t-1)*y.train.func.mse + Dist.L2(y.spline$X.fd, Time.t, y.train.fd) ) / t
    predict.y.spline<-eval.fd(Time.output, y.train.fd)
    y.train.real.mse<-((t-1)*y.train.real.mse + mean((predict.y.spline-y.data)^2) ) / t
    
    rule.ave<-((t-1)*rule.ave+n.clu)/t
    rule.monitor<-c(rule.monitor, n.clu)

    if(n.clu>rule.max){
      rule.max<-n.clu
    }

    ###################################################################################################
    #Get the new coming data.
    data.real<-data.X.Y$get_points(1)
    data.real<-as.numeric(data.real)
    null.x<-isempty(data.real)
    
    ######Comparision######################
    if (t==NN){
      y.train.func.mse.c<-0
      y.train.real.mse.c<-0
    }else if (t>NN){
      
      y.train.func.mse.c<-( (t-1-NN)*y.train.func.mse.c + Dist.L2(y.spline$X.fd, Time.t, y.train.fd) ) / (t-NN)
      predict.y.spline<-eval.fd(Time.output, y.train.fd)
      y.train.real.mse.c<-((t-1-NN)*y.train.real.mse.c + mean((predict.y.spline-y.data)^2) ) / (t-NN)
    }
    if (null.x==0){
      x.data<-data.real[1:ncol.input]
      y.data<-data.real[(ncol.input+1):(ncol.input+ncol.output)]
      t<-t+1
      Time.input<-Time.input.in#The time inputs for Spline fitting of the inputs data X.
      NA.ind<-which(is.na(x.data)==1)
      if(isempty(NA.ind)==0){
        Time.input<-Time.input[-NA.ind]
        x.data<-x.data[ -NA.ind]
      }
      Time.output<-Time.output.in#The time inputs for Spline fitting of the inputs data X.
      NA.ind<-which(is.na(y.data)==1)
      if(isempty(NA.ind)==0){
        Time.output<-Time.output[-NA.ind]
        y.data<-y.data[ -NA.ind]
      }
      x.spline<-spline.fit.vector( x.data , nbasis.input , orders , Time.input )
      y.spline<-spline.fit.vector( y.data , nbasis.output , orders , Time.output )
      
      x.para<-cbind(x.para, x.spline$X.fd$coefs)
      y.para<-cbind(y.para, y.spline$X.fd$coefs)
    }
    gc()
  }
  results<-list(y.train.func.mse=y.train.func.mse, y.test.func.mse=y.test.func.mse,
                y.train.real.mse=y.train.real.mse, y.test.real.mse=y.test.real.mse,
                y.train.func.mse.c=y.train.func.mse.c, y.test.func.mse.c=y.test.func.mse.c, 
                y.train.real.mse.c=y.train.real.mse.c, y.test.real.mse.c=y.test.real.mse.c,
                n.clu=n.clu, rule.ave=rule.ave, rule.max=rule.max, add.thr=add.thr,
                merge.thr=merge.thr, orders=orders, x.spline=x.spline, y.spline=y.spline, 
                x.para=x.para, y.para=y.para, rule.monitor=rule.monitor, center=center, 
                cc.fd=cc.fd, cons.para=cons.para, y.train.fd=y.train.fd)
  return(results)
}

