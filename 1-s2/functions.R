

## Data generation ---------------------------------------------------------

simulate_data<-function(case,n_obs=3000) {
  
  length_tot<-500
  grid_s<-grid_t<-seq(0,1,length.out = length_tot)
  
  
  
  # generate X --------------------------------------------------------------
  
  domain<-c(0,1)
  n_basis_x<-32     #random chosen between 10 and 50
  X_basis<-create.bspline.basis(domain,norder = 4,nbasis = n_basis_x)
  X_coef<-matrix(rnorm(n_obs*n_basis_x),nrow = n_basis_x,ncol = n_obs )
  X_fd<-fd(X_coef,X_basis)
  X<-eval.fd(grid_s,X_fd)
  
  # Generate ERROR ----------------------------------------------------------
  
  
  n_basis_eps<-20 #random chosen between 10 and 50
  eps_basis<-create.bspline.basis(domain,norder = 4,nbasis = n_basis_eps)
  eps_coef<-matrix(rnorm(n_obs*n_basis_eps),nrow = n_basis_eps,ncol = n_obs )
  eps_fd<-fd(eps_coef,eps_basis)
  Eps<-t(eval.fd(grid_t,eps_fd))
  
  
  # Define beta -----------------------------------------------------------
  
  if(case=="Scenario I"){
    cat("Scenario I")
    beta<-function(s,t){
      
      if(length(s)!=1){d=expand.grid(s,t)
      colnames(d)=c('s','t')
      
      z_matrix<-matrix(0,nrow=length(s),ncol = length(t),byrow=TRUE)
      z_matrix}
      else{
        z_matrix<-matrix(0,nrow=length(s),ncol = length(t),byrow=TRUE)
        z_matrix}
    }
    
    
    
  }
  if(case=="Scenario II"){
    cat("Scenario II")
    beta<-function(s,t){
      a=0.25
      b=0.25
      if(length(s)!=1){d=expand.grid(s,t)
      colnames(d)=c('s','t')
      z<- -(((d$s-0.5)/a)^2 + ((d$t-0.5)/b)^2) +1
      z[z<0]<-0
      z_matrix<-matrix(z,nrow=length(s))
      z_matrix}
      else{
        z<- -(((s)/a)^2 + ((d)/b)^2) + 1
        z[z<0]<-0
        z}
    }}
  if(case=="Scenario III"){
    cat("Scenario III")
    
    beta<-function(s,t){
      a<-0.05
      b<-0.5
      f_1<-function(s,t){b*(1-s)*sin(10*pi*(s-a-1+sqrt(1-(t-0.5)^2)))}
      f_2<-function(s,t){b*sin(10*pi*(s+a+1-sqrt(1-(t-0.5)^2)))}
      
      z<-matrix(0,length_tot,length_tot)
      for (ii in 1:length(grid_t)) {
        t<-grid_t[ii]
        s_0_1<-grid_s[grid_s>( a+1-sqrt(1-(t-0.5)^2))&grid_s<0.5]
        s_0_2<-grid_s[grid_s<( -a+sqrt(1-(t-0.5)^2))&grid_s>=0.5]
        s_n0_1<-grid_s[grid_s<=(a+1-sqrt(1-(t-0.5)^2))&grid_s<0.5]
        s_n0_2<-grid_s[grid_s>=(-a+sqrt(1-(t-0.5)^2))&grid_s>0.5]
        z_i<-c(f_1(s_n0_1,t),rep(0,length(s_0_1)),rep(0,length(s_0_2)),f_2(s_n0_2,t))
        z[ii,]=z_i
      }
      return(t(z))
    }
  }
  if(case=="Scenario IV"){
    cat("Scenario IV")
    beta<-function(s,t){
      a=0.5
      b=0.5
      c=0.5
      d=0.5
      f_1<-function(s,t){((t-0.5)/c)^3+((s-0.5)/d)^3+((t-0.5)/b)^2 + ((s-0.5)/a)^2+5}
      z<- outer(s,t,f_1)
      z
    }
  }
  
  G<-(1/length(grid_s))*t(eval.basis(grid_s,X_basis))%*%beta(grid_s,grid_t)
  
  Y_parz<-t(X_coef)%*%G
  
  signal_to_noise_ratio<-4
  ##for each observation SN=(sigma^2_signal/sigma^2_error)
  if(case=="Scenario I"){Y = Y_parz + Eps%*%diag(colVars(Eps)^(-1/2))}
  else{
    k <- sqrt((colVars(Y_parz)+max(colVars(Y_parz)))/(signal_to_noise_ratio*colVars(Eps)))
    Y = Y_parz + Eps%*%diag(k)
  }
  
  out<-list(X=X,
            Y=t(Y),
            X_fd=X_fd,
            Eps=Eps,
            beta_matrix=beta(grid_s,grid_t))
  
  return(out)
}

## TRU regression ----------------------------------

fr.tru.cv<-function(Y_fd,X_fd,K=10,nbasiss=seq(4,10,by=1),nbasist=seq(4,10,by=1),basis_type="bspline",ncores=1){##Cross-validation for the TRU regression
  
  X_fd_cen=center.fd(X_fd)
  Y_fd_cen=center.fd(Y_fd)
  domainx<-X_fd_cen$basis$rangeval
  domainy<-Y_fd_cen$basis$rangeval
  if(basis_type=="bspline"){
    create_fun=create.bspline.basis
    }
  else if(basis_type=="fourier"){
    create_fun=create.fourier.basis
  }
  
  n_basis_x_seq<-nbasiss
  n_basis_y_seq<-nbasist
  n_obs<-length(X_fd_cen$coefs[1,])
  
  comb_list<-expand.grid(n_basis_x_seq,n_basis_y_seq)
  inpr_vec<-list()
  parr_tru<-function(ii){
    
    basis<-as.numeric(comb_list[ii,])
    n_basis_x<-basis[1]
    n_basis_y<-basis[2]
    basis_x <- create_fun(domainx,nbasis = n_basis_x)
    basis_y <- create_fun(domainy,nbasis=n_basis_y)
    ran_seq<-sample(seq(1, n_obs), n_obs, replace=FALSE)
    split_vec<-split(ran_seq,cut(seq(1,n_obs),breaks=K,labels=FALSE))
    
    for(ll in 1:K){
      Y_i<-Y_fd_cen[split_vec[[ll]]]
      X_i<-X_fd_cen[split_vec[[ll]],]
      X_minus<-X_fd_cen[-split_vec[[ll]],]
      Y_minus<-Y_fd_cen[-split_vec[[ll]],]
      mod<-fregre.basis.fr(X_minus,Y_minus,basis.s =basis_x,basis.t = basis_y, lambda.s = 0,lambda.t = 0)
      Y_hat<-predict(mod,X_i)
      inpr_vec[[ll]]<-diag(inprod(Y_i-Y_hat,Y_i-Y_hat))
    }
    mean<-mean(unlist(inpr_vec))
    sd<-sd(unlist(inpr_vec))/sqrt(n_obs)
    out<- list(mean=mean,
               sd=sd)
    return(out)
    
  }
  
  vec_par<-mclapply(seq(1,length(comb_list[,1])),parr_tru,mc.cores = ncores)
  par<-sapply(vec_par,"[[",1)
  sd<-sapply(vec_par,"[[",2)
  n_basis_opt<-as.numeric(comb_list[max(which(par==min(par,na.rm = T))),])
  
  out<-list(par_opt=n_basis_opt,
            CV=par,
            CV_sd=sd,
            comb_list=comb_list,
            X=X_fd,
            Y=Y_fd,
            type="TRU")
  
  return(out)
  
}
fr.tru<-function(Y_fd,X_fd,nbasiss=5,nbasist=5,basis_type="bspline",...){##TRU regression for given values of nbasiss and nbasist
 
  n_obs<-length(X_fd$coefs[1,])
  
  if(basis_type=="bspline"){
    create_fun=create.bspline.basis
  }
  else if(basis_type=="fourier"){
    create_fun=create.fourier.basis
  }
  
  domainx<-X_fd$basis$rangeval
  domainy<-Y_fd$basis$rangeval
  basis_x <- create_fun(domainx,nbasis = nbasiss)
  basis_y <- create_fun(domainy,nbasis=nbasist)
  
  mod<-fregre.basis.fr(X_fd,Y_fd,basis.s =basis_x,basis.t = basis_y, lambda.s = 0,lambda.t = 0)
  
  B<-mod$beta.estbifd$coefs
  Beta_hat_fd<-mod$beta.estbifd
  cat("TRU:",nbasiss,nbasist,"     ")
  
  out<-list(B=B,
            Beta_hat_fd=Beta_hat_fd,
            nbasiss=nbasiss,
            nbasist=nbasist,
            mod=mod,
            type="TRU")
  
  return(out)
  
}

## SMOOTH regression -------------------------------
  
fr.usc.cv<-function(Y_fd,X_fd,basis_x,basis_y,K=10,lambdas_s=10^seq(5,15,by=1),lambdas_t=10^seq(5,15,by=1),ncores=1){##Cross-validation for the SMOOTH regression
    
    n_obs<-length(X_fd$coefs[1,])
    inpr_vec<-list()
    comb_list<-expand.grid(lambdas_s,lambdas_t)
    parr_smooth<-function(ii){
      
      lambdas<-as.numeric(comb_list[ii,])
      lambda_s<-lambdas[1]
      lambda_t<-lambdas[2]
      ran_seq<-sample(seq(1, n_obs), n_obs, replace=FALSE)
      split_vec<-split(ran_seq,cut(seq(1,n_obs),breaks=K,labels=FALSE))
      
      for(ll in 1:K){
        Y_i<-Y_fd[split_vec[[ll]]]
        X_i<-X_fd[split_vec[[ll]],]
        X_minus<-X_fd[-split_vec[[ll]],]
        Y_minus<-Y_fd[-split_vec[[ll]],]
        mod<-fregre.basis.fr(X_minus,Y_minus,basis.s =basis_x,basis.t = basis_y, Lfdobj.s = 2,Lfdobj.t = 2,lambda.s = lambda_s,lambda.t = lambda_t)
        Y_hat<-predict(mod,X_i)
        inpr_vec[[ll]]<-diag(inprod(Y_i-Y_hat,Y_i-Y_hat))
      }
      mean<-mean(unlist(inpr_vec))
      sd<-sd(unlist(inpr_vec))/sqrt(n_obs)
      out<- list(mena=mean,
                 sd=sd)
      
      return(out)
      
    }
    
    vec_par<-mclapply(seq(1,length(comb_list[,1])),parr_smooth,mc.cores = ncores)
    par<-sapply(vec_par,"[[",1)
    sd<-sapply(vec_par,"[[",2)
    l_opt<-as.numeric(comb_list[max(which(par==min(par))),])
    lambda_s_opt<-l_opt[1]
    lambda_t_opt<-l_opt[2]
    mod<-fregre.basis.fr(X_fd,Y_fd,basis.s =basis_x,basis.t = basis_y, Lfdobj.s = 2,Lfdobj.t = 2,lambda.s = lambda_s_opt,lambda.t = lambda_t_opt)
    B<-mod$beta.estbifd$coefs
    Beta_hat_fd<-mod$beta.estbifd
    
    out<-list(B=B,
              Beta_hat_fd=Beta_hat_fd,
              lambda_s_opt=lambda_s_opt,
              lambda_t_opt=lambda_t_opt,
              CV=par,
              CV_sd=sd,
              comb_list=comb_list,
              type="SMOOTH")
    return(out)
  }
fr.usc<-function(Y_fd,X_fd,basis_x,basis_y,K=10,lambdas_s=0,lambdas_t=0){##SMOOTH regression for given values of lambdas_s and lambdas_t
    
    n_obs<-length(X_fd$coefs[1,])
    mod<-fregre.basis.fr(X_fd,Y_fd,basis.s =basis_x,basis.t = basis_y, Lfdobj.s = 2,Lfdobj.t = 2,lambda.s = lambdas_s,lambda.t = lambdas_t)
    
    B<-mod$beta.estbifd$coefs
    Beta_hat_fd<-mod$beta.estbifd
    l_opt<-c(lambdas_s,lambdas_t)
    cat("SMOOTH:",l_opt,"     ")
    out<-list(B=B,
              Beta_hat_fd=Beta_hat_fd,
              lambda_s=lambdas_s,
              lambda_t=lambdas_t,
              mod=mod,
              type="SMOOTH")
    return(out)
  }

# PFFR regression ---------------------------------------------------------

PFFR<-function(Y_fd_tra,X_fd_tra){
  grid<-seq(0,1,length.out = 30)
  s<-attr(grid, "xindex")
  Y_eval<-t(eval.fd(grid,Y_fd_tra))
  X_eval<-t(eval.fd(grid,X_fd_tra))
  data=list()
  data$X=X_eval
  data$Y=Y_eval
  m2 <- refund::pffr(Y ~ -1+ff(X,xind = grid,splinepars = list(bs="ps", m=c(2, 2),k=15)),data = data,yind = grid,method="REML")
  coeff<-refund:::coef.pffr(m2,n2=500)
  eval_mat2=matrix(coeff$smterms[[1]]$value,500,500)
  out<-list(eval_mat=eval_mat2,
            mod=m2)
}  

## PCA regression ----------------------------------

fr.PCA.cv<-function(Y_fd,X_fd,K=10,ncores=1){##Cross-validation for the PCA regression
    
  X_fd_cen=center.fd(X_fd)
  Y_fd_cen=center.fd(Y_fd)
  X_mean<-mean.fd(X_fd)
  Y_mean<-mean.fd(Y_fd)
  n_obs<-length(X_fd$coefs[1,])
  domain_s=X_fd$basis$rangeval
  domain_t=Y_fd$basis$rangeval
  length_grid=500
  grid_s<-seq(domain_s[1],domain_s[2],length.out = length_grid)
  grid_t<-seq(domain_t[1],domain_t[2],length.out = length_grid)
  delta_t<-((domain_t[2]-domain_t[1])/length_grid)
  delta_s<-((domain_s[2]-domain_s[1])/length_grid)
  X_fd_eval<-t(eval.fd(grid_s,X_fd_cen))
  Y_fd_eval<-t(eval.fd(grid_t,Y_fd_cen))
  basis_int_s<-create.bspline.basis(domain_s,nbasis=length_grid,norder=1)
  basis_int_t<-create.bspline.basis(domain_t,nbasis=length_grid,norder=1)
  grid_list_s<-lapply(1:n_obs, function(ii)grid_s)
  grid_list_t<-lapply(1:n_obs, function(ii)grid_t)
  X_list<-split(t(X_fd_eval),rep(1:ncol(t(X_fd_eval)), each = length_grid))
  Y_list<-split(t(Y_fd_eval),rep(1:ncol(t(Y_fd_eval)), each = length_grid))
  
  
  PCA_x<-FPCA(X_list,grid_list_s,optns = list(useBinnedData="OFF"))
  PCA_y<-FPCA(Y_list,grid_list_t,optns = list(useBinnedData="OFF"))
  n_comp_x<-seq(3,PCA_x$selectK-2,by=2)
  n_comp_y<-seq(3,PCA_y$selectK-2,by=2)
  
  comb_list<-expand.grid(n_comp_x,n_comp_y)
  inpr_vec<-numeric()
  
  parr_PCA<-function(ii){
    
    n_components<-as.numeric(comb_list[ii,])
    n_comp_x<-n_components[1]
    n_comp_y<-n_components[2]
    
    ran_seq<-sample(seq(1, n_obs), n_obs, replace=FALSE)
    split_vec<-split(ran_seq,cut(seq(1,n_obs),breaks=K,labels=FALSE))
    CV_ll=list()
    
    for(ll in 1:K){
      Y_i<-Y_fd_cen[split_vec[[ll]]]
      X_i<-X_fd_eval[split_vec[[ll]],]
      X_minus<-X_fd_eval[-split_vec[[ll]],]
      Y_minus<-Y_fd_eval[-split_vec[[ll]],]
      X_fd_minus<-X_fd_cen[-split_vec[[ll]]]
      Y_fd_minus<-Y_fd_cen[-split_vec[[ll]]]
      
      grid_list_s<-lapply(1:ncol(t(X_minus)), function(ii)grid_s)
      grid_list_t<-lapply(1:ncol(t(Y_minus)), function(ii)grid_t)
      
      X_list<-split(t(X_minus),rep(1:ncol(t(X_minus)), each = nrow(t(X_minus))))
      Y_list<-split(t(Y_minus),rep(1:ncol(t(Y_minus)), each = nrow(t(Y_minus))))
      
      PCA_x<-FPCA(X_list,grid_list_s,optns = list(methodSelectK=n_comp_x,useBinnedData="OFF"))
      PCA_y<-FPCA(Y_list,grid_list_t,optns = list(methodSelectK=n_comp_y,useBinnedData="OFF"))
      
      var_xy<-var.fd(X_fd_minus,Y_fd_minus)
      var_xy_eval<-eval.bifd(grid_s,grid_t,var_xy)
      sigma_xy<-delta_s*delta_t*t(PCA_x$phi)%*%var_xy_eval%*%PCA_y$phi
      rho<-PCA_x$lambda
      beta_mat<-PCA_x$phi%*%diag(1/rho)%*%sigma_xy%*%t(PCA_y$phi)
      Y_hat_mat<-delta_s*X_i%*%beta_mat
      Y_hat<-fd(t(Y_hat_mat),basis_int_t)
      CV_ll[[ll]]<-diag(inprod(Y_i-Y_hat,Y_i-Y_hat))
    }
    mean<-mean(unlist(CV_ll))
    sd<-sd(unlist(CV_ll))/sqrt(n_obs)
    out<- list(mean=mean,
               sd=sd)
    return(out)
    
    
  }
  
  vec_par<-mclapply(seq(1,length(comb_list[,1])),parr_PCA,mc.cores = ncores)
  par<-sapply(vec_par,"[[",1)
  sd<-sapply(vec_par,"[[",2)
  
  n_comps_opt<-as.numeric(comb_list[min(which(par==min(par))),])
  n_comp_x<-n_comps_opt[1]
  n_comp_y<-n_comps_opt[2]
  
  cat("PCA:",n_comps_opt,"     ")
  grid_list_s<-lapply(1:n_obs, function(ii)grid_s)
  grid_list_t<-lapply(1:n_obs, function(ii)grid_t)
  X_list<-split(t(X_fd_eval),rep(1:ncol(t(X_fd_eval)), each = length_grid))
  Y_list<-split(t(Y_fd_eval),rep(1:ncol(t(Y_fd_eval)), each = length_grid))
  PCA_x<-FPCA(X_list,grid_list_s,optns = list(methodSelectK=n_comp_x,useBinnedData="OFF"))
  PCA_y<-FPCA(Y_list,grid_list_t,optns = list(methodSelectK=n_comp_y,useBinnedData="OFF"))
  var_xy<-var.fd(X_fd_cen,Y_fd_cen)
  var_xy_eval<-eval.bifd(grid_s,grid_t,var_xy)
  sigma_xy<-delta_s*delta_t*t(PCA_x$phi)%*%var_xy_eval%*%PCA_y$phi
  rho<-PCA_x$lambda
  
  beta_mat<-PCA_x$phi%*%diag(1/rho)%*%sigma_xy%*%t(PCA_y$phi)
  basis_int_s<-create.bspline.basis(domain_s,nbasis=length_grid,norder=1)
  basis_int_t<-create.bspline.basis(domain_t,nbasis=length_grid,norder=1)
  Beta_fd<-bifd(as.matrix(beta_mat),basis_int_s,basis_int_t)
  
  out<-list(Beta_hat_fd=Beta_fd,
            n_comp_opt_x=n_comp_x,
            n_comp_opt_y=n_comp_y,
            CV=par,
            CV_sd=sd,
            comb_list=comb_list,
            X=X_fd,
            Y=Y_fd,
            type="PCA")
  return(out)
  
}
fr.PCA<-function(Y_fd,X_fd,ncomps,ncompt){##PCA regression for given values of number of components ncomps and ncompt
  
  X_fd_cen=center.fd(X_fd)
  Y_fd_cen=center.fd(Y_fd)
  X_mean<-mean.fd(X_fd)
  Y_mean<-mean.fd(Y_fd)
  n_obs<-length(X_fd$coefs[1,])
  domain_s=X_fd$basis$rangeval
  domain_t=Y_fd$basis$rangeval
  length_grid=500
  grid_s<-seq(domain_s[1],domain_s[2],length.out = length_grid)
  grid_t<-seq(domain_t[1],domain_t[2],length.out = length_grid)
  delta_t<-((domain_t[2]-domain_t[1])/length_grid)
  delta_s<-((domain_s[2]-domain_s[1])/length_grid)
  X_fd_eval<-t(eval.fd(grid_s,X_fd_cen))
  Y_fd_eval<-t(eval.fd(grid_t,Y_fd_cen))
  basis_int_s<-create.bspline.basis(domain_s,nbasis=length_grid,norder=1)
  basis_int_t<-create.bspline.basis(domain_t,nbasis=length_grid,norder=1)
  grid_list_s<-lapply(1:n_obs, function(ii)grid_s)
  grid_list_t<-lapply(1:n_obs, function(ii)grid_t)
  
  n_comps<-c(ncomps,ncompt)
  n_comp_x<-n_comps[1]
  n_comp_y<-n_comps[2]
  
  cat("PCA:",n_comps,"     ")
  X_list<-split(t(X_fd_eval),rep(1:n_obs, each = length_grid))
  Y_list<-split(t(Y_fd_eval),rep(1:n_obs, each = length_grid))
  PCA_x<-FPCA(X_list,grid_list_s,optns = list(methodSelectK=n_comp_x,useBinnedData="OFF"))
  PCA_y<-FPCA(Y_list,grid_list_t,optns = list(methodSelectK=n_comp_y,useBinnedData="OFF"))
  
  var_xy<-var.fd(X_fd_cen,Y_fd_cen)
  var_xy_eval<-eval.bifd(grid_s,grid_t,var_xy)
  sigma_xy<-delta_s*delta_t*t(PCA_x$phi)%*%var_xy_eval%*%PCA_y$phi
  
  rho<-PCA_x$lambda
  
  B<-PCA_x$phi%*%diag(1/rho)%*%sigma_xy%*%t(PCA_y$phi)
  basis_int_s<-create.bspline.basis(domain_s,nbasis=length_grid,norder=1)
  basis_int_t<-create.bspline.basis(domain_t,nbasis=length_grid,norder=1)
  Beta<-bifd(as.matrix(B),basis_int_s,basis_int_t)
  X_mean_new<-inprod(X_mean,Beta$sbasis)
  alpha<-Y_mean-fd(t(X_mean_new%*%B),Beta$tbasis)
  
  
  out<-list(Beta_hat_fd=Beta,
            B=B,
            alpha=alpha,
            ncompx=n_comp_x,
            ncompy=n_comp_y,
            X=X_fd,
            Y=Y_fd,
            type="PCA")
  return(out)
  
}
 
## RIDGE regression ----------------------------------------------------

fr.ridge.cv<-function(Y_fd,X_fd,basiss,basist,K=10,alpha_seq=NA,ncores=1){##Cross-validation for the RIDGE regression
  
  X_fd_cen=center.fd(X_fd)
  Y_fd_cen=center.fd(Y_fd)
  X_mean<-mean.fd(X_fd)
  Y_mean<-mean.fd(Y_fd)
  n_obs<-length(X_fd$coefs[1,])
  domain_s=X_fd$basis$rangeval
  domain_t=Y_fd$basis$rangeval
  X_new<-inprod(X_fd_cen,basiss)
  Y_new<-inprod(Y_fd_cen,basist)
  W_s<-eval.penalty(basiss)
  W_t<-eval.penalty(basist)
  R_s<-eval.penalty(basiss,2)
  R_t<-eval.penalty(basist,2)
  
  W_t_inv<-solve(W_t)
  if(is.na(alpha_seq)[1]){
    alpha_seq<-seq(-5,5)
  }
  
  inpr_vec<-list()
  parr_ridge<-function(ii){
    
    alpha<-alpha_seq[ii]
    ran_seq<-sample(seq(1, n_obs), n_obs, replace=FALSE)
    split_vec<-split(ran_seq,cut(seq(1,n_obs),breaks=K,labels=FALSE))
    
    for(ll in 1:K){
      Y_i<-Y_fd_cen[split_vec[[ll]]]
      X_i<-X_new[split_vec[[ll]],]
      X_minus<-X_new[-split_vec[[ll]],]
      Y_minus<-Y_new[-split_vec[[ll]],]
      
      Beta_hat<-solve(t(X_minus)%*%X_minus+(10^alpha)*W_s)%*%t(X_minus)%*%Y_minus%*%W_t_inv
      Y_hat<-fd(as.matrix(t(X_i%*%Beta_hat)),basist)
      inpr_vec[[ll]]<-diag(inprod(Y_i-Y_hat,Y_i-Y_hat))
      
    }
    mean<-mean(unlist(inpr_vec))
    sd<-sd(unlist(inpr_vec))/sqrt(n_obs)
    out<- list(mean=mean,
               sd=sd)
    return(out)
  }
  
  vec_par<-mclapply(seq(1,length(alpha_seq)),parr_ridge,mc.cores = ncores)
  par<-sapply(vec_par,"[[",1)
  sd<-sapply(vec_par,"[[",2)
  
  l_opt<-alpha_seq[max(which(par==min(par)))]
  
  out<-list(par_opt=l_opt,
            CV=par,
            CV_sd=sd,
            comb_list=alpha_seq,
            X=X_fd,
            Y=Y_fd,
            type="RIDGE")
  return(out)
}
fr.ridge<-function(Y_fd,X_fd,basiss,basist,alpha_pen=10^0){##RIDGE regression for given values of alpha_pen
  
  X_fd_cen=center.fd(X_fd)
  Y_fd_cen=center.fd(Y_fd)
  X_mean<-mean.fd(X_fd)
  Y_mean<-mean.fd(Y_fd)
  n_obs<-length(X_fd$coefs[1,])
  domain_s=X_fd$basis$rangeval
  domain_t=Y_fd$basis$rangeval
  X_new<-inprod(X_fd_cen,basiss)
  Y_new<-inprod(Y_fd_cen,basist)
  W_s<-eval.penalty(basiss)
  W_t<-eval.penalty(basist)
  W_t_inv<-solve(W_t)
  cat("RIDGE:",alpha_pen,"     ")
  B<-solve(t(X_new)%*%X_new+(alpha_pen)*W_s)%*%t(X_new)%*%Y_new%*%W_t_inv
  Beta<-bifd(as.matrix(B),basiss,basist)
  X_mean_new<-inprod(X_mean,basiss)
  alpha<-Y_mean-fd(t(X_mean_new%*%B),basist)
  
  out<-list(Beta_hat_fd=Beta,
            B=B,
            alpha=alpha,
            alpha_pen=alpha_pen,
            X=X_fd,
            Y=Y_fd,
            type="RIDGE")
  
  return(out)
}


## ISE and PMSE and other functions------------------------------------------------------------

get_ISE<-function(beta_hat_fd,Beta_vero_fd,case){
  
  length_grid_int<-500
  delta<-1/length_grid_int
  grid_int<-seq(0,1,length.out = length_grid_int)
  delta_bifd<-difference_bifd(Beta_vero_fd,beta_hat_fd)
  eval_mat<-eval.bifd(grid_int,grid_int, delta_bifd)^2
  A<- eval.bifd(grid_int,grid_int,Beta_vero_fd)
  ind_0<-which(A==0,arr.ind = T)
  sum_0<-sum(eval_mat[ind_0])
  sum_1<-sum(eval_mat)-sum_0
  area_0<-length(which(A==0))/(length_grid_int^2)
  area_1<-1-area_0
  ISE_0<-(1/area_0)*delta*delta*sum_0
  ISE_1<-(1/area_1)*delta*delta*sum_1
  
  out<-list(ISE_0=ISE_0,
            ISE_1=ISE_1)
  return(out)
}
get_ISE_ref<-function(Beta_vero_fd,casei,eval_bhat){
  
  length_grid_int<-500
  delta<-1/length_grid_int
  grid_int<-seq(0,1,length.out = length_grid_int)
  A<- eval.bifd(grid_int,grid_int,Beta_vero_fd)
  eval_mat<-(A-eval_bhat)^2
  
  ind_0<-which(A==0,arr.ind = T)
  sum_0<-sum(eval_mat[ind_0])
  sum_1<-sum(eval_mat)-sum_0
  area_0<-length(which(A==0))/(length_grid_int^2)
  area_1<-1-area_0
  ISE_0<-(1/area_0)*delta*delta*sum_0
  ISE_1<-(1/area_1)*delta*delta*sum_1
  
  if(case=="Concurrent"){
    ISE_1<-(sum(diag(eval_mat))*delta*delta)/(delta^2*length_grid_int)
    ISE_tot<-sum(eval_mat)*delta*delta
    ISE_0<-(ISE_tot-sum(diag(eval_mat))*delta*delta)/(1-delta^2*length_grid_int)
  }
  out<-list(ISE_0=ISE_0,
            ISE_1=ISE_1)
  return(out)
}
difference_bifd<-function(bifd_1,bifd_2){
  grid_s<-grid_t<-seq(0,1,length.out = 500)
  X_1<-eval.bifd(grid_s,grid_t,bifd_1)
  X_2<-eval.bifd(grid_s,grid_t,bifd_2)
  diff<-(X_1-X_2)
  B_spline_0<-create.bspline.basis(domain,nbasis = length(grid_s),norder = 1)
  bifd(diff,B_spline_0,B_spline_0)
}
mean_rep<-function (fdobj,nobs){
  coef <- as.array(fdobj$coefs)
  coefd <- dim(coef)
  ndim <- length(coefd)
  basisobj <- fdobj$basis
  nbasis <- basisobj$nbasis
  coefmean <- apply(coef, 1, mean)
  mean_rep <- fd(matrix(rep(coefmean,nobs),coefd[1],nobs), basisobj)
  return(mean_rep)
}
get_PMSE<-function(Y_fd_test,X_fd_test,Beta,case_1=case){
  length_grid<-500
  grid<-seq(0,1,length.out = length_grid)
  delta<-1/length_grid
  X_fd_eval<-t(eval.fd(grid,X_fd_test))
  Y_fd_eval<-t(eval.fd(grid,Y_fd_test))
  Beta_mat<-eval.bifd(grid,grid,Beta)
  Y_hat<-delta*X_fd_eval%*%Beta_mat
  PMSE<-mean(delta*rowSums((Y_fd_eval-Y_hat)^2))
  
  return(PMSE)
}
get_PMSE_ref<-function(Y_fd_test,X_fd_test,eval_mat,case_1=case){
  
  length_grid<-500
  grid<-seq(0,1,length.out = length_grid)
  delta<-1/length_grid
  X_fd_eval<-t(eval.fd(grid,X_fd_test))
  Y_fd_eval<-t(eval.fd(grid,Y_fd_test))
  Beta_mat<-eval_mat
  if(case_1=="Concurrent"){
    Y_hat<-X_fd_eval%*%Beta_mat
  }
  else{
    Y_hat<-delta*X_fd_eval%*%Beta_mat
  }
  PMSE<-mean(delta*rowSums((Y_fd_eval-Y_hat)^2))
  
  return(PMSE)
}
