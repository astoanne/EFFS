#data
complete <- with(DTI, complete.cases(cca) & complete.cases(rcst))
cca <- DTI$cca[complete,]
rcst<-DTI$rcst[complete,]
###################################################################################
#Fit X by splines and compute the fitting MSPE. 
input<-cca
output<-rcst
NN<-50#When to compete
ntrain<-NN
ntest<-nrow(cca)-ntrain
train.index<-c(1:ntrain)
nbasis.input<-20
nbasis.output<-20

X.all = input
Y = output
Data.new<-list(N=ncol(cca), nsamples=nrow(cca), num_out=(nrow(cca)-50), Predictor=X.all, 
               Response=Y)
# Model Selection 
source('/Users/dongjiaoge/Desktop/TNNLS resubmit/RobustFLR-master/Model_select_FLR.R')

# max_number : the maximum number of FPCs to be estimated 
FLR_models = Model_select_FLR(Data.new, max_number=5)
print ("Model selected")

source('/Users/dongjiaoge/Desktop/TNNLS resubmit/RobustFLR-master/RobustFLR_functions.R')
max_number<-5
Sim_output<-Data.new
N = Sim_output$N 
nsamples  = Sim_output$nsamples
num_out = Sim_output$num_out

normal_samples = nsamples-num_out
outlier_samples = nsamples - normal_samples

X.all = Sim_output$Predictor
Y = Sim_output$Response
choice = 3 # B generated from unif(-choice,choice) 
out_var = 0.5 # Outliers from B'=B+R, where R~N(0,out_var)
alpha=0.3 #proportion of trimming 

breakdown_p = round(nsamples*(1-alpha))

# BIC and RBIC are estimated using classical FPCA  
BIC = matrix(nrow = max_number,ncol = max_number)
RBIC = matrix(nrow = max_number,ncol = max_number)

# BIC_rob and RBIC_rob are estimated using robust FPCA
BIC_rob = matrix(nrow = max_number,ncol = max_number)
RBIC_rob = matrix(nrow = max_number,ncol = max_number)




# number of eigenfunctions
num_comp_X = max_number+1
num_comp_Y = max_number+1

#### setup B-spline basis ####
N.X<-ncol(cca)
N.Y<-ncol(rcst)
gaittime.X <- seq(0,1,len=N.X)
gaittime.Y <- seq(0,1,len=N.Y)
gaitrange <- c(0,1)
nord=20
gaitbasis <- create.bspline.basis(c(0,1), nbasis = nord, norder = 2)




fdParobj = fdPar(gaitbasis,int2Lfd(2),lambda=0) #Applies smoothness to the functions
gaitfd.X <- smooth.basisPar(gaittime.X, t(X.all), gaitbasis, Lfdobj=NULL, lambda=0)$fd
gaitfd.Y <- smooth.basisPar(gaittime.Y, t(Y), gaitbasis, Lfdobj=NULL, lambda=0)$fd
# FPCA on Y
resp_fd<- gaitfd.Y[,2]
resp_pca = pca.fd(resp_fd,nharm=num_comp_Y,fdParobj)
resp_scores=resp_pca$scores

# FPCA components of Y
harmonics_tgt = resp_pca$harmonics
mu_y = resp_pca$meanfd
#print(sum(resp_pca$varprop[1:num_comp_Y]))

# FPCA on X
pcs = pca.fd(gaitfd.X[,1],nharm=num_comp_X,fdParobj)
harmonics_x = pcs$harmonics
pred_scores=pcs$scores
#print(sum(pcs$varprop[1:num_comp_X]))

coefX=pred_scores
coefY=resp_scores


#### Residuals using robust FPCA ####

robY = robust_FPCA_Bali_Bspline(Y[1:NN,], harmonics_tgt, gaitbasis, num_comp_Y, N.Y)

robust_FPCA_y = robY$harmonics
med_mu_y = robY$med_mu_y
rob_coefY = robY$rob_coef

# X
robX = robust_FPCA_Bali_Bspline(X.all, harmonics_x, gaitbasis, num_comp_X, N.X)

rob_coefX = robX$rob_coef

beta_est  = MLTS_iterated_reweighted(rob_coefX[1:NN,], rob_coefY,ntrain)$beta

Y_hat_rob = matrix(ncol= N.Y, nrow = ntest)
PE = rep(0,ntest)
rob_residual = matrix(ncol=ntest, nrow=N.Y)
rob_coefX.test<-rob_coefX[-(1:NN),]
for (dd in 1:ntest){
  
  Xtest <- rob_coefX[dd,]
  truth=Y[dd,] 
  rec_one_out= Pred_MFLR(Xtest, beta_est, robust_FPCA_y) + med_mu_y
  Y_hat_rob[dd,] = eval.fd(rec_one_out, gaittime.Y)
  rob_residual[,dd] = truth-eval.fd(rec_one_out, gaittime.Y)
}



##########################
#########################
########################













for(hhh in 2:(max_number+1)){
  for (ggg in 2:(max_number+1)){
    
    #### classical FLR ####
    
    # number of eigenfunctions
    num_comp_X = hhh
    num_comp_Y = ggg
    
    #### setup B-spline basis ####
    N.X<-ncol(cca)
    N.Y<-ncol(rcst)
    gaittime.X <- seq(0,1,len=N.X)
    gaittime.Y <- seq(0,1,len=N.Y)
    gaitrange <- c(0,1)
    nord=10
    gaitbasis <- create.bspline.basis(c(0,1), nbasis = nord, norder = 2)
    
    
  
    
    fdParobj = fdPar(gaitbasis,int2Lfd(2),lambda=0) #Applies smoothness to the functions
    gaitfd.X <- smooth.basisPar(gaittime.X, t(X.all), gaitbasis, Lfdobj=NULL, lambda=0)$fd
    gaitfd.Y <- smooth.basisPar(gaittime.Y, t(Y), gaitbasis, Lfdobj=NULL, lambda=0)$fd
    # FPCA on Y
    resp_fd<- gaitfd.Y[,2]
    resp_pca = pca.fd(resp_fd,nharm=num_comp_Y,fdParobj)
    resp_scores=resp_pca$scores
    
    # FPCA components of Y
    harmonics_tgt = resp_pca$harmonics
    mu_y = resp_pca$meanfd
    #print(sum(resp_pca$varprop[1:num_comp_Y]))
    
    # FPCA on X
    pcs = pca.fd(gaitfd.X[,1],nharm=num_comp_X,fdParobj)
    harmonics_x = pcs$harmonics
    pred_scores=pcs$scores
    #print(sum(pcs$varprop[1:num_comp_X]))
    
    coefX=pred_scores
    coefY=resp_scores
    
    # Estimate regression matrix
    classical_beta = MFLR(coefX, coefY) 
    
    
    
    #### residuals using classical FLR ####  
    Y_hat = matrix(ncol= N.Y, nrow = nsamples) # estimates 
    classical_residuals = matrix(ncol=nsamples, nrow=N.Y) #residuals
    
    for (dd in 1:nsamples){
      Xtest <- coefX[dd,]
      truth= Y[dd,] 
      rec_one_out= Pred_MFLR(Xtest, classical_beta, harmonics_tgt) + mu_y
      Y_hat[dd,] = eval.fd(rec_one_out, gaittime.Y)
      
      classical_residuals[,dd] = truth-eval.fd(rec_one_out, gaittime.Y)
    }
    
    func_data_class = fdata(t(classical_residuals),argvals=NULL,rangeval=NULL,names=NULL,fdata2d=FALSE)
    #plot(func_data_class, main='Residuals using classical FLR')
    #lines(func_data_class[(normal_samples+1):nsamples,], col=3)
    
    #### Classical BIC and RBIC ####
    
    dist_FLR = rep(0,nsamples)
    for (pp in 1:nsamples){
      dist_FLR[pp] = t(Y[pp,]-Y_hat[pp,])%*%(Y[pp,]-Y_hat[pp,])
    }
    
    scale_FLR = sd(dist_FLR)
    
    index = sort.int(dist_FLR, decreasing = FALSE, index.return = TRUE)
    num_trim = round(0.8*nsamples)
    new_ind = index$ix[1:num_trim]
    
    BIC[hhh-1,ggg-1] = (nsamples)*log( sum(dist_FLR)/(nsamples) ) + (num_comp_X*num_comp_Y+1)*log(nsamples)  
    RBIC[hhh-1,ggg-1] = (num_trim)*log( sum(dist_FLR[new_ind])/(num_trim) ) + (num_comp_X*num_comp_Y+1)*log(num_trim)
    
    
    #### Residuals using robust FPCA ####
    
    robY = robust_FPCA_Bali_Bspline(Y, harmonics_tgt, gaitbasis, num_comp_Y, N.Y)
    
    robust_FPCA_y = robY$harmonics
    med_mu_y = robY$med_mu_y
    rob_coefY = robY$rob_coef
    
    # X
    robX = robust_FPCA_Bali_Bspline(X.all, harmonics_x, gaitbasis, num_comp_X, N.X)
    
    rob_coefX = robX$rob_coef
    
    beta_est  = MLTS_iterated_reweighted(rob_coefX, rob_coefY,nsamples)$beta
    
    Y_hat_rob = matrix(ncol= N.Y, nrow = nsamples)
    PE = rep(0,nsamples)
    rob_residual = matrix(ncol=nsamples, nrow=N.Y)
    for (dd in 1:nsamples){
      Xtest <- rob_coefX[dd,]
      truth=Y[dd,] 
      rec_one_out= Pred_MFLR(Xtest, beta_est, robust_FPCA_y) + med_mu_y
      Y_hat_rob[dd,] = eval.fd(rec_one_out, gaittime.Y)
      rob_residual[,dd] = truth-eval.fd(rec_one_out, gaittime.Y)
    }
    
    
    func_data_rob = fdata(t(rob_residual),argvals=NULL,rangeval=NULL,names=NULL,fdata2d=FALSE)
    #plot(func_data_rob, main='Residuals using robust FLR', col=1)
    #lines(func_data_rob[(normal_samples+1):nsamples,], col=3)
    
    #### Robust BIC and RBIC ####
    
    dist_FLR = rep(0,nsamples)
    for (pp in 1:nsamples){
      dist_FLR[pp] = t(Y[pp,]-Y_hat_rob[pp,])%*%(Y[pp,]-Y_hat_rob[pp,])
    }
    
    scale_FLR_rob = sd(dist_FLR)
    
    index = sort.int(dist_FLR, decreasing = FALSE, index.return = TRUE)
    num_trim = round(0.8*nsamples)
    new_ind = index$ix[1:num_trim]
    
    BIC_rob[hhh-1,ggg-1] = (nsamples)*log( sum(dist_FLR)/(nsamples) ) + (num_comp_X*num_comp_Y+1)*log(nsamples)  
    RBIC_rob[hhh-1,ggg-1] = (0.8*nsamples)*log( sum(dist_FLR[new_ind])/(nsamples*0.8) ) + (num_comp_X*num_comp_Y+1)*log(nsamples*0.8)
    
  }
  #print ("Model iteration")
}




RBIC_comp_X = FLR_models$RBIC_comp_X
RBIC_comp_Y = FLR_models$RBIC_comp_Y

comp_X = FLR_models$comp_X
comp_Y = FLR_models$comp_Y

# Model fit 
RBIC_val = FLR_models$RBIC_val #model fit using robust model with RBIC
BIC_val = FLR_models$BIC_val #model fit using robust model with BIC
BIC_class_val = FLR_models$BIC_class_val #model fit using classical model with BIC

RBIC_vec[jj] = RBIC_val
BIC_vec[jj] = BIC_val
BIC_class_vec[jj] = BIC_class_val

#### Depth values ####
class_depth_val = class_depth(Output, comp_X,comp_Y)
rob_depth_val = rob_depth(Output, RBIC_comp_X,RBIC_comp_Y)

func_data = fdata(Y,argvals=NULL,rangeval=NULL,names=NULL,fdata2d=FALSE)
Y_depth = depth.mode(func_data)$dep


#### ROC curves ####


labels = c(rep(1, nsamples-num_out), rep(0,num_out))
pred_class <- prediction( class_depth_val, labels )
pred_rob <- prediction( rob_depth_val, labels )
pred_y <- prediction(Y_depth, labels )

perf_class = performance(pred_class, 'auc')
perf_rob = performance(pred_rob, 'auc')
perf_y = performance(pred_y, 'auc')

auc_class[jj] = as.numeric(perf_class@y.values)
auc_rob[jj] = as.numeric(perf_rob@y.values)
auc_y[jj] = as.numeric(perf_y@y.values)
