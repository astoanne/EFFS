

##Single run of the simulation study in Scenario II with sample size n=500.

library(slasso)
library(fda)
library(Rfast)
library(parallel)
library(fda.usc)
library(refund)
library(mgcv)
library(fdapace)
source("functions.R",encoding="UTF-8")



case<-"Scenario II"
n_obs_tra<-500
n_obs_test<-4000
n_obs<-n_obs_tra+n_obs_test

data<-simulate_data(case,n_obs=n_obs)
Beta_vero_mat<-data$beta_matrix
X<-data$X
Y<-data$Y

domain<-c(0,1)
length_grid<-500
norder<-4
grid<-seq(0,1,length.out = length_grid)


# Smoothing
n_basis_x<-min(80,length_grid)
n_basis_y<-min(80,length_grid)
breaks_x<-seq(0,1,length.out = (n_basis_x-2))
breaks_y<-seq(0,1,length.out = (n_basis_y-2))
basis_x <- create.bspline.basis(domain,breaks=breaks_x)
basis_y <- create.bspline.basis(domain,breaks=breaks_y)
X_fd <- smooth.basis(grid,X,basis_x)$fd
Y_fd <- smooth.basis(grid,Y,basis_y)$fd
inter_basis<-create.bspline.basis(domain,nbasis = length(grid),norder = 1)
Beta_vero_fd<-bifd(Beta_vero_mat,inter_basis,inter_basis)

# Basis s and t
n_basis_x<-min(60,length_grid)
n_basis_y<-min(60,length_grid)
breaks_x<-seq(0,1,length.out = (n_basis_x-2))
breaks_y<-seq(0,1,length.out = (n_basis_y-2))
basis_x <- create.bspline.basis(domain,breaks=breaks_x)
basis_y <- create.bspline.basis(domain,breaks=breaks_y)

# Matrices
W_X<-eval.penalty(basis_x)
W_X_sp<-sparseMatrix(which(W_X>0,arr.ind=T)[,1],which(W_X>0,arr.ind=T)[,2],x=W_X[which(W_X>0,arr.ind=T)])
W_Y<-eval.penalty(basis_y)
W_Y_sp<-sparseMatrix(which(W_Y>0,arr.ind=T)[,1],which(W_Y>0,arr.ind=T)[,2],x=W_Y[which(W_Y>0,arr.ind=T)])
W_XY<-inprod(basis_x,basis_y)
W_XY_sp<-sparseMatrix(which(W_XY>0,arr.ind=T)[,1],which(W_XY>0,arr.ind=T)[,2],x=W_XY[which(W_XY>0,arr.ind=T)])
R_X<-eval.penalty(basis_x,2)
R_X_sp<-sparseMatrix(which(R_X!=0,arr.ind=T)[,1],which(R_X!=0,arr.ind=T)[,2],x=R_X[which(R_X!=0,arr.ind=T)])
R_Y<-eval.penalty(basis_y,2)
R_Y_sp<-sparseMatrix(which(R_Y!=0,arr.ind=T)[,1],which(R_Y!=0,arr.ind=T)[,2],x=R_Y[which(R_Y!=0,arr.ind=T)])

# Training and Test set
X_fd_tra_nc<-X_fd[1:n_obs_tra]
Y_fd_tra_nc<-Y_fd[1:n_obs_tra]
X_fd_tra<-center.fd(X_fd_tra_nc)
Y_fd_tra<-center.fd(Y_fd_tra_nc)
X_fd_test<-X_fd[(n_obs_tra+1):(n_obs)]-mean_rep(X_fd_tra_nc,n_obs_test)
Y_fd_test<-Y_fd[(n_obs_tra+1):(n_obs)]-mean_rep(Y_fd_tra_nc,n_obs_test)
X_coef_tra<-t(X_fd_tra$coefs)
Y_coef_tra<-t(Y_fd_tra$coefs)
X_coef_test<-t(X_fd_test$coefs)
Y_coef_test<-t(Y_fd_test$coefs)   

ncores=detectCores()
# TRU regression ---------------------------------------------------------
mod_tru_cv<-fr.tru.cv(Y_fd_tra,X_fd_tra,K=10,nbasiss = seq(5,20,by=3),nbasist = seq(5,20,by=3),ncores=ncores)
mod_tru<-fr.tru(Y_fd_tra,X_fd_tra,nbasiss = mod_tru_cv$par_opt[1],mod_tru_cv$par_opt[1])
ISE_tru<-get_ISE(mod_tru$Beta_hat_fd,Beta_vero_fd,case)
PMSE_tru<-get_PMSE(Y_fd_test,X_fd_test,mod_tru$Beta_hat_fd)

# SMOOTH regression --------------------------------------------------------
mod_smooth_cv<-fr.usc.cv(Y_fd_tra,X_fd_tra,basis_x,basis_y,K=10,lambdas_s = 10^seq(-6,-2),lambdas_t = 10^seq(-6,-2),ncores=ncores)
mod_smooth<-fr.usc(Y_fd_tra,X_fd_tra,basis_x,basis_y,lambdas_s=mod_smooth_cv$lambda_s_opt,lambdas_t =mod_smooth_cv$lambda_t_opt)
ISE_smooth<-get_ISE(mod_smooth$Beta_hat_fd,Beta_vero_fd,case)
PMSE_smooth<-get_PMSE(Y_fd_test,X_fd_test,mod_smooth$Beta_hat_fd)


# PFFR regression ---------------------------------------------------------
mod_PFFR<-PFFR(Y_fd_tra,X_fd_tra)
ISE_PFFR<-get_ISE_ref(Beta_vero_fd,case,mod_PFFR$eval_mat)
PMSE_PFFR<-get_PMSE_ref(Y_fd_test,X_fd_test,mod_PFFR$eval_mat)

# PCA regression ----------------------------------------------------------
mod_PCA_cv<-fr.PCA.cv(Y_fd_tra,X_fd_tra,K=10,ncores=ncores)
mod_PCA<-fr.PCA(Y_fd_tra,X_fd_tra,ncomps = mod_PCA_cv$n_comp_opt_x,ncompt = mod_PCA_cv$n_comp_opt_y)
ISE_PCA<-get_ISE(mod_PCA$Beta_hat_fd,Beta_vero_fd,case)
PMSE_PCA<-get_PMSE(Y_fd_test,X_fd_test,mod_PCA$Beta_hat_fd)

# RIDGE regression --------------------------------------------------------
mod_ridge_cv<-fr.ridge.cv(Y_fd_tra,X_fd_tra,basis_x,basis_y,K = 10,alpha_seq = -4:4,ncores=ncores)
mod_ridge<-fr.ridge(Y_fd_tra,X_fd_tra,basis_x,basis_y,10^mod_ridge_cv$par_opt)
ISE_ridge<-get_ISE(mod_ridge$Beta_hat_fd,Beta_vero_fd,case)
PMSE_ridge<-get_PMSE(Y_fd_test,X_fd_test,mod_ridge$Beta_hat_fd)

# SLASSO regression ---------------------------------------------------------------------
mod_slasso_cv<-slasso::slasso.fr_cv(Y_fd_tra,X_fd_tra,basis_x,basis_y,lambda_L_vec =seq(-2,2,by=1),lambda_s_vec =seq(-6,-2),lambda_t_vec =seq(-6,-2),B0=mod_smooth$B,kss_rule_par=0.5,invisible=1,max_iterations=2000,ncores=ncores)
mod_slasso<-slasso::slasso.fr(Y_fd_tra,X_fd_tra,basis_x,basis_y,lambda_L = mod_slasso_cv$lambda_opt_vec[1],lambda_s = mod_slasso_cv$lambda_opt_vec[2],lambda_t = mod_slasso_cv$lambda_opt_vec[3],B0 = mod_smooth$B,invisible=1,max_iterations=10000)
ISE_slasso<-get_ISE(mod_slasso$Beta_hat_fd,Beta_vero_fd,case)
PMSE_slasso<-get_PMSE(Y_fd_test,X_fd_test,mod_slasso$Beta_hat_fd)


ISE= data.frame(matrix(unlist(ISE_tru),ncol=2,byrow = T))
names(ISE)[1]="ISE_0_tru"
names(ISE)[2]="ISE_1_tru"
a= data.frame(matrix(unlist(ISE_smooth),ncol=2,byrow = T))
names(a)[1]="ISE_0_smooth"
names(a)[2]="ISE_1_smooth"
ISE=cbind(ISE,a)
a= data.frame(matrix(unlist(ISE_PFFR),ncol=2,byrow = T))
names(a)[1]="ISE_0_PFFR"
names(a)[2]="ISE_1_PFFR"
ISE=cbind(ISE,a)
a= data.frame(matrix(unlist(ISE_PCA),ncol=2,byrow = T))
names(a)[1]="ISE_0_PCA"
names(a)[2]="ISE_1_PCA"
ISE=cbind(ISE,a)
a= data.frame(matrix(unlist(ISE_ridge),ncol=2,byrow = T))
names(a)[1]="ISE_0_ridge"
names(a)[2]="ISE_1_ridge"
ISE=cbind(ISE,a)
a= data.frame(matrix(unlist(ISE_slasso),ncol=2,byrow = T))
names(a)[1]="ISE_0_slasso_05"
names(a)[2]="ISE_1_slasso_05"
ISE=cbind(ISE,a)



PMSE=data.frame(PMSE_tru,PMSE_smooth,PMSE_PFFR,PMSE_PCA,PMSE_ridge,PMSE_slasso)


write.table(ISE,paste0("Chunk_",kkk,"_ISE.txt"))
write.table(PMSE,paste0("Chunk_",kkk,"_PMSE.txt"))
