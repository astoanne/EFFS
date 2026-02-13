rm(list=ls())
source("auxiliary_functions.R")

# Load ocean data (Please install "FRegSigCom" package using the provided .tar.gz file to load the Hawaii ocean data)
data(ocean)

# Number of simulations
nsim = 100

# Number of discrete time points
t_j = 101
# Equally spaced time points in the interval [0, 1]
tp_seq = seq(0, 1, length = t_j)
# Number of predictors
n_pred = 4
# Number of functions in training sample
n_train = 90
# Number of functions in test sample
n_test = 26


# Number of basis functions for predictor
nbasis_est_x = 20
# Number of basis functions for response
nbasis_est_y = 20

MSPE_pls_full = numeric()
MSPE_pls_selected = numeric()
MSPE_pls_main_full = numeric()

RMSPE_pls_full = numeric()
RMSPE_pls_selected = numeric()
RMSPE_pls_main_full = numeric()

MAPE_pls_full = numeric()
MAPE_pls_selected = numeric()
MAPE_pls_main_full = numeric()

Rsq_full_pls_fit = numeric()
Rsq_pls_selected_fit = numeric()
Rsq_full_pls_main_fit = numeric()

Rsq_full_pls_pred = numeric()
Rsq_pls_selected_pred = numeric()
Rsq_full_pls_main_pred = numeric()


for(sim in 1:nsim){
  
  train_index = sample(1: (n_train+n_test), n_train, replace = FALSE)
  Y=ocean[[1]]
  Y.train=Y[train_index,]
  Y.test=Y[-(train_index),]
  
  X.list=list()
  X.train.list=list()
  X.test.list=list()
  t.x.list=list()
  for(i in 1:4){
    X.list[[i]]=ocean[[i+1]]
    X.train.list[[i]]=X.list[[i]][train_index,]
    X.test.list[[i]]=X.list[[i]][-(train_index),]
    t.x.list[[i]]=seq(0,1, length.out = ncol(X.list[[i]]))
  }
  
  model_pls = pls_fun(response = Y.train, predictors_train = X.train.list,
                      predictors_test = X.test.list, nbf_x = nbasis_est_x, nbf_y = nbasis_est_y)
  

  MSPE_pls_full[sim] = mean((Y.test-model_pls$predictions_full_interaction)^2)
  MSPE_pls_selected[sim] = mean((Y.test-model_pls$predictions_selected_interaction)^2)
  MSPE_pls_main_full[sim] = mean((Y.test-model_pls$predictions_main_only)^2)
  
  RMSPE_pls_full[sim] = sqrt(mean(((Y.test - model_pls$predictions_full_interaction)/Y.test)^2))
  RMSPE_pls_selected[sim] = sqrt(mean(((Y.test - model_pls$predictions_selected_interaction)/Y.test)^2))
  RMSPE_pls_main_full[sim] = sqrt(mean(((Y.test - model_pls$predictions_main_only)/Y.test)^2))
  
  MAPE_pls_full[sim] = mean(abs((Y.test - model_pls$predictions_full_interaction)/Y.test))
  MAPE_pls_selected[sim] = mean(abs((Y.test - model_pls$predictions_selected_interaction)/Y.test))
  MAPE_pls_main_full[sim] = mean(abs((Y.test - model_pls$predictions_main_only)/Y.test))
  
  Rsq_full_pls_fit[sim] = 1-mean((Y.train - model_pls$fits_full_interaction)^2) / mean((Y.train - apply(Y.train, 2, mean))^2)
  Rsq_pls_selected_fit[sim] = 1-mean((Y.train - model_pls$fits_selected_interaction)^2) / mean((Y.train - apply(Y.train, 2, mean))^2)
  Rsq_full_pls_main_fit[sim] = 1-mean((Y.train - model_pls$fits_main_only)^2) / mean((Y.train - apply(Y.train, 2, mean))^2)
  
  Rsq_full_pls_pred[sim] = 1-mean((Y.test - model_pls$predictions_full_interaction)^2) / mean((Y.test - apply(Y.test, 2, mean))^2)
  Rsq_pls_selected_pred[sim] = 1-mean((Y.test - model_pls$predictions_selected_interaction)^2) / mean((Y.test - apply(Y.test, 2, mean))^2)
  Rsq_full_pls_main_pred[sim] = 1-mean((Y.test - model_pls$predictions_main_only)^2) / mean((Y.test - apply(Y.test, 2, mean))^2)
  
}

#ts.plot(t(Y.test))
#ts.plot(t(model_pls$predictions_full_interaction))
#ts.plot(t(model_pls$predictions_selected_interaction))
#ts.plot(t(model_pls$predictions_main_only))

mean(MSPE_pls_full);sd(MSPE_pls_full)
mean(MSPE_pls_selected);sd(MSPE_pls_selected)
mean(MSPE_pls_main_full);sd(MSPE_pls_main_full)

mean(RMSPE_pls_full);sd(RMSPE_pls_full)
mean(RMSPE_pls_selected);sd(RMSPE_pls_selected)
mean(RMSPE_pls_main_full);sd(RMSPE_pls_main_full)

mean(MAPE_pls_full);sd(MAPE_pls_full)
mean(MAPE_pls_selected);sd(MAPE_pls_selected)
mean(MAPE_pls_main_full);sd(MAPE_pls_main_full)

mean(Rsq_full_pls_fit); sd(Rsq_full_pls_fit)
mean(Rsq_full_pls_main_fit); sd(Rsq_full_pls_main_fit)
mean(Rsq_pls_selected_fit); sd(Rsq_pls_selected_fit)

mean(Rsq_full_pls_pred); sd(Rsq_full_pls_pred)
mean(Rsq_full_pls_main_pred); sd(Rsq_full_pls_main_pred)
mean(Rsq_pls_selected_pred); sd(Rsq_pls_selected_pred)

