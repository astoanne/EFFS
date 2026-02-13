# Packages
#install.packages(c("fda", "MASS", "combinat", "FRegSigCom", "plsdepot", "expm"))
library(fda)
library(MASS)
library(FRegSigCom)
library(combinat)
library(plsdepot)
library(expm)


### PLS function
# response: a matrix containing functional response variable
# predictors_train: a list containing functional predictor variables (training sample)
# predictors_test: a list containing functional predictor variables (test sample)
# nbf_x: number of basis functions to approximate the functional predictors
# nbf_y: number of basis functions to approximate the functional response

pls_fun = function(response, predictors_train, predictors_test, nbf_x, nbf_y, n.order){
  
  # Percentages of the data used to determine the optimum number of PLS components and significant variables
  forw_n_train = round(0.5 * nrow(response))
  forw_n_test = nrow(response) - forw_n_train
  
  # Number of predictors
  np = length(predictors_train)
  
  # Discrete time points
  mtp = seq(0, 1, length=ncol(response))
  
  ## Dongjiao added
  mtp.1 = seq(0, 1, length=ncol(predictors_train[[1]]))
  # B-spline basis for predictors
  B_spline_basis_x = create.bspline.basis(c(0,1), nbasis = nbf_x, norder = n.order)
  B_spline_basis_funs_x = eval.basis(mtp, B_spline_basis_x)
  
  # B-spline basis for response
  B_spline_basis_y = create.bspline.basis(c(0,1), nbasis = nbf_y, norder = n.order)
  B_spline_basis_fun_y = eval.basis(mtp, B_spline_basis_y)
  
  # Inner products
  Inner_prod_y = inprod(B_spline_basis_y, B_spline_basis_y)
  Inner_prod_x = inprod(B_spline_basis_x, B_spline_basis_x)
  
  # Square roots of inner products
  Inner_prod_y_sqrt = sqrtm(Inner_prod_y)
  Inner_prod_x_sqrt = sqrtm(Inner_prod_x)
  
  # Arguments for smoothing step
  smooth_arg_train = matrix(mtp, nrow = nrow(response), ncol = ncol(response), byrow=T)
  #Dongjiao changed this 
  smooth_arg_test = matrix(mtp.1, nrow = nrow(predictors_test[[1]]), ncol = ncol(predictors_test[[1]]), byrow=T)
  #Dongjiao changed this 
  smooth_arg_train_x = matrix(mtp.1, nrow = nrow(predictors_train[[1]]), ncol = ncol(predictors_train[[1]]), byrow=T)#dongjiao add this
  # Convert predictors and response into functional forms
  X_train_smooth = vector("list", length = np)
  X_test_smooth = vector("list", length = np)
  
  for(sm in 1:np){
    X_train_smooth[[sm]] = t(smooth.basis(argvals=t(smooth_arg_train_x), y=t(predictors_train[[sm]]), fdParobj=B_spline_basis_x)$fd$coefs)
    #X_train_smooth[[sm]] = t(smooth.basis(argvals=t(smooth_arg_train), y=t(predictors_train[[sm]]), fdParobj=B_spline_basis_x)$fd$coefs)#dongjiao comment this
    X_test_smooth[[sm]] = t(smooth.basis(argvals=t(smooth_arg_test), y=t(predictors_test[[sm]]), fdParobj=B_spline_basis_x)$fd$coefs)
  }
  
  Y_train_smooth = t(smooth.basis(argvals=t(smooth_arg_train), y=t(response), fdParobj=B_spline_basis_y)$fd$coefs)
  Y_train_reg = Y_train_smooth %*% Inner_prod_y_sqrt
  
  # Approximate response functions in the training sample 
  approx_response = (Y_train_smooth)%*%t(B_spline_basis_fun_y)
  
  ### Variable selection
  # Matrix for main effects
  
  main_effects_full_model = vector("list", length = np)
  for(main_full in 1:np)
    main_effects_full_model[[main_full]] = X_train_smooth[[main_full]] %*% Inner_prod_x_sqrt
  
  # Order of the importance of the main effects
  error_full_model = numeric()
  for(order in 1:np){
    Reg_mat_main_full_model = cbind(1,main_effects_full_model[[order]])
    pls_main_full_model = plsreg2(predictors = Reg_mat_main_full_model[,2:ncol(Reg_mat_main_full_model)],
                                        responses = Y_train_reg, comps = 4, crosval = TRUE)
    Coeff_pls_main_full_model = rbind(pls_main_full_model$reg.coefs[ncol(Reg_mat_main_full_model),], 
                                      pls_main_full_model$reg.coefs[1:(ncol(Reg_mat_main_full_model)-1),])
    
    predicted_curves_full_model = (Reg_mat_main_full_model %*% Coeff_pls_main_full_model) %*% solve(Inner_prod_y_sqrt) %*% t(B_spline_basis_fun_y)
    error_full_model[order] = mean((approx_response - predicted_curves_full_model)^2)
  }
  # Order of variables according to their MSE
  imp_var = order(error_full_model)
  
  # Forward procedure
  # Starting model
  Reg_mat_main_forw_start = cbind(1, main_effects_full_model[[imp_var[1]]])
  forward_error = min(error_full_model)
  
  np_forw = c(imp_var[1], rep(NA, (length(imp_var)-1)))
  selected_main = c(which.min(error_full_model))
  
  for(forw1 in 2:np){
    error_for_forward_selection = rbind(subset(imp_var, !(imp_var %in% selected_main)), NA)
    for(forw in 1:ncol(error_for_forward_selection)){
      Reg_mat_main_forw_model = cbind(Reg_mat_main_forw_start, main_effects_full_model[[error_for_forward_selection[1,forw]]])
      
      pls_main_forw_model = plsreg2(predictors = Reg_mat_main_forw_model[,2:ncol(Reg_mat_main_forw_model)],
                                    responses = Y_train_reg, comps = 4, crosval = TRUE)
      
      Coeff_pls_main_forw_model = rbind(pls_main_forw_model$reg.coefs[ncol(Reg_mat_main_forw_model),], 
                                        pls_main_forw_model$reg.coefs[1:(ncol(Reg_mat_main_forw_model)-1),])
      
      predicted_curves_forw_model = (Reg_mat_main_forw_model %*% Coeff_pls_main_forw_model) %*% solve(Inner_prod_y_sqrt) %*% t(B_spline_basis_fun_y)
      error_for_forward_selection[1,forw] = mean((approx_response - predicted_curves_forw_model)^2)#Dongjiao write this line
      #error_for_forward_selection[2,forw] = mean((approx_response - predicted_curves_forw_model)^2)#dongjiao comment this
    }
    
    err_next = error_for_forward_selection[1,][which.min(error_for_forward_selection[1,])]#Dongjiao write this line
    var_next = error_for_forward_selection[1,][which.min(error_for_forward_selection[1,])]#Dongjiao write this line
    #err_next = error_for_forward_selection[2,][which.min(error_for_forward_selection[2,])]#Dongjiao comment this
    #var_next = error_for_forward_selection[1,][which.min(error_for_forward_selection[2,])]#Dongjiao comment this line
    
    selected_main = c(selected_main,var_next)
    
    if(err_next < forward_error){
      Reg_mat_main_forw_start = cbind(Reg_mat_main_forw_start, main_effects_full_model[[var_next]])
      forward_error = err_next
      np_forw[forw1] = var_next
    }else if(err_next > forward_error){
      Reg_mat_main_forw_start = Reg_mat_main_forw_start
      forward_error = forward_error
      np_forw[forw1] = np_forw[forw1]
    }
  }
  
  np_forw = sort(subset(np_forw, !(np_forw %in% NA)))
  
  # Length of interaction terms
  len_int = length(np_forw)
  
  # A matrix for the quadratic and interaction terms
  interaction_matrix_full = matrix(0, len_int+len_int*(len_int-1)/2,2)
  kk = 1
  for(ik in 1:len_int){
    for(jk in ik:len_int){
      interaction_matrix_full[kk,1] = sort(np_forw)[ik];
      interaction_matrix_full[kk,2] = sort(np_forw)[jk];
      kk = kk+1;
    }
  }
  
  # Number of quadratic and interaction terms
  no_int = nrow(interaction_matrix_full)
  
  # Matrix for quadratic and interaction terms
  interaction_coefficient_full = vector("list",)
  for(it_full in 1:no_int){
    interaction_coefficient_full_i = matrix(NA, nrow = nrow(response) , ncol = nbf_x^2)
    for(it_full_i in 1:nrow(response)){
      interaction_coefficient_full_i[it_full_i,] = (Inner_prod_x_sqrt %x% Inner_prod_x_sqrt) %*% 
        (X_train_smooth[[interaction_matrix_full[it_full,1]]][it_full_i,] %x% X_train_smooth[[interaction_matrix_full[it_full,2]]][it_full_i,])
      
    }
    interaction_coefficient_full[[it_full]] = interaction_coefficient_full_i
  }
  
  # Order of interaction terms
  error_interaction_model = numeric()
  
  for(order_int in 1:no_int){
    Reg_mat_int_order = cbind(Reg_mat_main_forw_start, interaction_coefficient_full[[order_int]])
    
    model_int_order = plsreg2(predictors = Reg_mat_int_order[,2:ncol(Reg_mat_int_order)],
                              responses = Y_train_reg, comps = 8, crosval = TRUE)
    
    coef_int_order = rbind(model_int_order$reg.coefs[ncol(Reg_mat_int_order),], 
                           model_int_order$reg.coefs[1:(ncol(Reg_mat_int_order)-1),])
    
    pred_curves_int_order = (Reg_mat_int_order %*% coef_int_order) %*% solve(Inner_prod_y_sqrt) %*% t(B_spline_basis_fun_y)
    
    error_interaction_model[order_int] = mean((approx_response - pred_curves_int_order)^2)
  }
  
  int_var = order(error_interaction_model)
  
  # Forward procedure
  # Starting model
  Reg_mat_with_int_effs = Reg_mat_main_forw_start
  
  interaction_matrix_selected = rep(NA, nrow(interaction_matrix_full))
  selected_int = numeric()
  
  for(inter1 in 1: no_int){
    int_error_for_forward_selection = rbind(subset(int_var, !(int_var %in% selected_int)), NA)
    for(inter in 1:ncol(int_error_for_forward_selection)){
      Reg_mat_int_forw_model = cbind(Reg_mat_with_int_effs, interaction_coefficient_full[[int_error_for_forward_selection[1,inter]]])
      pls_interaction_forw_model = plsreg2(predictors = Reg_mat_int_forw_model[,2:ncol(Reg_mat_int_forw_model)],
                                           responses = Y_train_reg, comps = 8, crosval = TRUE)
      
      Coeff_pls_interaction_forw = rbind(pls_interaction_forw_model$reg.coefs[ncol(Reg_mat_int_forw_model),], 
                                         pls_interaction_forw_model$reg.coefs[1:(ncol(Reg_mat_int_forw_model)-1),])
      
      predicted_curves_interaction_forw =  (Reg_mat_int_forw_model %*% Coeff_pls_interaction_forw) %*% solve(Inner_prod_y_sqrt) %*% t(B_spline_basis_fun_y)
      
      int_error_for_forward_selection[2,inter] = mean((approx_response - predicted_curves_interaction_forw)^2)
    }
    
    err_int_next = int_error_for_forward_selection[2,][which.min(int_error_for_forward_selection[2,])]
    int_next = int_error_for_forward_selection[1,][which.min(int_error_for_forward_selection[2,])]
    selected_int = c(selected_int, int_next)
    
    if(err_int_next < forward_error){
      Reg_mat_with_int_effs = cbind(Reg_mat_with_int_effs, interaction_coefficient_full[[int_next]])
      forward_error = err_int_next
      interaction_matrix_selected[inter1] = int_next
    }else if(err_int_next > forward_error){
      Reg_mat_with_int_effs = Reg_mat_with_int_effs
      forward_error = forward_error
      interaction_matrix_selected[inter1] = interaction_matrix_selected[inter1]
    }
  }
  
  interaction_terms = subset(interaction_matrix_selected, !(interaction_matrix_selected %in% NA))
  if(is.logical(interaction_terms) == FALSE)
    interaction_matrix_selected = matrix(interaction_matrix_full[sort(interaction_terms),], ncol=2)
    
  
  ### Full interaction model
  # main effect terms
  
  full_model_main_effects = vector("list", )
  for(full_me in 1:np){
    full_model_main_effects[[full_me]] = X_train_smooth[[full_me]] %*% Inner_prod_x_sqrt
  }
  
  # Quadratic and interaction terms
  full_interaction_matrix = matrix(0, np+np*(np-1)/2,2);
  fk = 1
  for(fi in 1:np){
    for(fj in fi:np){
      full_interaction_matrix[fk,1] = (1:np)[fi];
      full_interaction_matrix[fk,2] = (1:np)[fj];
      fk = fk+1;
    }
  }
  
  full_model_interaction_coefficients = vector("list", )
  for(full_it in 1:nrow(full_interaction_matrix)){
    full_model_interaction_coefficients_i = matrix(NA, nrow = nrow(response) , ncol = nbf_x^2)
    for(full_iti in 1:nrow(response)){
      full_model_interaction_coefficients_i[full_iti,] = (Inner_prod_x_sqrt %x% Inner_prod_x_sqrt) %*% 
        (X_train_smooth[[full_interaction_matrix[full_it,1]]][full_iti,] %x% X_train_smooth[[full_interaction_matrix[full_it,2]]][full_iti,])
    }
    full_model_interaction_coefficients[[full_it]] = full_model_interaction_coefficients_i
  }
  
  # Regression matrix for full-interaction model
  Reg_mat_full = cbind(1, do.call(cbind, full_model_main_effects), do.call(cbind, full_model_interaction_coefficients))
  
  # Regression matrix only for main effects model
  Reg_mat_full_main_only = cbind(1, do.call(cbind, full_model_main_effects))
  
  ### Selected interaction model
  # Main effects
  selected_main_effects = np_forw
  selected_model_main_effects = vector("list", )
  for(selected_me in 1:length(selected_main_effects)){
    selected_model_main_effects[[selected_me]] = X_train_smooth[[np_forw[selected_me]]] %*% Inner_prod_x_sqrt 
  }
  
  # Quadratic and interaction terms
  if(is.logical(interaction_terms) == TRUE){
    Reg_mat_selected = cbind(1, do.call(cbind, selected_model_main_effects))
  }else if(is.logical(interaction_terms) == FALSE){
    selected_interaction_matrix = interaction_matrix_selected
    
    selected_model_interaction_coefficients = vector("list", )
    for(selected_it in 1:nrow(selected_interaction_matrix)){
      selected_model_interaction_coefficients_i = matrix(NA, nrow = nrow(response) , ncol = nbf_x^2)
      for(selected_iti in 1:nrow(response)){
        selected_model_interaction_coefficients_i[selected_iti,] = (Inner_prod_x_sqrt %x% Inner_prod_x_sqrt) %*% 
          (X_train_smooth[[selected_interaction_matrix[selected_it,1]]][selected_iti,] %x%
             X_train_smooth[[selected_interaction_matrix[selected_it,2]]][selected_iti,])
      }
      selected_model_interaction_coefficients[[selected_it]] = selected_model_interaction_coefficients_i
    }
    
    # Regression matrix for selected-interaction model
    Reg_mat_selected = cbind(1, do.call(cbind, selected_model_main_effects), do.call(cbind, selected_model_interaction_coefficients))
  }
  
  
  # Optimum number of PLS
  pls_comp_err_selected = numeric()
  for(icomp_selected in 1:10){
    comp_index_train_selected = sample(1: nrow(Reg_mat_selected), forw_n_train, replace=FALSE) 
    comp_index_test_selected = (1:nrow(Reg_mat_selected))[-comp_index_train_selected]
    
    model_pls_selected = plsreg2(predictors = Reg_mat_selected[comp_index_train_selected,][,2:ncol(Reg_mat_selected)], 
                                 responses = Y_train_reg[comp_index_train_selected,], comps = icomp_selected, crosval = TRUE)
    
    Coeff_pls_selected = rbind(model_pls_selected$reg.coefs[ncol(Reg_mat_selected),], 
                               model_pls_selected$reg.coefs[1:(ncol(Reg_mat_selected)-1),])
    
    sim_comp_err_selected1 =  (Reg_mat_selected[comp_index_test_selected,] %*% Coeff_pls_selected) %*% solve(Inner_prod_y_sqrt) %*% t(B_spline_basis_fun_y)
    
    pls_comp_err_selected[icomp_selected] = mean((approx_response[comp_index_test_selected,] - sim_comp_err_selected1)^2)
    
  }
  
  num_comp_selected = which.min(pls_comp_err_selected)
  
  ### Estimation
  # Full model
  full_model_pls = plsreg2(predictors = Reg_mat_full[,2:ncol(Reg_mat_full)], 
                           responses = Y_train_reg, comps = num_comp_selected, crosval = TRUE)
  full_model_coeff = rbind(full_model_pls$reg.coefs[ncol(Reg_mat_full),], 
                           full_model_pls$reg.coefs[1:(ncol(Reg_mat_full)-1),])
  
  # Selected model
  selected_model_pls = plsreg2(predictors = Reg_mat_selected[,2:ncol(Reg_mat_selected)], 
                               responses = Y_train_reg, comps = num_comp_selected, crosval = TRUE)
  
  selected_model_coeff = rbind(selected_model_pls$reg.coefs[ncol(Reg_mat_selected),], 
                               selected_model_pls$reg.coefs[1:(ncol(Reg_mat_selected)-1),])
  
  # Full model main effects only
  full_model_pls_main_only = plsreg2(predictors = Reg_mat_full_main_only[,2:ncol(Reg_mat_full_main_only)], 
                                     responses = Y_train_reg, comps = num_comp_selected, crosval = TRUE)
  
  full_model_coeff_main_only = rbind(full_model_pls_main_only$reg.coefs[ncol(Reg_mat_full_main_only),], 
                                     full_model_pls_main_only$reg.coefs[1:(ncol(Reg_mat_full_main_only)-1),])
  
  ### Testing phase
  # Full model
  # Main effects
  full_model_main_effects_test = vector("list", )
  for(full_me_test in 1:np){
    full_model_main_effects_test[[full_me_test]] = X_test_smooth[[full_me_test]] %*% Inner_prod_x_sqrt
  }
  
  # Quadratic and interaction effects
  full_model_interaction_coefficients_test = vector("list", )
  for(full_it_test in 1:nrow(full_interaction_matrix)){
    full_model_interaction_coefficients_test_i = matrix(NA, nrow = nrow(predictors_test[[1]]), ncol = nbf_x^2)
    for(full_iti_test in 1:nrow(predictors_test[[1]])){
      full_model_interaction_coefficients_test_i[full_iti_test,] = (Inner_prod_x_sqrt %x% Inner_prod_x_sqrt) %*% 
        (X_test_smooth[[full_interaction_matrix[full_it_test,1]]][full_iti_test,] %x%
           X_test_smooth[[full_interaction_matrix[full_it_test,2]]][full_iti_test,])
    }
    full_model_interaction_coefficients_test[[full_it_test]] = full_model_interaction_coefficients_test_i
  }
  
  # Regression matrix for full-interaction model
  Reg_mat_full_test = cbind(1, do.call(cbind, full_model_main_effects_test), do.call(cbind, full_model_interaction_coefficients_test))
  
  # Regression matrix only for main effects model
  Reg_mat_full_main_only_test = cbind(1, do.call(cbind, full_model_main_effects_test))
  
  
  # Selected model
  # Main effects
  selected_model_main_effects_test = vector("list", )
  for(selected_me_test in 1:length(selected_main_effects)){
    selected_model_main_effects_test[[selected_me_test]] = X_test_smooth[[np_forw[selected_me_test]]] %*% Inner_prod_x_sqrt
  }
  
  # Quadratic and interaction effects
  if(is.logical(interaction_terms) == TRUE){
    Reg_mat_selected_test = cbind(1, do.call(cbind, selected_model_main_effects_test))
  }else if(is.logical(interaction_terms) == FALSE){
    selected_model_interaction_coefficients_test = vector("list", )
    for(selected_it_test in 1:nrow(selected_interaction_matrix)){
      selected_model_interaction_coefficients_test_i = matrix(NA, nrow = nrow(predictors_test[[1]]), ncol = nbf_x^2)
      for(selected_iti_test in 1:nrow(predictors_test[[1]])){
        selected_model_interaction_coefficients_test_i[selected_iti_test,] = (Inner_prod_x_sqrt %x% Inner_prod_x_sqrt) %*% 
          (X_test_smooth[[selected_interaction_matrix[selected_it_test,1]]][selected_iti_test,] %x%
             X_test_smooth[[selected_interaction_matrix[selected_it_test,2]]][selected_iti_test,])
      }
      selected_model_interaction_coefficients_test[[selected_it_test]] = selected_model_interaction_coefficients_test_i
    }
    # Regression matrix for selected-interaction model
    Reg_mat_selected_test = cbind(1, do.call(cbind, selected_model_main_effects_test), do.call(cbind, selected_model_interaction_coefficients_test))
  }
  
  ### Fits
  # Full model with interaction
  fit_full_interaction = (Reg_mat_full %*% full_model_coeff) %*% solve(Inner_prod_y_sqrt)  %*% t(B_spline_basis_fun_y)
  
  # Full model main effects only
  fit_full_main = (Reg_mat_full_main_only %*% full_model_coeff_main_only) %*% solve(Inner_prod_y_sqrt)  %*% t(B_spline_basis_fun_y)
  
  # Selected model with interaction
  fit_selected_interaction = (Reg_mat_selected %*% selected_model_coeff) %*% solve(Inner_prod_y_sqrt)  %*% t(B_spline_basis_fun_y)
  
  
  ### Prediction
  # Full model with interaction
  pred_full_interaction = (Reg_mat_full_test %*% full_model_coeff) %*% solve(Inner_prod_y_sqrt) %*% t(B_spline_basis_fun_y)
  # Full model main effects only
  pred_full_main = (Reg_mat_full_main_only_test %*% full_model_coeff_main_only) %*% solve(Inner_prod_y_sqrt) %*% t(B_spline_basis_fun_y)
  # Selected model with interaction
  pred_selected_interaction = (Reg_mat_selected_test %*% selected_model_coeff) %*% solve(Inner_prod_y_sqrt) %*% t(B_spline_basis_fun_y)
  
  return(list("fits_full_interaction" = fit_full_interaction, "predictions_full_interaction" = pred_full_interaction,
              "fits_selected_interaction" = fit_selected_interaction, "predictions_selected_interaction" = pred_selected_interaction,
              "fits_main_only" = fit_full_main, "predictions_main_only" = pred_full_main))
}
