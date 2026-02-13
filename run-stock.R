## ====== 0) Dependencies ======
library(latex2exp)
library(rJava)
library(RMOA)
library(stream)
library(splines)
library(Matrix)
library(lars)
library(graphics)
library(fda)
library(ggplot2)
library(factoextra)
library(elasticnet)
library(SparseM)
library(cluster)
library(gmfd)
library(MASS)
library(pracma)
library(refund)
library(smooth)
library(FRegSigCom)
library(mgcv)
library(corpcor)
library(Rfast)
library(NbClust)
library(stats)
library(rbenchmark)
library(foreach)
library(doSNOW)
library(doParallel)
library(parallel)
library(SOAR)
library(compiler)

setwd("/Users/USER/Downloads/EFFS/organized_code")
source("main_effs.R") 

## ====== 1) Read and Shape Real Data (SINGLE INPUT, SINGLE OUTPUT) ======
t_data <- read.csv('real_data/stock.csv', header = FALSE)
t_data <- t_data[2:nrow(t_data), ]
len_data <- nrow(t_data)

# Reshape a single column time series into a "Samples x Grid Points" matrix and perform [0,1] normalization
get_format_matrix <- function(vec, n_rows = 22, n_cols = 457) {
  stopifnot(length(vec) >= n_rows * n_cols)
  mat <- matrix(as.numeric(vec[1:(n_rows * n_cols)]),
                nrow = n_rows, ncol = n_cols)
  mat <- apply(mat, 2, function(x)
    (x - min(x)) / (max(x) - min(x) + 1e-12))
  t(mat)  # Returns (Sample size = n_cols) x (Grid points = n_rows)
}

input_name  <- "AAPL"
output_name <- "NVDA"

column_indices <- c(
  AAPL  = 2,
  GOOGL = 3,
  MSFT  = 4,
  NVDA  = 5
)

# -------- Construct functional input/output --------
input.ori  <- get_format_matrix(t_data[[ column_indices[input_name] ]])
output.raw <- get_format_matrix(t_data[[ column_indices[output_name] ]])

stopifnot(nrow(input.ori) == nrow(output.raw))

n.row <- nrow(output.raw)

# -------- Log transform (reversible) --------
min_in  <- min(input.ori)
min_out <- min(output.raw)

input  <- if (min_in  < 0) log(input.ori  - min_in  + 1) else log(input.ori  + 1)
output <- if (min_out < 0) log(output.raw - min_out + 1) else log(output.raw + 1)

inv_out <- function(mat_log) {
  if (min_out < 0) exp(mat_log) - 1 + min_out else exp(mat_log) - 1
}

# -------- Dimensions and time grids --------
ncol.input  <- ncol(input)
ncol.output <- ncol(output)

Time.input  <- seq(0, 1, length.out = ncol.input)
Time.output <- seq(0, 1, length.out = ncol.output)

Time.t      <- seq(0, 1, length.out = 1000)
Time.t.len  <- length(Time.t)

# -------- Hyperparameters --------
orders <- c(3, 4, 5)
order.idx <- 2

nbasis.input  <- 10
nbasis.output <- 10

M.clu     <- 50
add.thr   <- 0.5
merge.thr <- 0.9

# -------- Train / test split --------
train_test_split <- 0.5
NN <- floor(train_test_split * n.row)

train_idx <- 1:NN
test_idx  <- setdiff(seq_len(n.row), train_idx)

cat("Running stock experiment:",
    input_name, "->", output_name, "\n")

## ====== 2) Parallel Initialization (Only for FFS) ======
n_cores <- max(1, detectCores()-1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, { source("main_effs.R") })

## ====== 3) Evaluation Tools ======
mse <- function(y, yhat) mean((as.vector(y) - as.vector(yhat))^2)

score_pair <- function(y_log, yhat_log) {
  tr_log  <- mse(y_log[train_idx, , drop=FALSE], yhat_log[train_idx, , drop=FALSE])
  te_log  <- mse(y_log[test_idx,  , drop=FALSE], yhat_log[test_idx,  , drop=FALSE])
  y_tr    <- inv_out(y_log[train_idx,,drop=FALSE])
  y_te    <- inv_out(y_log[test_idx, ,drop=FALSE])
  yh_tr   <- inv_out(yhat_log[train_idx,,drop=FALSE])
  yh_te   <- inv_out(yhat_log[test_idx, ,drop=FALSE])
  tr_real <- mse(y_tr, yh_tr)
  te_real <- mse(y_te, yh_te)
  c(train.func.mse=tr_log, test.func.mse=te_log,
    train.real.mse=tr_real, test.real.mse=te_real)
}

## ====== 4) Run FFS (Multiple pruning strategies) ======
strategies <- c("none","ema")
data.X.Y.data <- cbind(input, output)

# Store rule.monitor trajectories for each pruning strategy
rule_monitor_list <- list()

run_ffs_once <- function(ps) {
  re <- main.STR.FFS.visual(
    add.thr=add.thr, merge.thr=merge.thr,
    Time.t=Time.t, Time.t.len=Time.t.len,
    nbasis.input=nbasis.input, nbasis.output=nbasis.output,
    Time.input.in=Time.input, Time.output.in=Time.output,
    NN=NN, data.X.Y.d=data.X.Y.data,
    orders[order.idx], ncol.input, ncol.output, M.clu,
    pruning.strategy=ps
  )
  
  ## Save the evolution of rule counts for this strategy
  rule_monitor_list[[ps]] <<- re$rule.monitor
  
  data.frame(
    method=paste0("FFS-", ps),
    train.func.mse = re$y.train.func.mse,
    test.func.mse  = re$y.test.func.mse,
    train.real.mse = re$y.train.real.mse,
    test.real.mse  = re$y.test.real.mse,
    train.func.mse.c = re$y.train.func.mse.c,
    test.func.mse.c  = re$y.test.func.mse.c,
    train.real.mse.c = re$y.train.real.mse.c,
    test.real.mse.c  = re$y.test.real.mse.c,
    n.clu          = re$n.clu,
    rule.ave       = re$rule.ave,
    rule.max       = re$rule.max,
    note           = NA_character_,
    stringsAsFactors=FALSE
  )
}

ffs_list <- lapply(strategies, run_ffs_once)
tab_ffs  <- do.call(rbind, ffs_list)

## =========================
## 4b) Plot rule-number trajectories
##     Export PNGs automatically
## =========================
rm_none <- rule_monitor_list[["none"]]
rm_ema  <- rule_monitor_list[["ema"]]

any_diff <- any(rm_none != rm_ema, na.rm = TRUE)
sum_diff <- sum(rm_none != rm_ema, na.rm = TRUE)
cat("Any difference? ", any_diff, " | #positions with different rule counts: ", sum_diff, "\n")

y_range <- c(1, 10)

plot_rule_traj_base <- function(x, y, main_text, file_path, ylim_range) {
  png(file_path, width = 2000, height = 800, res = 200)  # Wide aspect ratio for readability
  on.exit(dev.off(), add = TRUE)
  plot(
    x, y,
    type = "s",
    xlab = "Data arrival index",
    ylab = "Rule numbers",
    ylim = ylim_range,
    main = main_text
  )
}

# Stock: No pruning trajectory
if (!is.null(rm_none)) {
  x1 <- seq_along(rm_none)
  plot_rule_traj_base(
    x = x1, y = rm_none,
    main_text = expression(
      alpha[add] == 0.5 ~ "," ~
        tau[merge] == 0.9),
    file_path = file.path("EEFFS_stock.png"),
    ylim_range = y_range
  )
}

# Stock: EMA pruning trajectory
if (!is.null(rm_ema)) {
  x2 <- seq_along(rm_ema)
  plot_rule_traj_base(
    x = x2, y = rm_ema,
    main_text = expression(
      alpha[add] == 0.5 ~ "," ~
        tau[merge] == 0.9 ~ "," ~
        tau[EMA] == 0.02 ~ "," ~
        beta == 0.95),
    file_path = file.path("EFFS_stock.png"),
    ylim_range = y_range
  )
}

## ====== 5) Baselines (Consistent splitting & evaluation) ======
safe <- function(expr) tryCatch(expr, error=function(e) e)

# 5.1 Interaction-PLS (Beyaztas & Shang)
# Collapse multi-channel X horizontally into a "single function" for pls_fun
run_pls_interaction <- function() {
  src <- safe(source("FLM_interaction-master/auxiliary_functions.R"))
  if (inherits(src, "error") && !exists("pls_fun")) {
    return(data.frame(
      method="PLS-interaction",
      train.func.mse=NA, test.func.mse=NA,
      train.real.mse=NA, test.real.mse=NA,
      train.func.mse.c = NA, test.func.mse.c = NA,
      train.real.mse.c = NA, test.real.mse.c = NA,
      n.clu=NA, rule.ave=NA, rule.max=NA,
      note=paste0("load err: ", src$message)
    ))
  }
  if (!exists("pls_fun")) {
    return(data.frame(
      method="PLS-interaction",
      train.func.mse=NA, test.func.mse=NA,
      train.real.mse=NA, test.real.mse=NA,
      train.func.mse.c = NA, test.func.mse.c = NA,
      train.real.mse.c = NA, test.real.mse.c = NA,
      n.clu=NA, rule.ave=NA, rule.max=NA,
      note="pls_fun not found"
    ))
  }
  
  # ---- Concatenate multiple channels into a single functional predictor (per author's demo) ----
  # ---- SINGLE INPUT: directly use input ----
  X.list0 <- list(input)
  
  # —— Training/Testing split
  X.train.list <- list(X.list0[[1]][train_idx, , drop = FALSE])
  X.test.list  <- list(X.list0[[1]][test_idx,  , drop = FALSE])
  
  # —— Response (log space)
  Y.train <- output[train_idx, , drop = FALSE]
  Y.test  <- output[test_idx,  , drop = FALSE]
  
  # ---- Basis safety valve: nbasis cannot exceed grid length ----
  grid_len_x <- ncol(X.train.list[[1]])
  grid_len_y <- ncol(Y.train)
  nbf_x_eff  <- min(nbasis.input,  grid_len_x - 1L)
  nbf_y_eff  <- min(nbasis.output, grid_len_y - 1L)
  if (nbf_x_eff < 2L) nbf_x_eff <- 2L
  if (nbf_y_eff < 2L) nbf_y_eff <- 2L
  norder_eff <- orders[order.idx]
  
  # ---- Training set prediction (in-sample) ----
  fit_tr <- safe(pls_fun(
    response          = Y.train,
    predictors_train = X.train.list,
    predictors_test  = X.train.list,
    nbf_x             = nbf_x_eff,
    nbf_y             = nbf_y_eff,
    n.order           = norder_eff
  ))
  if (inherits(fit_tr, "error")) {
    return(data.frame(
      method="PLS-interaction",
      train.func.mse=NA, test.func.mse=NA,
      train.real.mse=NA, test.real.mse=NA,
      train.func.mse.c = NA, test.func.mse.c = NA,
      train.real.mse.c = NA, test.real.mse.c = NA,
      n.clu=NA, rule.ave=NA, rule.max=NA,
      note=paste0("train predict err: ", fit_tr$message,
                  " | collapsed_to_single_X, nbf_x=", nbf_x_eff, ", nbf_y=", nbf_y_eff)
    ))
  }
  
  P.tr.full <- fit_tr$predictions_full_interaction
  P.tr.sel  <- fit_tr$predictions_selected_interaction
  P.tr.main <- fit_tr$predictions_main_only
  
  # ---- Test set prediction (out-of-sample) ----
  fit_te <- safe(pls_fun(
    response          = Y.train,
    predictors_train = X.train.list,
    predictors_test  = X.test.list,
    nbf_x             = nbf_x_eff,
    nbf_y             = nbf_y_eff,
    n.order           = norder_eff
  ))
  if (inherits(fit_te, "error")) {
    return(data.frame(
      method="PLS-interaction",
      train.func.mse=NA, test.func.mse=NA,
      train.real.mse=NA, test.real.mse=NA,
      train.func.mse.c = NA, test.func.mse.c = NA,
      train.real.mse.c = NA, test.real.mse.c = NA,
      n.clu=NA, rule.ave=NA, rule.max=NA,
      note=paste0("test predict err: ", fit_te$message,
                  " | collapsed_to_single_X, nbf_x=", nbf_x_eff, ", nbf_y=", nbf_y_eff)
    ))
  }
  
  P.te.full <- fit_te$predictions_full_interaction
  P.te.sel  <- fit_te$predictions_selected_interaction
  P.te.main <- fit_te$predictions_main_only
  
  eval_one <- function(name, Ptr, Pte) {
    tr_log <- mse(Y.train, Ptr)
    te_log <- mse(Y.test , Pte)
    if (is.function(inv_out)) {
      tr_real <- mse(inv_out(Y.train), inv_out(Ptr))
      te_real <- mse(inv_out(Y.test ), inv_out(Pte))
    } else {
      tr_real <- NA_real_; te_real <- NA_real_
    }
    data.frame(
      method=name,
      train.func.mse=tr_log, test.func.mse=te_log,
      train.real.mse=tr_real, test.real.mse=te_real,
      train.func.mse.c = NA, test.func.mse.c = NA,
      train.real.mse.c = NA, test.real.mse.c = NA,
      n.clu=NA, rule.ave=NA, rule.max=NA,
      note=paste0("collapsed_to_single_X, nbf_x=", nbf_x_eff, ", nbf_y=", nbf_y_eff),
      stringsAsFactors = FALSE
    )
  }
  
  rbind(
    eval_one("PLS-interaction(full)",      P.tr.full, P.te.full),
    eval_one("PLS-interaction(selected)", P.tr.sel , P.te.sel ),
    eval_one("PLS-main-only",              P.tr.main, P.te.main)
  )
}

# 5.2 fAPLS (Zhiyang Zhou): Collapsing multi-channel X horizontally
run_fapls <- function() {
  src <- safe(source("fAPLS-master/functions.R"))
  if (inherits(src,"error")) return(data.frame(
    method="fAPLS", train.func.mse=NA, test.func.mse=NA,
    train.real.mse=NA, test.real.mse=NA, train.func.mse.c = NA, test.func.mse.c = NA,
    train.real.mse.c = NA, test.real.mse.c = NA, n.clu=NA, rule.ave=NA, rule.max=NA,
    note=paste0("load err: ", src$message)
  ))
  Xc <- input   # Input is already a 5*48 concatenation
  Yc <- output
  dom.x <- seq(0,1,length.out=ncol(Xc))
  dom.y <- seq(0,1,length.out=ncol(Yc))
  ridge.par <- 0; K <- 10
  p.max <- safe(pUpper.compu(Xc, dom.x, .99, basis.name="bspline", K=K))
  if (inherits(p.max,"error")) return(data.frame(
    method="fAPLS", train.func.mse=NA, test.func.mse=NA, train.real.mse=NA, test.real.mse=NA,train.func.mse.c = NA, test.func.mse.c = NA,
    train.real.mse.c = NA, test.real.mse.c = NA, n.clu=NA, rule.ave=NA, rule.max=NA, note=paste0("pUpper err: ", p.max$message)
  ))
  
  tr <- safe(fAPLS(Xc[train_idx,], Yc[train_idx,], dom.x, dom.y,
                   Xc[train_idx,], Yc[train_idx,], ridge.par, p.max, tune='CV'))
  te <- safe(fAPLS(Xc[train_idx,], Yc[train_idx,], dom.x, dom.y,
                   Xc[test_idx, ], Yc[test_idx, ], ridge.par, p.max, tune='CV'))
  if (inherits(tr,"error") || inherits(te,"error")) return(data.frame(
    method="fAPLS", train.func.mse=NA, test.func.mse=NA, train.real.mse=NA, test.real.mse=NA,train.func.mse.c = NA, test.func.mse.c = NA,
    train.real.mse.c = NA, test.real.mse.c = NA,n.clu=NA, rule.ave=NA, rule.max=NA,
    note=paste("train/test err:",
               if(inherits(tr,"error")) tr$message else "",
               if(inherits(te,"error")) te$message else "")
  ))
  
  Ytr <- output[train_idx,,drop=FALSE]
  Yte <- output[test_idx, ,drop=FALSE]
  
  data.frame(
    method            = "fAPLS",
    train.func.mse    = mse(Ytr, tr$y.pred),
    test.func.mse     = mse(Yte, te$y.pred),
    train.real.mse    = mse(inv_out(Ytr), inv_out(tr$y.pred)),
    test.real.mse     = mse(inv_out(Yte), inv_out(te$y.pred)),
    train.func.mse.c  = NA_real_,
    test.func.mse.c   = NA_real_,
    train.real.mse.c  = NA_real_,
    test.real.mse.c   = NA_real_,
    n.clu             = NA_real_,
    rule.ave          = NA_real_,
    rule.max          = NA_real_,
    note              = NA_character_,
    check.names       = FALSE,
    stringsAsFactors  = FALSE
  )
}

# ---- Patch: missing helpers expected by Bernardi code (global scope) ----
if (!exists("vec", mode = "function", inherits = TRUE)) {
  vec <- function(x, byrow = FALSE) if (byrow) as.vector(t(x)) else as.vector(x)
}
if (!exists("invvec", mode = "function", inherits = TRUE)) {
  invvec <- function(v, nrow, ncol) matrix(v, nrow = nrow, ncol = ncol)
}

# 5.3 Smooth LASSO (Centofanti et al.) —— refund/fda implementation
run_slasso <- function() {
  if (!requireNamespace("slasso", quietly = TRUE)) {
    return(data.frame(
      method="SLASSO", train.func.mse=NA, test.func.mse=NA,
      train.real.mse=NA, test.real.mse=NA,train.func.mse.c = NA, test.func.mse.c = NA,
      train.real.mse.c = NA, test.real.mse.c = NA,
      n.clu=NA, rule.ave=NA, rule.max=NA,
      note="pkg slasso not installed"
    ))
  }
  
  X <- t(input)    # s_grid x n
  Y <- t(output)   # t_grid x n (Note: log space)
  domain <- c(0,1); n.basis <- 10; n.order <- 4
  grid.s <- seq(0,1,length.out=nrow(X))
  grid.t <- seq(0,1,length.out=nrow(Y))
  basis_s <- fda::create.bspline.basis(domain, nbasis=n.basis, norder=n.order)
  basis_t <- fda::create.bspline.basis(domain, nbasis=n.basis, norder=n.order)
  
  X_fd <- fda::smooth.basis(grid.s, X, basis_s)$fd
  Y_fd <- fda::smooth.basis(grid.t, Y, basis_t)$fd
  
  X_fd_tr <- X_fd[train_idx]
  Y_fd_tr <- Y_fd[train_idx]
  X_fd_te <- X_fd[test_idx]
  Y_fd_te <- Y_fd[test_idx]
  
  # Cross-validation to select lambda
  cv <- try(slasso::slasso.fr_cv(
    Y_fd = Y_fd_tr, X_fd = X_fd_tr,
    basis_s = basis_s, basis_t = basis_t,
    K = 5,
    lambda_L_vec = 10^seq(-2, 2, by=1),
    lambda_s_vec  = 10^seq(-6, -2, by=1),
    lambda_t_vec  = 10^seq(-6, -2, by=1),
    ncores = 1, invisible = 1, max_iterations = 2000
  ), silent = TRUE)
  
  if (inherits(cv, "try-error")) {
    return(data.frame(
      method="SLASSO", train.func.mse=NA, test.func.mse=NA,
      train.real.mse=NA, test.real.mse=NA,train.func.mse.c = NA, test.func.mse.c = NA,
      train.real.mse.c = NA, test.real.mse.c = NA,
      n.clu=NA, rule.ave=NA, rule.max=NA,
      note=paste0("CV err: ", as.character(cv))
    ))
  }
  
  lam <- cv$lambda_opt_vec
  fit <- try(slasso::slasso.fr(
    Y_fd = Y_fd_tr, X_fd = X_fd_tr,
    basis_s = basis_s, basis_t = basis_t,
    lambda_L = lam[1], lambda_s = lam[2], lambda_t = lam[3],
    B0 = NULL, invisible = 1, max_iterations = 10000
  ), silent = TRUE)
  
  if (inherits(fit, "try-error")) {
    return(data.frame(
      method="SLASSO", train.func.mse=NA, test.func.mse=NA,
      train.real.mse=NA, test.real.mse=NA,train.func.mse.c = NA, test.func.mse.c = NA,
      train.real.mse.c = NA, test.real.mse.c = NA,
      n.clu=NA, rule.ave=NA, rule.max=NA,
      note=paste0("fit err: ", as.character(fit))
    ))
  }
  
  # Prediction in log space: y_hat(t) = alpha(t) + ∫ x(s) * beta(s,t) ds
  delta <- (max(grid.s)-min(grid.s)) / (length(grid.s)-1)
  Bmat  <- fda::eval.bifd(grid.s, grid.t, fit$Beta_hat_fd)  # |S| x |T|
  alpha_t <- as.numeric(fda::eval.fd(grid.t, fit$alpha))     # |T|
  
  Xtr_eval <- t(fda::eval.fd(grid.s, X_fd_tr))               # n_tr x |S|
  Xte_eval <- t(fda::eval.fd(grid.s, X_fd_te))               # n_te x |S|
  
  Yhat_tr_log <- Xtr_eval %*% Bmat * delta + matrix(alpha_t, nrow(Xtr_eval), length(alpha_t), byrow=TRUE)
  Yhat_te_log <- Xte_eval %*% Bmat * delta + matrix(alpha_t, nrow(Xte_eval), length(alpha_t), byrow=TRUE)
  
  # Align with evaluation protocol
  Ytr_log <- output[train_idx,,drop=FALSE]
  Yte_log <- output[test_idx, ,drop=FALSE]
  
  tr_log  <- mse(Ytr_log, Yhat_tr_log)
  te_log  <- mse(Yte_log, Yhat_te_log)
  
  Ytr_real_hat <- inv_out(Yhat_tr_log)
  Yte_real_hat <- inv_out(Yhat_te_log)
  Ytr_real     <- inv_out(Ytr_log)
  Yte_real     <- inv_out(Yte_log)
  
  tr_real <- mse(Ytr_real, Ytr_real_hat)
  te_real <- mse(Yte_real, Yte_real_hat)
  
  data.frame(method="SLASSO",
             train.func.mse=tr_log, test.func.mse=te_log,
             train.real.mse=tr_real, test.real.mse=te_real,
             train.func.mse.c = NA, test.func.mse.c = NA,
             train.real.mse.c = NA, test.real.mse.c = NA,
             n.clu=NA, rule.ave=NA, rule.max=NA, note=NA)
}

## ====== 6) Execute Baselines ======
tab_pls     <- run_pls_interaction()
tab_fapls   <- run_fapls()
tab_slasso  <- run_slasso()

## ====== 7) Summary and Persistence ======
# --- Align all result tables by column names and order ---
ensure_cols <- function(df) {
  need <- c(
    "method",
    "train.func.mse","test.func.mse","train.real.mse","test.real.mse",
    "train.func.mse.c","test.func.mse.c","train.real.mse.c","test.real.mse.c",
    "n.clu","rule.ave","rule.max","note"
  )
  for (nm in need) if (!nm %in% names(df)) df[[nm]] <- NA
  df[need]
}

tab_ffs     <- ensure_cols(tab_ffs)
tab_pls     <- ensure_cols(tab_pls)
tab_fapls   <- ensure_cols(tab_fapls)
tab_slasso  <- ensure_cols(tab_slasso)

tab_all <- rbind(
  transform(tab_ffs,   method=as.character(method)),
  transform(tab_pls,   method=as.character(method)),
  transform(tab_fapls, method=as.character(method)),
  transform(tab_slasso,method=as.character(method))
)

tab_all$train.func.mse <- as.numeric(tab_all$train.func.mse)
tab_all$test.func.mse  <- as.numeric(tab_all$test.func.mse)
tab_all$train.real.mse <- as.numeric(tab_all$train.real.mse)
tab_all$test.real.mse  <- as.numeric(tab_all$test.real.mse)
tab_all$train.func.mse.c <- as.numeric(tab_all$train.func.mse.c)
tab_all$test.func.mse.c  <- as.numeric(tab_all$test.func.mse.c)
tab_all$train.real.mse.c <- as.numeric(tab_all$train.real.mse.c)
tab_all$test.real.mse.c  <- as.numeric(tab_all$test.real.mse.c)
tab_all$rule.ave       <- as.numeric(tab_all$rule.ave)
tab_all$rule.max       <- as.numeric(tab_all$rule.max)
tab_all$n.clu          <- as.numeric(tab_all$n.clu)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
csv_name  <- paste0("FoF_stock_compare_run_", train_test_split,"_m",merge.thr,"_a",add.thr, "_t",timestamp, ".csv")
write.csv(tab_all, csv_name, row.names = FALSE)

cat("==== Summary (rounded) ====\n")
print(within(tab_all, {
  train.func.mse <- round(train.func.mse, 6)
  test.func.mse  <- round(test.func.mse,  6)
  train.real.mse <- round(train.real.mse, 6)
  test.real.mse  <- round(test.real.mse,  6)
  rule.ave       <- round(rule.ave, 3)
}), row.names=FALSE)
cat("Saved to:", csv_name, "\n")

## Cleanup Parallel Cluster
stopCluster(cl)