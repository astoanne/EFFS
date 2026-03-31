## ====== 0) ä¾èµ– ======
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

resolve_candidate_dirs <- function() {
  dirs <- character()
  add_dir <- function(path_value) {
    if (is.null(path_value) || length(path_value) == 0) {
      return()
    }
    path_value <- path_value[1]
    if (is.na(path_value) || !nzchar(path_value)) {
      return()
    }
    dirs <<- c(
      dirs,
      dirname(normalizePath(path_value, winslash = "/", mustWork = FALSE))
    )
  }

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_match <- grep(file_arg, args)
  if (length(file_match) > 0) {
    add_dir(sub(file_arg, "", args[file_match[1]]))
  }

  source_calls <- Filter(function(call_obj) {
    is.call(call_obj) && identical(call_obj[[1]], as.name("source"))
  }, sys.calls())
  for (call_obj in source_calls) {
    if (length(call_obj) >= 2 && is.character(call_obj[[2]])) {
      add_dir(as.character(call_obj[[2]]))
    }
  }

  for (env in sys.frames()) {
    if (!is.null(env$ofile)) {
      add_dir(env$ofile)
    }
    if (!is.null(env$srcfile) && !is.null(env$srcfile$filename)) {
      add_dir(env$srcfile$filename)
    }
  }
  cwd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  known_runner_paths <- c(
    file.path(cwd, "run_stock_paper.R"),
    file.path(cwd, "run_agriculture_paper.R"),
    file.path(cwd, "organized_code2", "run_stock_paper.R"),
    file.path(cwd, "organized_code2", "run_agriculture_paper.R"),
    file.path(cwd, "EFFS", "organized_code2", "run_stock_paper.R"),
    file.path(cwd, "EFFS", "organized_code2", "run_agriculture_paper.R")
  )
  for (path_value in known_runner_paths) {
    if (file.exists(path_value)) {
      add_dir(path_value)
    }
  }
  dirs <- c(
    dirs,
    cwd,
    dirname(cwd),
    normalizePath(file.path(cwd, "organized_code2"), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(cwd, "..", "organized_code2"), winslash = "/", mustWork = FALSE)
  )
  unique(dirs[nzchar(dirs)])
}

find_existing_file <- function(filename, dirs) {
  for (dir_path in dirs) {
    candidate <- file.path(dir_path, filename)
    if (file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  NA_character_
}

find_existing_dir <- function(dirname_hint, dirs) {
  for (dir_path in dirs) {
    candidate <- file.path(dir_path, dirname_hint)
    if (dir.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  NA_character_
}

candidate_dirs <- resolve_candidate_dirs()
main_file <- find_existing_file("main_STR_FFS_more_pruning.R", candidate_dirs)
if (is.na(main_file)) {
  stop(
    "Could not locate main_STR_FFS_more_pruning.R. Checked:\n",
    paste(file.path(candidate_dirs, "main_STR_FFS_more_pruning.R"), collapse = "\n")
  )
}
script_dir <- dirname(main_file)
project_root <- script_dir
data_dir <- find_existing_dir("real_data", c(script_dir, candidate_dirs))
flm_dir <- find_existing_dir("FLM_interaction-master", c(script_dir, candidate_dirs))
fapls_dir <- find_existing_dir("fAPLS-master", c(script_dir, candidate_dirs))
if (is.na(data_dir) || is.na(flm_dir) || is.na(fapls_dir)) {
  stop(
    "organized_code2 is missing required local folders.\n",
    "real_data: ", data_dir, "\n",
    "FLM_interaction-master: ", flm_dir, "\n",
    "fAPLS-master: ", fapls_dir, "\n",
    "Make sure this folder includes its bundled datasets and helper code."
  )
}
output_dir <- file.path(script_dir, "outputs", "agriculture")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message("Using EEFFS core: ", main_file)
message("Using data directory: ", data_dir)
message("Using PLS helper directory: ", flm_dir)
message("Using fAPLS helper directory: ", fapls_dir)
message("Switching working directory to: ", project_root)
message("Writing outputs to: ", output_dir)
setwd(project_root)
source(main_file)

## ====== 1) è¯»å–ä¸Žæˆå½¢ï¼ˆç»Ÿä¸€å£å¾„, SINGLE INPUT) ======
t_data <- read.csv(file.path(data_dir, "north_dakota_agricultural.csv"), header = FALSE)
t_data <- t_data[6:nrow(t_data), ]
len_data <- nrow(t_data)

get_format_matrix <- function(data) {
  mat <- matrix(as.numeric(data), nrow = 365, ncol = 91)
  mat <- t(mat)
  mat <- apply(mat, 2, function(x) {
    (x - min(x)) / (max(x) - min(x))
  })
  return(mat)
}

# -------- Choose ONE input and ONE output --------
input_name  <- "solar_radiation"
output_name <- "bare_soil_temp"

column_indices <- c(
  wind_speed       = 12,
  wind_chill       = 16,
  air_temp         = 8,
  solar_radiation  = 14,
  bare_soil_temp   = 10
)

# -------- Construct functional input/output --------
input.ori  <- get_format_matrix(t_data[[ column_indices[input_name] ]])
output.raw <- get_format_matrix(t_data[[ column_indices[output_name] ]])

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
add.thr  <- 0.5
merge.thr<- 0.9

# -------- Train / test split --------
train_test_split <- 0.5
NN <- floor(train_test_split * n.row)

train_idx <- 1:NN
test_idx  <- setdiff(seq_len(n.row), train_idx)

cat("Running SISO experiment:", input_name, "->", output_name, "\n")

## ====== 2) å¹¶è¡Œï¼ˆä»…ç”¨äºŽFFSï¼‰ ======
n_cores <- max(1, detectCores()-1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterExport(cl, "main_file")
clusterEvalQ(cl, { source(main_file) })

## ====== 3) ç»Ÿä¸€è¯„æµ‹å·¥å…· ======
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

## ====== 4) è·‘ FFSï¼ˆå¤šç§ pruningï¼‰ ======
# strategies <- c("none","usage","age","importance","utility","ema")
strategies <- c("none","ema")
data.X.Y.data <- cbind(input, output)

# å­˜æ¯ç§ pruning ç­–ç•¥çš„ rule.monitor è½¨è¿¹
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
  
  ## ä¿å­˜è¯¥ç­–ç•¥ä¸‹çš„è§„åˆ™æ•°æ¼”åŒ–
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

## ====== 4b) ç”»è§„åˆ™æ•°éšæ ·æœ¬æ•°å˜åŒ–ï¼ˆnone vs emaï¼Œå„ä¸€å¼ å›¾ï¼‰ ======
rm_none <- rule_monitor_list[["none"]]
rm_ema  <- rule_monitor_list[["ema"]]

any_diff <- any(rm_none != rm_ema, na.rm = TRUE)
sum_diff <- sum(rm_none != rm_ema, na.rm = TRUE)
cat("Any difference? ", any_diff, " | #positions with different rule counts: ", sum_diff, "\n")


# ä¸ºäº†å¯æ¯”æ€§ï¼Œå…ˆç®—ä¸€ä¸ªç»Ÿä¸€çš„ y è½´èŒƒå›´
# y_range <- range(c(rm_none, rm_ema), na.rm = TRUE)
y_range <- c(1, 10)
plot_rule_traj_base <- function(x, y, main_text, file_path, ylim_range) {
  png(file_path, width = 2000, height = 800, res = 200)  # wide, readable
  on.exit(dev.off(), add = TRUE)
  plot(
    x, y,
    type = "s",
    xlab = "Data arrival index",
    # xlab = "Normalized data arrival",
    ylab = "Rule numbers",
    ylim = ylim_range,
    main = main_text
  )
}

#  no pruning
# if (!is.null(rm_none)) {
#   # x1 <- seq_along(rm_none)
#   x1 <- seq(0, 1, length.out = length(rm_none))
#   plot_rule_traj_base(
#     x = x1, y = rm_none,
#     main_text = paste0("Îµ_a = ", add.thr, ", Îµ_m = ", merge.thr, " (no pruning)"),
#     file_path = file.path("noprune_agriculture_old.png"),
#     ylim_range = y_range
#   )
# }
if (!is.null(rm_none)) {
  x1 <- seq_along(rm_none)
  plot_rule_traj_base(
    x = x1, y = rm_none,
    main_text = expression(
      alpha[add] == 0.5 ~ "," ~
        tau[merge] == 0.9),
    file_path = file.path(output_dir, "noprune_agriculture.png"),
    ylim_range = y_range
  )
}

# EMA pruning
if (!is.null(rm_ema)) {
  x2 <- seq_along(rm_ema)
  plot_rule_traj_base(
    x = x2, y = rm_ema,
    main_text = expression(
      alpha[add] == 0.5 ~ "," ~
        tau[merge] == 0.9 ~ "," ~
        tau[EMA] == 0.02 ~ "," ~
        beta == 0.95),
    file_path = file.path(output_dir, "prune_agriculture.png"),
    ylim_range = y_range
  )
}

## ====== 5) Baselinesï¼ˆç»Ÿä¸€åˆ‡åˆ† & å£å¾„ï¼‰ ======
safe <- function(expr) tryCatch(expr, error=function(e) e)

# 5.1 Interaction-PLSï¼ˆBeyaztas & Shangï¼‰â€”â€”æŠŠå¤šé€šé“ X æ¨ªå‘æ‹¼æˆâ€œå•å‡½æ•°â€å–‚ pls_fun
run_pls_interaction <- function() {
  src <- safe(source(file.path(flm_dir, "auxiliary_functions.R")))
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
  
  # ---- æŠŠå¤šé€šé“æŒ‰åˆ—æ‹¼æŽ¥æˆä¸€ä¸ªåŠŸèƒ½è‡ªå˜é‡ï¼ˆä¸Žä½œè€…demoä¸€è‡´ï¼šX.list[[1]]ï¼‰----
  # ---- SINGLE INPUT: directly use input ----
  X.list0 <- list(input)
  
  
  # â€”â€” è®­ç»ƒ/æµ‹è¯•åˆ‡åˆ†
  X.train.list <- list(X.list0[[1]][train_idx, , drop = FALSE])
  X.test.list  <- list(X.list0[[1]][test_idx,  , drop = FALSE])
  
  # â€”â€” å“åº”ï¼ˆlog ç©ºé—´ï¼‰
  Y.train <- output[train_idx, , drop = FALSE]
  Y.test  <- output[test_idx,  , drop = FALSE]
  
  # ---- åŸºæ•°å®‰å…¨é˜€ï¼šæ ·æ¡åŸºæ•°ä¸èƒ½è¶…è¿‡ç½‘æ ¼é•¿åº¦ ----
  grid_len_x <- ncol(X.train.list[[1]])
  grid_len_y <- ncol(Y.train)
  nbf_x_eff  <- min(nbasis.input,  grid_len_x - 1L)
  nbf_y_eff  <- min(nbasis.output, grid_len_y - 1L)
  if (nbf_x_eff < 2L) nbf_x_eff <- 2L
  if (nbf_y_eff < 2L) nbf_y_eff <- 2L
  norder_eff <- orders[order.idx]
  
  # ---- è®­ç»ƒé›†é¢„æµ‹ï¼ˆin-sampleï¼‰----
  fit_tr <- safe(pls_fun(
    response         = Y.train,
    predictors_train = X.train.list,
    predictors_test  = X.train.list,
    nbf_x            = nbf_x_eff,
    nbf_y            = nbf_y_eff,
    n.order          = norder_eff
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
  
  # ---- æµ‹è¯•é›†é¢„æµ‹ï¼ˆout-of-sampleï¼‰----
  fit_te <- safe(pls_fun(
    response         = Y.train,
    predictors_train = X.train.list,
    predictors_test  = X.test.list,
    nbf_x            = nbf_x_eff,
    nbf_y            = nbf_y_eff,
    n.order          = norder_eff
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
    eval_one("PLS-interaction(full)",     P.tr.full, P.te.full),
    eval_one("PLS-interaction(selected)", P.tr.sel , P.te.sel ),
    eval_one("PLS-main-only",             P.tr.main, P.te.main)
  )
}


# 5.2 fAPLSï¼ˆZhiyang Zhouï¼‰ï¼šæŠŠå¤šé€šé“ X æ¨ªå‘æ‹¼æŽ¥æˆâ€œå•å‡½æ•°â€
run_fapls <- function() {
  src <- safe(source(file.path(fapls_dir, "functions.R")))
  if (inherits(src,"error")) return(data.frame(
    method="fAPLS", train.func.mse=NA, test.func.mse=NA,
    train.real.mse=NA, test.real.mse=NA, train.func.mse.c = NA, test.func.mse.c = NA,
    train.real.mse.c = NA, test.real.mse.c = NA, n.clu=NA, rule.ave=NA, rule.max=NA,
    note=paste0("load err: ", src$message)
  ))
  Xc <- input   # å·²ç»æ˜¯ 5*48 çš„æ‹¼æŽ¥
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

# 5.3 Smooth LASSOï¼ˆCentofanti et al.ï¼‰â€”â€” refund/fda æ–¹å¼
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
  Y <- t(output)   # t_grid x n   (æ³¨æ„ï¼šè¿™é‡Œæ˜¯logç©ºé—´)
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
  
  # äº¤å‰éªŒè¯é€‰æ‹© lambdaï¼ˆç”¨å®˜æ–¹å¯¼å‡ºçš„ slasso.fr_cvï¼‰
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
  
  # åœ¨ log ç©ºé—´åšé¢„æµ‹ï¼š y_hat(t) = alpha(t) + âˆ« x(s) * beta(s,t) ds
  delta <- (max(grid.s)-min(grid.s)) / (length(grid.s)-1)
  Bmat  <- fda::eval.bifd(grid.s, grid.t, fit$Beta_hat_fd)  # |S| x |T|
  alpha_t <- as.numeric(fda::eval.fd(grid.t, fit$alpha))     # |T|
  
  Xtr_eval <- t(fda::eval.fd(grid.s, X_fd_tr))               # n_tr x |S|
  Xte_eval <- t(fda::eval.fd(grid.s, X_fd_te))               # n_te x |S|
  
  Yhat_tr_log <- Xtr_eval %*% Bmat * delta + matrix(alpha_t, nrow(Xtr_eval), length(alpha_t), byrow=TRUE)
  Yhat_te_log <- Xte_eval %*% Bmat * delta + matrix(alpha_t, nrow(Xte_eval), length(alpha_t), byrow=TRUE)
  
  # å’Œä½ çš„è¯„æµ‹å£å¾„å¯¹é½
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

# 5.4 OVGLASSOï¼ˆBernardi et al.ï¼‰ï¼šä¹ŸæŠŠå¤šé€šé“æ‹¼æˆå•å‡½æ•°
# 5.4 OVGLASSOï¼ˆBernardi et al.ï¼‰ï¼šæŠŠå¤šé€šé“æ‹¼æˆå•å‡½æ•° â€”â€” è‡ªé€‚åº”é€šé“ä¸Žç½‘æ ¼é•¿åº¦
run_ovglasso <- function() {
  need <- c("21266123/functions.R","21266123/utilfuns.R","21266123/linreg_LASSO.R",
            "21266123/f2freg_LASSO.R","21266123/f2freg_RIDGE.R","21266123/linreg_RIDGE.R",
            "21266123/f2freg_OVGLASSO.R")
  for (f in need) if (!file.exists(f)) return(data.frame(
    method="OVGLASSO", train.func.mse=NA, test.func.mse=NA,
    train.real.mse=NA, test.real.mse=NA, n.clu=NA, rule.ave=NA, rule.max=NA,
    note=paste0("missing file: ", f)))
  for (f in need) source(f)
  
  Xc <- input   # n Ã— (n_channels * L)
  Yc <- output  # n Ã— r
  n      <- nrow(Xc)
  P_full <- ncol(Xc)
  r      <- ncol(Yc)
  
  # ---- ä»Žå…¨å±€å¯¹è±¡è‡ªåŠ¨æŽ¨æ–­é€šé“æ•°å’Œæ¯é€šé“é•¿åº¦ ----
  # è‹¥å·²æœ‰ ncol.in.eachï¼ˆç”±ä½ å‰é¢ä»£ç ç®—å‡ºï¼‰ï¼Œç”¨å®ƒï¼›å¦åˆ™å›žé€€åˆ°ç­‰åˆ†çŒœæµ‹
  if (exists("ncol.in.each")) {
    stopifnot(length(ncol.in.each) >= 1)
    L <- as.integer(ncol.in.each[1])
    if (!all(ncol.in.each == L))
      stop("å„é€šé“é•¿åº¦ä¸ä¸€è‡´ï¼Œå½“å‰å®žçŽ°è¦æ±‚ç›¸åŒé•¿åº¦ã€‚")
    n_channels <- as.integer(P_full / L)
    if (n_channels * L != P_full)
      stop("P_full ä¸Ž ncol.in.each ä¸æ•´é™¤ï¼Œæ£€æŸ¥ input æž„é€ ã€‚")
  } else {
    # å…œåº•ï¼šå°è¯•ç”¨ r çš„å› æ•°åŽ»çŒœæµ‹ï¼ˆä¸å¸¸ç”¨è·¯å¾„ï¼‰
    facs <- which(P_full %% (1:200) == 0)
    # é€‰ä¸€ä¸ªä¸Ž r æŽ¥è¿‘çš„ Lï¼ˆå…œåº•ï¼Œä¸èµ°è¿™ï¼‰
    L <- (1:200)[facs][which.min(abs((1:200)[facs] - r))]
    n_channels <- P_full / L
  }
  
  train_idx <- 1:NN
  test_idx  <- setdiff(seq_len(n), train_idx)
  
  x_train <- Xc[train_idx, , drop=FALSE]  # (n_tr Ã— P_full)
  y_train <- Yc[train_idx, , drop=FALSE]  # (n_tr Ã— r)
  x_test  <- Xc[test_idx,  , drop=FALSE]
  y_test  <- Yc[test_idx,  , drop=FALSE]
  
  # ---- åŸºå‡½æ•°ï¼ˆå¯¹ t æ–¹å‘ q ä¸ªæ ·æ¡ï¼›å¯¹æ¯ä¸ªé€šé“çš„ s æ–¹å‘ p_per_ch ä¸ªæ ·æ¡ï¼‰----
  q <- 20L
  grid_t <- seq(0, 1, length.out = r)
  phi <- fda::eval.basis(grid_t, fda::create.bspline.basis(c(0,1), nbasis=q, norder=4))  # r Ã— q
  
  p_per_ch <- 4L
  p <- as.integer(n_channels * p_per_ch)
  grid_s <- seq(0, 1, length.out = L)
  psi_ch <- fda::eval.basis(grid_s, fda::create.bspline.basis(c(0,1), nbasis=p_per_ch, norder=4)) # L Ã— p_per_ch
  # Psi ç»´åº¦ï¼š(n_channels*L) Ã— (n_channels*p_per_ch) = P_full Ã— p
  Psi <- as.matrix(Matrix::bdiag(replicate(n_channels, psi_ch, simplify = FALSE)))
  
  # ---------- å¼ºåˆ¶ä¸º double çŸ©é˜µï¼Œé¿å… Matrix ç±»æ–¹æ³•åˆ†æ´¾ ----------
  to_dmat <- function(M) { M <- as.matrix(M); storage.mode(M) <- "double"; M }
  x_train <- to_dmat(x_train); x_test <- to_dmat(x_test)
  y_train <- to_dmat(y_train); y_test <- to_dmat(y_test)
  Psi     <- to_dmat(Psi);     phi    <- to_dmat(phi)
  
  # ---------- é¢„æŠ•å½±ï¼ˆC ç³»åˆ—é…ç½®è¦ç”¨ï¼‰ ----------
  Xs_train <- x_train %*% Psi        # n_tr Ã— p
  Xs_test  <- x_test  %*% Psi        # n_te Ã— p
  I_p      <- diag(1.0, nrow = p, ncol = p)
  
  ## ---- lambda åºåˆ—ï¼ˆä½ åŽŸæœ¬çš„ï¼‰----
  lam.ret <- try(getlambda_OVGLASSO(
    y = matrix(as.vector(t(y_train)), ncol = 1),
    X = base::kronecker(phi, x_train %*% Psi),  # ç”¨æŠ•å½±åŽçš„ Xs ä¼° Î»
    p = p, q = q, splOrd = 4, method = "ogl&1",
    lambda.min = 1e-5, maxl = 30
  ), silent = TRUE)
  lam <- if (inherits(lam.ret, "try-error")) exp(seq(log(1e-4), log(1e+1), length.out = 30)) else lam.ret$lambda.seq * 40
  
  ## ---- ç»´åº¦æ–­è¨€ï¼ˆå…³é”®å››æ¡ï¼‰----
  stopifnot(
    ncol(x_train) == ncol(t(Psi)),     # P_full == P_full
    nrow(t(Psi)) == p,                 # rows of mH = p (X-basis)
    nrow(t(phi))  == q,                # rows of mTheta = q (Y-basis)
    ncol(t(phi))  == r                 # cols of mTheta = r
  )
  
  ## ---- æ­£ç¡®çš„ API æ˜ å°„ï¼ˆæ¥è‡ªæºç ï¼štcrossprod ä¸Ž kronecker çš„ç»´åº¦çº¦æŸï¼‰----
  fit <- try(f2freg_OVGLASSO(
    mY = as.matrix(y_train),
    mX = as.matrix(x_train),
    mTheta = t(as.matrix(phi)),   # q Ã— r
    mH = t(as.matrix(Psi)),       # p Ã— P_full
    vRegP_INIT = NULL,
    standardize.data = TRUE,      # â† switch to TRUE
    splOrd = 4, method = "MM-CD", group.all = TRUE,
    lambda = lam, nlambda = NULL,
    init.rnd = TRUE, beta.threshold.toler = 1e-12, toler = 1e-2,
    maxiter = 15000, cReFr = 100, print.out = FALSE
  ), silent = TRUE)
  
  
  if (inherits(fit, "try-error")) return(data.frame(
    method="OVGLASSO", train.func.mse=NA, test.func.mse=NA,
    train.real.mse=NA, test.real.mse=NA, n.clu=NA, rule.ave=NA, rule.max=NA,
    note=paste0("fit err: ", as.character(fit))
  ))
  
  ## ---- é‡å»º Beta(s,t) å¹¶é¢„æµ‹ï¼ˆå…¼å®¹ 3D/2D coef.pathï¼‰----
  cp <- fit$coef.path
  nd <- length(dim(cp))
  
  if (nd == 3L) {
    # Expected shape: (nlambda, L, M) == (nlambda, p, q)
    maxl <- dim(cp)[1L]
    beta.arr <- array(NA_real_, dim = c(maxl, P_full, r))
    for (k in seq_len(maxl)) {
      C_pk <- cp[k, , , drop = FALSE]        # 1 Ã— p Ã— q
      C_pk <- matrix(C_pk, nrow = p, ncol = q)
      beta.arr[k, , ] <- Psi %*% C_pk %*% t(phi)  # P_full Ã— r
    }
  } else if (nd == 2L) {
    # Fallback: rows correspond to lambda; each row is vec(C_pk) of length p*q
    maxl <- nrow(cp)
    stopifnot(ncol(cp) == p * q)
    beta.arr <- array(NA_real_, dim = c(maxl, P_full, r))
    for (k in seq_len(maxl)) {
      C_pk <- matrix(cp[k, ], nrow = p, ncol = q, byrow = FALSE)
      beta.arr[k, , ] <- Psi %*% C_pk %*% t(phi)
    }
  } else {
    stop("Unexpected coef.path dimensions: ", paste(dim(cp), collapse = "Ã—"))
  }
  
  # ç»Ÿä¸€åœ¨åŽŸç½‘æ ¼ä¸Šé¢„æµ‹
  mse_mat <- matrix(NA_real_, nrow = maxl, ncol = 2)
  yhat_tr_list <- vector("list", maxl)
  yhat_te_list <- vector("list", maxl)
  for (k in seq_len(maxl)) {
    Bk <- beta.arr[k, , ]
    yhat_tr <- x_train %*% Bk
    yhat_te <- x_test  %*% Bk
    yhat_tr_list[[k]] <- yhat_tr
    yhat_te_list[[k]] <- yhat_te
    mse_mat[k, 1] <- mean((as.vector(y_train) - as.vector(yhat_tr))^2)
    mse_mat[k, 2] <- mean((as.vector(y_test)  - as.vector(yhat_te))^2)
  }
  best <- which.min(mse_mat[, 2])
  yhat_tr_log <- yhat_tr_list[[best]]
  yhat_te_log <- yhat_te_list[[best]]
  
  tr_log  <- mean((as.vector(y_train) - as.vector(yhat_tr_log))^2)
  te_log  <- mean((as.vector(y_test)  - as.vector(yhat_te_log))^2)
  tr_real <- mean((as.vector(inv_out(y_train)) - as.vector(inv_out(yhat_tr_log)))^2)
  te_real <- mean((as.vector(inv_out(y_test))  - as.vector(inv_out(yhat_te_log)))^2)
  
  data.frame(method="OVGLASSO(Bernardi)",
             train.func.mse=tr_log, test.func.mse=te_log,
             train.real.mse=tr_real, test.real.mse=te_real,
             n.clu=NA, rule.ave=NA, rule.max=NA, note=NA)
  
}




## ====== 6) æ‰§è¡Œ baselines ======
tab_pls     <- run_pls_interaction()
tab_fapls   <- run_fapls()
tab_slasso  <- run_slasso()

## ====== 7) æ±‡æ€»ä¸Žä¿å­˜ ======
# --- å¯¹é½æ‰€æœ‰ç»“æžœè¡¨çš„åˆ—åä¸Žé¡ºåº ---
ensure_cols <- function(df) {
  need <- c(
    "method",
    "train.func.mse","test.func.mse","train.real.mse","test.real.mse",
    "train.func.mse.c","test.func.mse.c","train.real.mse.c","test.real.mse.c",
    "n.clu","rule.ave","rule.max","note"
  )
  # ç¼ºå•¥è¡¥å•¥ï¼ˆç”¨ NAï¼‰ï¼Œå¤šä½™çš„å…ˆä¿ç•™ï¼Œæœ€åŽå†è£å‰ªåˆ° need
  for (nm in need) if (!nm %in% names(df)) df[[nm]] <- NA
  # åªä¿ç•™å¹¶æŒ‰ç»Ÿä¸€é¡ºåºè¿”å›ž
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

csv_name  <- file.path(output_dir, "agriculture_paper_results.csv")
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

## æ¸…ç†å¹¶è¡Œ
stopCluster(cl)
