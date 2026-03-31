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
  known_script_paths <- c(
    file.path(cwd, "plot_stock_rule_curves.R"),
    file.path(cwd, "organized_code2", "plot_stock_rule_curves.R"),
    file.path(cwd, "EFFS", "organized_code2", "plot_stock_rule_curves.R"),
    file.path(cwd, "run_stock_paper.R"),
    file.path(cwd, "organized_code2", "run_stock_paper.R")
  )
  for (path_value in known_script_paths) {
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
if (is.na(data_dir)) {
  stop(
    "organized_code2 is missing the bundled real_data folder.\n",
    "real_data: ", data_dir
  )
}

output_dir <- file.path(script_dir, "outputs", "stock")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Using EEFFS core: ", main_file)
message("Using data directory: ", data_dir)
message("Saving figures to: ", output_dir)

setwd(project_root)
source(main_file)

get_format_matrix <- function(vec, n_rows = 22, n_cols = 457) {
  stopifnot(length(vec) >= n_rows * n_cols)
  mat <- matrix(as.numeric(vec[1:(n_rows * n_cols)]),
                nrow = n_rows, ncol = n_cols)
  mat <- apply(mat, 2, function(x) {
    (x - min(x)) / (max(x) - min(x) + 1e-12)
  })
  t(mat)
}

evaluate_rule_curves <- function(coef_matrix, fd_template, eval_grid) {
  coef_matrix <- as.matrix(coef_matrix)
  curves <- sapply(seq_len(nrow(coef_matrix)), function(i) {
    fd_i <- fd_template
    fd_i$coefs <- as.matrix(coef_matrix[i, ])
    as.numeric(eval.fd(eval_grid, fd_i))
  })
  if (is.null(dim(curves))) {
    curves <- matrix(curves, ncol = 1)
  }
  curves
}

rule_index_text <- function(n_rules) {
  if (n_rules <= 1) {
    "i = 1"
  } else {
    paste0("i = 1, ..., ", n_rules)
  }
}

center_legend_expr <- function(n_rules) {
  as.expression(lapply(seq_len(n_rules), function(i) {
    bquote(c[.(i)](s))
  }))
}

consequent_legend_expr <- function(n_rules) {
  as.expression(lapply(seq_len(n_rules), function(i) {
    bquote(y[c]^{"(" * .(i) * ")"} * group("(", t, ")"))
  }))
}

center_title_expr <- function(n_rules) {
  bquote("If " ~ x(s) ~ " is " ~ c[k](s) ~ " (" * k == 1 * "," ~ ldots * "," ~ .(n_rules) * ")")
}

consequent_title_expr <- function(n_rules) {
  bquote("Then " ~ y(t) ~ " is " ~ y[c]^{"(" * k * ")"} * group("(", t, ")") ~ " (" * k == 1 * "," ~ ldots * "," ~ .(n_rules) * ")")
}

plot_rule_snapshot <- function(center_coefs, cons_coefs, input_fd_template,
                               output_fd_template, input_grid, output_grid,
                               output_file, snapshot_label,
                               step_label = NA_integer_) {
  centers <- evaluate_rule_curves(center_coefs, input_fd_template, input_grid)
  consequents <- evaluate_rule_curves(cons_coefs, output_fd_template, output_grid)
  n_rules <- ncol(centers)
  
  colors <- c("black", "#f1a3b5", "#cfe8c3", "#7fc9ff", "#fdb462", "#8da0cb")
  ltys <- c(1, 2, 3, 4, 5, 6)
  lwds <- c(2.4, 1.8, 1.8, 1.8, 1.8, 1.8)
  colors <- rep(colors, length.out = n_rules)
  ltys <- rep(ltys, length.out = n_rules)
  lwds <- rep(lwds, length.out = n_rules)
  
  top_legend <- center_legend_expr(n_rules)
  bottom_legend <- consequent_legend_expr(n_rules)
  top_title <- center_title_expr(n_rules)
  bottom_title <- consequent_title_expr(n_rules)
  outer_title <- sprintf("EEFFS stock: %s", snapshot_label)
  if (!is.na(step_label)) {
    outer_title <- sprintf("%s (step %s)", outer_title, step_label)
  }
  
  png(filename = output_file, width = 1200, height = 1500, res = 220)
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)
  
  par(mfrow = c(2, 1), mar = c(5, 5, 4.5, 2), oma = c(0, 0, 2.5, 5), xpd = NA)
  
  matplot(input_grid, centers, type = "l",
          col = colors, lty = ltys, lwd = lwds,
          xlab = "s", ylab = "centers",
          main = top_title,
          xaxt = "n", yaxt = "n")
  axis(1)
  axis(2, las = 0)
  box()
  legend("topright", inset = c(-0.22, 0.02), legend = top_legend,
         col = colors, lty = ltys, lwd = lwds, bty = "o", cex = 0.9)
  
  matplot(output_grid, consequents, type = "l",
          col = colors, lty = ltys, lwd = lwds,
          xlab = "t", ylab = "consequents",
          main = bottom_title,
          xaxt = "n", yaxt = "n")
  axis(1)
  axis(2, las = 0)
  box()
  legend("topright", inset = c(-0.22, 0.02), legend = bottom_legend,
         col = colors, lty = ltys, lwd = lwds, bty = "o", cex = 0.9)
  
  mtext(outer_title, outer = TRUE, cex = 1.0, font = 2, line = 0.4)
}

plot_rule_snapshot_compact <- function(center_coefs, cons_coefs, input_fd_template,
                                       output_fd_template, input_grid, output_grid,
                                       output_file) {
  centers <- evaluate_rule_curves(center_coefs, input_fd_template, input_grid)
  consequents <- evaluate_rule_curves(cons_coefs, output_fd_template, output_grid)
  n_rules <- ncol(centers)
  
  colors <- c("black", "#f1a3b5", "#cfe8c3", "#7fc9ff", "#fdb462", "#8da0cb")
  ltys <- c(1, 2, 3, 4, 5, 6)
  lwds <- c(2.0, 1.5, 1.5, 1.5, 1.5, 1.5)
  colors <- rep(colors, length.out = n_rules)
  ltys <- rep(ltys, length.out = n_rules)
  lwds <- rep(lwds, length.out = n_rules)
  
  png(filename = output_file, width = 1800, height = 640, res = 220)
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)
  
  par(mfrow = c(1, 2), mar = c(5.2, 4.0, 2.8, 0.8), oma = c(0.2, 0, 1.4, 0))
  
  matplot(input_grid, centers, type = "l",
          col = colors, lty = ltys, lwd = lwds,
          xlab = "s", ylab = "centers",
          main = center_title_expr(n_rules))
  box()
  legend("topright", legend = center_legend_expr(n_rules),
         col = colors, lty = ltys, lwd = lwds, bty = "o", cex = 0.70)
  
  matplot(output_grid, consequents, type = "l",
          col = colors, lty = ltys, lwd = lwds,
          xlab = "t", ylab = "consequents",
          main = consequent_title_expr(n_rules),
          yaxt = "n")
  axis(2, las = 0)
  box()
  legend("topright", legend = consequent_legend_expr(n_rules),
         col = colors, lty = ltys, lwd = lwds, bty = "o", cex = 0.70)
  
  mtext("EEFFS stock: peak snapshot", outer = TRUE, cex = 0.95, font = 2, line = 0.2)
}

t_data <- read.csv(file.path(data_dir, "stock.csv"), header = FALSE)
t_data <- t_data[2:nrow(t_data), ]

input_name <- "AAPL"
output_name <- "NVDA"

column_indices <- c(
  AAPL = 2,
  GOOGL = 3,
  MSFT = 4,
  NVDA = 5
)

input.ori <- get_format_matrix(t_data[[column_indices[input_name]]])
output.raw <- get_format_matrix(t_data[[column_indices[output_name]]])

min_in <- min(input.ori)
min_out <- min(output.raw)
input <- if (min_in < 0) log(input.ori - min_in + 1) else log(input.ori + 1)
output <- if (min_out < 0) log(output.raw - min_out + 1) else log(output.raw + 1)

n.row <- nrow(output)
ncol.input <- ncol(input)
ncol.output <- ncol(output)

Time.input <- seq(0, 1, length.out = ncol.input)
Time.output <- seq(0, 1, length.out = ncol.output)
Time.t <- seq(0, 1, length.out = 1000)
Time.t.len <- length(Time.t)

orders <- c(3, 4, 5)
order.idx <- 2
nbasis.input <- 10
nbasis.output <- 10
M.clu <- 50
add.thr <- 0.5
merge.thr <- 0.9
train_test_split <- 0.5
NN <- floor(train_test_split * n.row)

message("Running stock EEFFS visualization: ", input_name, " -> ", output_name)

result <- main.STR.FFS.visual(
  add.thr = add.thr,
  merge.thr = merge.thr,
  Time.t = Time.t,
  Time.t.len = Time.t.len,
  nbasis.input = nbasis.input,
  nbasis.output = nbasis.output,
  Time.input.in = Time.input,
  Time.output.in = Time.output,
  NN = NN,
  data.X.Y.d = cbind(input, output),
  orders = orders[order.idx],
  ncol.input = ncol.input,
  ncol.output = ncol.output,
  M.clu = M.clu,
  pruning.strategy = "ema"
)

final_center <- result$center[seq_len(result$n.clu), , drop = FALSE]
final_cons <- result$cons.para[seq_len(result$n.clu), , drop = FALSE]

final_file <- file.path(output_dir, "stock_rule_curves_final.png")
plot_rule_snapshot(
  center_coefs = final_center,
  cons_coefs = final_cons,
  input_fd_template = result$x.spline$X.fd,
  output_fd_template = result$y.spline$X.fd,
  input_grid = Time.input,
  output_grid = Time.output,
  output_file = final_file,
  snapshot_label = sprintf("final rule base (%d rules)", result$n.clu)
)

max_file <- file.path(output_dir, "stock_rule_curves_max_rules.png")
plot_rule_snapshot(
  center_coefs = result$center.max.snapshot,
  cons_coefs = result$cons.para.max.snapshot,
  input_fd_template = result$x.spline$X.fd,
  output_fd_template = result$y.spline$X.fd,
  input_grid = Time.input,
  output_grid = Time.output,
  output_file = max_file,
  snapshot_label = sprintf("peak snapshot (%d rules)", result$n.clu.max.snapshot),
  step_label = result$max.rule.time
)

max_compact_file <- file.path(output_dir, "stock_rule_curves_max_rules_compact.png")
plot_rule_snapshot_compact(
  center_coefs = result$center.max.snapshot,
  cons_coefs = result$cons.para.max.snapshot,
  input_fd_template = result$x.spline$X.fd,
  output_fd_template = result$y.spline$X.fd,
  input_grid = Time.input,
  output_grid = Time.output,
  output_file = max_compact_file
)

metadata_file <- file.path(output_dir, "stock_rule_curve_plots.csv")
write.csv(
  data.frame(
    snapshot = c("final", "max_rules", "max_rules_compact"),
    n_rules = c(result$n.clu, result$n.clu.max.snapshot, result$n.clu.max.snapshot),
    stream_step = c(NA_integer_, result$max.rule.time, result$max.rule.time),
    file = c(final_file, max_file, max_compact_file),
    stringsAsFactors = FALSE
  ),
  metadata_file,
  row.names = FALSE
)

cat("Saved plot:", final_file, "\n")
cat("Saved plot:", max_file, "\n")
cat("Saved plot:", max_compact_file, "\n")
cat("Saved metadata:", metadata_file, "\n")
