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
    file.path(cwd, "validate_paper_results.R"),
    file.path(cwd, "run_stock_paper.R"),
    file.path(cwd, "run_agriculture_paper.R"),
    file.path(cwd, "organized_code2", "validate_paper_results.R"),
    file.path(cwd, "organized_code2", "run_stock_paper.R"),
    file.path(cwd, "organized_code2", "run_agriculture_paper.R"),
    file.path(cwd, "EFFS", "organized_code2", "validate_paper_results.R"),
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
script_dir <- find_existing_dir("reference_outputs", candidate_dirs)
if (is.na(script_dir)) {
  stop("Could not locate organized_code2/reference_outputs from the current session.")
}
script_dir <- dirname(script_dir)
reference_dir <- file.path(script_dir, "reference_outputs")
generated_dir <- file.path(script_dir, "outputs")

read_case <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

paper_table_from_reference <- function(df, case_name) {
  if (case_name == "stock") {
    out <- data.frame(
      Method = c("EEFFS", "EFFNS", "PLS-int", "PLS-main", "fAPLS", "SLASSO"),
      train = c(
        df[df$method == "FFS-ema", "train.real.mse.c"],
        df[df$method == "FFS-none", "train.real.mse.c"],
        df[df$method == "PLS-interaction(full)", "train.real.mse"],
        df[df$method == "PLS-main-only", "train.real.mse"],
        df[df$method == "fAPLS", "train.real.mse"],
        df[df$method == "SLASSO", "train.real.mse"]
      ),
      test = c(
        df[df$method == "FFS-ema", "test.real.mse.c"],
        df[df$method == "FFS-none", "test.real.mse.c"],
        df[df$method == "PLS-interaction(full)", "test.real.mse"],
        df[df$method == "PLS-main-only", "test.real.mse"],
        df[df$method == "fAPLS", "test.real.mse"],
        df[df$method == "SLASSO", "test.real.mse"]
      ),
      stringsAsFactors = FALSE
    )
  } else {
    out <- data.frame(
      Method = c("EEFFS", "EFFNS", "PLS-int", "PLS-main", "fAPLS", "SLASSO"),
      train = c(
        df[df$method == "FFS-ema", "train.real.mse.c"],
        df[df$method == "FFS-none", "train.real.mse.c"],
        df[df$method == "PLS-interaction(full)", "train.real.mse"],
        df[df$method == "PLS-main-only", "train.real.mse"],
        df[df$method == "fAPLS", "train.real.mse"],
        df[df$method == "SLASSO", "train.real.mse"]
      ),
      test = c(
        df[df$method == "FFS-ema", "test.real.mse.c"],
        df[df$method == "FFS-none", "test.real.mse.c"],
        df[df$method == "PLS-interaction(full)", "test.real.mse"],
        df[df$method == "PLS-main-only", "test.real.mse"],
        df[df$method == "fAPLS", "test.real.mse"],
        df[df$method == "SLASSO", "test.real.mse"]
      ),
      stringsAsFactors = FALSE
    )
  }

  out$train <- round(100 * out$train, 4)
  out$test <- round(100 * out$test, 4)
  out
}

rule_table_from_reference <- function(df) {
  out <- data.frame(
    Method = c("EEFFS", "EFFNS"),
    Final = c(df[df$method == "FFS-ema", "n.clu"], df[df$method == "FFS-none", "n.clu"]),
    Avg = c(df[df$method == "FFS-ema", "rule.ave"], df[df$method == "FFS-none", "rule.ave"]),
    Max = c(df[df$method == "FFS-ema", "rule.max"], df[df$method == "FFS-none", "rule.max"]),
    stringsAsFactors = FALSE
  )
  out$Avg <- round(out$Avg, 4)
  out
}

compare_to_reference <- function(case_name, generated_path, reference_path) {
  cat("\n==", toupper(case_name), "==\n")
  ref <- read_case(reference_path)
  cat("Reference paper table values:\n")
  print(paper_table_from_reference(ref, case_name), row.names = FALSE)
  cat("Reference rule statistics:\n")
  print(rule_table_from_reference(ref), row.names = FALSE)

  if (!file.exists(generated_path)) {
    cat("Generated output not found:", generated_path, "\n")
    return(invisible(NULL))
  }

  gen <- read_case(generated_path)
  merged <- merge(ref, gen, by = "method", suffixes = c(".ref", ".gen"))
  numeric_cols <- c(
    "train.func.mse", "test.func.mse", "train.real.mse", "test.real.mse",
    "train.func.mse.c", "test.func.mse.c", "train.real.mse.c", "test.real.mse.c",
    "n.clu", "rule.ave", "rule.max"
  )
  numeric_cols <- numeric_cols[numeric_cols %in% names(ref) & numeric_cols %in% names(gen)]

  diffs <- lapply(numeric_cols, function(col) {
    abs(merged[[paste0(col, ".ref")]] - merged[[paste0(col, ".gen")]])
  })
  names(diffs) <- numeric_cols
  max_diffs <- vapply(diffs, max, numeric(1), na.rm = TRUE)
  max_diffs[!is.finite(max_diffs)] <- NA_real_

  cat("Max absolute difference vs reference by column:\n")
  print(data.frame(column = names(max_diffs), max_abs_diff = max_diffs), row.names = FALSE)
}

compare_to_reference(
  case_name = "stock",
  generated_path = file.path(generated_dir, "stock", "stock_paper_results.csv"),
  reference_path = file.path(reference_dir, "stock_paper_reference.csv")
)

compare_to_reference(
  case_name = "agriculture",
  generated_path = file.path(generated_dir, "agriculture", "agriculture_paper_results.csv"),
  reference_path = file.path(reference_dir, "agriculture_paper_reference.csv")
)
