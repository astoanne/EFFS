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
    file.path(cwd, "EFFS", "organized_code2", "run_agriculture_paper.R"),
    file.path(cwd, "run_all_paper.R"),
    file.path(cwd, "organized_code2", "run_all_paper.R"),
    file.path(cwd, "EFFS", "organized_code2", "run_all_paper.R")
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

candidate_dirs <- resolve_candidate_dirs()
stock_script <- find_existing_file("run_stock_paper.R", candidate_dirs)
agriculture_script <- find_existing_file("run_agriculture_paper.R", candidate_dirs)
if (is.na(stock_script) || is.na(agriculture_script)) {
  stop("Could not locate organized_code2 runner scripts.")
}

source(stock_script)
source(agriculture_script)
