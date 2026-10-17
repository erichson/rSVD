#!/usr/bin/env Rscript

# Reproduce a comparison using two separately installed rsvd versions:
#   Rscript tools/benchmark-lowrank.R --library /path/to/old-library \
#     --output old.csv --reps 3
#   Rscript tools/benchmark-lowrank.R --library /path/to/new-library \
#     --output new.csv --reps 3
# Run versions sequentially with the same R, Matrix, BLAS and thread settings.
# For threaded BLAS, set its thread limit before starting R; for example,
# OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript ...
#
# Fixed fixtures have an exact prescribed low rank, except the small full-rank
# ID case. Fixture creation and residual checks are outside all measurements.
# Each case has one warmup, at least three timed repetitions, and one separate
# allocation profile. GC and RNG setup are outside elapsed-time measurements.
#
# Rprofmem records R-managed allocations. Its largest and cumulative allocation
# sizes are NOT peak resident memory, retained memory, or all native allocations.
# Raw profiles are saved beside the output for inspection. These are developer
# benchmarks, not timing assertions in the package's unit tests.

benchmark_options <- function(arguments) {
  usage <- paste(
    "Usage: Rscript tools/benchmark-lowrank.R --library PATH --output FILE.csv",
    "       [--reps 3] [--cases comma,separated,case,names] [--profile-dir PATH]",
    "Cases: qb_sparse_tall,qb_sparse_wide,rsvd_dense_tall,rpca_no_scores,id_full_rank,",
    "       rpca_sparse_tall,rpca_sparse_wide,rpca_sparse_no_scores,rpca_sparse_predict",
    sep = "\n")
  if (length(arguments) == 0L || "--help" %in% arguments) {
    cat(usage, "\n")
    return(NULL)
  }
  if (length(arguments) %% 2L != 0L) stop(usage, call. = FALSE)
  options <- list(reps = "3", cases = NULL, `profile-dir` = NULL)
  allowed <- c("library", "output", "reps", "cases", "profile-dir")
  for (index in seq.int(1L, length(arguments), by = 2L)) {
    key <- sub("^--", "", arguments[index])
    if (!startsWith(arguments[index], "--") || !key %in% allowed) {
      stop("Unknown option: ", arguments[index], "\n", usage, call. = FALSE)
    }
    options[[key]] <- arguments[index + 1L]
  }
  if (is.null(options$library) || is.null(options$output)) {
    stop(usage, call. = FALSE)
  }
  repetitions <- suppressWarnings(as.numeric(options$reps))
  if (length(repetitions) != 1L || !is.finite(repetitions) ||
      repetitions < 3 || repetitions != floor(repetitions)) {
    stop("--reps must be an integer of at least 3.", call. = FALSE)
  }
  options$reps <- as.integer(repetitions)
  options$library <- normalizePath(options$library, mustWork = TRUE)
  dir.create(dirname(options$output), recursive = TRUE, showWarnings = FALSE)
  options$output <- file.path(normalizePath(dirname(options$output), mustWork = TRUE),
                              basename(options$output))
  if (is.null(options[["profile-dir"]])) {
    options[["profile-dir"]] <- paste0(tools::file_path_sans_ext(options$output),
                                       "-profiles")
  }
  dir.create(options[["profile-dir"]], recursive = TRUE, showWarnings = FALSE)
  options[["profile-dir"]] <- normalizePath(options[["profile-dir"]],
                                            mustWork = TRUE)
  options
}

benchmark_dense_fixture <- function(m, n, d) {
  basis <- function(size) {
    out <- sqrt(2 / size) * cos(pi * outer(seq_len(size) - 0.5,
                                          seq_along(d) - 1) / size)
    out[, 1L] <- out[, 1L] / sqrt(2)
    out
  }
  basis(m) %*% (d * t(basis(n)))
}

benchmark_sparse_fixture <- function() {
  blocks <- lapply(seq_len(20L), function(block) {
    left <- sin(seq_len(200L) * (block + 1) / 19) +
      cos(seq_len(200L) * (block + 2) / 23)
    right <- cos(seq_len(100L) * (block + 3) / 17) +
      sin(seq_len(100L) * (block + 1) / 29)
    left <- left / sqrt(sum(left^2))
    right <- right / sqrt(sum(right^2))
    Matrix::Matrix((21 - block) * tcrossprod(left, right), sparse = TRUE)
  })
  Matrix::bdiag(blocks)
}

benchmark_case <- function(name) {
  if (name %in% c("qb_sparse_tall", "qb_sparse_wide")) {
    A <- benchmark_sparse_fixture()
    if (name == "qb_sparse_wide") A <- Matrix::t(A)
    return(list(A = A, algorithm = "rqb", k = 20L, p = 10L, q = 2L,
                rand = TRUE, retx = NA,
                run = function(A) rsvd::rqb(A, k = 20, p = 10, q = 2)))
  }
  if (name == "rsvd_dense_tall") {
    return(list(A = benchmark_dense_fixture(2200L, 1000L, 20:1),
                algorithm = "rsvd", k = 20L, p = 10L, q = 2L,
                rand = TRUE, retx = NA,
                run = function(A) rsvd::rsvd(A, k = 20, p = 10, q = 2)))
  }
  if (name == "rpca_no_scores") {
    return(list(A = benchmark_dense_fixture(6000L, 600L, 40:1),
                algorithm = "rpca", k = 40L, p = 10L, q = 2L,
                rand = TRUE, retx = FALSE,
                run = function(A) rsvd::rpca(A, k = 40, p = 10, q = 2,
                                              center = FALSE, scale = FALSE,
                                              retx = FALSE)))
  }
  if (name %in% c("rpca_sparse_tall", "rpca_sparse_wide",
                   "rpca_sparse_no_scores", "rpca_sparse_predict")) {
    A <- benchmark_sparse_fixture()
    if (name == "rpca_sparse_wide") A <- Matrix::t(A)
    scores <- name != "rpca_sparse_no_scores"
    if (name == "rpca_sparse_predict") {
      # Fit outside the measured prediction call.
      set.seed(731)
      model <- rsvd::rpca(A, k = 20, p = 10, q = 2, retx = FALSE)
      return(list(A = A, algorithm = "predict.rpca", k = 20L, p = 10L,
                  q = 2L, rand = TRUE, retx = TRUE, center = TRUE, scale = TRUE,
                  model = model, run = function(A) stats::predict(model, A)))
    }
    return(list(A = A, algorithm = "rpca", k = 20L, p = 10L, q = 2L,
                rand = TRUE, retx = scores, center = TRUE, scale = TRUE,
                run = function(A) rsvd::rpca(A, k = 20, p = 10, q = 2,
                                              retx = scores)))
  }
  if (name == "id_full_rank") {
    return(list(A = benchmark_dense_fixture(600L, 250L,
                                            seq(2, 1, length.out = 250)),
                algorithm = "rid", k = 250L, p = NA_integer_, q = NA_integer_,
                rand = FALSE, retx = NA,
                run = function(A) rsvd::rid(A, k = ncol(A), rand = FALSE)))
  }
  stop("Unknown benchmark case: ", name, call. = FALSE)
}

benchmark_residual <- function(A, reconstructed_rows, center = FALSE, scale = FALSE) {
  squared_error <- squared_norm <- 0
  for (first in seq.int(1L, nrow(A), by = 256L)) {
    rows <- seq.int(first, min(first + 255L, nrow(A)))
    reference <- as.matrix(A[rows, , drop = FALSE])
    if (!is.logical(center)) reference <- sweep(reference, 2, center, "-")
    if (!is.logical(scale)) reference <- sweep(reference, 2, scale, "/")
    difference <- reference - reconstructed_rows(rows)
    squared_error <- squared_error + sum(Mod(difference)^2)
    squared_norm <- squared_norm + sum(Mod(reference)^2)
  }
  sqrt(squared_error / squared_norm)
}

benchmark_orthogonality <- function(Q) {
  sqrt(sum((crossprod(Q) - diag(ncol(Q)))^2))
}

benchmark_validate <- function(case, fit) {
  A <- case$A
  interpolation_error <- NA_real_
  center <- scale <- FALSE
  if (case$algorithm == "rqb") {
    reconstructed_rows <- function(rows) fit$Q[rows, , drop = FALSE] %*% fit$B
    orthogonality <- benchmark_orthogonality(fit$Q)
    width <- ncol(fit$Q)
  } else if (case$algorithm == "rsvd") {
    weighted_vectors <- fit$d * t(fit$v)
    reconstructed_rows <- function(rows) {
      fit$u[rows, , drop = FALSE] %*% weighted_vectors
    }
    orthogonality <- max(benchmark_orthogonality(fit$u),
                          benchmark_orthogonality(fit$v))
    width <- length(fit$d)
  } else if (case$algorithm %in% c("rpca", "predict.rpca")) {
    model <- if (case$algorithm == "rpca") fit else case$model
    if (case$algorithm == "rpca") {
      if (case$retx) stopifnot(is.matrix(fit$x)) else stopifnot(is.null(fit$x))
      scores <- fit$x
    } else {
      stopifnot(is.matrix(fit), identical(dim(fit), c(nrow(A), case$k)))
      scores <- fit
    }
    center <- model$center
    scale <- model$scale
    transposed_rotation <- t(model$rotation)
    reconstructed_rows <- function(rows) {
      projected <- if (is.null(scores)) {
        stats::predict(model, A[rows, , drop = FALSE])
      } else scores[rows, , drop = FALSE]
      projected %*% transposed_rotation
    }
    orthogonality <- benchmark_orthogonality(model$rotation)
    width <- ncol(model$rotation)
  } else {
    reconstructed_rows <- function(rows) fit$C[rows, , drop = FALSE] %*% fit$Z
    interpolation_error <- max(abs(fit$Z[, fit$idx, drop = FALSE] - diag(case$k)))
    orthogonality <- NA_real_
    width <- nrow(fit$Z)
  }
  residual <- benchmark_residual(A, reconstructed_rows, center, scale)
  if (!is.finite(residual) || residual > 1e-8 ||
      (!is.na(orthogonality) && orthogonality > 1e-8) ||
      (!is.na(interpolation_error) && interpolation_error > 1e-8)) {
    stop("Numerical validation failed for ", case$algorithm, call. = FALSE)
  }
  list(relative_frobenius_residual = residual,
       orthogonality_frobenius_error = orthogonality,
       interpolation_max_error = interpolation_error,
       output_width = width)
}

benchmark_allocations <- function(case, seed, path) {
  gc()
  set.seed(seed)
  utils::Rprofmem(path)
  on.exit(utils::Rprofmem(NULL), add = TRUE)
  fit <- case$run(case$A)
  utils::Rprofmem(NULL)
  rm(fit)
  lines <- readLines(path, warn = FALSE)
  records <- lines[grepl("^[0-9]+[[:space:]]*:", lines)]
  sizes <- as.numeric(sub("^([0-9]+).*", "\\1", records))
  list(r_allocation_count = length(sizes),
       r_allocation_total_bytes = sum(sizes),
       r_allocation_largest_bytes = if (length(sizes)) max(sizes) else 0)
}

benchmark_main <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  options <- benchmark_options(arguments)
  if (is.null(options)) return(invisible(NULL))
  if (!capabilities("profmem")) {
    stop("This R build does not support Rprofmem allocation profiling.", call. = FALSE)
  }
  .libPaths(c(options$library, .libPaths()))
  loadNamespace("rsvd", lib.loc = options$library)
  loadNamespace("Matrix")
  cases <- c("qb_sparse_tall", "qb_sparse_wide", "rsvd_dense_tall",
             "rpca_no_scores", "id_full_rank", "rpca_sparse_tall",
             "rpca_sparse_wide", "rpca_sparse_no_scores", "rpca_sparse_predict")
  if (!is.null(options$cases)) {
    selected <- strsplit(options$cases, ",", fixed = TRUE)[[1L]]
    if (!length(selected) || any(!selected %in% cases)) {
      stop("Unknown or empty --cases selection.", call. = FALSE)
    }
    cases <- selected
  }
  metadata <- list(
    rsvd_version = as.character(utils::packageVersion("rsvd")),
    Matrix_version = as.character(utils::packageVersion("Matrix")),
    R_version = as.character(getRversion()), R_platform = R.version$platform,
    BLAS = unname(extSoftVersion()["BLAS"]),
    OPENBLAS_NUM_THREADS = Sys.getenv("OPENBLAS_NUM_THREADS", unset = ""),
    OMP_NUM_THREADS = Sys.getenv("OMP_NUM_THREADS", unset = ""),
    MKL_NUM_THREADS = Sys.getenv("MKL_NUM_THREADS", unset = ""),
    repetitions = options$reps, warmup_runs = 1L, allocation_profile_runs = 1L,
    memory_measure = "Rprofmem allocations; not peak RSS or all native allocations"
  )
  results <- list()
  for (name in cases) {
    cat("Benchmarking ", name, " with rsvd ", metadata$rsvd_version, "\n", sep = "")
    case <- benchmark_case(name)
    seed <- 731L
    set.seed(seed)
    warmup <- case$run(case$A)
    rm(warmup)
    elapsed <- numeric(options$reps)
    fit <- NULL
    for (repetition in seq_len(options$reps)) {
      fit <- NULL
      gc()
      set.seed(seed)
      elapsed[repetition] <- system.time(fit <- case$run(case$A),
                                          gcFirst = FALSE)[["elapsed"]]
    }
    numerical <- benchmark_validate(case, fit)
    fit <- NULL
    profile_path <- file.path(options[["profile-dir"]], paste0(name, ".Rprofmem"))
    allocations <- benchmark_allocations(case, seed, profile_path)
    timing <- list(elapsed_median_s = stats::median(elapsed),
                   elapsed_min_s = min(elapsed), elapsed_max_s = max(elapsed))
    samples <- as.list(elapsed)
    names(samples) <- paste0("elapsed_rep", seq_along(elapsed), "_s")
    row <- c(list(case = name, algorithm = case$algorithm,
                  input_class = class(case$A)[1L],
                  m = nrow(case$A), n = ncol(case$A), k = case$k,
                  p = case$p, q = case$q, rand = case$rand, retx = case$retx,
                  center = if (!is.null(case$center)) case$center else
                    if (case$algorithm == "rpca") FALSE else NA,
                  scale = if (!is.null(case$scale)) case$scale else
                    if (case$algorithm == "rpca") FALSE else NA,
                  seed = seed, input_nonzeros = Matrix::nnzero(case$A),
                  input_object_bytes = as.numeric(object.size(case$A))),
             metadata, timing, samples, allocations, numerical,
             list(allocation_profile = profile_path))
    results[[name]] <- as.data.frame(row, stringsAsFactors = FALSE)
    utils::write.csv(do.call(rbind, results), options$output, row.names = FALSE)
    cat(sprintf("  median %.3f s; largest R allocation %.2f MiB; total %.2f MiB; residual %.3g\n",
                timing$elapsed_median_s,
                allocations$r_allocation_largest_bytes / 1024^2,
                allocations$r_allocation_total_bytes / 1024^2,
                numerical$relative_frobenius_residual))
    rm(case)
    gc()
  }
  cat("Saved ", options$output, "\n", sep = "")
  invisible(do.call(rbind, results))
}

if (sys.nframe() == 0L) benchmark_main()
