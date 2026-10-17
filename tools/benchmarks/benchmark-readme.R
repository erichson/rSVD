#!/usr/bin/env Rscript

# Reproducible comparisons for README.md. Requires only rsvd and Matrix.
# Run numerical jobs serially and set BLAS thread limits BEFORE starting R.
# Example:
# OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 Rscript \
#   tools/benchmarks/benchmark-readme.R --library candidate-library \
#   --output benchmark-results --suite all --reps 5 --source-ref COMMIT
#
# References request the same rank and returned factors. base::svd computes
# the complete singular spectrum even when nu/nv are small; it is not a
# specialized truncated-SVD solver. Sparse-to-dense conversion belongs to
# the measured deterministic call. ID/CUR references are the package's
# deterministic LAPACK pivoted-QR paths, not functions supplied by base R.
#
# Fixture creation, RNG setup, explicit GC, and numerical validation are
# outside timing. Method positions rotate across repetitions. Each method
# has a warmup and a separate Rprofmem run. R allocations are not peak RSS,
# live memory, or a complete accounting of native-library allocations.

readme_options <- function(args) {
  usage <- paste(
    "Rscript tools/benchmarks/benchmark-readme.R --library PATH --output DIR",
    "  [--suite all|tables|scaling|smoke] [--reps 5] [--cases name,name]",
    "  [--source-ref COMMIT] [--machine DESCRIPTION]", sep = "\n")
  if (!length(args) || "--help" %in% args) {
    cat(usage, "\n")
    return(NULL)
  }
  if (length(args) %% 2L) stop(usage, call. = FALSE)
  options <- list(suite = "all", reps = "5", cases = NULL,
                  `source-ref` = "unspecified", machine = "unspecified")
  allowed <- c("library", "output", names(options))
  for (i in seq.int(1L, length(args), by = 2L)) {
    key <- sub("^--", "", args[i])
    if (!startsWith(args[i], "--") || !key %in% allowed) {
      stop("Unknown argument: ", args[i], "\n", usage, call. = FALSE)
    }
    options[[key]] <- args[i + 1L]
  }
  if (is.null(options$library) || is.null(options$output) ||
      !options$suite %in% c("all", "tables", "scaling", "smoke")) {
    stop(usage, call. = FALSE)
  }
  reps <- suppressWarnings(as.numeric(options$reps))
  if (length(reps) != 1L || !is.finite(reps) || reps < 3 || reps != floor(reps)) {
    stop("--reps must be an integer of at least three.", call. = FALSE)
  }
  options$reps <- as.integer(reps)
  options$library <- normalizePath(options$library, mustWork = TRUE)
  dir.create(options$output, recursive = TRUE, showWarnings = FALSE)
  options$output <- normalizePath(options$output, mustWork = TRUE)
  if (file.exists(file.path(options$output, "summary.csv"))) {
    stop("Use a new output directory; existing measurements will not be overwritten.",
         call. = FALSE)
  }
  dir.create(file.path(options$output, "profiles"), showWarnings = FALSE)
  options
}

readme_cases <- function() {
  case <- function(id, family, storage, m, n, k = 40L, table = FALSE,
                   scaling = FALSE, description = "") {
    data.frame(case = id, family = family, storage = storage, m = m, n = n,
               k = k, p = 10L, q = if (family %in% c("id", "cur")) 0L else 2L,
               table = table, scaling = scaling, description = description,
               stringsAsFactors = FALSE)
  }
  out <- list()
  for (m in c(10000L, 15000L, 20000L)) {
    out[[length(out) + 1L]] <- case(
      paste0("svd_tall_", m), "svd", "dense", m, 3L * m / 10L,
      table = m %in% c(10000L, 20000L), scaling = TRUE,
      description = "Rank-80 geometric signal plus dense Gaussian noise")
  }
  out[[length(out) + 1L]] <- case(
    "pca_dense", "pca", "dense", 15000L, 4500L, table = TRUE,
    description = "Geometric signal, dense noise, and nonzero column offsets")
  for (m in c(10000L, 15000L, 20000L)) {
    out[[length(out) + 1L]] <- case(
      paste0("pca_sparse_scale_", m), "pca", "sparse", m, 3L * m / 10L,
      table = m %in% c(10000L, 20000L), scaling = TRUE,
      description = "Twenty noisy rank-two blocks; 5 percent density")
  }
  for (family in c("id", "cur")) {
    for (storage in c("dense", "sparse")) {
      out[[length(out) + 1L]] <- case(
        paste(family, storage, sep = "_"), family, storage, 15000L, 4500L,
        table = TRUE,
        description = if (storage == "dense")
          "Rank-80 geometric signal plus dense Gaussian noise" else
          "Twenty noisy rank-two blocks; 5 percent density")
    }
  }
  do.call(rbind, out)
}

readme_dense_fixture <- function(m, n, seed = 20261016L, offsets = FALSE,
                                 block_columns = 128L) {
  set.seed(seed)
  rank <- min(80L, min(m, n) - 1L)
  left <- qr.Q(qr(matrix(stats::rnorm(m * rank), m, rank)))
  right <- qr.Q(qr(matrix(stats::rnorm(n * rank), n, rank)))
  values <- 10^seq(0, -2, length.out = rank)
  A <- left %*% (values * t(right))
  # The noise makes the matrix full rank; its RMS Frobenius norm is 0.01.
  # Column blocks preserve the full-matrix RNG ordering while avoiding a
  # second full-size noise matrix and full-size preprocessing temporaries.
  # Fixture generation is outside all timing and allocation profiles.
  noise_sd <- 0.01 / sqrt(m * n)
  for (start in seq.int(1L, n, by = block_columns)) {
    columns <- seq.int(start, min(n, start + block_columns - 1L))
    block <- A[, columns, drop = FALSE] +
      matrix(stats::rnorm(m * length(columns)), m, length(columns)) * noise_sd
    if (offsets) block <- sweep(block, 2L, 0.2 * cos(columns / 11), "+")
    A[, columns] <- block
  }
  A
}

readme_sparse_fixture <- function(m, n, seed = 20261016L) {
  blocks <- 20L
  stopifnot(m %% blocks == 0L, n %% blocks == 0L)
  bm <- m / blocks
  bn <- n / blocks
  stopifnot(bm >= 2L, bn >= 2L)
  set.seed(seed)
  pieces <- lapply(seq_len(blocks), function(index) {
    slope <- seq(-1, 1, length.out = bm)
    slope <- slope / sqrt(sum(slope^2))
    left <- cbind(rep(1 / sqrt(bm), bm), slope)
    right <- qr.Q(qr(matrix(stats::rnorm(bn * 2L), bn, 2L)))
    values <- c(1 + index / blocks, 0.4 + 0.1 * index / blocks)
    A <- left %*% (values * t(right))
    # Noise stays within the blocks, preserving a precisely stated density.
    A <- A + matrix(stats::rnorm(bm * bn), bm, bn) *
      (0.01 / sqrt(m * n / blocks))
    Matrix::Matrix(A, sparse = TRUE)
  })
  Matrix::bdiag(pieces)
}

readme_fixture <- function(spec) {
  if (spec$storage == "sparse") {
    return(readme_sparse_fixture(spec$m, spec$n))
  }
  readme_dense_fixture(spec$m, spec$n, offsets = spec$family == "pca")
}

readme_nonzeros <- function(A, block_columns = 128L) {
  if (inherits(A, "sparseMatrix")) return(Matrix::nnzero(A))
  count <- 0
  for (start in seq.int(1L, ncol(A), by = block_columns)) {
    columns <- seq.int(start, min(ncol(A), start + block_columns - 1L))
    count <- count + sum(A[, columns, drop = FALSE] != 0)
  }
  count
}

readme_methods <- function(spec) {
  k <- spec$k
  p <- spec$p
  q <- spec$q
  if (spec$family == "svd") {
    return(list(
      base_svd = function(A) base::svd(as.matrix(A), nu = k, nv = k),
      rsvd = function(A) rsvd::rsvd(A, k = k, p = p, q = q)))
  }
  if (spec$family == "pca") {
    methods <- list(
      prcomp = function(A) stats::prcomp(as.matrix(A), rank. = k, retx = TRUE,
                                        center = TRUE, scale. = TRUE),
      rpca = function(A) rsvd::rpca(A, k = k, retx = TRUE, center = TRUE,
                                   scale = TRUE, p = p, q = q))
    if (spec$storage == "sparse") {
      methods$rpca_dense <- function(A) rsvd::rpca(as.matrix(A), k = k,
                                                retx = TRUE, center = TRUE,
                                                scale = TRUE, p = p, q = q)
    }
    return(methods)
  }
  if (spec$family == "id") {
    return(list(
      id_deterministic = function(A) rsvd::rid(A, k = k, rand = FALSE),
      rid = function(A) rsvd::rid(A, k = k, p = p, q = q, rand = TRUE)))
  }
  if (spec$family == "cur") {
    return(list(
      cur_deterministic = function(A) rsvd::rcur(A, k = k, rand = FALSE),
      rcur = function(A) rsvd::rcur(A, k = k, p = p, q = q, rand = TRUE)))
  }
  stop("Unknown family")
}

readme_orthogonality <- function(Q) {
  sqrt(sum((crossprod(Q) - diag(ncol(Q)))^2))
}

readme_validate <- function(spec, A, fit, oracle = NULL, block_rows = 128L) {
  k <- spec$k
  center <- scale <- NULL
  orthogonality <- NA_real_
  spectrum <- NULL
  if (spec$family == "svd") {
    stopifnot(identical(dim(fit$u), c(nrow(A), k)),
              identical(dim(fit$v), c(ncol(A), k)),
              length(fit$d) >= k)
    spectrum <- fit$d[seq_len(k)]
    weighted <- sweep(fit$v, 2L, spectrum, "*")
    reconstruct <- function(rows) tcrossprod(fit$u[rows, , drop = FALSE], weighted)
    orthogonality <- max(readme_orthogonality(fit$u),
                          readme_orthogonality(fit$v))
  } else if (spec$family == "pca") {
    stopifnot(identical(dim(fit$x), c(nrow(A), k)),
              identical(dim(fit$rotation), c(ncol(A), k)))
    center <- fit$center
    scale <- fit$scale
    spectrum <- fit$sdev[seq_len(k)] * sqrt(nrow(A) - 1)
    reconstruct <- function(rows) tcrossprod(fit$x[rows, , drop = FALSE],
                                             fit$rotation)
    orthogonality <- readme_orthogonality(fit$rotation)
  } else if (spec$family == "id") {
    stopifnot(identical(dim(fit$C), c(nrow(A), k)),
              identical(dim(fit$Z), c(k, ncol(A))))
    stopifnot(max(abs(fit$Z[, fit$idx, drop = FALSE] - diag(k))) < 1e-10)
    reconstruct <- function(rows) fit$C[rows, , drop = FALSE] %*% fit$Z
  } else {
    stopifnot(identical(dim(fit$C), c(nrow(A), k)),
              identical(dim(fit$R), c(k, ncol(A))),
              identical(dim(fit$U), c(k, k)))
    reconstruct <- function(rows) (fit$C[rows, , drop = FALSE] %*% fit$U) %*% fit$R
  }
  error_sq <- input_sq <- 0
  for (start in seq.int(1L, nrow(A), by = block_rows)) {
    rows <- seq.int(start, min(nrow(A), start + block_rows - 1L))
    block <- as.matrix(A[rows, , drop = FALSE])
    if (!is.null(center) && !identical(center, FALSE)) {
      block <- sweep(block, 2L, center, "-")
    }
    if (!is.null(scale) && !identical(scale, FALSE)) {
      block <- sweep(block, 2L, scale, "/")
    }
    input_sq <- input_sq + sum(block^2)
    error_sq <- error_sq + sum((block - reconstruct(rows))^2)
  }
  relative_error <- sqrt(error_sq / input_sq)
  if (!is.finite(relative_error) || relative_error > 1 ||
      (!is.na(orthogonality) && orthogonality > 1e-8)) {
    stop("Numerical validation failed for ", spec$case, call. = FALSE)
  }
  ratio <- spectrum_error <- NA_real_
  if (!is.null(oracle)) {
    tail <- sqrt(sum(oracle[-seq_len(k)]^2))
    ratio <- sqrt(error_sq) / tail
    spectrum_error <- sqrt(sum((spectrum - oracle[seq_len(k)])^2) /
                             sum(oracle[seq_len(k)]^2))
  }
  list(relative_frobenius_error = relative_error,
       error_over_optimal_rank_k = ratio,
       retained_spectrum_relative_error = spectrum_error,
       orthogonality_error = orthogonality)
}

readme_profile <- function(run, A, seed, filename) {
  gc()
  set.seed(seed)
  utils::Rprofmem(filename)
  on.exit(utils::Rprofmem(NULL), add = TRUE)
  fit <- run(A)
  utils::Rprofmem(NULL)
  fit <- NULL
  lines <- readLines(filename, warn = FALSE)
  records <- lines[grepl("^[0-9]+[[:space:]]*:", lines)]
  sizes <- as.numeric(sub("^([0-9]+).*", "\\1", records))
  list(r_allocation_count = length(sizes),
       r_allocation_total_bytes = sum(sizes),
       r_allocation_largest_bytes = if (length(sizes)) max(sizes) else 0)
}

readme_write <- function(rows, filename) {
  utils::write.csv(do.call(rbind, rows), filename, row.names = FALSE, na = "")
}

readme_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  options <- readme_options(args)
  if (is.null(options)) return(invisible(NULL))
  if (!capabilities("profmem")) stop("This benchmark needs Rprofmem support.")
  .libPaths(c(options$library, .libPaths()))
  loadNamespace("rsvd", lib.loc = options$library)
  loadNamespace("Matrix")
  metadata <- c(
    timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    source_ref = options[["source-ref"]], machine = options$machine,
    R_version = as.character(getRversion()), R_platform = R.version$platform,
    rsvd_version = as.character(utils::packageVersion("rsvd")),
    Matrix_version = as.character(utils::packageVersion("Matrix")),
    BLAS = unname(extSoftVersion()["BLAS"]),
    OPENBLAS_NUM_THREADS = Sys.getenv("OPENBLAS_NUM_THREADS"),
    OMP_NUM_THREADS = Sys.getenv("OMP_NUM_THREADS"),
    MKL_NUM_THREADS = Sys.getenv("MKL_NUM_THREADS"),
    repetitions = options$reps, warmups_per_method = 1L,
    fixture_seed = 20261016L, run_seed_start = 8100L,
    memory_measure = "Rprofmem allocations; excludes setup; not peak RSS",
    timing_measure = "Elapsed time; excludes setup, explicit GC, RNG setup and validation")
  utils::write.csv(data.frame(key = names(metadata), value = unname(metadata)),
                   file.path(options$output, "metadata.csv"), row.names = FALSE)
  writeLines(trimws(capture.output(sessionInfo()), which = "right"),
             file.path(options$output, "session-info.txt"))
  cases <- readme_cases()
  if (!is.null(options$cases)) {
    selected <- strsplit(options$cases, ",", fixed = TRUE)[[1L]]
    if (!length(selected) || any(!selected %in% cases$case)) stop("Unknown case.")
    cases <- cases[match(selected, cases$case), , drop = FALSE]
  } else if (options$suite == "tables") {
    cases <- cases[cases$table, , drop = FALSE]
  } else if (options$suite == "scaling") {
    cases <- cases[cases$scaling, , drop = FALSE]
  } else if (options$suite == "smoke") {
    cases <- cases[match(c("svd_tall_10000", "pca_dense", "pca_sparse_scale_10000",
                           "id_dense", "cur_sparse"), cases$case), , drop = FALSE]
    cases$m <- 240L
    cases$n <- 120L
    cases$k <- 20L
  }
  summaries <- timings <- fixtures <- list()
  for (case_index in seq_len(nrow(cases))) {
    spec <- cases[case_index, , drop = FALSE]
    spec$k <- as.integer(spec$k)
    cat("Case ", case_index, "/", nrow(cases), ": ", spec$case, "\n", sep = "")
    A <- readme_fixture(spec)
    spec$m <- nrow(A)
    spec$n <- ncol(A)
    methods <- readme_methods(spec)
    oracle <- NULL
    for (method in names(methods)) {
      gc()
      set.seed(8100L)
      fit <- methods[[method]](A)
      if (method == "base_svd") oracle <- fit$d
      if (method == "prcomp") oracle <- fit$sdev * sqrt(nrow(A) - 1L)
      fit <- NULL
    }
    if (!is.null(oracle)) {
      tail_energy <- sum(oracle[-seq_len(spec$k)]^2)
      stopifnot(length(oracle) == min(dim(A)), all(is.finite(oracle)),
                is.finite(tail_energy), tail_energy > 0)
    }
    input_nonzeros <- readme_nonzeros(A)
    fixture <- c(as.list(spec), list(
      input_class = class(A)[1L], input_nonzeros = input_nonzeros,
      density = input_nonzeros / (nrow(A) * ncol(A)),
      input_object_bytes = as.numeric(object.size(A))))
    fixtures[[length(fixtures) + 1L]] <- as.data.frame(fixture)
    case_times <- list()
    for (rep in seq_len(options$reps)) {
      order <- names(methods)
      order <- order[(seq_along(order) + rep - 2L) %% length(order) + 1L]
      for (method in order) {
        fit <- NULL
        gc()
        seed <- 8100L + rep
        set.seed(seed)
        elapsed <- system.time(fit <- methods[[method]](A), gcFirst = FALSE)[["elapsed"]]
        numerical <- readme_validate(spec, A, fit, oracle)
        fit <- NULL
        row <- as.data.frame(c(list(case = spec$case, method = method,
                                    repetition = rep, seed = seed,
                                    elapsed_s = elapsed), numerical))
        case_times[[length(case_times) + 1L]] <- row
        timings[[length(timings) + 1L]] <- row
        cat(sprintf("  %s rep %d: %.3f s, relative error %.3g\n",
                    method, rep, elapsed, numerical$relative_frobenius_error))
      }
    }
    case_times <- do.call(rbind, case_times)
    for (method in names(methods)) {
      rows <- case_times[case_times$method == method, , drop = FALSE]
      profile <- readme_profile(methods[[method]], A, seed = 8101L,
        filename = file.path(options$output, "profiles",
                             paste0(spec$case, "-", method, ".Rprofmem")))
      all_na <- function(x, fun) if (all(is.na(x))) NA_real_ else fun(x, na.rm = TRUE)
      row <- c(fixture, list(method = method, repetitions = options$reps,
        elapsed_median_s = stats::median(rows$elapsed_s),
        elapsed_q25_s = unname(stats::quantile(rows$elapsed_s, 0.25)),
        elapsed_q75_s = unname(stats::quantile(rows$elapsed_s, 0.75)),
        elapsed_min_s = min(rows$elapsed_s), elapsed_max_s = max(rows$elapsed_s),
        relative_error_median = stats::median(rows$relative_frobenius_error),
        relative_error_max = max(rows$relative_frobenius_error),
        error_over_optimal_median = all_na(rows$error_over_optimal_rank_k, stats::median),
        error_over_optimal_max = all_na(rows$error_over_optimal_rank_k, max),
        spectrum_relative_error_max = all_na(rows$retained_spectrum_relative_error, max),
        orthogonality_error_max = all_na(rows$orthogonality_error, max)), profile)
      summaries[[length(summaries) + 1L]] <- as.data.frame(row)
    }
    readme_write(summaries, file.path(options$output, "summary.csv"))
    readme_write(timings, file.path(options$output, "timings.csv"))
    readme_write(fixtures, file.path(options$output, "fixtures.csv"))
    rm(A, methods, oracle)
    gc()
  }
  cat("Saved measurements to ", options$output, "\n", sep = "")
  invisible(do.call(rbind, summaries))
}

if (sys.nframe() == 0L) readme_main()
