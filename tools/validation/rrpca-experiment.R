#!/usr/bin/env Rscript

# Run from the package root, with other benchmarks idle:
# Rscript tools/validation/rrpca-experiment.R LIBRARY OUTPUT.csv [REFERENCE.R]
#
# LIBRARY contains an installed rsvd package supplying the unchanged rsvd core.
# The candidate robust-PCA function is sourced from R/rrpca.R. The reference
# defaults to R/rrpca.R at commit ad8f05a (requires git history); an explicitly
# supplied REFERENCE.R makes this script usable without a git checkout.
# One warmup, three timed runs, and a separate Rprofmem run use the same seed.
# Fixture creation, output validation, and SVD-call instrumentation are untimed.
# Before/after stopping criteria differ: these are default-call cost comparisons,
# not a benchmark claiming equally certified solutions or an RPCA speedup.
# Rprofmem reports R allocations, not peak RSS or all native allocations.

arguments <- commandArgs(trailingOnly = TRUE)
if (!length(arguments) %in% c(2L, 3L)) {
  stop("Usage: Rscript tools/validation/rrpca-experiment.R LIBRARY OUTPUT.csv [REFERENCE.R]")
}
.libPaths(c(normalizePath(arguments[1L], mustWork = TRUE), .libPaths()))
invisible(loadNamespace("rsvd"))
invisible(loadNamespace("Matrix"))
candidate <- new.env(parent = asNamespace("rsvd"))
sys.source("R/rrpca.R", envir = candidate)
reference <- new.env(parent = asNamespace("rsvd"))
reference_id <- "ad8f05a08491f4a75f12e8d7f14a524ce025d80d"
if (length(arguments) == 3L) {
  sys.source(arguments[3L], envir = reference)
  reference_id <- paste0("file-md5:", unname(tools::md5sum(arguments[3L])))
} else {
  reference_text <- system2("git", c("show", paste0(reference_id, ":R/rrpca.R")),
                            stdout = TRUE)
  if (!is.null(attr(reference_text, "status"))) {
    stop("Cannot read the reference commit; supply its R/rrpca.R as REFERENCE.R.")
  }
  eval(parse(text = reference_text), envir = reference)
}

fixture <- function(m, n, d, number) {
  basis <- function(size) {
    x <- sqrt(2 / size) * cos(pi * outer(seq_len(size) - 0.5,
                                         seq_along(d) - 1) / size)
    x[, 1L] <- x[, 1L] / sqrt(2)
    x
  }
  L <- basis(m) %*% (d * t(basis(n)))
  S <- matrix(0, m, n)
  set.seed(372)
  entries <- sample(length(S), number)
  S[entries] <- 20 * sample(c(-1, 1), number, replace = TRUE)
  list(A = L + S, L = L, S = S, rank = length(d), nonzeros = number)
}

inspect_run <- function(run, A, ...) {
  monitor <- new.env(parent = environment(run))
  monitor$powers <- numeric()
  monitor$full_calls <- 0L
  monitor$rsvd <- function(A, k, ..., q = 2) {
    monitor$powers <- c(monitor$powers, q)
    rsvd::rsvd(A, k = k, q = q, ...)
  }
  monitor$svd <- function(...) {
    monitor$full_calls <- monitor$full_calls + 1L
    base::svd(...)
  }
  environment(run) <- monitor
  set.seed(918)
  progress <- capture.output(fit <- run(A, trace = TRUE, ...))
  stationarity <- progress[grepl("stationarity error =", progress, fixed = TRUE)]
  stationarity <- if (length(stationarity)) {
    as.numeric(sub(".*stationarity error = ", "", tail(stationarity, 1L)))
  } else NA_real_
  list(fit = fit, stationarity = stationarity,
       # The first rsvd call estimates the initial spectral norm.
       random_calls = length(monitor$powers) - 1L,
       full_calls = monitor$full_calls, max_q = max(monitor$powers))
}

output <- arguments[2L]
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
profile_directory <- paste0(tools::file_path_sans_ext(output), "-profiles")
dir.create(profile_directory, recursive = TRUE, showWarnings = FALSE)
fixtures <- list(small = fixture(80L, 60L, c(30, 15, 8), 50L),
                 medium = fixture(600L, 400L, c(100, 75, 50, 30, 15), 2400L),
                 tall = fixture(1200L, 120L, c(500, 100, 50), 1440L))
rows <- list()
for (case in names(fixtures)) {
  example <- fixtures[[case]]
  for (method in c("before", "after")) {
    run <- if (method == "before") reference$rrpca.default else candidate$rrpca.default
    set.seed(918)
    invisible(run(example$A))
    elapsed <- numeric(3L)
    for (repetition in seq_len(3L)) {
      gc()
      set.seed(918)
      elapsed[repetition] <- system.time(fit <- run(example$A),
                                         gcFirst = FALSE)[["elapsed"]]
    }
    profile <- file.path(profile_directory, paste0(case, "-", method, ".Rprofmem"))
    gc()
    set.seed(918)
    utils::Rprofmem(profile)
    fit <- run(example$A)
    utils::Rprofmem(NULL)
    lines <- readLines(profile, warn = FALSE)
    allocations <- lines[grepl("^[0-9]+[[:space:]]*:", lines)]
    sizes <- as.numeric(sub("^([0-9]+).*", "\\1", allocations))
    inspected <- inspect_run(run, example$A)
    row <- data.frame(
      case = case, method = method, reference_revision = reference_id,
      candidate_file_md5 = unname(tools::md5sum("R/rrpca.R")),
      R_version = as.character(getRversion()),
      Matrix_version = as.character(utils::packageVersion("Matrix")),
      BLAS = basename(unname(extSoftVersion()["BLAS"])),
      m = nrow(example$A), n = ncol(example$A), expected_rank = example$rank,
      corrupted_entries = example$nonzeros, maxiter = 50L, tol = 1e-5,
      initial_p = 10L, initial_q = 2L, repetitions = 3L,
      median_seconds = median(elapsed), min_seconds = min(elapsed),
      max_seconds = max(elapsed), iterations = length(fit$err),
      final_relative_residual = tail(fit$err, 1L),
      final_stationarity_diagnostic = inspected$stationarity,
      L_relative_error = norm(fit$L - example$L, "F") / norm(example$L, "F"),
      S_relative_error = norm(fit$S - example$S, "F") / norm(example$S, "F"),
      recovered_rank = sum(svd(fit$L, nu = 0, nv = 0)$d > 1e-4),
      randomized_svd_calls = inspected$random_calls,
      exact_svd_calls = inspected$full_calls,
      maximum_power_iterations = inspected$max_q,
      r_allocation_total_bytes = sum(sizes),
      r_allocation_largest_bytes = max(sizes)
    )
    stopifnot(row$final_relative_residual <= row$tol,
              row$L_relative_error < 1e-4, row$S_relative_error < 1e-4,
              row$recovered_rank == row$expected_rank)
    if (method == "after") stopifnot(row$final_stationarity_diagnostic <= row$tol)
    rows[[length(rows) + 1L]] <- row
    utils::write.csv(do.call(rbind, rows), output, row.names = FALSE)
    cat(sprintf("%s / %s: %.3f s, %d iterations, %d exact SVD calls, max q=%d\n",
                case, method, row$median_seconds, row$iterations,
                row$exact_svd_calls, row$maximum_power_iterations))
  }
}
cat("Saved ", output, "\n", sep = "")

# Untimed examples with analytically known optima make the corrected stopping
# defect independently reproducible. The high-rank shrinkage example checks a
# single proximal update, so it is intentionally limited to one iteration.
rank_twenty <- matrix(0, 80, 60)
diag(rank_twenty)[seq_len(20L)] <- 10
regressions <- list(
  constant_small = list(A = matrix(1, 2, 2), lambda = 2^-0.5),
  constant_thin = list(A = matrix(1, 500, 4), lambda = 500^-0.5),
  diagonal_lowrank = list(A = diag(10), lambda = 2),
  diagonal_sparse = list(A = diag(10), lambda = 0.5, sparse = TRUE),
  rank_twenty = list(A = rank_twenty, lambda = 3, maxiter = 1L)
)
regression_rows <- list()
for (case in names(regressions)) {
  example <- regressions[[case]]
  expected_L <- if (isTRUE(example$sparse)) example$A * 0 else example$A
  expected_S <- example$A - expected_L
  optimum <- sum(svd(expected_L, nu = 0, nv = 0)$d) +
    example$lambda * sum(abs(expected_S))
  maxiter <- if (is.null(example$maxiter)) 100L else example$maxiter
  for (method in c("before", "after")) {
    run <- if (method == "before") reference$rrpca.default else candidate$rrpca.default
    result <- inspect_run(run, example$A, lambda = example$lambda,
                           maxiter = maxiter, tol = 1e-9)
    fit <- result$fit
    row <- data.frame(
      case = case, method = method, lambda = example$lambda, maxiter = maxiter,
      iterations = length(fit$err), final_relative_residual = tail(fit$err, 1L),
      final_stationarity_diagnostic = result$stationarity,
      L_error_over_input_norm = norm(fit$L - expected_L, "F") / norm(example$A, "F"),
      S_error_over_input_norm = norm(fit$S - expected_S, "F") / norm(example$A, "F"),
      objective = sum(svd(fit$L, nu = 0, nv = 0)$d) + example$lambda * sum(abs(fit$S)),
      known_optimum = optimum
    )
    if (method == "after") {
      stopifnot(row$L_error_over_input_norm < 1e-8,
                row$S_error_over_input_norm < 1e-8)
    }
    regression_rows[[length(regression_rows) + 1L]] <- row
  }
}
regression_output <- paste0(tools::file_path_sans_ext(output), "-regressions.csv")
utils::write.csv(do.call(rbind, regression_rows), regression_output, row.names = FALSE)
cat("Saved ", regression_output, "\n", sep = "")
