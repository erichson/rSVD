#!/usr/bin/env Rscript

# Developer experiment, not a proposed public API. Real matrices only, and ell
# must be a multiple of zeta: the exact SparseStack recipe in Definition 1.7
# of https://arxiv.org/abs/2508.21189. Each row has one independent signed entry
# in each of zeta equal-width blocks; this is not SparseIID or one CountSketch.
# No claim that constant zeta=4 gives a general high-probability guarantee.
# Rscript tools/validation/sparsestack-experiment.R LIBRARY OUTPUT_DIRECTORY [accuracy|timing]
# The default mode checks accuracy. Timing measures only sketch construction and
# the initial A %*% Omega product, not the full decomposition. Run timed modes
# sequentially with all other numerical benchmarks/tests idle. Allocation data
# describe R-managed allocations, not peak RSS or all native allocations.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(2L, 3L))
mode <- if (length(args) == 3L) args[3L] else "accuracy"
stopifnot(mode %in% c("accuracy", "timing"))
.libPaths(c(normalizePath(args[1L], mustWork = TRUE), .libPaths()))
invisible(loadNamespace("rsvd"))
invisible(loadNamespace("Matrix"))
output <- args[2L]
dir.create(output, recursive = TRUE, showWarnings = FALSE)

sparse_stack <- function(n, ell, zeta = 4L) {
  stopifnot(ell >= zeta, ell %% zeta == 0L)
  block_width <- ell %/% zeta
  rows <- rep.int(seq_len(n), zeta)
  columns <- rep.int(seq.int(0L, ell - block_width, by = block_width),
                     rep.int(n, zeta)) +
    sample.int(block_width, n * zeta, replace = TRUE)
  signs <- (2 * sample.int(2L, n * zeta, replace = TRUE) - 3) / sqrt(zeta)
  Matrix::sparseMatrix(i = rows, j = columns, x = signs, dims = c(n, ell))
}

draw <- function(n, ell, method) {
  if (method == "gaussian") matrix(rnorm(n * ell), n, ell) else sparse_stack(n, ell)
}

sketch <- function(A, ell, method) as.matrix(A %*% draw(ncol(A), ell, method))

prototype_svd <- function(A, k, ell, q, method) {
  flipped <- nrow(A) < ncol(A)
  if (flipped) A <- Matrix::t(A)
  Y <- sketch(A, ell, method)
  for (iteration in seq_len(q)) {
    Q <- qr.Q(qr(Y))
    Z <- qr.Q(qr(as.matrix(Matrix::crossprod(A, Q))))
    Y <- as.matrix(A %*% Z)
  }
  Q <- qr.Q(qr(Y))
  B <- as.matrix(Matrix::crossprod(Q, A))
  fit <- svd(B, nu = k, nv = k)
  fit$d <- fit$d[seq_len(k)]
  fit$u <- Q %*% fit$u
  if (flipped) {
    temporary <- fit$u
    fit$u <- fit$v
    fit$v <- temporary
  }
  fit
}

basis <- function(size, rank, coherent = FALSE) {
  if (coherent) return(diag(size)[, seq_len(rank), drop = FALSE])
  Q <- sqrt(2 / size) * cos(pi * outer(seq_len(size)-0.5,
                                      seq_len(rank)-1) / size)
  Q[, 1L] <- Q[, 1L] / sqrt(2)
  Q
}

fixture <- function(m, n, d, coherent = FALSE) {
  basis(m, length(d), coherent) %*%
    (d * t(basis(n, length(d), coherent)))
}

frobenius <- function(A) sqrt(sum(Mod(A)^2))
residual <- function(A, fit) frobenius(as.matrix(A) - fit$u %*% (fit$d * t(fit$v)))

run_accuracy <- function() {
  # Validate the experimental driver's Gaussian branch against the package.
  for (A in list(fixture(160, 100, 10:1), fixture(100, 160, 10:1),
                  Matrix::Diagonal(100, c(10:1, rep(0, 90))))) {
    for (q in c(0L, 2L)) {
      set.seed(591)
      expected <- rsvd::rsvd(A, k = 10, p = 10, q = q)
      set.seed(591)
      actual <- prototype_svd(A, 10L, 20L, q, "gaussian")
      stopifnot(max(abs(actual$d - expected$d)) < 1e-10,
                frobenius(actual$u %*% t(actual$u) -
                            expected$u %*% t(expected$u)) < 1e-9)
    }
  }

  accuracy <- list()
  for (k in c(10L, 22L)) {
    ell <- k + 10L
    for (coherent in c(FALSE, TRUE)) {
      for (spectrum in c("exact", "geometric", "slow", "gap")) {
        d <- switch(spectrum, exact = seq(k, 1),
                     geometric = 10 * 0.85^(0:99),
                     slow = 1 / sqrt(seq_len(100)),
                     gap = c(rep(1, k), rep(0.01, 100-k)))
        A <- fixture(160, 100, d, coherent)
        optimal <- if (length(d) > k) sqrt(sum(d[-seq_len(k)]^2)) else 0
        for (q in c(0L, 2L)) {
          for (method in c("gaussian", "sparsestack4")) {
            for (seed in seq_len(20L)) {
              set.seed(seed)
              fit <- prototype_svd(A, k, ell, q, method)
              error <- residual(A, fit)
              accuracy[[length(accuracy) + 1L]] <- data.frame(
                m = nrow(A), n = ncol(A), k = k, ell = ell,
                coherence = if (coherent) "coordinate" else "smooth",
                spectrum = spectrum, q = q, method = method, seed = seed,
                relative_residual = error / frobenius(A),
                error_over_optimal = if (optimal > 0) error / optimal else NA_real_)
            }
          }
        }
      }
    }
  }
  accuracy <- do.call(rbind, accuracy)
  write.csv(accuracy, file.path(output, "accuracy.csv"), row.names = FALSE)
  print(aggregate(relative_residual ~ k + coherence + spectrum + q + method,
                    accuracy, max))

  # An edge case: width four, rank two. A sparse sign sketch becomes a scaled
  # dense sign matrix when each block has width one, with a nonzero probability
  # of exactly colliding rows. Power iteration may rescue some failures through
  # the numerical completion of a deficient intermediate basis; record q=0/2.
  edge_A <- basis(20, 2) %*% (c(2, 1) * t(basis(4, 2, coherent = TRUE)))
  edge_rows <- list()
  for (method in c("gaussian", "sparsestack4")) {
    for (q in c(0L, 2L)) {
      for (seed in seq_len(500L)) {
        set.seed(seed)
        fit <- prototype_svd(edge_A, 2L, 4L, q, method)
        edge_rows[[length(edge_rows) + 1L]] <- data.frame(
          method = method, q = q, seed = seed,
          relative_residual = residual(edge_A, fit) / frobenius(edge_A))
      }
    }
  }
  edges <- do.call(rbind, edge_rows)
  write.csv(edges, file.path(output, "small-width.csv"), row.names = FALSE)
  print(aggregate(relative_residual > 1e-9 ~ method + q, edges, sum))

  # Keep compact summaries in the repository; full per-seed data are generated
  # here so every reported maximum and failure count can be reproduced.
  metric_summary <- function(x) {
    if (all(is.na(x))) return(c(minimum = NA, median = NA, maximum = NA))
    c(minimum = min(x), median = median(x), maximum = max(x))
  }
  keys <- c("m", "n", "k", "ell", "coherence", "spectrum", "q", "method")
  groups <- split(accuracy, interaction(accuracy[keys], drop = TRUE))
  summary <- do.call(rbind, lapply(groups, function(group) {
    relative <- metric_summary(group$relative_residual)
    optimal <- metric_summary(group$error_over_optimal)
    data.frame(group[1L, keys], seeds = nrow(group), seed_first = min(group$seed),
                 seed_last = max(group$seed),
                 relative_error_min = relative[1L], relative_error_median = relative[2L],
                 relative_error_max = relative[3L],
                 error_over_optimal_min = optimal[1L],
                 error_over_optimal_median = optimal[2L],
                 error_over_optimal_max = optimal[3L],
                 R_version = as.character(getRversion()),
                 Matrix_version = as.character(utils::packageVersion("Matrix")))
  }))
  rownames(summary) <- NULL
  write.csv(summary, file.path(output, "sparsestack-accuracy-summary.csv"), row.names = FALSE)
  edge_groups <- split(edges, interaction(edges[c("method", "q")], drop = TRUE))
  edge_summary <- do.call(rbind, lapply(edge_groups, function(group) {
    data.frame(method = group$method[1L], q = group$q[1L], m = 20L, n = 4L,
                 k = 2L, ell = 4L, seeds = nrow(group), seed_first = min(group$seed),
                 seed_last = max(group$seed), failure_threshold = 1e-9,
                 failures = sum(group$relative_residual > 1e-9),
                 relative_error_max = max(group$relative_residual),
                 R_version = as.character(getRversion()),
                 Matrix_version = as.character(utils::packageVersion("Matrix")))
  }))
  write.csv(edge_summary, file.path(output, "sparsestack-small-width-summary.csv"),
              row.names = FALSE)
  complex_operations <- list(
    complex_sparse_construction = function() {
      Matrix::sparseMatrix(i = 1L, j = 1L, x = 1+1i, dims = c(3L, 2L))
    },
    complex_dense_by_real_sparse_product = function() {
      matrix(1:12+1i, 4, 3) %*% Matrix::sparseMatrix(
        i = 1:3, j = c(1, 2, 1), x = c(1, -1, 1), dims = c(3L, 2L))
    })
  complex_results <- data.frame(
    operation = names(complex_operations),
    result = vapply(complex_operations, function(operation) {
      tryCatch({ operation(); "supported" }, error = function(error) conditionMessage(error))
    }, character(1)),
    R_version = as.character(getRversion()),
    Matrix_version = as.character(utils::packageVersion("Matrix")))
  write.csv(complex_results, file.path(output, "sparsestack-complex-support.csv"),
              row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output, "session-info.txt"))
  invisible(summary)
}

sparse_fixture <- function() {
  Matrix::bdiag(lapply(seq_len(40L), function(block) {
    left <- sin(seq_len(50L) * (block+1) / 19)
    right <- cos(seq_len(25L) * (block+2) / 13)
    left <- left / sqrt(sum(left^2))
    right <- right / sqrt(sum(right^2))
    Matrix::Matrix((41-block) * tcrossprod(left, right), sparse = TRUE)
  }))
}

measure <- function(run, profile) {
  set.seed(592)
  warmup <- run()
  rm(warmup)
  elapsed <- numeric(3L)
  for (repetition in seq_len(3L)) {
    gc()
    set.seed(592)
    elapsed[repetition] <- system.time(result <- run(), gcFirst = FALSE)[["elapsed"]]
    rm(result)
  }
  gc()
  set.seed(592)
  utils::Rprofmem(profile)
  result <- run()
  utils::Rprofmem(NULL)
  records <- readLines(profile, warn = FALSE)
  records <- records[grepl("^[0-9]+[[:space:]]*:", records)]
  sizes <- as.numeric(sub("^([0-9]+).*", "\\1", records))
  list(median_seconds = median(elapsed), min_seconds = min(elapsed),
       max_seconds = max(elapsed),
       r_allocation_total_bytes = sum(sizes),
       r_allocation_largest_bytes = max(sizes))
}

run_timing <- function() {
  timings <- list()
  for (case in c("dense_tall", "dense_wide", "sparse_tall", "sparse_wide")) {
    A <- if (grepl("dense", case)) fixture(1800, 1000, 40:1) else sparse_fixture()
    if (grepl("wide", case)) A <- Matrix::t(A)
    for (method in c("gaussian", "sparsestack4")) {
      for (stage in c("sample", "sample_and_AX")) {
        run <- switch(stage,
          sample = function() draw(ncol(A), 80L, method),
          sample_and_AX = function() sketch(A, 80L, method))
        metrics <- measure(run, file.path(output, paste0(case, "-", method, "-", stage,
                                                         ".Rprofmem")))
        row <- as.data.frame(c(list(case = case, method = method, stage = stage,
                                     m = nrow(A), n = ncol(A), k = 40L, ell = 80L,
                                     R_version = as.character(getRversion()),
                                     Matrix_version = as.character(utils::packageVersion("Matrix")),
                                     BLAS = basename(unname(extSoftVersion()["BLAS"])),
                                     seed = 592L, repetitions = 3L), metrics))
        timings[[length(timings) + 1L]] <- row
        write.csv(do.call(rbind, timings), file.path(output, "timings.csv"),
                    row.names = FALSE)
        cat(sprintf("%s / %s / %s: %.3f s, %.3f MiB total allocations\n",
                    case, method, stage, metrics$median_seconds,
                    metrics$r_allocation_total_bytes / 1024^2))
      }
    }
  }
  writeLines(capture.output(sessionInfo()), file.path(output, "session-info.txt"))
}

if (mode == "accuracy") run_accuracy() else run_timing()
