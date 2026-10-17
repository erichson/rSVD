#!/usr/bin/env Rscript

# Run from the package root:
# Rscript tools/validation/qr-experiment.R /path/to/installed/rsvd-library output.csv [large]
#
# Compares the current source implementation with the full dense base QR used
# at revision ad8f05a. Each path has one warmup, three timed runs, and a separate
# Rprofmem profile. All fixtures and numerical checks are outside measurements.
# Rprofmem reports R allocations, not peak RSS or all native allocations.
# This is a developer experiment; elapsed times are not unit-test assertions.
# Sparse fixtures default to 2000x1000 and 1000x2000. The optional "large"
# argument uses 4000x2000 and 2000x4000 and makes the full-QR reference slower.

arguments <- commandArgs(trailingOnly = TRUE)
if (!length(arguments) %in% c(2L, 3L)) {
  stop("Usage: Rscript tools/validation/qr-experiment.R LIBRARY OUTPUT.csv [large]")
}
large <- length(arguments) == 3L && identical(arguments[3L], "large")
if (length(arguments) == 3L && !large) stop("Unknown third argument")
block_rows <- if (large) 200L else 100L
block_columns <- if (large) 100L else 50L
.libPaths(c(normalizePath(arguments[1L], mustWork = TRUE), .libPaths()))
invisible(loadNamespace("rsvd"))
invisible(loadNamespace("Matrix"))
candidate <- new.env(parent = asNamespace("rsvd"))
sys.source("R/rqb.R", envir = candidate)

reference_qb <- function(A, k) {
  A <- as.matrix(A)
  identity <- diag(if (is.complex(A)) 1+0i else 1,
                   nrow = nrow(A), ncol = k)
  Q <- qr.qy(qr(A), identity)
  B <- if (is.complex(A)) crossprod(Conj(Q), A) else crossprod(Q, A)
  list(Q = Q, B = B)
}

dense_fixture <- function() {
  m <- 2500L
  n <- 1000L
  k <- 20L
  U <- sqrt(2 / m) * cos(pi * outer(seq_len(m) - 0.5,
                                    seq_len(k) - 1) / m)
  U[, 1L] <- U[, 1L] / sqrt(2)
  coefficients <- cbind(diag(k), matrix(sin(seq_len(k * (n-k))), k, n-k))
  U %*% (seq_len(k) * coefficients)
}

sparse_fixture <- function() {
  blocks <- lapply(seq_len(20L), function(block) {
    left <- sin(seq_len(block_rows) * (block + 1) / 19) +
      cos(seq_len(block_rows) * (block + 2) / 23)
    right <- cos(seq_len(block_columns) * (block + 3) / 17) +
      sin(seq_len(block_columns) * (block + 1) / 29)
    left <- left / sqrt(sum(left^2))
    right <- right / sqrt(sum(right^2))
    Matrix::Matrix((21-block) * tcrossprod(left, right), sparse = TRUE)
  })
  A <- Matrix::bdiag(blocks)
  leading <- seq.int(1L, ncol(A), by = block_columns)
  # One independent column from each rank-one block comes first. This is the
  # documented case where the partial QR can retain the full-QR basis exactly.
  A[, c(leading, setdiff(seq_len(ncol(A)), leading)), drop = FALSE]
}

frobenius <- function(A) sqrt(sum(Mod(A)^2))
relative_residual <- function(A, fit) {
  error <- total <- 0
  for (first in seq.int(1L, nrow(A), by = 256L)) {
    rows <- seq.int(first, min(first + 255L, nrow(A)))
    reference <- as.matrix(A[rows, , drop = FALSE])
    difference <- reference - fit$Q[rows, , drop = FALSE] %*% fit$B
    error <- error + sum(Mod(difference)^2)
    total <- total + sum(Mod(reference)^2)
  }
  sqrt(error / total)
}

output <- arguments[2L]
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
profile_directory <- paste0(tools::file_path_sans_ext(output), "-profiles")
dir.create(profile_directory, recursive = TRUE, showWarnings = FALSE)
rows <- list()
for (case in c("dense_tall", "sparse_tall", "sparse_wide")) {
  A <- if (case == "dense_tall") dense_fixture() else sparse_fixture()
  if (case == "sparse_wide") {
    # After transposition, reorder columns to expose independent directions.
    A <- Matrix::t(A)
    leading <- seq.int(1L, ncol(A), by = block_rows)
    A <- A[, c(leading, setdiff(seq_len(ncol(A)), leading)), drop = FALSE]
  }
  expected <- reference_qb(A, 20L)
  for (method in c("full_base_qr", "partial_base_qr")) {
    run <- if (method == "full_base_qr") reference_qb else {
      function(A, k) candidate$rqb.default(A, k = k, rand = FALSE)
    }
    warmup <- run(A, 20L)
    stopifnot(frobenius(warmup$Q - expected$Q) < 1e-10,
              frobenius(warmup$B - expected$B) < 1e-9)
    rm(warmup)
    elapsed <- numeric(3L)
    for (repetition in seq_len(3L)) {
      gc()
      elapsed[repetition] <- system.time(fit <- run(A, 20L),
                                         gcFirst = FALSE)[["elapsed"]]
      rm(fit)
    }
    profile <- file.path(profile_directory, paste0(case, "-", method, ".Rprofmem"))
    gc()
    utils::Rprofmem(profile)
    fit <- run(A, 20L)
    utils::Rprofmem(NULL)
    lines <- readLines(profile, warn = FALSE)
    allocations <- lines[grepl("^[0-9]+[[:space:]]*:", lines)]
    sizes <- as.numeric(sub("^([0-9]+).*", "\\1", allocations))
    row <- data.frame(case = case, method = method,
                      m = nrow(A), n = ncol(A), k = 20L,
                      input_class = class(A)[1L],
                      R_version = as.character(getRversion()),
                      Matrix_version = as.character(utils::packageVersion("Matrix")),
                      BLAS = basename(unname(extSoftVersion()["BLAS"])),
                      median_seconds = median(elapsed),
                      min_seconds = min(elapsed), max_seconds = max(elapsed),
                      repetitions = 3L,
                      r_allocation_total_bytes = sum(sizes),
                      r_allocation_largest_bytes = max(sizes),
                      relative_frobenius_residual = relative_residual(A, fit),
                      q_difference_frobenius = frobenius(fit$Q - expected$Q),
                      b_difference_frobenius = frobenius(fit$B - expected$B))
    stopifnot(row$relative_frobenius_residual < 1e-10)
    rows[[length(rows) + 1L]] <- row
    utils::write.csv(do.call(rbind, rows), output, row.names = FALSE)
    cat(sprintf("%s / %s: %.3f s, %.3f MiB total R allocations, %.3f MiB largest\n",
                case, method, row$median_seconds,
                row$r_allocation_total_bytes / 1024^2,
                row$r_allocation_largest_bytes / 1024^2))
    rm(fit)
  }
  rm(A, expected)
}
cat("Saved ", output, "\n", sep = "")
