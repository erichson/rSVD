#!/usr/bin/env Rscript
# Compute the actual image approximations used by render-tiger.py.
# This is an illustration, not a timing benchmark. All reconstruction errors
# are measured on the unmodified numeric matrices, before display clipping.
#
# From the repository root, after installing the candidate package:
#   Rscript tools/benchmarks/illustrate-tiger.R \
#     --results tools/benchmarks/tiger-results --images /tmp/rsvd-tiger
# Optional: --library /path/to/private/R/library --repo /path/to/rSVD

args <- commandArgs(trailingOnly = TRUE)
options <- list(results = "tools/benchmarks/tiger-results",
                images = file.path(tempdir(), "rsvd-tiger"),
                library = NULL, repo = ".")
if (length(args) %% 2L != 0L) {
  stop("Supply options as --name value pairs.")
}
if (length(args)) {
  for (i in seq.int(1L, length(args), by = 2L)) {
    key <- sub("^--", "", args[i])
    if (!startsWith(args[i], "--") || !key %in% names(options)) {
      stop("Unknown option: ", args[i])
    }
    options[key] <- args[i + 1L]
  }
}
if (!is.null(options$library)) .libPaths(c(options$library, .libPaths()))
if (!requireNamespace("rsvd", quietly = TRUE)) stop("Install rsvd first.")
if (!requireNamespace("Matrix", quietly = TRUE)) stop("Install Matrix first.")
if (getRversion() < "4.4.0" || packageVersion("Matrix") < "1.7-0") {
  stop("This illustration uses the supported R >= 4.4 / Matrix >= 1.7-0 runtime.")
}

dir.create(options$results, recursive = TRUE, showWarnings = FALSE)
dir.create(options$images, recursive = TRUE, showWarnings = FALSE)
repo <- normalizePath(options$repo, mustWork = TRUE)
dataset_file <- file.path(repo, "data", "tiger.RData")
dataset <- new.env(parent = emptyenv())
load(dataset_file, envir = dataset)
A <- dataset$tiger
stopifnot(is.matrix(A), is.numeric(A), identical(dim(A), c(1600L, 1200L)),
          all(is.finite(A)), min(A) >= 0, max(A) <= 1)
m <- nrow(A)
n <- ncol(A)
original_norm <- norm(A, "F")
seed <- 42L
p <- 10L
q <- 2L
ranks <- c(10L, 50L, 150L)

write_matrix <- function(A, filename) {
  target <- file.path(options$images, filename)
  connection <- file(target, open = "wb")
  on.exit(close(connection))
  writeBin(as.double(A), connection, size = 8L, endian = "little")
  invisible(target)
}

write_matrix(A, "original-f64.bin")
rows <- list(data.frame(
  panel = "original", m = m, n = n, k = NA_integer_,
  seed = NA_integer_, p = NA_integer_, q = NA_integer_,
  relative_frobenius_error = 0, factor_entries = NA_real_,
  original_entries = as.double(m) * n, factor_entry_fraction = 1,
  original_entries_per_factor_entry = 1,
  min_value = min(A), max_value = max(A),
  fraction_outside_display_limits = 0,
  numeric_file = "original-f64.bin"
))

for (k in ranks) {
  # Independent, reproducible fits. In particular, k = 150 exactly follows
  # the README's set.seed(42); rsvd(tiger, k = 150) example.
  set.seed(seed)
  fit <- rsvd::rsvd(A, k = k, p = p, q = q)
  stopifnot(identical(dim(fit$u), c(m, k)),
            identical(dim(fit$v), c(n, k)), length(fit$d) == k,
            all(is.finite(fit$u)), all(is.finite(fit$d)),
            all(is.finite(fit$v)))
  approximation <- tcrossprod(sweep(fit$u, 2L, fit$d, "*"), fit$v)
  stopifnot(identical(dim(approximation), dim(A)),
            all(is.finite(approximation)))
  residual <- norm(A - approximation, "F") / original_norm
  factor_entries <- as.double(k) * (m + n + 1)
  filename <- sprintf("rank-%03d-f64.bin", k)
  write_matrix(approximation, filename)
  rows[[length(rows) + 1L]] <- data.frame(
    panel = paste0("rank_", k), m = m, n = n, k = k,
    seed = seed, p = p, q = q,
    relative_frobenius_error = residual,
    factor_entries = factor_entries,
    original_entries = as.double(m) * n,
    factor_entry_fraction = factor_entries / (as.double(m) * n),
    original_entries_per_factor_entry = as.double(m) * n / factor_entries,
    min_value = min(approximation), max_value = max(approximation),
    fraction_outside_display_limits = mean(approximation < 0 | approximation > 1),
    numeric_file = filename
  )
  message(sprintf("k = %d: error %.9f; factor entries %.9f%%", k,
                  residual, 100 * factor_entries / (as.double(m) * n)))
}
summary <- do.call(rbind, rows)
write.csv(summary, file.path(options$results, "summary.csv"),
          row.names = FALSE, na = "")

git <- function(...) {
  output <- tryCatch(system2("git", c("-C", shQuote(repo), ...),
                            stdout = TRUE, stderr = FALSE),
                     error = function(e) "unavailable")
  paste(output, collapse = "\n")
}
source_files <- c("R/rsvd.R", "R/rqb.R", "R/wrapper_function.R")
source_hashes <- tools::md5sum(file.path(repo, source_files))
metadata <- c(
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  package = as.character(packageVersion("rsvd")),
  R = R.version.string,
  Matrix = as.character(packageVersion("Matrix")),
  BLAS = unname(extSoftVersion()["BLAS"]),
  RNG = paste(RNGkind(), collapse = "; "),
  seed = as.character(seed),
  ranks = paste(ranks, collapse = "; "),
  p = as.character(p), q = as.character(q),
  source_commit = git("rev-parse", "HEAD"),
  dataset = "data/tiger.RData",
  dataset_md5 = unname(tools::md5sum(dataset_file)),
  setNames(as.character(source_hashes), paste0("md5_", source_files)),
  numeric_format = "little-endian IEEE-754 float64; R column-major order",
  display_orientation = "transpose, with origin at the lower left; equivalent to image(tiger)",
  display_limits = "0 to 1 for every panel; values outside this interval clipped only for display",
  residual = "norm(tiger - reconstruction, 'F') / norm(tiger, 'F'), before clipping",
  factor_entry_fraction = "k * (m + n + 1) / (m * n); numeric entry count, not image-file bytes",
  image_credit = "S. Taheri (2006), edited by Fir0002; bundled grayscale conversion; rsvd reconstructions",
  image_source = "https://commons.wikimedia.org/wiki/File:Sibirischer_tiger_de_edit02.jpg",
  image_license = "CC BY-SA 2.5",
  image_license_url = "https://creativecommons.org/licenses/by-sa/2.5/",
  image_license_accessed = "2026-10-17"
)
write.csv(data.frame(key = names(metadata), value = unname(metadata)),
          file.path(options$results, "metadata.csv"), row.names = FALSE)
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
           file.path(options$results, "session-info.txt"))
message("Numeric image intermediates: ", normalizePath(options$images))
