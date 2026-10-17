# Explicitly execute downstream testthat files, including packages without a
# connected tests/testthat.R runner. Original downloaded sources are unchanged.
# Usage: Rscript release-downstream-testthat.R WORK RSVD_LIB LABEL PACKAGE
#          [--fix-spca-name] [--filter=EXPRESSION]
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 4L)
work <- normalizePath(args[[1L]])
rsvd_library <- normalizePath(args[[2L]])
label <- args[[3L]]
package <- args[[4L]]
.libPaths(c(rsvd_library, file.path(work, "dependency-library"), .libPaths()))
stopifnot(find.package("rsvd") == file.path(rsvd_library, "rsvd"))
path <- file.path(work, "sources", package, "tests", "testthat")
stopifnot(dir.exists(path))
output <- file.path(work, paste0("explicit-tests-", label), package)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
renamed <- "--fix-spca-name" %in% args
if (renamed) {
  stopifnot(package == "sparsepca")
  original <- path
  path <- file.path(output, "compatibility-copy")
  dir.create(path, showWarnings = FALSE)
  edits <- list()
  for (file in list.files(original, pattern = "\\.[Rr]$", full.names = TRUE)) {
    lines <- readLines(file, warn = FALSE)
    changed <- grepl("library(spca)", lines, fixed = TRUE)
    writeLines(gsub("library(spca)", "library(sparsepca)", lines, fixed = TRUE),
               file.path(path, basename(file)))
    edits[[basename(file)]] <- data.frame(file = basename(file),
      substitutions = sum(changed), change = "library(spca) -> library(sparsepca)")
  }
  write.csv(do.call(rbind, edits), file.path(output, "harness-edits.csv"), row.names = FALSE)
}
filter_argument <- args[startsWith(args, "--filter=")]
filter <- if (length(filter_argument)) sub("^--filter=", "", filter_argument[[1L]]) else NULL
set.seed(217031L)
transcript <- capture.output({
  result <- tryCatch(testthat::test_dir(path, package = package, filter = filter,
                      reporter = "summary", stop_on_failure = FALSE,
                      stop_on_warning = FALSE), error = identity)
})
writeLines(transcript, file.path(output, "test-transcript.txt"))
cat(paste(transcript, collapse = "\n"), "\n")
if (inherits(result, "error")) {
  totals <- data.frame(blocks = NA_integer_, assertions = NA_integer_,
                       failed = NA_integer_, errors = NA_integer_,
                       warnings = NA_integer_, skipped = NA_integer_,
                       passed = NA_integer_, outer_error = conditionMessage(result))
} else {
  frame <- as.data.frame(result)
  frame$result <- NULL
  write.csv(frame, file.path(output, "test-results.csv"), row.names = FALSE)
  saveRDS(result, file.path(output, "test-results.rds"))
  totals <- data.frame(blocks = nrow(frame), assertions = sum(frame$nb),
    failed = sum(frame$failed), errors = sum(frame$error),
    warnings = sum(frame$warning), skipped = sum(frame$skipped),
    passed = sum(frame$passed), outer_error = "")
}
summary <- cbind(data.frame(package = package,
  version = as.character(packageVersion(package)),
  rsvd_version = as.character(packageVersion("rsvd")),
  compatibility_name_fix = renamed,
  filter = if (is.null(filter)) "all" else filter,
  files_skipped_outside_tests = sum(grepl("(code run outside of `test_that()`)",
                                         transcript, fixed = TRUE))), totals)
write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
if (nzchar(totals$outer_error) || is.na(totals$failed) || totals$failed || totals$errors)
  quit(status = 1L)
