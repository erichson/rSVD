# Collect the executed checks and adapter/testthat evidence without calling
# unexecuted or skipped tests successes. Run only after the check drivers finish.
# Usage: Rscript release-downstream-summarize.R WORK OUTPUT_DIR
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
work <- normalizePath(args[[1L]])
output <- normalizePath(args[[2L]], mustWork = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
read_csv <- function(file) read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)

check_files <- Sys.glob(file.path(work, "checks-*", "summary.csv"))
# These preliminary checks preceded the final numerical core fix.
check_files <- check_files[!grepl("checks-candidate-light/", check_files, fixed = TRUE)]
checks <- lapply(check_files, function(file) {
  rows <- read_csv(file)
  label <- sub("^checks-", "", basename(dirname(file)))
  rows$run <- label
  rows$check_completed_utc <- NA_character_
  rows$connected_test_summary <- ""
  rows$connected_test_passed <- NA_integer_
  rows$connected_test_failed <- NA_integer_
  rows$connected_test_skipped <- NA_integer_
  for (i in seq_len(nrow(rows))) {
    root <- file.path(dirname(file), rows$package[[i]],
                       paste0(rows$package[[i]], ".Rcheck"))
    check_log <- file.path(root, "00check.log")
    if (file.exists(check_log))
      rows$check_completed_utc[[i]] <- format(file.info(check_log)$mtime,
                                             tz = "UTC", usetz = TRUE)
    logs <- list.files(file.path(root, "tests"), pattern = "\\.Rout(\\.fail)?$",
                       full.names = TRUE)
    lines <- unlist(lapply(logs, readLines, warn = FALSE), use.names = FALSE)
    summaries <- unique(lines[grepl("^\\[ FAIL [0-9]+ \\| WARN [0-9]+", lines)])
    rows$connected_test_summary[[i]] <- paste(summaries, collapse = " | ")
    if (length(summaries)) {
      number <- function(field)
        sum(as.integer(sub(paste0(".*", field, " ([0-9]+).*"), "\\1", summaries)))
      rows$connected_test_passed[[i]] <- number("PASS")
      rows$connected_test_failed[[i]] <- number("FAIL")
      rows$connected_test_skipped[[i]] <- number("SKIP")
    }
  }
  rows
})
if (length(checks)) {
  checks <- do.call(rbind, checks)
  write.csv(checks, file.path(output, "release-downstream-check-results.csv"), row.names = FALSE)
  # Select the latest completed candidate check for each package. A driver
  # summary's mtime cannot be used: earlier failures may be followed by a
  # successful independent rerun before a long multi-package driver ends.
  final <- checks[checks$rsvd_version == "1.0.7", ]
  final <- final[order(final$package, final$check_completed_utc, decreasing = TRUE), ]
  final <- final[!duplicated(final$package), ]
  final <- final[order(final$relationship, final$package), ]
  write.csv(final, file.path(output, "release-downstream-check-summary.csv"), row.names = FALSE)
}

probe_files <- file.path(work, c("probes-candidate", "probes-baseline"), "summary.csv")
stopifnot(all(file.exists(probe_files)))
probes <- do.call(rbind, lapply(probe_files, read_csv))
write.csv(probes, file.path(output, "release-downstream-probe-results.csv"), row.names = FALSE)
candidate <- readRDS(file.path(work, "probes-candidate", "outputs.rds"))
baseline <- readRDS(file.path(work, "probes-baseline", "outputs.rds"))
stopifnot(setequal(names(candidate), names(baseline)))
comparisons <- do.call(rbind, lapply(names(candidate), function(name) {
  result <- all.equal(candidate[[name]], baseline[[name]],
                       tolerance = 1e-8, check.attributes = TRUE)
  data.frame(probe = name, consistent = isTRUE(result),
               detail = paste(result, collapse = " | "))
}))
write.csv(comparisons, file.path(output, "release-downstream-baseline-comparison.csv"),
            row.names = FALSE)

test_files <- Sys.glob(file.path(work, "explicit-tests-*", "*", "summary.csv"))
tests <- lapply(test_files, function(file) {
  rows <- read_csv(file)
  rows$run <- sub("^explicit-tests-", "", basename(dirname(dirname(file))))
  rows
})
if (length(tests)) write.csv(do.call(rbind, tests),
  file.path(output, "release-downstream-testthat-results.csv"), row.names = FALSE)

file.copy(file.path(work, "call-sites.csv"),
          file.path(output, "release-downstream-call-sites.csv"), overwrite = TRUE)
file.copy(file.path(work, "reverse-dependency-inventory.csv"),
          file.path(output, "release-reverse-dependency-inventory.csv"), overwrite = TRUE)
file.copy(file.path(work, "candidate-source-manifest.json"),
          file.path(output, "release-downstream-source-manifest.json"), overwrite = TRUE)
if (file.exists(file.path(work, "bioconductor-install-status.csv")))
  file.copy(file.path(work, "bioconductor-install-status.csv"),
            file.path(output, "release-bioconductor-dependency-versions.csv"), overwrite = TRUE)
environment_file <- file.path(work, "installed-packages.csv")
if (file.exists(environment_file)) {
  environment <- read_csv(environment_file)
  # installed.packages() follows .libPaths(); retain the first available copy.
  environment <- environment[!duplicated(environment$Package), c("Package", "Version", "Built")]
  environment <- environment[order(environment$Package), ]
  write.csv(environment, file.path(output, "release-downstream-environment.csv"), row.names = FALSE)
}
cat(nrow(probes), "probe executions and", nrow(comparisons), "baseline comparisons collected.\n")
