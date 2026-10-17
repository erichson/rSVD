# Run source-package checks with one selected rsvd installation.
# Vignettes and the PDF manual are excluded; missing optional dependencies are
# reported by R CMD check rather than forcing an early dependency error.
# Usage: Rscript release-downstream-check.R WORK_DIR RSVD_LIBRARY LABEL PKG...
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 4L)
work <- normalizePath(args[[1L]])
rsvd_library <- normalizePath(args[[2L]])
label <- args[[3L]]
packages <- args[-seq_len(3L)]
.libPaths(c(rsvd_library, file.path(work, "dependency-library"), .libPaths()))
stopifnot(find.package("rsvd") == file.path(rsvd_library, "rsvd"))
# Minimal runtime containers may have shared libraries without development
# linker names. When present, use private aliases rather than changing R or
# system libraries. On a normal R development installation this is unnecessary.
private_link <- file.path(work, "native", "link")
if (dir.exists(private_link)) {
  paths <- c(private_link, Sys.getenv("LIBRARY_PATH"))
  Sys.setenv(LIBRARY_PATH = paste(paths[nzchar(paths)], collapse = .Platform$path.sep))
}
inventory <- read.csv(file.path(work, "reverse-dependency-inventory.csv"))
if (identical(packages, "ALL")) packages <- inventory$package
checks <- file.path(work, paste0("checks-", label))
dir.create(checks, recursive = TRUE, showWarnings = FALSE)
writeLines(as.character(Sys.getpid()), file.path(checks, "driver.pid"))
writeLines(capture.output(sessionInfo()), file.path(checks, "session-info.txt"))
results <- list()
for (pkg in packages) {
  if (file.exists(file.path(work, "PAUSE_CHECKS")))
    stop("Checks paused at the caller's request before ", pkg)
  item <- inventory[match(pkg, inventory$package), ]
  archive <- file.path(work, "downloads", basename(item$source_url))
  stopifnot(file.exists(archive))
  output <- file.path(checks, pkg)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  log <- file.path(output, "command.log")
  old <- setwd(output)
  start <- Sys.time()
  status <- system2(file.path(R.home("bin"), "R"),
    c("CMD", "check", "--no-manual", "--no-vignettes", shQuote(archive)),
    stdout = log, stderr = log,
    env = c(paste0("R_LIBS=", shQuote(paste(.libPaths(), collapse = .Platform$path.sep))),
            "_R_CHECK_FORCE_SUGGESTS_=false", "_R_CHECK_CRAN_INCOMING_=false",
            "_R_CHECK_TESTS_NLINES_=0"), timeout = 1200)
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  setwd(old)
  check_log <- file.path(output, paste0(pkg, ".Rcheck"), "00check.log")
  lines <- if (file.exists(check_log)) readLines(check_log, warn = FALSE) else
    readLines(log, warn = FALSE)
  count <- function(pattern) sum(grepl(pattern, lines))
  result <- data.frame(
    package = pkg, version = item$version, relationship = item$relationship,
    rsvd_version = as.character(packageVersion("rsvd")),
    r_version = as.character(getRversion()),
    matrix_version = as.character(packageVersion("Matrix")),
    exit_status = status,
    errors = count("^\\* .* ERROR$"), warnings = count("^\\* .* WARNING$"),
    notes = count("^\\* .* NOTE$"),
    tests = paste(lines[grepl("^\\* checking tests", lines)], collapse = " | "),
    elapsed_seconds = elapsed,
    check_status = paste(lines[grepl("^Status:", lines)], collapse = " | "),
    stringsAsFactors = FALSE)
  results[[pkg]] <- result
  write.csv(do.call(rbind, results), file.path(checks, "summary.csv"), row.names = FALSE)
  cat(pkg, item$version, result$check_status, "exit", status, "\n")
  flush.console()
}
