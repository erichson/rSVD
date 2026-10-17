# Install reverse dependencies in an isolated library.  The selected rsvd
# library must be first; its package must not be replaced by CRAN's release.
# Usage: Rscript release-downstream-install.R WORK_DIR RSVD_LIBRARY [REPOSITORY]
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2L)
work <- normalizePath(args[[1L]])
rsvd_library <- normalizePath(args[[2L]])
library_path <- file.path(work, "dependency-library")
dir.create(library_path, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(rsvd_library, library_path, .libPaths()))
stopifnot(find.package("rsvd") == file.path(rsvd_library, "rsvd"))
repository <- if (length(args) >= 3L) args[[3L]] else
  "https://packagemanager.posit.co/cran/__linux__/noble/latest"
options(repos = c(CRAN = repository), timeout = 180)
# Posit selects its Linux binaries using R's HTTP user agent. Pass the header
# explicitly: current download.file() does not read an HTTPUserAgent option.
headers <- character()
if (grepl("packagemanager.posit.co/cran/__linux__", repository, fixed = TRUE))
  headers <- c("User-Agent" = sprintf("R (%s.%s.0 %s %s %s)",
    R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1L]][[1L]],
    R.version$platform, R.version$arch, R.version$os))
inventory <- read.csv(file.path(work, "reverse-dependency-inventory.csv"))
# These small optional packages are used by unguarded examples or connected
# tests in the reviewed sources. Installing every Suggests closure would add
# unrelated services and large frameworks, so retain the observed minimum.
check_suggests <- c("RcppML", "bench", "reformulas", "GPArotation", "tm",
                    "SnowballC", "wordcloud")
requested <- unique(c(inventory$package, check_suggests))
needed <- requested[!requested %in% installed.packages()[, "Package"]]
selected_versions <- vapply(inventory$package, function(package)
  tryCatch(as.character(packageVersion(package)), error = function(e) NA_character_),
  character(1))
needed <- unique(c(needed, inventory$package[
  is.na(selected_versions) | selected_versions != inventory$version]))
if (length(needed)) install.packages(
  needed, lib = library_path, dependencies = c("Depends", "Imports", "LinkingTo"),
  Ncpus = 1L, INSTALL_opts = "--no-lock", headers = headers)
installed <- as.data.frame(installed.packages()[, c("Package", "Version", "LibPath", "Built")])
write.csv(installed, file.path(work, "installed-packages.csv"), row.names = FALSE)
status <- data.frame(package = inventory$package, expected_version = inventory$version,
                     installed = FALSE, actual_version = NA_character_,
                     load_error = NA_character_)
for (i in seq_len(nrow(status))) {
  pkg <- status$package[[i]]
  tryCatch({
    ns <- loadNamespace(pkg)
    status$installed[[i]] <- TRUE
    status$actual_version[[i]] <- as.character(getNamespaceVersion(ns))
  }, error = function(e) status$load_error[[i]] <<- conditionMessage(e))
}
write.csv(status, file.path(work, "install-status.csv"), row.names = FALSE)
print(status, row.names = FALSE)
if (!all(status$installed & !is.na(status$actual_version) &
         status$actual_version == status$expected_version))
  stop("Downstream installations do not match the metadata snapshot; see install-status.csv")
