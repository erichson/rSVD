# Install a checksum-verified official BiocSingular runtime dependency closure.
# Download and install serially so every source archive is verified before use.
# Usage: Rscript release-downstream-bioconductor.R WORK RSVD_LIBRARY [BIOC_REPO]
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2L)
work <- normalizePath(args[[1L]])
rsvd_library <- normalizePath(args[[2L]])
repository <- if (length(args) >= 3L) sub("/$", "", args[[3L]]) else
  "https://bioconductor.org/packages/3.23/bioc"
library <- file.path(work, "dependency-library")
.libPaths(c(rsvd_library, library, .libPaths()))
stopifnot(find.package("rsvd") == file.path(rsvd_library, "rsvd"))
options(timeout = 180)
metadata <- file.path(work, "Bioc-PACKAGES.gz")
if (!file.exists(metadata))
  download.file(paste0(repository, "/src/contrib/PACKAGES.gz"), metadata,
                mode = "wb", quiet = TRUE)
db <- read.dcf(gzfile(metadata))
rownames(db) <- db[, "Package"]
dep_names <- function(value) {
  if (is.na(value) || !nzchar(value)) return(character())
  trimws(sub("\\s*\\(.*$", "", strsplit(value, ",", fixed = TRUE)[[1L]]))
}
ordered <- character()
visit <- function(package) {
  if (package %in% ordered) return(invisible(NULL))
  dependencies <- unique(unlist(lapply(c("Depends", "Imports", "LinkingTo"),
    function(field) dep_names(db[package, field]))))
  for (dependency in intersect(dependencies, rownames(db))) visit(dependency)
  ordered <<- c(ordered, package)
}
visit("BiocSingular")
records <- list()
for (package in ordered) {
  row <- db[package, ]
  archive <- file.path(work, "downloads", paste0(package, "_", row[["Version"]], ".tar.gz"))
  url <- paste0(repository, "/src/contrib/", basename(archive))
  valid <- function() file.exists(archive) &&
    identical(unname(tools::md5sum(archive)), row[["MD5sum"]])
  if (!valid()) download.file(url, archive, mode = "wb", quiet = TRUE)
  stopifnot(valid())
  # Parsing the archive list catches corrupt downloads before installation.
  stopifnot(paste0(package, "/DESCRIPTION") %in% untar(archive, list = TRUE))
  installed <- file.exists(file.path(library, package, "DESCRIPTION")) &&
    identical(as.character(packageVersion(package, lib.loc = library)), row[["Version"]])
  if (!installed) install.packages(archive, lib = library, repos = NULL,
                                   INSTALL_opts = "--no-lock", Ncpus = 1L)
  stopifnot(identical(as.character(packageVersion(package, lib.loc = library)),
                      row[["Version"]]), requireNamespace(package, quietly = TRUE))
  records[[package]] <- data.frame(package = package, version = row[["Version"]],
    md5 = row[["MD5sum"]], source_url = url, installed = TRUE)
  write.csv(do.call(rbind, records), file.path(work, "bioconductor-install-status.csv"),
            row.names = FALSE)
  cat(package, row[["Version"]], "source checksum verified; loaded\n")
}
