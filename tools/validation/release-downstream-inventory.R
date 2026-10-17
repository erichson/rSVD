# Snapshot CRAN's direct reverse dependencies and their exact source archives.
# Usage: Rscript release-downstream-inventory.R WORK_DIR [CRAN_REPOSITORY]
# Generated sources and libraries belong outside the package source tree.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L)
work <- normalizePath(args[[1L]], mustWork = FALSE)
cran <- if (length(args) >= 2L) sub("/$", "", args[[2L]]) else
  "https://cran.r-project.org"
dir.create(work, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(work, "downloads"), showWarnings = FALSE)
dir.create(file.path(work, "sources"), showWarnings = FALSE)
metadata <- file.path(work, "PACKAGES.gz")
if (!file.exists(metadata))
  download.file(paste0(cran, "/src/contrib/PACKAGES.gz"), metadata,
                mode = "wb", quiet = TRUE)
db <- read.dcf(gzfile(metadata))
fields <- c("Depends", "Imports", "LinkingTo", "Suggests", "Enhances")
dep_names <- function(value) {
  if (is.na(value) || !nzchar(value)) return(character())
  trimws(sub("\\s*\\(.*$", "", strsplit(value, ",", fixed = TRUE)[[1L]]))
}
rows <- list()
for (i in seq_len(nrow(db))) {
  direct <- fields[vapply(fields, function(field)
    "rsvd" %in% dep_names(db[i, field]), logical(1))]
  if (length(direct)) rows[[length(rows) + 1L]] <- data.frame(
    package = db[i, "Package"], version = db[i, "Version"],
    relationship = paste(direct, collapse = ";"),
    needs_compilation = db[i, "NeedsCompilation"],
    published = db[i, "Published"], md5 = db[i, "MD5sum"],
    source_url = paste0(cran, "/src/contrib/", db[i, "Package"], "_",
                        db[i, "Version"], ".tar.gz"),
    stringsAsFactors = FALSE)
}
inventory <- do.call(rbind, rows)
write.csv(inventory, file.path(work, "reverse-dependency-inventory.csv"),
          row.names = FALSE)
writeLines(c(paste("Metadata downloaded UTC:",
                   format(file.info(metadata)$mtime, tz = "UTC", usetz = TRUE)),
             paste("Inventory generated UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
             paste("Metadata:", paste0(cran, "/src/contrib/PACKAGES.gz")),
             paste("Metadata MD5:", tools::md5sum(metadata)),
             paste("Published rsvd:", db[db[, "Package"] == "rsvd", "Version"])),
           file.path(work, "snapshot.txt"))
saveRDS(db, file.path(work, "PACKAGES.rds"))
for (i in seq_len(nrow(inventory))) {
  item <- inventory[i, ]
  dest <- file.path(work, "downloads", basename(item$source_url))
  if (!file.exists(dest))
    download.file(item$source_url, dest, mode = "wb", quiet = TRUE)
  stopifnot(unname(tools::md5sum(dest)) == item$md5)
  untar(dest, exdir = file.path(work, "sources"))
  cat(item$package, item$version, item$relationship, "verified\n")
}

# Use R's parser so comments, examples in documentation, and the word "rsvd"
# inside diagnostic strings do not count as production call sites.
calls <- list()
public_functions <- c("rsvd", "rpca", "rqb", "rrpca", "rid", "rcur", "id", "cur")
for (i in seq_len(nrow(inventory))) {
  item <- inventory[i, ]
  source_dir <- file.path(work, "sources", item$package, "R")
  for (file in list.files(source_dir, pattern = "\\.[Rr]$", full.names = TRUE)) {
    parsed <- parse(file, keep.source = TRUE)
    tokens <- getParseData(parsed)
    selected <- tokens[tokens$token == "SYMBOL_FUNCTION_CALL" &
                       tokens$text %in% public_functions, , drop = FALSE]
    lines <- readLines(file, warn = FALSE)
    for (j in seq_len(nrow(selected))) calls[[length(calls) + 1L]] <- data.frame(
      package = item$package, version = item$version,
      file = paste(item$package, "R", basename(file), sep = "/"),
      line = selected$line1[[j]], function_name = selected$text[[j]],
      code = trimws(lines[[selected$line1[[j]]]]), stringsAsFactors = FALSE)
  }
}
write.csv(do.call(rbind, calls), file.path(work, "call-sites.csv"), row.names = FALSE)
