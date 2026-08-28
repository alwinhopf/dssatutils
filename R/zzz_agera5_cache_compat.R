# AgERA5 cache compatibility helpers
#
# ecmwfr inspects the returned CDS media type and may normalize a requested
# target such as `foo.zip.partial` to `foo.zip.zip`. Older versions of the
# AgERA5 downloader used that `.partial` suffix, then looked only for the
# original temporary name. A successful CDS transfer could therefore be
# treated as missing and requested again.
#
# This file is deliberately collated after weather_agera5.R (no explicit
# Collate field is used by the package), so the corrected downloader below
# replaces the legacy implementation while remaining backward-compatible with
# already-downloaded `.zip.zip` artifacts.

.agera5_valid_zip <- function(path) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0) {
    return(FALSE)
  }
  listing <- try(utils::unzip(path, list = TRUE), silent = TRUE)
  !inherits(listing, "try-error") && nrow(listing) > 0 &&
    any(grepl("\\.nc$", listing$Name, ignore.case = TRUE) & listing$Length > 0)
}

.agera5_recover_cache_zip <- function(zip_dest) {
  # Canonical cache already exists and is healthy.
  if (.agera5_valid_zip(zip_dest)) return(zip_dest)

  # Remove a corrupt canonical cache before promoting a valid legacy artifact.
  if (file.exists(zip_dest)) unlink(zip_dest)

  # `.zip.zip` is the filename observed when ecmwfr normalized a requested
  # `*.zip.partial` target. The other names cover related interrupted/older
  # variants and are cheap to check.
  candidates <- unique(c(
    paste0(zip_dest, ".zip"),
    paste0(zip_dest, ".partial"),
    paste0(zip_dest, ".partial.zip")
  ))
  candidates <- candidates[vapply(candidates, .agera5_valid_zip, logical(1))]
  if (!length(candidates)) return(NA_character_)

  source <- candidates[[1]]
  moved <- file.rename(source, zip_dest)
  if (!moved) {
    moved <- file.copy(source, zip_dest, overwrite = TRUE)
    if (moved) unlink(source)
  }
  if (!moved || !.agera5_valid_zip(zip_dest)) return(NA_character_)

  message(sprintf("  AgERA5 recovered legacy cache artifact: %s -> %s",
                  basename(source), basename(zip_dest)))
  zip_dest
}

.agera5_download_job <- function(job) {
  # First salvage any successful download that ecmwfr previously named
  # `*.zip.zip`, so reruns do not resubmit an already-completed CDS job.
  .agera5_recover_cache_zip(job$zip_dest)

  data_files <- .agera5_data_files(job$nc_dest, job$zip_dest, job$unzip_dir)
  if (length(data_files)) {
    return(list(ok = TRUE, cached = TRUE, job = job, data_files = data_files,
                message = sprintf("  AgERA5 cache hit (%s)", job$tag)))
  }

  .agera5_ensure_ecmwfr_key(quiet = TRUE)

  # ecmwfr already stages downloads internally before moving them into place.
  # Request the final `.zip` filename directly. Supplying `*.zip.partial` is
  # what caused ecmwfr to rewrite the target to `*.zip.zip`.
  if (file.exists(job$zip_dest)) unlink(job$zip_dest)
  req <- list(dataset_short_name = "sis-agrometeorological-indicators",
              variable = job$spec$var, year = as.character(job$yr),
              month = sprintf("%02d", 1:12), day = sprintf("%02d", 1:31),
              area = job$area, version = "2_0",
              target = basename(job$zip_dest))
  if (!is.na(job$spec$sel_kind)) req[[job$spec$sel_kind]] <- job$spec$sel

  err <- NULL
  returned <- tryCatch(
    ecmwfr::wf_request(request = req, path = job$cache_dir),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  # Current ecmwfr should now leave the canonical target in place. Still be
  # defensive about future/media-type filename normalization and about legacy
  # behavior by checking both the returned path and known duplicate suffixes.
  if (is.character(returned) && length(returned)) {
    returned <- returned[nzchar(returned)]
    for (candidate in returned) {
      candidate_path <- if (file.exists(candidate)) candidate else file.path(job$cache_dir, basename(candidate))
      if (.agera5_valid_zip(candidate_path) && normalizePath(candidate_path, mustWork = FALSE) !=
          normalizePath(job$zip_dest, mustWork = FALSE)) {
        if (file.exists(job$zip_dest)) unlink(job$zip_dest)
        moved <- file.rename(candidate_path, job$zip_dest)
        if (!moved) {
          moved <- file.copy(candidate_path, job$zip_dest, overwrite = TRUE)
          if (moved) unlink(candidate_path)
        }
        if (moved) break
      }
    }
  }
  .agera5_recover_cache_zip(job$zip_dest)

  data_files <- .agera5_data_files(job$nc_dest, job$zip_dest, job$unzip_dir)
  if (length(data_files)) {
    return(list(ok = TRUE, cached = FALSE, job = job, data_files = data_files,
                message = sprintf("  AgERA5 downloaded (%s)", job$tag)))
  }

  list(ok = FALSE, cached = FALSE, job = job, data_files = character(),
       message = sprintf("  AgERA5 download failed (%s): %s",
                         job$tag, if (is.null(err)) "no data file returned" else err))
}
