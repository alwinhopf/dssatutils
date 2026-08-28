# Final AgERA5 download hardening.
#
# ecmwfr 2.x normalizes a target's extension to the extension of the returned
# CDS download URL.  Therefore targets such as `foo.zip.partial` and
# `foo.csv.partial` become `foo.zip.zip` / `foo.csv.csv`.  The canonical target
# must be requested directly.  These overrides are intentionally self-contained
# so `.agera5_download_job` also works when exported to a PSOCK worker: the
# worker does not need additional cache-helper functions exported separately.

.agera5_download_job <- function(job) {
  valid_zip <- function(path) {
    if (!file.exists(path)) return(FALSE)
    info <- file.info(path)
    if (is.na(info$size) || info$size <= 0) return(FALSE)
    listing <- try(utils::unzip(path, list = TRUE), silent = TRUE)
    !inherits(listing, "try-error") && nrow(listing) > 0 &&
      any(grepl("\\.nc$", listing$Name, ignore.case = TRUE) & listing$Length > 0)
  }

  promote <- function(source, dest) {
    if (!valid_zip(source)) return(FALSE)
    if (normalizePath(source, mustWork = FALSE) == normalizePath(dest, mustWork = FALSE)) {
      return(TRUE)
    }
    if (file.exists(dest)) unlink(dest)
    moved <- suppressWarnings(file.rename(source, dest))
    if (!moved) {
      moved <- file.copy(source, dest, overwrite = TRUE)
      if (moved) unlink(source)
    }
    isTRUE(moved) && valid_zip(dest)
  }

  # Recover successful transfers made by the legacy `*.zip.partial` request.
  if (!valid_zip(job$zip_dest)) {
    if (file.exists(job$zip_dest)) unlink(job$zip_dest)
    legacy <- unique(c(
      paste0(job$zip_dest, ".zip"),          # observed: foo.zip.zip
      paste0(job$zip_dest, ".partial"),
      paste0(job$zip_dest, ".partial.zip")
    ))
    for (candidate in legacy) {
      if (promote(candidate, job$zip_dest)) {
        message(sprintf("  AgERA5 recovered legacy cache artifact: %s -> %s",
                        basename(candidate), basename(job$zip_dest)))
        break
      }
    }
  }

  data_files <- .agera5_data_files(job$nc_dest, job$zip_dest, job$unzip_dir)
  if (length(data_files)) {
    return(list(ok = TRUE, cached = TRUE, job = job, data_files = data_files,
                message = sprintf("  AgERA5 cache hit (%s)", job$tag)))
  }

  .agera5_ensure_ecmwfr_key(quiet = TRUE)

  # Ask ecmwfr for the canonical ZIP path.  Its own downloader already stages
  # to a temporary file before moving the completed transfer into place.
  if (file.exists(job$zip_dest)) unlink(job$zip_dest)
  req <- list(
    dataset_short_name = "sis-agrometeorological-indicators",
    variable = job$spec$var,
    year = as.character(job$yr),
    month = sprintf("%02d", 1:12),
    day = sprintf("%02d", 1:31),
    area = job$area,
    version = "2_0",
    target = basename(job$zip_dest)
  )
  if (!is.na(job$spec$sel_kind)) req[[job$spec$sel_kind]] <- job$spec$sel

  err <- NULL
  returned <- tryCatch(
    ecmwfr::wf_request(request = req, path = job$cache_dir),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  # wf_request() returns the actual downloaded filename.  Prefer that path and
  # normalize it into our canonical cache name if ecmwfr changed it anyway.
  if (is.character(returned) && length(returned)) {
    for (candidate in returned[nzchar(returned)]) {
      candidate_path <- if (file.exists(candidate)) {
        candidate
      } else {
        file.path(job$cache_dir, basename(candidate))
      }
      if (promote(candidate_path, job$zip_dest)) break
    }
  }

  # Defensive compatibility with old/new ecmwfr extension normalization.
  if (!valid_zip(job$zip_dest)) {
    for (candidate in c(paste0(job$zip_dest, ".zip"),
                        paste0(job$zip_dest, ".partial.zip"),
                        paste0(job$zip_dest, ".partial"))) {
      if (promote(candidate, job$zip_dest)) break
    }
  }

  data_files <- .agera5_data_files(job$nc_dest, job$zip_dest, job$unzip_dir)
  if (length(data_files)) {
    return(list(ok = TRUE, cached = FALSE, job = job, data_files = data_files,
                message = sprintf("  AgERA5 downloaded (%s)", job$tag)))
  }

  list(ok = FALSE, cached = FALSE, job = job, data_files = character(),
       message = sprintf("  AgERA5 download failed (%s): %s",
                         job$tag,
                         if (is.null(err)) "no valid NetCDF archive returned" else err))
}

.agera5_download_timeseries_job <- function(job) {
  data_format <- tolower(job$data_format)
  if (data_format != "csv") {
    stop("AgERA5 time-series backend currently supports data_format='csv'.", call. = FALSE)
  }

  bounds <- .agera5_date_bounds_for_year(job$year)
  if (is.null(bounds)) return(NULL)
  dest <- .agera5_timeseries_cache_path(job$cache_dir, job$year, job$area, data_format)

  valid_csv <- function(path) {
    if (!file.exists(path)) return(FALSE)
    info <- file.info(path)
    if (is.na(info$size) || info$size <= 0) return(FALSE)
    x <- try(utils::read.csv(path, nrows = 5, check.names = FALSE), silent = TRUE)
    if (inherits(x, "try-error") || !nrow(x)) return(FALSE)
    nms <- tolower(names(x))
    all(c("valid_time", "latitude", "longitude") %in% nms)
  }

  promote_csv <- function(source, dest) {
    if (!valid_csv(source)) return(FALSE)
    if (normalizePath(source, mustWork = FALSE) == normalizePath(dest, mustWork = FALSE)) {
      return(TRUE)
    }
    if (file.exists(dest)) unlink(dest)
    moved <- suppressWarnings(file.rename(source, dest))
    if (!moved) {
      moved <- file.copy(source, dest, overwrite = TRUE)
      if (moved) unlink(source)
    }
    isTRUE(moved) && valid_csv(dest)
  }

  if (valid_csv(dest)) return(dest)
  if (file.exists(dest)) unlink(dest)

  # Recover the historical `foo.csv.csv` / partial variants before resubmitting.
  for (candidate in c(paste0(dest, ".csv"),
                      paste0(dest, ".partial.csv"),
                      paste0(dest, ".partial"))) {
    if (promote_csv(candidate, dest)) {
      message(sprintf("  AgERA5 recovered legacy time-series cache: %s -> %s",
                      basename(candidate), basename(dest)))
      return(dest)
    }
  }

  req <- list(
    dataset_short_name = "sis-agrometeorological-indicators-timeseries",
    variable = vapply(.agera5_timeseries_vars, `[[`, character(1), "var"),
    date = unname(bounds),
    data_format = data_format,
    area = as.numeric(job$area),
    target = basename(dest)
  )

  err <- NULL
  returned <- tryCatch(
    ecmwfr::wf_request(request = req, path = job$cache_dir),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  if (is.character(returned) && length(returned)) {
    for (candidate in returned[nzchar(returned)]) {
      candidate_path <- if (file.exists(candidate)) {
        candidate
      } else {
        file.path(job$cache_dir, basename(candidate))
      }
      if (promote_csv(candidate_path, dest)) return(dest)
    }
  }

  for (candidate in c(paste0(dest, ".csv"), paste0(dest, ".partial.csv"))) {
    if (promote_csv(candidate, dest)) return(dest)
  }

  message(sprintf("  AgERA5 time-series download failed (%d, area=%s): %s",
                  job$year, paste(job$area, collapse = ","),
                  if (is.null(err)) "no valid CSV returned" else err))
  NULL
}
