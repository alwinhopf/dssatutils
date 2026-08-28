# Final AgERA5 download hardening.
#
# ecmwfr 2.x normalizes a target's extension to the extension of the returned
# CDS download URL.  Therefore targets such as `foo.zip.partial` and
# `foo.csv.partial` become `foo.zip.zip` / `foo.csv.csv`.  The canonical target
# must be requested directly.  These overrides are intentionally self-contained
# so `.agera5_download_job` also works when exported to a PSOCK worker: the
# worker does not need additional cache-helper functions exported separately.

.agera5_cds_job_url <- function(message) {
  hit <- regexpr(
    "https://cds\\.climate\\.copernicus\\.eu/api/retrieve/v1/jobs/[0-9A-Fa-f-]{36}",
    message,
    perl = TRUE
  )
  if (hit[[1]] < 0L) return(NULL)
  regmatches(message, hit)[[1]]
}

.agera5_cds_retry_delay <- function(message, attempt,
                                    base_seconds = 2,
                                    max_seconds = 60) {
  wait_hit <- regexec("Please wait[[:space:]]+([0-9]+)[[:space:]]+seconds?", message,
                      ignore.case = TRUE, perl = TRUE)
  wait_parts <- regmatches(message, wait_hit)[[1]]
  server_wait <- if (length(wait_parts) >= 2L) suppressWarnings(as.numeric(wait_parts[[2]])) else 0
  if (!is.finite(server_wait)) server_wait <- 0
  exponential <- min(as.numeric(max_seconds),
                     as.numeric(base_seconds) * 2^(max(1L, as.integer(attempt)) - 1L))
  max(1, server_wait, exponential)
}

.agera5_cds_transient_error <- function(message) {
  grepl(
    paste(c(
      "(^|[^0-9])429([^0-9]|$)", "rate[ -]?limit", "too many requests",
      "(^|[^0-9])50[234]([^0-9]|$)", "temporar", "timed?[ -]?out",
      "timeout", "connection reset", "connection error", "unavailable"
    ), collapse = "|"),
    message,
    ignore.case = TRUE,
    perl = TRUE
  )
}

.agera5_wf_request_with_retry <- function(request, path, target,
                                          request_fn = NULL,
                                          resume_fn = NULL,
                                          sleep_fn = Sys.sleep,
                                          max_attempts = NULL) {
  if (is.null(request_fn)) {
    request_fn <- function(request, path) {
      # A one-minute status interval materially reduces pressure on the CDS job
      # endpoint compared with ecmwfr's default 30-second polling interval.
      ecmwfr::wf_request(request = request, path = path, retry = 60, verbose = FALSE)
    }
  }
  if (is.null(resume_fn)) {
    resume_fn <- function(url, path, target) {
      ecmwfr::wf_transfer(url = url, path = path, filename = target, verbose = FALSE)
    }
  }

  if (is.null(max_attempts)) {
    max_attempts <- suppressWarnings(as.integer(Sys.getenv("AGERA5_CDS_MAX_ATTEMPTS", "8")))
  }
  if (is.na(max_attempts) || max_attempts < 1L) max_attempts <- 8L
  base_seconds <- suppressWarnings(as.numeric(Sys.getenv("AGERA5_CDS_RETRY_BASE_SECONDS", "2")))
  max_seconds <- suppressWarnings(as.numeric(Sys.getenv("AGERA5_CDS_RETRY_MAX_SECONDS", "60")))
  if (!is.finite(base_seconds) || base_seconds <= 0) base_seconds <- 2
  if (!is.finite(max_seconds) || max_seconds < base_seconds) max_seconds <- max(60, base_seconds)

  job_url <- NULL
  last_error <- NULL
  for (attempt in seq_len(max_attempts)) {
    err <- NULL
    returned <- tryCatch(
      if (is.null(job_url)) {
        request_fn(request, path)
      } else {
        resume_fn(job_url, path, target)
      },
      error = function(e) {
        err <<- conditionMessage(e)
        NULL
      }
    )
    if (is.null(err)) {
      return(list(value = returned, error = NULL, job_url = job_url,
                  attempts = attempt))
    }

    last_error <- err
    discovered_url <- .agera5_cds_job_url(err)
    if (!is.null(discovered_url)) job_url <- discovered_url
    retryable <- !is.null(job_url) || .agera5_cds_transient_error(err)
    if (!retryable || attempt >= max_attempts) break

    delay <- .agera5_cds_retry_delay(err, attempt, base_seconds, max_seconds)
    message(sprintf(
      "  AgERA5 CDS request throttled/transient failure; retrying %s in %.0f s (attempt %d/%d).",
      if (is.null(job_url)) "submission" else "existing job",
      delay, attempt + 1L, max_attempts
    ))
    sleep_fn(delay)
  }

  list(value = NULL, error = last_error, job_url = job_url,
       attempts = max_attempts)
}

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
  req <- .agera5_gridded_request(job, basename(job$zip_dest))

  transfer <- .agera5_wf_request_with_retry(
    request = req,
    path = job$cache_dir,
    target = basename(job$zip_dest)
  )
  returned <- transfer$value
  err <- transfer$error

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

  req <- .agera5_timeseries_request(job, basename(dest))

  transfer <- .agera5_wf_request_with_retry(
    request = req,
    path = job$cache_dir,
    target = basename(dest)
  )
  returned <- transfer$value
  err <- transfer$error

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
