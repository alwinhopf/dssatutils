# General DSSAT .WTH weather repair utilities.
#
# Short gaps and invalid temperature pairs are fixed after weather retrieval,
# not inside individual providers, so optional behavior applies to every weather
# source.

.weather_repair_default_vars <- c("SRAD", "TMAX", "TMIN", "RAIN", "TDEW", "RH2M", "WIND")

.weather_repair_find_header <- function(lines) {
  hits <- grep("^\\s*@\\s*DATE\\s+SRAD\\s+TMAX\\s+TMIN\\s+RAIN\\s+TDEW\\s+RH2M\\s+WIND\\s*$",
               lines, ignore.case = TRUE)
  if (length(hits)) hits[1] else NA_integer_
}

.weather_repair_date_label <- function(date_value) {
  d <- .weather_repair_date_from_code(date_value)
  if (!is.na(d)) return(as.character(d))
  as.character(date_value)
}

.weather_repair_date_from_code <- function(date_value) {
  txt <- as.character(date_value)
  year <- suppressWarnings(as.integer(substr(txt, 1, nchar(txt) - 3)))
  doy <- suppressWarnings(as.integer(substr(txt, nchar(txt) - 2, nchar(txt))))
  if (is.na(year) || is.na(doy)) return(as.Date(NA))
  if (year < 100L) year <- if (year >= 80L) 1900L + year else 2000L + year
  out <- as.Date(doy - 1L, origin = sprintf("%04d-01-01", year))
  if (is.na(out) || as.integer(format(out, "%Y")) != year) return(as.Date(NA))
  out
}

.weather_repair_code_from_date <- function(date_value) {
  if (!inherits(date_value, "Date")) date_value <- as.Date(date_value, origin = "1970-01-01")
  sprintf("%04d%03d",
          as.integer(format(date_value, "%Y")),
          as.integer(format(date_value, "%j")))
}

.weather_repair_runs <- function(mask) {
  idx <- which(mask)
  if (!length(idx)) return(data.frame(start = integer(), end = integer(), length = integer()))
  breaks <- c(0, which(diff(idx) > 1), length(idx))
  data.frame(
    start = idx[breaks[-length(breaks)] + 1],
    end = idx[breaks[-1]],
    stringsAsFactors = FALSE
  ) |>
    transform(length = end - start + 1L)
}

.weather_repair_log_lines <- function(path, lines) {
  if (is.null(path) || !nzchar(path) || !length(lines)) return(invisible())
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "a", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con = con, sep = "\n", useBytes = TRUE)
}

.weather_repair_parse_wth <- function(wth_file, log_file = NULL, issue = "WEATHER_QA") {
  lines <- readLines(wth_file, warn = FALSE)
  header_idx <- .weather_repair_find_header(lines)
  id <- tools::file_path_sans_ext(basename(wth_file))
  if (is.na(header_idx) || header_idx >= length(lines)) {
    .weather_repair_log_lines(log_file, sprintf(
      "%s file=%s id=%s issue=%s status=skipped reason=no_DSSAT_weather_header",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"), basename(wth_file), id, issue
    ))
    return(list(lines = lines, header_idx = header_idx, id = id,
                dat = data.frame(), status = "skipped_no_header"))
  }
  data_lines <- lines[(header_idx + 1L):length(lines)]
  data_lines <- data_lines[grepl("^\\s*\\d{5,7}\\s+", data_lines)]
  if (!length(data_lines)) {
    .weather_repair_log_lines(log_file, sprintf(
      "%s file=%s id=%s issue=%s status=skipped reason=no_weather_rows",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"), basename(wth_file), id, issue
    ))
    return(list(lines = lines, header_idx = header_idx, id = id,
                dat = data.frame(), status = "skipped_no_rows"))
  }
  dat <- utils::read.table(
    text = paste(data_lines, collapse = "\n"),
    col.names = c("DATE", .weather_repair_default_vars),
    colClasses = c("character", rep("numeric", length(.weather_repair_default_vars))),
    na.strings = c("-99", "-99.0", "-99.00", "-99.000"),
    stringsAsFactors = FALSE
  )
  dat$..DATE_OBJ <- as.Date(vapply(dat$DATE, function(x) {
    as.character(.weather_repair_date_from_code(x))
  }, character(1)))
  list(lines = lines, header_idx = header_idx, id = id, dat = dat, status = "ok")
}

.weather_repair_write_wth <- function(wth_file, lines, header_idx, dat) {
  write_dat <- dat
  write_dat$..DATE_OBJ <- NULL
  write_dat[.weather_repair_default_vars] <- lapply(write_dat[.weather_repair_default_vars], function(x) {
    x[is.na(x) | !is.finite(x)] <- -99
    x[x >= 9999.95 | x <= -999.95] <- -99
    x
  })
  formatted <- sprintf("%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f",
                       write_dat$DATE, write_dat$SRAD, write_dat$TMAX, write_dat$TMIN,
                       write_dat$RAIN, write_dat$TDEW, write_dat$RH2M, write_dat$WIND)
  formatted <- gsub("-99.0", "  -99", formatted, fixed = TRUE)
  writeLines(c(lines[seq_len(header_idx)], formatted), wth_file, useBytes = TRUE)
}

#' Repair short missing-value gaps in one DSSAT weather file
#'
#' Missing values encoded as `-99`/`-99.0` are repaired only when the contiguous
#' gap is no longer than `max_gap_days` and the full two-day window before and
#' after the gap is present. Each corrected value is filled with the mean of
#' those surrounding days for the same weather variable.
#'
#' @param wth_file Path to one `.WTH` file.
#' @param max_gap_days Maximum contiguous missing run length to repair.
#' @param window_days Number of non-missing days before and after the gap to use.
#' @param variables DSSAT weather variables to inspect.
#' @param log_file Optional log file receiving correction notes.
#' @param dry_run If TRUE, report repairs without writing the `.WTH`.
#' @return A data.frame summary with one row per inspected file/variable.
repair_weather_file_missing_values <- function(wth_file,
                                               max_gap_days = 3L,
                                               window_days = 2L,
                                               variables = .weather_repair_default_vars,
                                               log_file = NULL,
                                               dry_run = FALSE) {
  stopifnot(length(wth_file) == 1L)
  if (!file.exists(wth_file)) stop(sprintf("Weather file not found: %s", wth_file), call. = FALSE)
  max_gap_days <- as.integer(max_gap_days)
  window_days <- as.integer(window_days)
  if (is.na(max_gap_days) || max_gap_days < 1L) stop("max_gap_days must be >= 1", call. = FALSE)
  if (is.na(window_days) || window_days < 1L) stop("window_days must be >= 1", call. = FALSE)

  lines <- readLines(wth_file, warn = FALSE)
  header_idx <- .weather_repair_find_header(lines)
  id <- tools::file_path_sans_ext(basename(wth_file))
  if (is.na(header_idx) || header_idx >= length(lines)) {
    msg <- sprintf("%s file=%s id=%s status=skipped reason=no_DSSAT_weather_header",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S"), basename(wth_file), id)
    .weather_repair_log_lines(log_file, msg)
    return(data.frame(file = wth_file, id = id, variable = NA_character_,
                      repaired_count = 0L, unrepaired_count = 0L,
                      status = "skipped_no_header", stringsAsFactors = FALSE))
  }

  data_lines <- lines[(header_idx + 1L):length(lines)]
  data_lines <- data_lines[grepl("^\\s*\\d{5,7}\\s+", data_lines)]
  if (!length(data_lines)) {
    msg <- sprintf("%s file=%s id=%s status=skipped reason=no_weather_rows",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S"), basename(wth_file), id)
    .weather_repair_log_lines(log_file, msg)
    return(data.frame(file = wth_file, id = id, variable = NA_character_,
                      repaired_count = 0L, unrepaired_count = 0L,
                      status = "skipped_no_rows", stringsAsFactors = FALSE))
  }

  dat <- utils::read.table(
    text = paste(data_lines, collapse = "\n"),
    col.names = c("DATE", .weather_repair_default_vars),
    colClasses = c("character", rep("numeric", length(.weather_repair_default_vars))),
    na.strings = c("-99", "-99.0", "-99.00", "-99.000"),
    stringsAsFactors = FALSE
  )
  variables <- intersect(toupper(as.character(variables)), .weather_repair_default_vars)
  if (!length(variables)) variables <- .weather_repair_default_vars

  log_lines <- character()
  summary_rows <- list()
  original <- dat

  for (var in variables) {
    before <- original[[var]]
    missing_before <- is.na(before) | !is.finite(before)
    runs <- .weather_repair_runs(missing_before)
    repaired_count <- 0L
    unrepaired_count <- 0L
    repaired_runs <- 0L

    if (nrow(runs)) {
      for (rr in seq_len(nrow(runs))) {
        s <- runs$start[rr]; e <- runs$end[rr]; len <- runs$length[rr]
        if (len <= max_gap_days) {
          neighbor_idx <- c((s - window_days):(s - 1L), (e + 1L):(e + window_days))
          in_bounds <- all(neighbor_idx >= 1L & neighbor_idx <= nrow(original))
          neighbor_vals <- if (in_bounds) before[neighbor_idx] else numeric()
          usable <- in_bounds && all(is.finite(neighbor_vals)) && !any(is.na(neighbor_vals))
          if (usable) {
            fill_value <- mean(neighbor_vals)
            dat[[var]][s:e] <- fill_value
            repaired_count <- repaired_count + len
            repaired_runs <- repaired_runs + 1L
            log_lines <- c(log_lines, sprintf(
              "%s file=%s id=%s variable=%s status=repaired dates=%s..%s gap_days=%d fill_value=%.4f method=mean_%d_days_before_after neighbor_dates=%s..%s;%s..%s",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), basename(wth_file), id, var,
              .weather_repair_date_label(dat$DATE[s]), .weather_repair_date_label(dat$DATE[e]),
              len, fill_value, window_days,
              .weather_repair_date_label(dat$DATE[neighbor_idx[1]]),
              .weather_repair_date_label(dat$DATE[neighbor_idx[window_days]]),
              .weather_repair_date_label(dat$DATE[neighbor_idx[window_days + 1L]]),
              .weather_repair_date_label(dat$DATE[neighbor_idx[length(neighbor_idx)]])
            ))
          } else {
            unrepaired_count <- unrepaired_count + len
            log_lines <- c(log_lines, sprintf(
              "%s file=%s id=%s variable=%s status=unrepaired dates=%s..%s gap_days=%d reason=insufficient_%d_day_neighbors",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), basename(wth_file), id, var,
              .weather_repair_date_label(dat$DATE[s]), .weather_repair_date_label(dat$DATE[e]),
              len, window_days
            ))
          }
        } else {
          unrepaired_count <- unrepaired_count + len
          log_lines <- c(log_lines, sprintf(
            "%s file=%s id=%s variable=%s status=unrepaired dates=%s..%s gap_days=%d reason=gap_exceeds_max_%d_days",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), basename(wth_file), id, var,
            .weather_repair_date_label(dat$DATE[s]), .weather_repair_date_label(dat$DATE[e]),
            len, max_gap_days
          ))
        }
      }
    }

    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      file = wth_file, id = id, variable = var,
      repaired_count = repaired_count, unrepaired_count = unrepaired_count,
      repaired_runs = repaired_runs,
      status = if (repaired_count > 0L) "repaired" else if (unrepaired_count > 0L) "unrepaired_missing" else "unchanged",
      stringsAsFactors = FALSE
    )
  }

  .weather_repair_log_lines(log_file, log_lines)

  if (!dry_run && any(vapply(summary_rows, function(x) x$repaired_count > 0L, logical(1)))) {
    write_dat <- dat
    write_dat[.weather_repair_default_vars] <- lapply(write_dat[.weather_repair_default_vars], function(x) {
      x[is.na(x) | !is.finite(x)] <- -99
      x[x >= 9999.95 | x <= -999.95] <- -99
      x
    })
    formatted <- sprintf("%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f",
                         write_dat$DATE, write_dat$SRAD, write_dat$TMAX, write_dat$TMIN,
                         write_dat$RAIN, write_dat$TDEW, write_dat$RH2M, write_dat$WIND)
    formatted <- gsub("-99.0", "  -99", formatted, fixed = TRUE)
    writeLines(c(lines[seq_len(header_idx)], formatted), wth_file, useBytes = TRUE)
  }

  do.call(rbind, summary_rows)
}

#' Repair short missing-value gaps in DSSAT weather files
#'
#' Applies `repair_weather_file_missing_values()` to all `.WTH` files in a
#' weather directory, optionally restricted to point IDs.
#'
#' @param weather_dir Directory containing DSSAT `.WTH` files.
#' @param ids Optional point IDs to restrict processing.
#' @inheritParams repair_weather_file_missing_values
#' @return A data.frame summary with repair counts per file/variable.
repair_weather_missing_values <- function(weather_dir,
                                          ids = NULL,
                                          max_gap_days = 3L,
                                          window_days = 2L,
                                          variables = .weather_repair_default_vars,
                                          log_file = file.path(weather_dir, "weather_repair.log"),
                                          dry_run = FALSE) {
  if (!dir.exists(weather_dir)) stop(sprintf("Weather directory not found: %s", weather_dir), call. = FALSE)
  files <- list.files(weather_dir, pattern = "\\.WTH$", full.names = TRUE)
  if (!is.null(ids) && length(ids)) {
    wanted <- sprintf("%s.WTH", as.character(ids))
    files <- files[basename(files) %in% wanted]
  }
  if (!length(files)) {
    msg <- sprintf("%s weather_dir=%s status=skipped reason=no_weather_files",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S"), weather_dir)
    .weather_repair_log_lines(log_file, msg)
    return(data.frame(file = character(), id = character(), variable = character(),
                      repaired_count = integer(), unrepaired_count = integer(),
                      repaired_runs = integer(), status = character(),
                      stringsAsFactors = FALSE))
  }

  header <- sprintf(
    "%s weather_dir=%s status=started files=%d max_gap_days=%d window_days=%d variables=%s dry_run=%s",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), weather_dir, length(files),
    as.integer(max_gap_days), as.integer(window_days), paste(variables, collapse = ","), dry_run)
  .weather_repair_log_lines(log_file, c("", header))

  out <- lapply(files, repair_weather_file_missing_values,
                max_gap_days = max_gap_days,
                window_days = window_days,
                variables = variables,
                log_file = log_file,
                dry_run = dry_run)
  summary <- do.call(rbind, out)
  footer <- sprintf(
    "%s weather_dir=%s status=finished repaired_values=%d unrepaired_values=%d log_file=%s",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), weather_dir,
    sum(summary$repaired_count, na.rm = TRUE),
    sum(summary$unrepaired_count, na.rm = TRUE),
    log_file)
  .weather_repair_log_lines(log_file, footer)
  summary
}

#' Repair short Tmax/Tmin inversion runs in one DSSAT weather file
#'
#' Rows where both temperatures are present but `TMIN > TMAX` are repaired only
#' when the contiguous inversion run is no longer than `max_gap_days` and the
#' full `window_days` before and after the run has valid, non-inverted
#' temperatures. `TMAX` and `TMIN` are replaced independently using the mean of
#' their respective neighboring values.
#'
#' @inheritParams repair_weather_file_missing_values
#' @return A data.frame summary with one row for the inspected file.
repair_weather_file_temperature_inversions <- function(wth_file,
                                                       max_gap_days = 3L,
                                                       window_days = 2L,
                                                       log_file = NULL,
                                                       dry_run = FALSE) {
  stopifnot(length(wth_file) == 1L)
  if (!file.exists(wth_file)) stop(sprintf("Weather file not found: %s", wth_file), call. = FALSE)
  max_gap_days <- as.integer(max_gap_days)
  window_days <- as.integer(window_days)
  if (is.na(max_gap_days) || max_gap_days < 1L) stop("max_gap_days must be >= 1", call. = FALSE)
  if (is.na(window_days) || window_days < 1L) stop("window_days must be >= 1", call. = FALSE)

  lines <- readLines(wth_file, warn = FALSE)
  header_idx <- .weather_repair_find_header(lines)
  id <- tools::file_path_sans_ext(basename(wth_file))
  if (is.na(header_idx) || header_idx >= length(lines)) {
    msg <- sprintf("%s file=%s id=%s issue=TMIN_GT_TMAX status=skipped reason=no_DSSAT_weather_header",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S"), basename(wth_file), id)
    .weather_repair_log_lines(log_file, msg)
    return(data.frame(file = wth_file, id = id, issue = "TMIN_GT_TMAX",
                      repaired_count = 0L, unrepaired_count = 0L,
                      repaired_runs = 0L, status = "skipped_no_header",
                      stringsAsFactors = FALSE))
  }

  data_lines <- lines[(header_idx + 1L):length(lines)]
  data_lines <- data_lines[grepl("^\\s*\\d{5,7}\\s+", data_lines)]
  if (!length(data_lines)) {
    msg <- sprintf("%s file=%s id=%s issue=TMIN_GT_TMAX status=skipped reason=no_weather_rows",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S"), basename(wth_file), id)
    .weather_repair_log_lines(log_file, msg)
    return(data.frame(file = wth_file, id = id, issue = "TMIN_GT_TMAX",
                      repaired_count = 0L, unrepaired_count = 0L,
                      repaired_runs = 0L, status = "skipped_no_rows",
                      stringsAsFactors = FALSE))
  }

  dat <- utils::read.table(
    text = paste(data_lines, collapse = "\n"),
    col.names = c("DATE", .weather_repair_default_vars),
    colClasses = c("character", rep("numeric", length(.weather_repair_default_vars))),
    na.strings = c("-99", "-99.0", "-99.00", "-99.000"),
    stringsAsFactors = FALSE
  )
  original <- dat
  inversion <- is.finite(original$TMAX) & is.finite(original$TMIN) &
    !is.na(original$TMAX) & !is.na(original$TMIN) &
    original$TMIN > original$TMAX
  runs <- .weather_repair_runs(inversion)

  log_lines <- character()
  repaired_count <- 0L
  unrepaired_count <- 0L
  repaired_runs <- 0L

  if (nrow(runs)) {
    for (rr in seq_len(nrow(runs))) {
      s <- runs$start[rr]; e <- runs$end[rr]; len <- runs$length[rr]
      if (len <= max_gap_days) {
        neighbor_idx <- c((s - window_days):(s - 1L), (e + 1L):(e + window_days))
        in_bounds <- all(neighbor_idx >= 1L & neighbor_idx <= nrow(original))
        neighbor_tmax <- if (in_bounds) original$TMAX[neighbor_idx] else numeric()
        neighbor_tmin <- if (in_bounds) original$TMIN[neighbor_idx] else numeric()
        neighbor_ok <- in_bounds &&
          all(is.finite(neighbor_tmax)) && all(is.finite(neighbor_tmin)) &&
          !any(is.na(neighbor_tmax)) && !any(is.na(neighbor_tmin)) &&
          all(neighbor_tmin <= neighbor_tmax)
        if (neighbor_ok) {
          fill_tmax <- mean(neighbor_tmax)
          fill_tmin <- mean(neighbor_tmin)
          dat$TMAX[s:e] <- fill_tmax
          dat$TMIN[s:e] <- fill_tmin
          repaired_count <- repaired_count + len
          repaired_runs <- repaired_runs + 1L
          log_lines <- c(log_lines, sprintf(
            "%s file=%s id=%s issue=TMIN_GT_TMAX status=repaired dates=%s..%s gap_days=%d fill_TMAX=%.4f fill_TMIN=%.4f method=mean_%d_days_before_after neighbor_dates=%s..%s;%s..%s",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), basename(wth_file), id,
            .weather_repair_date_label(dat$DATE[s]), .weather_repair_date_label(dat$DATE[e]),
            len, fill_tmax, fill_tmin, window_days,
            .weather_repair_date_label(dat$DATE[neighbor_idx[1]]),
            .weather_repair_date_label(dat$DATE[neighbor_idx[window_days]]),
            .weather_repair_date_label(dat$DATE[neighbor_idx[window_days + 1L]]),
            .weather_repair_date_label(dat$DATE[neighbor_idx[length(neighbor_idx)]])
          ))
        } else {
          unrepaired_count <- unrepaired_count + len
          log_lines <- c(log_lines, sprintf(
            "%s file=%s id=%s issue=TMIN_GT_TMAX status=unrepaired dates=%s..%s gap_days=%d reason=insufficient_%d_day_valid_temperature_neighbors",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), basename(wth_file), id,
            .weather_repair_date_label(dat$DATE[s]), .weather_repair_date_label(dat$DATE[e]),
            len, window_days
          ))
        }
      } else {
        unrepaired_count <- unrepaired_count + len
        log_lines <- c(log_lines, sprintf(
          "%s file=%s id=%s issue=TMIN_GT_TMAX status=unrepaired dates=%s..%s gap_days=%d reason=gap_exceeds_max_%d_days",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S"), basename(wth_file), id,
          .weather_repair_date_label(dat$DATE[s]), .weather_repair_date_label(dat$DATE[e]),
          len, max_gap_days
        ))
      }
    }
  }

  .weather_repair_log_lines(log_file, log_lines)

  if (!dry_run && repaired_count > 0L) {
    write_dat <- dat
    write_dat[.weather_repair_default_vars] <- lapply(write_dat[.weather_repair_default_vars], function(x) {
      x[is.na(x) | !is.finite(x)] <- -99
      x[x >= 9999.95 | x <= -999.95] <- -99
      x
    })
    formatted <- sprintf("%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f",
                         write_dat$DATE, write_dat$SRAD, write_dat$TMAX, write_dat$TMIN,
                         write_dat$RAIN, write_dat$TDEW, write_dat$RH2M, write_dat$WIND)
    formatted <- gsub("-99.0", "  -99", formatted, fixed = TRUE)
    writeLines(c(lines[seq_len(header_idx)], formatted), wth_file, useBytes = TRUE)
  }

  data.frame(
    file = wth_file, id = id, issue = "TMIN_GT_TMAX",
    repaired_count = repaired_count, unrepaired_count = unrepaired_count,
    repaired_runs = repaired_runs,
    status = if (repaired_count > 0L) "repaired" else if (unrepaired_count > 0L) "unrepaired_temperature_inversion" else "unchanged",
    stringsAsFactors = FALSE
  )
}

#' Repair short Tmax/Tmin inversion runs in DSSAT weather files
#'
#' Applies `repair_weather_file_temperature_inversions()` to all `.WTH` files in
#' a weather directory, optionally restricted to point IDs.
#'
#' @param weather_dir Directory containing DSSAT `.WTH` files.
#' @param ids Optional point IDs to restrict processing.
#' @inheritParams repair_weather_file_temperature_inversions
#' @return A data.frame summary with repair counts per file.
repair_weather_temperature_inversions <- function(weather_dir,
                                                  ids = NULL,
                                                  max_gap_days = 3L,
                                                  window_days = 2L,
                                                  log_file = file.path(weather_dir, "weather_repair.log"),
                                                  dry_run = FALSE) {
  if (!dir.exists(weather_dir)) stop(sprintf("Weather directory not found: %s", weather_dir), call. = FALSE)
  files <- list.files(weather_dir, pattern = "\\.WTH$", full.names = TRUE)
  if (!is.null(ids) && length(ids)) {
    wanted <- sprintf("%s.WTH", as.character(ids))
    files <- files[basename(files) %in% wanted]
  }
  if (!length(files)) {
    msg <- sprintf("%s weather_dir=%s issue=TMIN_GT_TMAX status=skipped reason=no_weather_files",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S"), weather_dir)
    .weather_repair_log_lines(log_file, msg)
    return(data.frame(file = character(), id = character(), issue = character(),
                      repaired_count = integer(), unrepaired_count = integer(),
                      repaired_runs = integer(), status = character(),
                      stringsAsFactors = FALSE))
  }

  header <- sprintf(
    "%s weather_dir=%s issue=TMIN_GT_TMAX status=started files=%d max_gap_days=%d window_days=%d dry_run=%s",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), weather_dir, length(files),
    as.integer(max_gap_days), as.integer(window_days), dry_run)
  .weather_repair_log_lines(log_file, c("", header))

  out <- lapply(files, repair_weather_file_temperature_inversions,
                max_gap_days = max_gap_days,
                window_days = window_days,
                log_file = log_file,
                dry_run = dry_run)
  summary <- do.call(rbind, out)
  footer <- sprintf(
    "%s weather_dir=%s issue=TMIN_GT_TMAX status=finished repaired_values=%d unrepaired_values=%d log_file=%s",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), weather_dir,
    sum(summary$repaired_count, na.rm = TRUE),
    sum(summary$unrepaired_count, na.rm = TRUE),
    log_file)
  .weather_repair_log_lines(log_file, footer)
  summary
}

#' Repair short missing DATE-row gaps in one DSSAT weather file
#'
#' Missing calendar rows are repaired only when the contiguous missing-date run
#' is no longer than `max_gap_days` and the full `window_days` before and after
#' the run has finite values. Inserted rows receive the mean of neighboring
#' values for each requested weather variable.
#'
#' @inheritParams repair_weather_file_missing_values
#' @return A data.frame summary with one row for the inspected file.
repair_weather_file_date_gaps <- function(wth_file,
                                          max_gap_days = 3L,
                                          window_days = 2L,
                                          variables = .weather_repair_default_vars,
                                          log_file = NULL,
                                          dry_run = FALSE) {
  stopifnot(length(wth_file) == 1L)
  if (!file.exists(wth_file)) stop(sprintf("Weather file not found: %s", wth_file), call. = FALSE)
  max_gap_days <- as.integer(max_gap_days)
  window_days <- as.integer(window_days)
  if (is.na(max_gap_days) || max_gap_days < 1L) stop("max_gap_days must be >= 1", call. = FALSE)
  if (is.na(window_days) || window_days < 1L) stop("window_days must be >= 1", call. = FALSE)

  parsed <- .weather_repair_parse_wth(wth_file, log_file, issue = "DATE_GAP")
  id <- parsed$id
  if (parsed$status != "ok") {
    return(data.frame(file = wth_file, id = id, issue = "DATE_GAP",
                      repaired_count = 0L, unrepaired_count = 0L,
                      repaired_runs = 0L, status = parsed$status,
                      stringsAsFactors = FALSE))
  }
  dat <- parsed$dat
  if (any(is.na(dat$..DATE_OBJ))) {
    bad <- sum(is.na(dat$..DATE_OBJ))
    .weather_repair_log_lines(log_file, sprintf(
      "%s file=%s id=%s issue=DATE_GAP status=unrepaired reason=unparseable_date_codes count=%d",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"), basename(wth_file), id, bad
    ))
    return(data.frame(file = wth_file, id = id, issue = "DATE_GAP",
                      repaired_count = 0L, unrepaired_count = bad,
                      repaired_runs = 0L, status = "unrepaired_unparseable_dates",
                      stringsAsFactors = FALSE))
  }
  dup <- duplicated(dat$..DATE_OBJ) | duplicated(dat$..DATE_OBJ, fromLast = TRUE)
  if (any(dup)) {
    bad <- sum(dup)
    .weather_repair_log_lines(log_file, sprintf(
      "%s file=%s id=%s issue=DATE_GAP status=unrepaired reason=duplicate_dates count=%d",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"), basename(wth_file), id, bad
    ))
    return(data.frame(file = wth_file, id = id, issue = "DATE_GAP",
                      repaired_count = 0L, unrepaired_count = bad,
                      repaired_runs = 0L, status = "unrepaired_duplicate_dates",
                      stringsAsFactors = FALSE))
  }

  variables <- intersect(toupper(as.character(variables)), .weather_repair_default_vars)
  if (!length(variables)) variables <- .weather_repair_default_vars
  original_dates <- dat$..DATE_OBJ
  dat <- dat[order(dat$..DATE_OBJ), , drop = FALSE]
  row.names(dat) <- NULL
  sorted_rows <- !identical(as.numeric(original_dates), as.numeric(dat$..DATE_OBJ))

  expected <- seq(min(dat$..DATE_OBJ), max(dat$..DATE_OBJ), by = "day")
  missing <- expected[!(expected %in% dat$..DATE_OBJ)]
  runs <- .weather_repair_runs(rep(FALSE, length(expected)))
  if (length(missing)) {
    miss_mask <- expected %in% missing
    runs <- .weather_repair_runs(miss_mask)
    runs$start_date <- expected[runs$start]
    runs$end_date <- expected[runs$end]
  }

  log_lines <- character()
  repaired_count <- 0L
  unrepaired_count <- 0L
  repaired_runs <- 0L
  new_rows <- list()

  if (nrow(runs)) {
    date_key <- as.character(dat$..DATE_OBJ)
    for (rr in seq_len(nrow(runs))) {
      s_date <- runs$start_date[rr]; e_date <- runs$end_date[rr]; len <- runs$length[rr]
      if (len <= max_gap_days) {
        neighbor_dates <- c(seq(s_date - window_days, s_date - 1, by = "day"),
                            seq(e_date + 1, e_date + window_days, by = "day"))
        idx <- match(as.character(neighbor_dates), date_key)
        usable <- !any(is.na(idx))
        if (usable) {
          usable <- all(vapply(variables, function(v) {
            vals <- dat[[v]][idx]
            all(is.finite(vals)) && !any(is.na(vals))
          }, logical(1)))
        }
        if (usable) {
          fill_values <- vapply(variables, function(v) mean(dat[[v]][idx]), numeric(1))
          for (d in seq(s_date, e_date, by = "day")) {
            d <- as.Date(d, origin = "1970-01-01")
            row <- as.list(stats::setNames(rep(NA_real_, length(.weather_repair_default_vars)), .weather_repair_default_vars))
            row$DATE <- .weather_repair_code_from_date(d)
            row$..DATE_OBJ <- d
            for (v in names(fill_values)) row[[v]] <- fill_values[[v]]
            new_rows[[length(new_rows) + 1L]] <- as.data.frame(row, stringsAsFactors = FALSE)
          }
          repaired_count <- repaired_count + len
          repaired_runs <- repaired_runs + 1L
          fill_txt <- paste(sprintf("%s=%.4f", names(fill_values), fill_values), collapse = ",")
          log_lines <- c(log_lines, sprintf(
            "%s file=%s id=%s issue=DATE_GAP status=repaired dates=%s..%s gap_days=%d fill_values=%s method=mean_%d_days_before_after neighbor_dates=%s..%s;%s..%s",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), basename(wth_file), id,
            s_date, e_date, len, fill_txt, window_days,
            neighbor_dates[1], neighbor_dates[window_days],
            neighbor_dates[window_days + 1L], neighbor_dates[length(neighbor_dates)]
          ))
        } else {
          unrepaired_count <- unrepaired_count + len
          log_lines <- c(log_lines, sprintf(
            "%s file=%s id=%s issue=DATE_GAP status=unrepaired dates=%s..%s gap_days=%d reason=insufficient_%d_day_neighbors",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), basename(wth_file), id,
            s_date, e_date, len, window_days
          ))
        }
      } else {
        unrepaired_count <- unrepaired_count + len
        log_lines <- c(log_lines, sprintf(
          "%s file=%s id=%s issue=DATE_GAP status=unrepaired dates=%s..%s gap_days=%d reason=gap_exceeds_max_%d_days",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S"), basename(wth_file), id,
          s_date, e_date, len, max_gap_days
        ))
      }
    }
  }

  .weather_repair_log_lines(log_file, log_lines)
  if (!dry_run && (length(new_rows) || sorted_rows)) {
    if (length(new_rows)) dat <- rbind(dat, do.call(rbind, new_rows))
    dat <- dat[order(dat$..DATE_OBJ), , drop = FALSE]
    row.names(dat) <- NULL
    .weather_repair_write_wth(wth_file, parsed$lines, parsed$header_idx, dat)
  }

  data.frame(
    file = wth_file, id = id, issue = "DATE_GAP",
    repaired_count = repaired_count, unrepaired_count = unrepaired_count,
    repaired_runs = repaired_runs,
    status = if (repaired_count > 0L) "repaired" else if (sorted_rows) "sorted" else if (unrepaired_count > 0L) "unrepaired_date_gap" else "unchanged",
    stringsAsFactors = FALSE
  )
}

#' Repair short missing DATE-row gaps in DSSAT weather files
#'
#' Applies `repair_weather_file_date_gaps()` to all `.WTH` files in a weather
#' directory, optionally restricted to point IDs.
#'
#' @param weather_dir Directory containing DSSAT `.WTH` files.
#' @param ids Optional point IDs to restrict processing.
#' @inheritParams repair_weather_file_date_gaps
#' @return A data.frame summary with repair counts per file.
repair_weather_date_gaps <- function(weather_dir,
                                     ids = NULL,
                                     max_gap_days = 3L,
                                     window_days = 2L,
                                     variables = .weather_repair_default_vars,
                                     log_file = file.path(weather_dir, "weather_repair.log"),
                                     dry_run = FALSE) {
  if (!dir.exists(weather_dir)) stop(sprintf("Weather directory not found: %s", weather_dir), call. = FALSE)
  files <- list.files(weather_dir, pattern = "\\.WTH$", full.names = TRUE)
  if (!is.null(ids) && length(ids)) {
    wanted <- sprintf("%s.WTH", as.character(ids))
    files <- files[basename(files) %in% wanted]
  }
  if (!length(files)) {
    msg <- sprintf("%s weather_dir=%s issue=DATE_GAP status=skipped reason=no_weather_files",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S"), weather_dir)
    .weather_repair_log_lines(log_file, msg)
    return(data.frame(file = character(), id = character(), issue = character(),
                      repaired_count = integer(), unrepaired_count = integer(),
                      repaired_runs = integer(), status = character(),
                      stringsAsFactors = FALSE))
  }
  .weather_repair_log_lines(log_file, c("", sprintf(
    "%s weather_dir=%s issue=DATE_GAP status=started files=%d max_gap_days=%d window_days=%d variables=%s dry_run=%s",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), weather_dir, length(files),
    as.integer(max_gap_days), as.integer(window_days), paste(variables, collapse = ","), dry_run
  )))
  out <- lapply(files, repair_weather_file_date_gaps,
                max_gap_days = max_gap_days,
                window_days = window_days,
                variables = variables,
                log_file = log_file,
                dry_run = dry_run)
  summary <- do.call(rbind, out)
  .weather_repair_log_lines(log_file, sprintf(
    "%s weather_dir=%s issue=DATE_GAP status=finished repaired_values=%d unrepaired_values=%d log_file=%s",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), weather_dir,
    sum(summary$repaired_count, na.rm = TRUE),
    sum(summary$unrepaired_count, na.rm = TRUE),
    log_file
  ))
  summary
}

.weather_quality_issue <- function(file, id, issue, severity, count,
                                   first_date = "", last_date = "", details = "") {
  if (is.na(count) || count <= 0L) return(NULL)
  data.frame(file = file, id = id, issue = issue, severity = severity,
             count = as.integer(count), first_date = as.character(first_date),
             last_date = as.character(last_date), details = details,
             stringsAsFactors = FALSE)
}

.weather_flatline_issues <- function(dat, var, min_days, file, id) {
  x <- dat[[var]]
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    if (is.na(x[i]) || !is.finite(x[i])) {
      i <- i + 1L
      next
    }
    j <- i
    while (j + 1L <= length(x) && is.finite(x[j + 1L]) && !is.na(x[j + 1L]) && x[j + 1L] == x[i]) {
      j <- j + 1L
    }
    len <- j - i + 1L
    if (len >= min_days) {
      out[[length(out) + 1L]] <- data.frame(
        start = i, end = j, length = len, value = x[i],
        first_date = dat$..DATE_OBJ[i], last_date = dat$..DATE_OBJ[j]
      )
    }
    i <- j + 1L
  }
  if (!length(out)) return(NULL)
  runs <- do.call(rbind, out)
  .weather_quality_issue(
    file, id, paste0(var, "_flatline"), "info", sum(runs$length),
    runs$first_date[1], runs$last_date[nrow(runs)],
    sprintf("min_days=%d; examples=%s", min_days,
            paste(sprintf("%.2fx%dd", runs$value, runs$length), collapse = ","))
  )
}

#' Audit one DSSAT weather file for suspicious values
#'
#' This is flag-only QA: it does not modify the `.WTH` file.
#'
#' @param wth_file Path to one `.WTH` file.
#' @param flatline_days Minimum identical-value run length to flag.
#' @param log_file Optional log file receiving audit notes.
#' @return A data.frame of findings.
audit_weather_file_quality <- function(wth_file,
                                       flatline_days = 10L,
                                       log_file = NULL) {
  parsed <- .weather_repair_parse_wth(wth_file, log_file, issue = "WEATHER_QA")
  id <- parsed$id
  file <- wth_file
  findings <- list()
  add <- function(x) if (!is.null(x)) findings[[length(findings) + 1L]] <<- x
  if (parsed$status != "ok") {
    add(.weather_quality_issue(file, id, parsed$status, "error", 1L, details = parsed$status))
    return(do.call(rbind, findings))
  }
  dat <- parsed$dat
  bad_dates <- is.na(dat$..DATE_OBJ)
  add(.weather_quality_issue(file, id, "unparseable_date_codes", "error", sum(bad_dates)))
  if (any(bad_dates)) return(do.call(rbind, findings))

  dup <- duplicated(dat$..DATE_OBJ) | duplicated(dat$..DATE_OBJ, fromLast = TRUE)
  add(.weather_quality_issue(file, id, "duplicate_dates", "error", sum(dup)))
  add(.weather_quality_issue(file, id, "out_of_order_dates", "warning",
                             as.integer(any(diff(dat$..DATE_OBJ) < 0))))
  expected <- seq(min(dat$..DATE_OBJ), max(dat$..DATE_OBJ), by = "day")
  missing <- expected[!(expected %in% dat$..DATE_OBJ)]
  add(.weather_quality_issue(file, id, "missing_date_rows", "error", length(missing),
                             if (length(missing)) missing[1] else "",
                             if (length(missing)) missing[length(missing)] else ""))

  for (var in .weather_repair_default_vars) {
    idx <- which(is.na(dat[[var]]) | !is.finite(dat[[var]]))
    add(.weather_quality_issue(file, id, paste0(var, "_missing_values"), "warning", length(idx),
                               if (length(idx)) dat$..DATE_OBJ[idx[1]] else "",
                               if (length(idx)) dat$..DATE_OBJ[idx[length(idx)]] else ""))
  }

  checks <- list(
    list("tmin_gt_tmax", "error", is.finite(dat$TMAX) & is.finite(dat$TMIN) & dat$TMIN > dat$TMAX, ""),
    list("tmax_out_of_range", "warning", is.finite(dat$TMAX) & (dat$TMAX < -60 | dat$TMAX > 60), "bounds=-60..60C"),
    list("tmin_out_of_range", "warning", is.finite(dat$TMIN) & (dat$TMIN < -70 | dat$TMIN > 50), "bounds=-70..50C"),
    list("diurnal_range_extreme", "warning", is.finite(dat$TMAX) & is.finite(dat$TMIN) & (dat$TMAX - dat$TMIN > 45), "TMAX-TMIN>45C"),
    list("rain_negative", "error", is.finite(dat$RAIN) & dat$RAIN < 0, ""),
    list("rain_extreme", "warning", is.finite(dat$RAIN) & dat$RAIN > 500, "RAIN>500mm"),
    list("srad_out_of_range", "warning", is.finite(dat$SRAD) & (dat$SRAD < 0 | dat$SRAD > 40), "bounds=0..40MJ/m2/day"),
    list("rh2m_out_of_range", "warning", is.finite(dat$RH2M) & (dat$RH2M < 0 | dat$RH2M > 100), "bounds=0..100%"),
    list("wind_out_of_range", "warning", is.finite(dat$WIND) & (dat$WIND < 0 | dat$WIND > 75), "bounds=0..75m/s"),
    list("tdew_gt_tmax", "warning", is.finite(dat$TDEW) & is.finite(dat$TMAX) & dat$TDEW > dat$TMAX, "")
  )
  for (chk in checks) {
    idx <- which(chk[[3]])
    add(.weather_quality_issue(file, id, chk[[1]], chk[[2]], length(idx),
                               if (length(idx)) dat$..DATE_OBJ[idx[1]] else "",
                               if (length(idx)) dat$..DATE_OBJ[idx[length(idx)]] else "",
                               chk[[4]]))
  }
  for (var in .weather_repair_default_vars) {
    add(.weather_flatline_issues(dat, var, as.integer(flatline_days), file, id))
  }

  if (!length(findings)) {
    return(data.frame(file = character(), id = character(), issue = character(),
                      severity = character(), count = integer(),
                      first_date = character(), last_date = character(),
                      details = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, findings)
}

#' Audit DSSAT weather files for suspicious values
#'
#' Writes a flag-only CSV summary and appends a brief note to the repair log.
#'
#' @param weather_dir Directory containing DSSAT `.WTH` files.
#' @param ids Optional point IDs to restrict processing.
#' @param audit_csv Output CSV path for findings.
#' @inheritParams audit_weather_file_quality
#' @return A data.frame of findings.
audit_weather_quality <- function(weather_dir,
                                  ids = NULL,
                                  audit_csv = file.path(weather_dir, "weather_quality_audit.csv"),
                                  flatline_days = 10L,
                                  log_file = file.path(weather_dir, "weather_repair.log")) {
  if (!dir.exists(weather_dir)) stop(sprintf("Weather directory not found: %s", weather_dir), call. = FALSE)
  files <- list.files(weather_dir, pattern = "\\.WTH$", full.names = TRUE)
  if (!is.null(ids) && length(ids)) {
    wanted <- sprintf("%s.WTH", as.character(ids))
    files <- files[basename(files) %in% wanted]
  }
  out <- lapply(files, audit_weather_file_quality,
                flatline_days = flatline_days,
                log_file = log_file)
  summary <- if (length(out)) do.call(rbind, out) else data.frame()
  if (is.null(summary) || !nrow(summary)) {
    summary <- data.frame(file = character(), id = character(), issue = character(),
                          severity = character(), count = integer(),
                          first_date = character(), last_date = character(),
                          details = character(), stringsAsFactors = FALSE)
  }
  dir.create(dirname(audit_csv), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(summary, audit_csv, row.names = FALSE)
  .weather_repair_log_lines(log_file, sprintf(
    "%s weather_dir=%s issue=WEATHER_QA status=finished files=%d findings=%d audit_csv=%s",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), weather_dir, length(files), nrow(summary), audit_csv
  ))
  summary
}
