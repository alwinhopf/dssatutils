# Validate a DSSAT weather file, including fixed-width rows whose adjacent
# negative values are not separated by whitespace.

.parse_wth_data_line <- function(line) {
  fixed <- NULL
  if (nchar(line, type = "chars") >= 49L) {
    fixed <- c(trimws(substr(line, 1L, 7L)), vapply(seq.int(8L, 44L, by = 6L),
      function(start) trimws(substr(line, start, start + 5L)), character(1)))
  }
  fields <- if (!is.null(fixed) && grepl("^[0-9]{5,7}$", fixed[1]) &&
                all(nzchar(fixed[-1]))) fixed else strsplit(trimws(line), "\\s+")[[1]]
  if (length(fields) != 8L || !grepl("^[0-9]{5,7}$", fields[1])) return(NULL)
  values <- suppressWarnings(as.numeric(fields[-1]))
  if (any(!is.finite(values))) return(NULL)
  list(code = fields[1], values = values)
}

.wth_code_to_date <- function(code) {
  code <- as.character(code)
  if (nchar(code) == 5L) {
    yy <- as.integer(substr(code, 1L, 2L))
    year <- if (yy < 80L) 2000L + yy else 1900L + yy
    doy <- as.integer(substr(code, 3L, 5L))
  } else if (nchar(code) == 7L) {
    year <- as.integer(substr(code, 1L, 4L))
    doy <- as.integer(substr(code, 5L, 7L))
  } else return(as.Date(NA))
  if (!is.finite(year) || !is.finite(doy) || doy < 1L || doy > 366L) return(as.Date(NA))
  date <- as.Date(doy - 1L, origin = sprintf("%04d-01-01", year))
  if (as.integer(format(date, "%Y")) != year) as.Date(NA) else date
}

#' Validate a DSSAT weather file
#'
#' Accepts the DSSAT fixed-width daily layout as well as whitespace-delimited
#' rows. Dates must be unique and consecutive. Year bounds require January 1
#' through December 31; explicit `start_date` / `end_date` override those bounds.
#' Set `required_columns` to reject DSSAT `-99` missing markers in forcing
#' variables that a particular model configuration requires.
#'
#' @export
is_wth_valid <- function(path, end_year = NULL, required_columns = NULL, start_year = NULL,
                         start_date = NULL, end_date = NULL) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0) return(FALSE)
  tryCatch({
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    data_lines <- grep("^\\s*[0-9]{5,7}", lines, value = TRUE)
    if (!length(data_lines)) return(FALSE)
    # Parse columns in vectors, rather than running trimws/date conversion for
    # each of the millions of daily rows in a regional cache. Preserve the
    # scalar parser as a fallback for whitespace-delimited, nonstandard rows.
    fields <- vapply(c(1L, seq.int(8L, 44L, by = 6L)), function(start) {
      trimws(substr(data_lines, start, if (start == 1L) 7L else start + 5L))
    }, character(length(data_lines)))
    fields <- matrix(fields, ncol = 8L)
    fixed <- nchar(data_lines) >= 49L & grepl("^[0-9]{5,7}$", fields[, 1L]) &
      rowSums(fields[, -1L, drop = FALSE] == "") == 0L
    codes <- fields[, 1L]
    weather <- matrix(suppressWarnings(as.numeric(fields[, -1L])), ncol = 7L)
    for (i in which(!fixed)) {
      row <- .parse_wth_data_line(data_lines[i])
      if (is.null(row)) return(FALSE)
      codes[i] <- row$code
      weather[i, ] <- row$values
    }
    if (any(!is.finite(weather))) return(FALSE)
    colnames(weather) <- c("SRAD", "TMAX", "TMIN", "RAIN", "TDEW", "RH2M", "WIND")
    observed <- function(x) abs(x + 99) > 1e-6
    within <- function(x, lower, upper) {
      keep <- observed(x)
      !any(x[keep] < lower | x[keep] > upper)
    }
    if (!within(weather[, "SRAD"], 0, 60) ||
        !within(weather[, "TMAX"], -90, 70) ||
        !within(weather[, "TMIN"], -90, 70) ||
        !within(weather[, "RAIN"], 0, 2000) ||
        !within(weather[, "TDEW"], -100, 70) ||
        !within(weather[, "RH2M"], 0, 100) ||
        !within(weather[, "WIND"], 0, 100)) return(FALSE)
    comparable <- observed(weather[, "TMAX"]) & observed(weather[, "TMIN"])
    if (any(weather[comparable, "TMAX"] < weather[comparable, "TMIN"])) return(FALSE)
    required_columns <- intersect(toupper(as.character(required_columns)), colnames(weather))
    if (length(required_columns) &&
        any(!observed(weather[, required_columns, drop = FALSE]))) return(FALSE)
    short <- nchar(codes) == 5L
    if (any(!nchar(codes) %in% c(5L, 7L))) return(FALSE)
    years <- as.integer(ifelse(short, substr(codes, 1L, 2L), substr(codes, 1L, 4L)))
    years[short] <- years[short] + ifelse(years[short] < 80L, 2000L, 1900L)
    days <- as.integer(ifelse(short, substr(codes, 3L, 5L), substr(codes, 5L, 7L)))
    leap <- years %% 4L == 0L & (years %% 100L != 0L | years %% 400L == 0L)
    if (anyNA(c(years, days)) || any(days < 1L | days > 365L + leap)) return(FALSE)
    unique_years <- unique(years)
    starts <- as.Date(sprintf("%04d-01-01", unique_years))
    dates <- starts[match(years, unique_years)] + days - 1L
    if (any(is.na(dates)) || anyDuplicated(dates) ||
        (length(dates) > 1L && any(diff(dates) != 1))) return(FALSE)
    requested_start <- if (!is.null(start_date)) as.Date(start_date) else if (!is.null(start_year)) {
      as.Date(sprintf("%04d-01-01", as.integer(start_year)))
    } else NULL
    requested_end <- if (!is.null(end_date)) as.Date(end_date) else if (!is.null(end_year)) {
      as.Date(sprintf("%04d-12-31", as.integer(end_year)))
    } else NULL
    if (anyNA(c(requested_start, requested_end))) return(FALSE)
    if (!is.null(requested_start) && !is.null(requested_end) && requested_start > requested_end) return(FALSE)
    if (!is.null(requested_start) && dates[1L] > requested_start) return(FALSE)
    if (!is.null(requested_end) && tail(dates, 1L) < requested_end) return(FALSE)
    TRUE
  }, error = function(e) FALSE)
}
