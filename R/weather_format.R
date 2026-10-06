# DSSAT PARSE_HEADERS skips the separator before each daily numeric field.
# Preserve that blank and fit the numeric token into the remaining five columns.
.format_wth_value <- function(x, decimals = 1L) {
  values <- suppressWarnings(as.numeric(x))
  precision <- rep_len(as.integer(decimals), length(values))
  vapply(seq_along(values), function(i) {
    value <- values[i]
    if (!is.finite(value) || value == -99) return("   -99")
    for (digits in seq.int(precision[i], 0L)) {
      token <- sprintf(paste0("%.", digits, "f"), value)
      if (nchar(token) <= 5L) return(sprintf("%6s", token))
    }
    stop(sprintf("DSSAT daily value cannot fit a five-character token: %s", value), call. = FALSE)
  }, character(1), USE.NAMES = FALSE)
}

.weather_wind_run <- function(x) {
  # Provider m/s -> DSSAT km/day; do not multiply missing sentinels.
  x <- suppressWarnings(as.numeric(x))
  good <- is.finite(x) & x >= 0
  x[good] <- x[good] * 86.4
  x[!good] <- -99
  x
}
