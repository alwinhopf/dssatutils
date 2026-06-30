# Provider-side weather missing-value normalization.
#
# Some APIs can return true NA/NaN/Inf values rather than their documented
# numeric missing sentinel. DSSAT weather files must remain numeric, so provider
# writers normalize those values to -99 and log the correction beside the
# post-download weather repair logs.

.weather_log_lines <- function(path, lines) {
  if (is.null(path) || !nzchar(path) || !length(lines)) return(invisible())
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "a", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con = con, sep = "\n", useBytes = TRUE)
}

.normalize_weather_missing_values <- function(weather_data,
                                              point_id,
                                              output_file,
                                              output_dir,
                                              source_label,
                                              variables = c("SRAD", "TMAX", "TMIN", "RAIN", "TDEW", "RH2M", "WIND")) {
  variables <- intersect(variables, names(weather_data))
  if (!length(variables)) return(weather_data)

  log_lines <- character()
  for (var in variables) {
    x <- weather_data[[var]]
    bad <- is.na(x) | !is.finite(x)
    if (!any(bad)) next

    dates <- as.character(weather_data$DATE[bad])
    weather_data[[var]][bad] <- -99
    log_lines <- c(log_lines, sprintf(
      "%s file=%s id=%s variable=%s status=normalized_missing dates=%s..%s count=%d fill_value=-99 method=source_missing_to_DSSAT_missing source=%s",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      basename(output_file),
      point_id,
      var,
      dates[1],
      dates[length(dates)],
      length(dates),
      source_label
    ))
  }

  if (length(log_lines)) {
    .weather_log_lines(file.path(output_dir, "weather_repair.log"), log_lines)
  }
  weather_data
}
