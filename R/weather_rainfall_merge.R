# Shared helper: replace a weather frame's RAIN column from daily rainfall.

merge_rainfall_into_weather <- function(weather_data, rainfall,
                                        date_col = "DATE", rain_col = "RAIN") {
  if (is.null(rainfall) || length(rainfall) == 0) {
    return(list(weather_data = weather_data, n_replaced = 0L))
  }
  if (is.null(names(rainfall))) {
    stop("rainfall must be a named vector keyed by DSSAT DATE codes (YYYYDOY).",
         call. = FALSE)
  }
  idx <- match(as.character(weather_data[[date_col]]), names(rainfall))
  hit <- !is.na(idx)
  if (any(hit)) {
    weather_data[[rain_col]] <- as.numeric(weather_data[[rain_col]])
    weather_data[[rain_col]][hit] <- as.numeric(rainfall[idx[hit]])
  }
  list(weather_data = weather_data, n_replaced = sum(hit))
}
