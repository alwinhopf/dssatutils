# Shared solar radiation helpers (extraterrestrial radiation and Bristow-Campbell estimation).
# Twin of python/cropmodel_data/weather_solar.py

extraterrestrial_radiation <- function(lat_deg, doy) {
  # Daily extraterrestrial radiation Ra (MJ/m²/day) and sunset hour angle ws (rad).
  # Implements FAO-56 Eq. 21 (Allen et al. 1998).
  phi <- lat_deg * pi / 180
  dr <- 1 + 0.033 * cos(2 * pi / 365 * doy)
  decl <- 0.409 * sin(2 * pi / 365 * doy - 1.39)
  arg <- pmin(pmax(-tan(phi) * tan(decl), -1), 1)
  ws <- acos(arg)
  Gsc <- 0.0820
  Ra <- (24 * 60 / pi) * Gsc * dr * (ws * sin(phi) * sin(decl) +
                                     cos(phi) * cos(decl) * sin(ws))
  list(Ra = pmax(Ra, 0), ws = ws)
}

estimate_srad_bristow_campbell <- function(dates, tmax, tmin, lat_deg, a = 0.75, c = 2.4) {
  # Estimate daily solar radiation (MJ/m²/day) from temperature range.
  # Method: Bristow, K. L., & Campbell, G. S. (1984). Agricultural and Forest Meteorology 31(2): 159-166.
  tmax <- as.numeric(tmax)
  tmin <- as.numeric(tmin)
  n <- length(tmax)
  if (n == 0) return(numeric())

  # Parse dates to get DOY, calendar month, and continuity check
  if (inherits(dates, "Date") || inherits(dates, "POSIXt")) {
    dts <- as.Date(dates)
  } else {
    s <- as.character(dates)
    if (length(s) > 0 && nchar(s[1]) == 7 && grepl("^[0-9]{7}$", s[1])) {
      dts <- as.Date(s, format = "%Y%j")
    } else {
      dts <- as.Date(s)
    }
  }
  doy <- as.integer(format(dts, "%j"))
  months <- as.integer(format(dts, "%m"))

  # Diurnal temperature range: Delta_T(t) = Tmax(t) - [Tmin(t) + Tmin(t+1)] / 2
  # For terminal day or across date gaps (non-consecutive days), Delta_T(t) = Tmax(t) - Tmin(t).
  tmin_next <- c(tmin[-1], tmin[n])
  if (n > 1) {
    consecutive <- c(as.numeric(diff(dts)) == 1, FALSE)
  } else {
    consecutive <- FALSE
  }

  use_tmin_next <- consecutive & is.finite(tmin_next)
  tmin_effective <- ifelse(use_tmin_next, (tmin + tmin_next) / 2, tmin)
  delta_t <- tmax - tmin_effective
  delta_t <- ifelse(is.finite(delta_t), pmax(delta_t, 0), NA_real_)

  # Monthly mean diurnal range Delta_T_bar_m across valid days
  valid_dt <- delta_t[is.finite(delta_t) & delta_t > 0]
  fallback_dt_bar <- if (length(valid_dt) > 0) mean(valid_dt) else 10.0
  if (!is.finite(fallback_dt_bar) || fallback_dt_bar <= 0) fallback_dt_bar <- 10.0

  dt_bar_m <- rep(fallback_dt_bar, 12)
  for (m in 1:12) {
    m_mask <- (months == m) & is.finite(delta_t) & (delta_t > 0)
    if (any(m_mask)) {
      dt_bar_m[m] <- mean(delta_t[m_mask])
    }
  }

  dt_bar <- dt_bar_m[months]
  b <- 0.036 * exp(-0.154 * dt_bar)

  # Atmospheric transmittance Tt = A * [1 - exp(-B * Delta_T^C)]
  tt <- a * (1 - exp(-b * (delta_t ^ c)))

  # Extraterrestrial radiation Ra
  e <- extraterrestrial_radiation(lat_deg, doy)

  # Solar radiation Rs = Tt * Ra clipped to [0, 0.80 * Ra]
  rs <- tt * e$Ra
  rs <- ifelse(e$Ra <= 0, 0, rs)
  rs <- pmin(pmax(rs, 0), 0.80 * e$Ra)

  # Missing inputs produce missing output
  invalid <- !is.finite(tmax) | !is.finite(tmin)
  rs[invalid] <- NA_real_
  rs
}
