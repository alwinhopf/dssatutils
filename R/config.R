# Shared package configuration helpers.

.dssatutils_config_cache <- new.env(parent = emptyenv())

.dssatutils_config_candidates <- function() {
  candidates <- c(
    system.file("config.yml", package = "dssatutils"),
    file.path(getwd(), "config.yml"),
    file.path(getwd(), "config.yaml"),
    Sys.getenv("DSSATUTILS_CONFIG", unset = "")
  )
  unique(candidates[nzchar(candidates)])
}

.dssatutils_deep_merge <- function(base, override) {
  if (!is.list(base)) base <- list()
  if (!is.list(override)) return(base)
  for (name in names(override)) {
    if (is.list(base[[name]]) && is.list(override[[name]])) {
      base[[name]] <- .dssatutils_deep_merge(base[[name]], override[[name]])
    } else {
      base[[name]] <- override[[name]]
    }
  }
  base
}

.dssatutils_load_config <- function(refresh = FALSE) {
  if (!refresh && exists("config", envir = .dssatutils_config_cache, inherits = FALSE)) {
    return(get("config", envir = .dssatutils_config_cache))
  }
  cfg <- list()
  for (candidate in .dssatutils_config_candidates()) {
    if (!file.exists(candidate)) next
    if (requireNamespace("yaml", quietly = TRUE)) {
      next_cfg <- yaml::read_yaml(candidate)
      if (is.null(next_cfg)) next_cfg <- list()
      cfg <- .dssatutils_deep_merge(cfg, next_cfg)
      next
    }
    warning("Package 'yaml' is not installed; dssatutils config.yml was not loaded.")
    break
  }
  assign("config", cfg, envir = .dssatutils_config_cache)
  cfg
}

.dssatutils_config_get <- function(path, default = NULL, refresh = FALSE) {
  cfg <- .dssatutils_load_config(refresh = refresh)
  value <- cfg
  for (part in strsplit(path, ".", fixed = TRUE)[[1]]) {
    if (!is.list(value) || is.null(value[[part]])) return(default)
    value <- value[[part]]
  }
  if (is.null(value)) default else value
}

.dssatutils_config_bool <- function(path, default = FALSE) {
  value <- .dssatutils_config_get(path, default)
  if (is.logical(value)) return(isTRUE(value))
  if (is.numeric(value)) return(!is.na(value) && value != 0)
  if (is.character(value)) return(tolower(value) %in% c("true", "t", "yes", "y", "1"))
  isTRUE(default)
}

.dssatutils_config_number <- function(path, default) {
  value <- suppressWarnings(as.numeric(.dssatutils_config_get(path, default)))
  if (is.na(value)) default else value
}
