# Shared credential helpers for providers that need Copernicus CDS access.

.dssatutils_cds_default_url <- function() {
  if (exists(".dssatutils_config_get", mode = "function")) {
    .dssatutils_config_get(
      "cds.default_url",
      "https://cds.climate.copernicus.eu/api"
    )
  } else {
    "https://cds.climate.copernicus.eu/api"
  }
}

.dssatutils_cds_rc_candidates <- function() {
  candidates <- c(
    Sys.getenv("CDSAPI_RC", unset = ""),
    path.expand("~/.cdsapirc"),
    file.path(Sys.getenv("USERPROFILE", unset = ""), ".cdsapirc"),
    file.path(Sys.getenv("HOME", unset = ""), ".cdsapirc")
  )
  unique(candidates[nzchar(candidates)])
}

.dssatutils_read_cdsapirc <- function(paths = .dssatutils_cds_rc_candidates()) {
  for (path in paths) {
    if (!file.exists(path)) next
    lines <- readLines(path, warn = FALSE)
    key_line <- grep("^\\s*key\\s*:", lines, value = TRUE)
    if (!length(key_line)) next
    url_line <- grep("^\\s*url\\s*:", lines, value = TRUE)
    token <- trimws(sub("^\\s*key\\s*:\\s*", "", key_line[[1]]))
    url <- if (length(url_line)) trimws(sub("^\\s*url\\s*:\\s*", "", url_line[[1]])) else ""
    if (nzchar(token)) {
      return(list(
        key = token,
        url = if (nzchar(url)) url else .dssatutils_cds_default_url(),
        path = path
      ))
    }
  }
  NULL
}

.dssatutils_prompt_secret <- function(prompt) {
  if (requireNamespace("getPass", quietly = TRUE)) {
    return(getPass::getPass(prompt))
  }
  readline(prompt)
}

#' Configure Copernicus CDS credentials for dssatutils
#'
#' CDS-backed sources such as AgERA5 and ERA5-Land need a Copernicus Climate
#' Data Store Personal Access Token. This helper accepts an explicit token, uses
#' `CDSAPI_KEY`/`CDSAPI_URL`, imports an existing `.cdsapirc`, or prompts in an
#' interactive session. When `ecmwfr` is installed it also stores the token in
#' the local keyring entry used by R download functions.
#'
#' @param token CDS Personal Access Token. If omitted, env vars, `.cdsapirc`, or
#'   an interactive prompt are used.
#' @param url CDS API URL.
#' @param user `ecmwfr` keyring user entry.
#' @param rc_path Path to write/read `.cdsapirc`.
#' @param write_cdsapirc Write a Python/cdsapi-compatible `.cdsapirc`.
#' @param set_ecmwfr_key Store the token for `ecmwfr::wf_request()`.
#' @param overwrite Replace an existing `.cdsapirc` when writing.
#' @param prompt Ask for a token if none is found and the session is interactive.
#' @param quiet Suppress setup messages.
#' @return Invisibly returns credential metadata without echoing the token.
#' @export
setup_cds_credentials <- function(token = NULL,
                                  url = NULL,
                                  user = "ecmwfr",
                                  rc_path = NULL,
                                  write_cdsapirc = TRUE,
                                  set_ecmwfr_key = TRUE,
                                  overwrite = FALSE,
                                  prompt = interactive(),
                                  quiet = FALSE) {
  if (is.null(token) || !nzchar(token)) token <- Sys.getenv("CDSAPI_KEY", unset = "")
  default_url <- .dssatutils_cds_default_url()
  if (is.null(url) || !nzchar(url)) url <- default_url
  env_url <- Sys.getenv("CDSAPI_URL", unset = "")
  if (nzchar(env_url) && identical(url, default_url)) url <- env_url

  rc <- NULL
  if (!nzchar(token) && !isTRUE(overwrite)) {
    rc <- .dssatutils_read_cdsapirc(if (is.null(rc_path)) .dssatutils_cds_rc_candidates() else rc_path)
    if (!is.null(rc)) {
      token <- rc$key
      url <- rc$url
      rc_path <- rc$path
    }
  }

  if (!nzchar(token) && isTRUE(prompt) && interactive()) {
    message(
      "Copernicus CDS credentials are required. Create a Personal Access Token at ",
      "https://cds.climate.copernicus.eu/how-to-api"
    )
    token <- .dssatutils_prompt_secret("Enter Copernicus CDS Personal Access Token: ")
  }

  if (!nzchar(token)) {
    stop(
      "Copernicus CDS credentials were not found. Set CDSAPI_KEY/CDSAPI_URL, ",
      "create ~/.cdsapirc, or run setup_cds_credentials(token = '<PAT>').",
      call. = FALSE
    )
  }

  if (is.null(rc_path) || !nzchar(rc_path)) {
    rc_env <- Sys.getenv("CDSAPI_RC", unset = "")
    rc_path <- if (nzchar(rc_env)) rc_env else path.expand("~/.cdsapirc")
  }

  if (isTRUE(write_cdsapirc) && (isTRUE(overwrite) || !file.exists(rc_path))) {
    dir.create(dirname(rc_path), recursive = TRUE, showWarnings = FALSE)
    writeLines(c(sprintf("url: %s", url), sprintf("key: %s", token)), rc_path, useBytes = TRUE)
    try(Sys.chmod(rc_path, mode = "0600"), silent = TRUE)
    if (!quiet) message("  Copernicus CDS credentials written to ", rc_path)
  }

  Sys.setenv(CDSAPI_URL = url, CDSAPI_KEY = token, CDSAPI_RC = rc_path)

  if (isTRUE(set_ecmwfr_key) && requireNamespace("ecmwfr", quietly = TRUE)) {
    tryCatch(
      ecmwfr::wf_set_key(key = token, user = user),
      error = function(e) {
        if (!quiet) {
          message("  Could not store CDS token for ecmwfr keyring entry '", user, "': ",
                  conditionMessage(e))
        }
      }
    )
  }

  invisible(list(url = url, path = rc_path, user = user, has_key = TRUE))
}

.dssatutils_ensure_cds_credentials <- function(user = "ecmwfr",
                                               prompt = interactive(),
                                               quiet = FALSE,
                                               require_ecmwfr = FALSE) {
  if (isTRUE(require_ecmwfr) && !requireNamespace("ecmwfr", quietly = TRUE)) {
    stop("This CDS source needs the 'ecmwfr' package. install.packages('ecmwfr')",
         call. = FALSE)
  }

  if (nzchar(Sys.getenv("CDSAPI_KEY", unset = ""))) {
    setup_cds_credentials(
      user = user,
      write_cdsapirc = FALSE,
      set_ecmwfr_key = requireNamespace("ecmwfr", quietly = TRUE),
      prompt = FALSE,
      quiet = TRUE
    )
    return(invisible(TRUE))
  }

  if (requireNamespace("ecmwfr", quietly = TRUE)) {
    existing <- tryCatch(ecmwfr::wf_get_key(user = user), error = function(e) "")
    if (nzchar(existing)) return(invisible(TRUE))
  }

  rc <- .dssatutils_read_cdsapirc()
  if (!is.null(rc)) {
    setup_cds_credentials(
      token = rc$key,
      url = rc$url,
      user = user,
      rc_path = rc$path,
      write_cdsapirc = FALSE,
      set_ecmwfr_key = requireNamespace("ecmwfr", quietly = TRUE),
      prompt = FALSE,
      quiet = quiet
    )
    return(invisible(TRUE))
  }

  setup_cds_credentials(user = user, prompt = prompt, quiet = quiet)
  invisible(TRUE)
}

#' Save the CDS API token in the local keyring for ecmwfr
#'
#' Backwards-compatible alias for older ERA5-Land scripts. Prefer
#' `setup_cds_credentials()` for new code.
#'
#' @param token Personal API token from the CDS profile page.
#' @param user Keyring entry name. Defaults to "ecmwfr".
#' @export
era5land_set_cds_key <- function(token, user = "ecmwfr") {
  setup_cds_credentials(token = token, user = user, write_cdsapirc = TRUE)
  invisible(user)
}
