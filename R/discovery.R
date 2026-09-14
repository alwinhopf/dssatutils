# Centralized executable discovery helpers for DSSAT workflows.

.is_executable_file <- function(path) {
  if (is.null(path) || !nzchar(path)) return(FALSE)
  if (!file.exists(path) || dir.exists(path)) return(FALSE)
  if (.Platform$OS.type == "windows") return(TRUE)
  file.access(path, 1L) == 0L
}

.resolve_candidate <- function(candidate) {
  if (is.null(candidate) || !length(candidate)) return(NULL)
  cand <- trimws(as.character(candidate)[1])
  if (!nzchar(cand)) return(NULL)

  # Check if candidate has directory components or exists directly
  has_dir <- grepl("[/\\\\]", cand)
  expanded <- path.expand(cand)

  if (has_dir || file.exists(expanded)) {
    if (.is_executable_file(expanded)) {
      return(normalizePath(expanded, winslash = "/", mustWork = FALSE))
    }
    if (.Platform$OS.type == "windows") {
      for (ext in c(".exe", ".cmd", ".bat")) {
        with_ext <- paste0(expanded, ext)
        if (.is_executable_file(with_ext)) {
          return(normalizePath(with_ext, winslash = "/", mustWork = FALSE))
        }
      }
    }
  }

  # Check via Sys.which
  found <- Sys.which(cand)
  if (nzchar(found) && .is_executable_file(found)) {
    return(normalizePath(found, winslash = "/", mustWork = FALSE))
  }

  if (.Platform$OS.type == "windows" && !grepl("\\.(exe|cmd|bat)$", cand, ignore.case = TRUE)) {
    for (ext in c(".exe", ".cmd", ".bat")) {
      found <- Sys.which(paste0(cand, ext))
      if (nzchar(found) && .is_executable_file(found)) {
        return(normalizePath(found, winslash = "/", mustWork = FALSE))
      }
    }
  }

  NULL
}

#' Find Rscript executable
#'
#' @param explicit Optional explicit path to check first.
#' @return Absolute path to Rscript if found, NULL otherwise.
#' @export
find_rscript <- function(explicit = NULL) {
  if (!is.null(explicit) && nzchar(trimws(as.character(explicit)[1]))) {
    return(.resolve_candidate(explicit))
  }

  env_val <- Sys.getenv("RSCRIPT", unset = "")
  if (nzchar(env_val)) {
    resolved <- .resolve_candidate(env_val)
    if (!is.null(resolved)) return(resolved)
  }

  r_home_bin <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  if (.is_executable_file(r_home_bin)) {
    return(normalizePath(r_home_bin, winslash = "/", mustWork = FALSE))
  }

  for (cand in c("Rscript", "Rscript.exe")) {
    resolved <- .resolve_candidate(cand)
    if (!is.null(resolved)) return(resolved)
  }

  candidates <- if (Sys.info()["sysname"] == "Darwin") {
    c("/usr/local/bin/Rscript", "/opt/homebrew/bin/Rscript",
      "/Library/Frameworks/R.framework/Resources/bin/Rscript")
  } else if (.Platform$OS.type == "windows") {
    c("C:/Program Files/R/R-4.3.0/bin/Rscript.exe",
      "C:/Program Files/R/R-4.3.1/bin/Rscript.exe",
      "C:/Program Files/R/R-4.3.2/bin/Rscript.exe",
      "C:/Program Files/R/R-4.3.3/bin/Rscript.exe")
  } else {
    c("/usr/bin/Rscript", "/usr/local/bin/Rscript")
  }

  for (loc in candidates) {
    if (.is_executable_file(loc)) {
      return(normalizePath(loc, winslash = "/", mustWork = FALSE))
    }
  }

  NULL
}

#' Find DSSAT CSM executable
#'
#' @param explicit Optional explicit path to check first.
#' @return Absolute path to DSSAT executable if found, NULL otherwise.
#' @export
find_dssat <- function(explicit = NULL) {
  if (!is.null(explicit) && nzchar(trimws(as.character(explicit)[1]))) {
    return(.resolve_candidate(explicit))
  }

  env_val <- Sys.getenv("DSSAT_EXE", unset = "")
  if (nzchar(env_val)) {
    resolved <- .resolve_candidate(env_val)
    if (!is.null(resolved)) return(resolved)
  }

  binary_names <- c("dscsm048", "dscsm048.exe", "dscsm047", "dscsm047.exe",
                    "dscsm046", "dscsm046.exe", "dscsm045", "dscsm045.exe",
                    "dscsm", "dscsm.exe")

  for (env_key in c("DSSAT_HOME", "DSSAT_DIR")) {
    env_dir <- Sys.getenv(env_key, unset = "")
    if (nzchar(env_dir) && dir.exists(env_dir)) {
      for (name in binary_names) {
        cand <- file.path(env_dir, name)
        if (.is_executable_file(cand)) {
          return(normalizePath(cand, winslash = "/", mustWork = FALSE))
        }
      }
    }
  }

  for (name in binary_names) {
    resolved <- .resolve_candidate(name)
    if (!is.null(resolved)) return(resolved)
  }

  search_dirs <- if (.Platform$OS.type == "windows") {
    c("C:/DSSAT48", "C:/DSSAT47", "C:/DSSAT46", "C:/dssat48")
  } else {
    c("/opt/dssat48", "/usr/local/dssat48", path.expand("~/DSSAT48"), path.expand("~/dssat48"))
  }

  for (sdir in search_dirs) {
    if (dir.exists(sdir)) {
      for (name in binary_names) {
        cand <- file.path(sdir, name)
        if (.is_executable_file(cand)) {
          return(normalizePath(cand, winslash = "/", mustWork = FALSE))
        }
      }
    }
  }

  NULL
}

#' Find MPI runner executable
#'
#' @param explicit Optional explicit path to check first.
#' @return Absolute path to MPI runner if found, NULL otherwise.
#' @export
find_mpi_runner <- function(explicit = NULL) {
  if (!is.null(explicit) && nzchar(trimws(as.character(explicit)[1]))) {
    return(.resolve_candidate(explicit))
  }

  for (env_key in c("MPIEXEC", "MPIRUN")) {
    env_val <- Sys.getenv(env_key, unset = "")
    if (nzchar(env_val)) {
      resolved <- .resolve_candidate(env_val)
      if (!is.null(resolved)) return(resolved)
    }
  }

  for (cand in c("mpiexec", "mpiexec.exe", "mpirun", "mpirun.exe", "srun")) {
    resolved <- .resolve_candidate(cand)
    if (!is.null(resolved)) return(resolved)
  }

  NULL
}

#' Find Python executable
#'
#' @param explicit Optional explicit path to check first.
#' @return Absolute path to Python executable if found, NULL otherwise.
#' @export
find_python <- function(explicit = NULL) {
  if (!is.null(explicit) && nzchar(trimws(as.character(explicit)[1]))) {
    return(.resolve_candidate(explicit))
  }

  for (env_key in c("PYTHON", "RETICULATE_PYTHON")) {
    env_val <- Sys.getenv(env_key, unset = "")
    if (nzchar(env_val)) {
      resolved <- .resolve_candidate(env_val)
      if (!is.null(resolved)) return(resolved)
    }
  }

  for (cand in c("python3", "python3.exe", "python", "python.exe")) {
    resolved <- .resolve_candidate(cand)
    if (!is.null(resolved)) return(resolved)
  }

  NULL
}
