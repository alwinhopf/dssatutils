# R-only platform policy; Python providers use spawn-safe process executors.
.soil_use_fork <- function(n_cores, os_type = .Platform$OS.type,
                           sysname = Sys.info()[["sysname"]]) {
  n_cores > 1L && os_type != "windows" && sysname != "Darwin"
}
