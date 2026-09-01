# Fixed-column soil preflight. Formatting only, not agronomic plausibility.
# DSSAT's -99 missing-value sentinel is allowed.
soil_file_issue <- function(path) {
  if (!file.exists(path)) return("SOIL.SOL is missing")
  lines <- tryCatch(readLines(path, warn = FALSE, encoding = "UTF-8"), error = function(e) character())
  if (!length(lines)) return("SOIL.SOL is empty or unreadable")
  errors <- grep("^\\*SOIL ERROR|^Source missing|^No Soil ID", lines, value = TRUE)
  if (length(errors)) return(paste(trimws(errors), collapse = " | "))
  headers <- grep("^@\\s+SLB\\b", lines)
  if (!length(headers)) return("SOIL.SOL has no @ SLB layer table")
  for (h in headers) {
    matches <- gregexpr("\\S+", lines[h])[[1]]
    ends <- (matches + attr(matches, "match.length") - 1L)[-1]
    fields <- regmatches(lines[h], list(matches))[[1]][-1]
    starts <- c(1L, head(ends, -1L) + 2L)
    depths <- numeric()
    for (idx in seq_len(length(lines) - h) + h) {
      row <- lines[idx]
      if (grepl("^[@*]", row)) break
      if (!nzchar(trimws(row)) || grepl("^\\s*!", row)) next
      token <- strsplit(trimws(row), "\\s+")[[1]][1]
      depth <- suppressWarnings(as.numeric(substr(row, 1, ends[1])))
      token_depth <- suppressWarnings(as.numeric(token))
      values <- suppressWarnings(as.numeric(substring(row, starts, ends)))
      separators <- substring(row, head(ends, -1L) + 1L, head(ends, -1L) + 1L)
      if (!is.finite(depth) || !is.finite(token_depth) || depth != token_depth ||
          depth <= 0 || depth != floor(depth) || any(!is.finite(values[fields != "SLMH"])) ||
          any(nzchar(separators) & !grepl("^\\s$", separators))) {
        return("SOIL.SOL has invalid fixed-width layer fields; regenerate the derived SOL with the corrected writer")
      }
      depths <- c(depths, depth)
    }
    if (!length(depths)) return("SOIL.SOL has no parseable SLB layer depths")
    if (length(depths) > 19) return(sprintf("SOIL.SOL has %d layers; DSSAT accepts at most 19", length(depths)))
    if (any(diff(depths) <= 0)) return(paste0("SOIL.SOL layer depths are not strictly increasing: ", paste(depths, collapse = ",")))
  }
  NULL
}

# Reformat a saved mapping into a NEW directory: no API calls, no recalculation
# of hydraulic properties, and no overwriting existing caches (even empty dirs).
rebuild_soil_files_from_mapping <- function(mapping_csv, output_dir, soil_source) {
  source <- toupper(soil_source)
  if (!source %in% c("SSURGO", "GNATSGO")) stop("soil_source must be SSURGO or GNATSGO")
  if (file.exists(output_dir)) stop("output_dir must be a new directory; existing caches are never overwritten")
  data <- utils::read.csv(mapping_csv, colClasses = c(ID = "character"), stringsAsFactors = FALSE)
  required <- c("ID", "latitude", "longitude", "SLLL", "SDUL", "SSAT", "bulk_density", "om_pct", "clay_pct", "silt_pct")
  if (!nrow(data) || !all(required %in% names(data))) stop("soil mapping is empty or lacks required profile columns")
  if (anyNA(data$ID) || any(!grepl("^[A-Za-z0-9_-]{1,10}$", data$ID))) stop("soil mapping contains invalid IDs")
  if (!any(c("depth_range", "depth_bottom") %in% names(data))) stop("soil mapping lacks layer depths")
  if (!"depth_range" %in% names(data)) data$depth_range <- paste0("0-", data$depth_bottom, "cm")
  dir.create(output_dir, recursive = TRUE)
  writer <- if (source == "SSURGO") format_dssat_soil_single else format_dssat_soil_gnatsgo
  ids <- unique(data$ID)
  paths <- file.path(output_dir, paste0(ids, ".SOL"))
  for (i in seq_along(ids)) {
    writer(data[data$ID == ids[i], , drop = FALSE], output_dir)
    issue <- soil_file_issue(paths[i])
    if (!is.null(issue)) stop(sprintf("%s: %s", ids[i], issue))
  }
  data.frame(ID = ids, path = paths, stringsAsFactors = FALSE)
}
