library(testthat)
library(dssatutils)

test_that('mixed POLARIS batch preserves its usable point', {
  rows <- read.csv(test_path('..','fixtures','polaris_mixed_profiles.csv'), colClasses=c(ID='character'))
  local_mocked_bindings(fetch_polaris=function(...) rows, .package='dssatutils')
  g <- sf::st_as_sf(data.frame(ID=c('00000001','00000002'), x=c(-84.6,-84.61),y=c(30.55,30.55)), coords=c('x','y'),crs=4326)
  out <- tempfile(); dir.create(out); on.exit(unlink(out, recursive=TRUE))
  process_soils_polaris(g,file.path(out,'map.csv'),out,'ID')
  expect_null(soil_file_issue(file.path(out,'00000001.SOL')))
  expect_false(file.exists(file.path(out,'00000002.SOL')))
  expect_match(paste(readLines(file.path(out,'soil_processing_errors.log')),collapse=' '),'00000002')
})

test_that('connection retries stop with a distinguishable failure', {
  calls <- 0L; delays <- numeric()
  op <- function() { calls <<- calls+1L; stop('Could not resolve hostname [example.org]') }
  expect_error(dssatutils:::.provider_retry(op,sleep=function(x) delays <<- c(delays,x)), class='dssat_connectivity_error')
  expect_equal(calls,3L); expect_equal(delays,c(5,10))
  calls <- 0L
  expect_error(dssatutils:::.provider_retry(function(){calls <<- calls+1L;stop('invalid response schema')}), 'invalid response schema')
  expect_equal(calls,1L)
})

test_that('empty WCS response is retryable, not no coverage', {
  local_mocked_bindings(mukey.wcs=function(...) NULL, .package='soilDB')
  result <- dssatutils:::gnatsgo_mukey_at_point(30.55,-84.6,max_retries=1)
  expect_false(result$ok)
  expect_match(result$error,'WCS extraction failed')
})
