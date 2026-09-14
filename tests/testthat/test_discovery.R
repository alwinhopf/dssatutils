test_that("discovery helpers return expected paths or NULL", {
  expect_null(find_rscript("non_existent_binary_xyz_12345"))
  expect_null(find_dssat("non_existent_binary_xyz_12345"))
  expect_null(find_mpi_runner("non_existent_binary_xyz_12345"))

  # find_rscript should find the running Rscript
  rscript <- find_rscript()
  if (nzchar(Sys.which("Rscript"))) {
    expect_false(is.null(rscript))
    expect_true(file.exists(rscript))
  }

  # Environment variable overrides
  temp_file <- tempfile("fake_dssat")
  writeLines("#!/bin/sh\nexit 0", temp_file)
  Sys.chmod(temp_file, "0755")
  on.exit(unlink(temp_file), add = TRUE)

  withr_available <- requireNamespace("withr", quietly = TRUE)
  if (withr_available) {
    withr::with_envvar(c(DSSAT_EXE = temp_file), {
      found <- find_dssat()
      expect_equal(normalizePath(found, winslash = "/", mustWork = FALSE),
                   normalizePath(temp_file, winslash = "/", mustWork = FALSE))
    })
  }
})
