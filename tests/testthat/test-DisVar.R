# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(DisVar)

# tests/testthat/test-DisVar.R
context("DisVar")

test_that("DisVar processes VCF file correctly", {
  # Create a temporary directory to store the test output file
  temp_output_dir <- tempfile()
  dir.create(temp_output_dir)

  # Run DisVar function on the test VCF file
  vcf_file <- "tests/testthat/test_file.vcf"  # Update the file path to the test VCF file
  result <- DisVar(vcf_file)

  # Load the generated 'aligned_df' from the test output file
  output_file <- file.path(temp_output_dir, "test_file_DisVar.txt")
  expected <- read.table(output_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Clean up temporary directory
  unlink(temp_output_dir, recursive = TRUE)

  # Perform the test assertions
  expect_equal(nrow(result), nrow(expected))
  expect_identical(colnames(result), colnames(expected))
  expect_identical(result$Disease, expected$Disease)
  # Add more specific comparisons for other columns as needed
})
