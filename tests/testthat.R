# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

# testthat harness entry point. `R CMD check` (and `devtools::test()`) sources
# this single file; everything under tests/testthat/test-*.R is discovered and
# run by the `test_check()` call below.
#
# `library(pkg)` attaches a package's exported names to the search path (like a
# wildcard import). testthat is the test framework; razer is the package under
# test.
library(testthat)
library(razer)

# `test_check("razer")` finds, runs, and reports every test-*.R file in this
# directory against the installed `razer` package. This is the line that wires
# the suite into `R CMD check`.
test_check("razer")
