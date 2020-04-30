# Tests for the example dataset
#
# This file is commented to describe the structure of the testing functions but
# these comments can probably be removed in the future.
#
# We also have a few more tests than are really necessary as examples.
#
# Tests can be run using devtools::test() but are also run automatically
# whenever you do devtools::check().
#
# If you want to which parts of your code need tests try using the
# covr::report() function.

# We can run code outside a test if we need to. For example sometimes we get the
# output of a function here and run multiple tests on it.
# output_to_test <- function_to_test(...)

# Each test starts with the `test_that()` function. You can have as many one
# these as you like but each one should only test on thing.
#
# The first argument is the description of the test. This is usually read as a
# sentence for example "test that something works". The second argument is a
# chunk of code containing the actual test.
#
# This first test checks that the example dataset is a data.frame.
test_that("example data is a data.frame", {
    # Inside the test we have a series of expectations which compare some output
    # to what we expect it to be.
    #
    # Expections all have this form:
    #
    # expect_TYPE_OF_EXPECTATION(thing_to_check, expected_value)

    # We can run code outside an expection. Again often it is neater to get some
    # output first then check it.
    example_data <- speckle_example_data()

    # First we check that the example data is a list
    expect_type(example_data, "list")

    # Then we can check that example data is the list class
    expect_s3_class(example_data, "data.frame")
})

# Here is a second test that checks the column names are correct
test_that("example data has correct column names", {
    col_names <- colnames(speckle_example_data())

    expect_equal(col_names, c("clusters", "samples", "group"))
})

# We can also check the size of the dataset
test_that("example data has correct dimensions", {
    dims <- dim(speckle_example_data())

    expect_equal(dims[1], 4600)
    expect_equal(dims[2], 3)
})

# Often it is really helpful to check that things break when they are supposed,
# for example if someone provides incorrect arguments.
#
# Here is an example of checking for an error.
test_that("example data function doesn't need arguments", {
    expect_error(speckle_example_data("unused argument"),
                 # The second argument to this expectation is a regex that
                 # matches the error message
                 "unused argument")
})

# There are many other expect functions. It is worth taking a quick look through
# them before writing tests try and use the simplest one.

# There are some other handy functions in the testthat package. One set is the
# skip functions that can be used to skip tests for some reason. One
# particularly useful one is skip_if_not_installed which can be used to skip
# tests that need suggested packages which might not be installed:
# test_that("only runs if suggested packages is installed", {
#     skip_if_not_installed("suggestd_pkg")
#     expect_something(...)
# })
#
# There are also skip_on_bioc and skip_on_cran if you don't want to run tests
# on their build systems (for example if they are really long).
