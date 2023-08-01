library('move2')

test_data <- test_data("input2.rds") #file must be move2!

test_that("happy filter path", {
  actual <- rFunction(data = test_data, duration = 172800, radius = 30000, annot = FALSE)
  expect_equal(nrow(actual), 6994)
})

test_that("happy annot path", {
  actual <- rFunction(data = test_data, duration = 172800, radius = 30000, annot = TRUE)
  expect_equal(nrow(actual), 7770)
})

test_that("happy one stopover", {
  actual <- rFunction(data = test_data, duration = 172800, radius = 200, annot = FALSE)
  expect_equal(length(unique(mt_track_id(actual))),1)
})

test_that("no stopovers null", {
  actual <- rFunction(data = test_data, duration = 172800, radius = 100, annot = FALSE)
  expect_null(actual)
})

#note that artefacts are not correctly written in the test.. this test passes if the write.csv lines are commented out