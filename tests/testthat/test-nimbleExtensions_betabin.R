library(extraDistr)     # Import extra distribution library for comparisson to NIMBLE extension functions

## 1. ------ TEST THE BETA-BINOMIAL DENSITY FUNCTION ------
test_that("beta-binomial density function gives correct values", {
  # Create a matrix of test values to use in the unit test
  testGrid <- as.matrix(do.call(rbind, lapply(X = 1:5 * 2, FUN = function(curSize, inGrid) {
    # Create a data.frame of a number of different sizes and values to calculate the density of
    do.call(rbind, lapply(X = 0:curSize, FUN = function(curK, curSize, inGrid) {
      cbind(inGrid, data.frame(
        size = rep(curSize, nrow(inGrid)),
        K = rep(curK, nrow(inGrid))
      ))
    }, curSize = curSize, inGrid = inGrid))
  },
    # Create a data.frame of different alpha and beta values for the beta-binomial
    inGrid = expand.grid(alpha = seq(0.1, 20.0, length.out = 6), beta = seq(0.1, 20.0, length.out = 6))
  )))
  expect_equal(
    # Call the dbetabin function in this package
    apply(X = testGrid, FUN = function(curRow) {
      dbetabin(curRow[4], curRow[1], curRow[2], curRow[3], 1)
    }, MARGIN = 1),
    # Compare the results with the output of the dbbinom function of the 'extraDistr' package
    dbbinom(testGrid[, "K"], testGrid[, "size"], testGrid[, "alpha"], testGrid[, "beta"], TRUE)
  )
})
## 2. ------ TEST THE BETA-BINOMIAL SIMULATION FUNCTION ------
test_that("beta-binomial simulation function gives samples with the correct distributional properties", {
  # Create a matrix of test values to use in the unit test
  testGrid <- as.matrix(expand.grid(alpha = seq(0.1, 20.0, length.out = 6), beta = seq(0.1, 20.0, length.out = 6), size = 1:5 * 2))
  # Simulate 1000 draws from the beta-binomial distribution for each parameter combination and test it against the expected
  # probabilities from a beta-binomial
  chiSqP <- apply(X = testGrid, FUN = function(curRow) {
    # Generate the random draws
    sampleSize <- 1000
    testSamples <- as.integer(replicate(sampleSize, rbetabin(1, curRow[1], curRow[2], curRow[3])))
    # Count the number of draws for each output value
    testCounts <- sapply(X = 0:curRow[3], FUN = function(curVal, testSamples) {
      sum(ifelse(curVal == testSamples, 1, 0))
    }, testSamples = testSamples)
    # Calculate the probability of each output value according to the extraDistr implementation
    calcProbs <- dbbinom(0:curRow[3], curRow[3], curRow[1], curRow[2], FALSE)
    # Calculate the chi-squared probabilities of observing the drawn values
    chiSqP <- 1.0
    if(length(testCounts) > 1) {
      chiSqP <- chisq.test(testCounts, p = calcProbs)$p.value
    }
    chiSqP
  }, MARGIN = 1)
  # By chance roughly 5% of the chi-squared probabilities will be < 0.05. If more than 10% are then it is likely that there is a problems with
  # simulation function
  expect_lt(sum(chiSqP < 0.05) / length(chiSqP), 0.1)
})
