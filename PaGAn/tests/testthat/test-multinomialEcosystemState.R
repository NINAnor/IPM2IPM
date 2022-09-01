library(nimble)

test_that("that it is possible to fit the multinomial ecosystem state model to simulated data (with Gaussian error)", {
  # Set the MCMC parameters
  mcmcIters <- 10000
  mcmcBurnin <- 5000
  mcmcChains <- 4
  mcmcThin <- 1
  # The number of data points to simulate
  numDataPoints <- 500
  padZero <- function(inVec, size) {
    outVec <- rep(0.0, size)
    outVec[1:min(size, length(inVec))] <- inVec[1:min(size, length(inVec))]
    outVec
  }
  # Create a data frame with test data covariates
  testFrame <- data.frame(
    # Create a NULL response variable
    respVariable = rep(NA, numDataPoints),
    # Create two scalar covariates
    covA = seq(0.0, 1.0, length.out = numDataPoints),
    covB = padZero(c(seq(0.0, 1.0, length.out = floor(numDataPoints / 2)), seq(1.0, 0.0, length.out = floor(numDataPoints / 2))), numDataPoints),
    # Create a factor covariate
    covC = as.factor(c(rep("levelA", floor(numDataPoints / 3)), rep("levelB", floor(numDataPoints / 3)), rep("levelC", numDataPoints - 2 * floor(numDataPoints / 3))))
  )
  # Set formulae controlling the behaviour of the ecosystem multistate model
  ecoValFormula <- respVariable ~ covB
  ecoPrecFormula <- ~ covC
  ecoProbFormula <- ~ covA
  # testModelSpec <- modelSpecificationMultinomialEcosystemState(ecoValFormula, ecoProbFormula, ecoPrecFormula, testFrame, 3)
  # Set the coefficient values for the simulation
  coefficientValues <- list(
    # Set the coefficient values for the ecosystem state value model
    "intercept_stateVal" = c(-1.0, 0.5, 0.5),
    "covB_stateVal" = c(-0.05, 0.05, 0.1),
    # Set the coefficient values for the ecosystem state probability model
    "intercept_stateProb" = c(NA, -1.0, 1.0),
    "covA_stateProb" = c(NA, 1.0, -1.0),
    # Set the coefficient values for the ecosystem state precision model
    "intercept_statePrec" = c(0.5, 1.0, 1.5),
    "covClevelB_statePrec" = c(0.1, -0.1, 0.5),
    "covClevelC_statePrec" = c(0.4, -0.4, 0.5)
  )
  # Simulate a set of values
  simOutput <- simulateMultinomialEcosystemState(1, ecoValFormula, ecoProbFormula, ecoPrecFormula, testFrame, 3, gaussian, coefficientValues)
  testFrame$respVariable <- simOutput$respVariable[, 1]
  fittedModel <- fitMultinomialEcosystemState(ecoValFormula, ecoProbFormula, ecoPrecFormula, testFrame, 3, gaussian, mcmcIters, mcmcBurnin, mcmcChains, mcmcThin)
})
