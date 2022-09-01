library(nimble)
library(DHARMa)
library(ggplot2)

source(paste(Sys.getenv("WORKSPACE_BASELOC"), "Work", "EcosystemResilience", "PaGAn", "R", "mcmcInternals.R", sep = "/"))
source(paste(Sys.getenv("WORKSPACE_BASELOC"), "Work", "EcosystemResilience", "PaGAn", "R", "glmNIMBLE.R", sep = "/"))

## 1. ------ TEST THE GLMNIMBLE MODELLING FUNCTION ------
### 1.1. ==== Gaussian regression test ====
test_that("Gaussian regression model can be adequately parameterised from simulated dataset", {
  # Regression coefficients to simulate data
  trueRegCoeffs <- c(4.0, 1.0, -1.0, 6.0)
  precisionVal <- 1.0
  # Simulate matrix of regression covariates
  numDataPoints <- 5000
  regCovariates <- data.frame(
    A = seq(0.0, 1.0, length.out = numDataPoints),
    B = seq(0.0, 1.0, length.out = 20)[(1:numDataPoints - 1) %% 20 + 1],
    C = factor(c(rep("levelOne", floor(numDataPoints / 2)), rep("levelTwo", numDataPoints - floor(numDataPoints / 2)))),
    testOffset = runif(numDataPoints, -5.0, 5.0)
  )
  # Set the MCMC control parameters
  mcmcControl <- list(
    numRuns = 10000,
    numChains = 4,
    numBurnIn = 5000,
    thinDensity = 1,
    predictThinDensity = 1
  )
  # Create a set of response values
  meanValues <- applyInverseLink(regCovariates$A * trueRegCoeffs[1] + regCovariates$B * trueRegCoeffs[2] + ifelse(as.integer(regCovariates$C) == 1, 0.0, 1.0) * trueRegCoeffs[3] + trueRegCoeffs[4] + regCovariates$testOffset, "identity")
  respValues <- simulateFromErrorFamily(meanValues, precisionVal, "gaussian")
  # Add the response variable to the data frame
  inputData <- cbind(data.frame(respVar = respValues), regCovariates)
  # Run the regression model
  modelOutput <- glmNIMBLE(respVar ~ A + B + C + offset(testOffset), inputData, gaussian, "none", "gaussianTest", mcmcControl)
  # Test the residuals of the model
  residualTests <- testResiduals(modelOutput$DHARMaResiduals)
  # Ensure that the residuals are not significantly different from what is expected
  expect_gt(residualTests$uniformity$p.value, 0.05)
  expect_gt(residualTests$outliers$p.value, 0.05)
})
### 1.2. ==== Gamma regression test ====
test_that("Gamma regression model can be adequately parameterised from simulated dataset", {
  # Regression coefficients to simulate data
  trueRegCoeffs <- c(-0.5, 1.0, -1.0, 4.0)
  scaleVal <- 1.0
  # Simulate matrix of regression covariates
  numDataPoints <- 5000
  regCovariates <- data.frame(
    A = seq(0.0, 1.0, length.out = numDataPoints),
    B = seq(0.0, 1.0, length.out = 20)[(1:numDataPoints - 1) %% 20 + 1],
    C = factor(c(rep("levelOne", floor(numDataPoints / 2)), rep("levelTwo", numDataPoints - floor(numDataPoints / 2)))),
    testOffset = runif(numDataPoints, -2.0, 2.0)
  )
  # Set the MCMC control parameters
  mcmcControl <- list(
    numRuns = 10000,
    numChains = 4,
    numBurnIn = 5000,
    thinDensity = 1,
    predictThinDensity = 1
  )
  # Create a set of response values
  meanValues <- applyInverseLink(regCovariates$A * trueRegCoeffs[1] + regCovariates$B * trueRegCoeffs[2] + ifelse(as.integer(regCovariates$C) == 1, 0.0, 1.0) * trueRegCoeffs[3] + trueRegCoeffs[4] + regCovariates$testOffset, "log")
  respValues <- simulateFromErrorFamily(meanValues, scaleVal, "gamma")
  # Add the response variable to the data frame
  inputData <- cbind(data.frame(respVar = respValues), regCovariates)
  # Run the regression model
  modelOutput <- glmNIMBLE(respVar ~ A + B + C + offset(testOffset), inputData, list(family = "gamma", link = "log"), "none", "gammaTest", mcmcControl)
  # Test the residuals of the model
  residualTests <- testResiduals(modelOutput$DHARMaResiduals)
  # Ensure that the residuals are not significantly different from what is expected
  expect_gt(residualTests$uniformity$p.value, 0.05)
  expect_gt(residualTests$outliers$p.value, 0.05)
})
### 1.3. ==== Beta regression test ====
test_that("Beta regression model can be adequately parameterised from simulated dataset", {
  # Regression coefficients to simulate data
  trueRegCoeffs <- c(1.0, -1.0, -1.5, 0.0)
  scaleVal <- 1.0
  # Simulate matrix of regression covariates
  numDataPoints <- 5000
  regCovariates <- data.frame(
    A = seq(0.0, 1.0, length.out = numDataPoints),
    B = seq(0.0, 1.0, length.out = 20)[(1:numDataPoints - 1) %% 20 + 1],
    C = factor(c(rep("levelOne", floor(numDataPoints / 2)), rep("levelTwo", numDataPoints - floor(numDataPoints / 2)))),
    testOffset = runif(numDataPoints, -0.5, 0.5)
  )
  # Set the MCMC control parameters
  mcmcControl <- list(
    numRuns = 10000,
    numChains = 4,
    numBurnIn = 5000,
    thinDensity = 1,
    predictThinDensity = 1
  )
  # Create a set of response values
  meanValues <- applyInverseLink(regCovariates$A * trueRegCoeffs[1] + regCovariates$B * trueRegCoeffs[2] + ifelse(as.integer(regCovariates$C) == 1, 0.0, 1.0) * trueRegCoeffs[3] + trueRegCoeffs[4] + regCovariates$testOffset, "logit")
  respValues <- simulateFromErrorFamily(meanValues, scaleVal, "beta")
  # Add the response variable to the data frame
  inputData <- cbind(data.frame(respVar = respValues), regCovariates)
  # Run the regression model
  modelOutput <- glmNIMBLE(respVar ~ A + B + C + offset(testOffset) - 1, inputData, list(family = "beta", link = "logit"), "none", "betaTest", mcmcControl)
  # Test the residuals of the model
  residualTests <- testResiduals(modelOutput$DHARMaResiduals)
  # Ensure that the residuals are not significantly different from what is expected
  expect_gt(residualTests$uniformity$p.value, 0.05)
  expect_gt(residualTests$outliers$p.value, 0.05)
})
### 1.4. ==== Poisson regression test ====
test_that("Poisson regression model can be adequately parameterised from simulated dataset", {
  # Regression coefficients to simulate data
  trueRegCoeffs <- c(1.0, -2.0, 3.5, 5.0)
  scaleVal <- NULL
  # Simulate matrix of regression covariates
  numDataPoints <- 5000
  regCovariates <- data.frame(
    A = seq(0.0, 1.0, length.out = numDataPoints),
    B = seq(0.0, 1.0, length.out = 20)[(1:numDataPoints - 1) %% 20 + 1],
    C = factor(c(rep("levelOne", floor(numDataPoints / 2)), rep("levelTwo", numDataPoints - floor(numDataPoints / 2)))),
    testOffset = runif(numDataPoints, -0.5, 0.5)
  )
  # Set the MCMC control parameters
  mcmcControl <- list(
    numRuns = 10000,
    numChains = 4,
    numBurnIn = 5000,
    thinDensity = 1,
    predictThinDensity = 1
  )
  # Create a set of response values
  meanValues <- applyInverseLink(regCovariates$A * trueRegCoeffs[1] + regCovariates$B * trueRegCoeffs[2] + ifelse(as.integer(regCovariates$C) == 1, 0.0, 1.0) * trueRegCoeffs[3] + trueRegCoeffs[4] + regCovariates$testOffset, "log")
  respValues <- simulateFromErrorFamily(meanValues, scaleVal, "poisson")
  # Add the response variable to the data frame
  inputData <- cbind(data.frame(respVar = respValues), regCovariates)
  # Run the regression model
  modelOutput <- glmNIMBLE(respVar ~ A + B + C + offset(testOffset), inputData, list(family = "poisson", link = "log"), "none", "poissonTest", mcmcControl)
  # Test the residuals of the model
  residualTests <- testResiduals(modelOutput$DHARMaResiduals)
  # Ensure that the residuals are not significantly different from what is expected
  expect_gt(residualTests$uniformity$p.value, 0.05)
  expect_gt(residualTests$outliers$p.value, 0.05)
})
### 1.5. ==== Binomial regression test ====
test_that("Binomial model can be adequately parameterised from simulated dataset", {
  # Regression coefficients to simulate data
  trueRegCoeffs <- c(4.0, -2.0, -2.0, -0.25)
  numTrials <- 20
  # Simulate matrix of regression covariates
  numDataPoints <- 5000
  regCovariates <- data.frame(
    A = seq(0.0, 1.0, length.out = numDataPoints),
    ASquared = seq(0.0, 1.0, length.out = numDataPoints) * seq(0.0, 1.0, length.out = numDataPoints),
    C = factor(c(rep("levelOne", floor(numDataPoints / 2)), rep("levelTwo", numDataPoints - floor(numDataPoints / 2)))),
    testOffset = runif(numDataPoints, -0.5, 0.5)
  )
  # Set the MCMC control parameters
  mcmcControl <- list(
    numRuns = 10000,
    numChains = 4,
    numBurnIn = 5000,
    thinDensity = 1,
    predictThinDensity = 1
  )
  # Create a set of response values
  meanValues <- applyInverseLink(regCovariates$A * trueRegCoeffs[1] + regCovariates$ASquared * trueRegCoeffs[2] + ifelse(as.integer(regCovariates$C) == 1, 0.0, 1.0) * trueRegCoeffs[3] + trueRegCoeffs[4] + regCovariates$testOffset, "cloglog")
  respValues <- simulateFromErrorFamily(meanValues, numTrials, "binomial")
  # Add the response variable to the data frame
  inputData <- cbind(data.frame(respVar = respValues), regCovariates)
  # Run the regression model
  modelOutput <- glmNIMBLE(cbind(respValues, numTrials - respValues) ~ A + I(A^2) + C + offset(testOffset), inputData, list(family = "binomial", link = "cloglog"), "none", "binomialTest", mcmcControl)
  # Test the residuals of the model
  residualTests <- testResiduals(modelOutput$DHARMaResiduals)
  # Ensure that the residuals are not significantly different from what is expected
  expect_gt(residualTests$uniformity$p.value, 0.05)
  expect_gt(residualTests$outliers$p.value, 0.05)
})
### 1.6. ==== Negative-binomial regression test ====
test_that("Negative-binomial model can be adequately parameterised from simulated dataset", {
  # Regression coefficients to simulate data
  trueRegCoeffs <- c(1.0, -2.0, 3.5, 5.0)
  scaleVal <- 1.0
  # Simulate matrix of regression covariates
  numDataPoints <- 5000
  regCovariates <- data.frame(
    A = seq(0.0, 1.0, length.out = numDataPoints),
    B = seq(0.0, 1.0, length.out = 20)[(1:numDataPoints - 1) %% 20 + 1],
    C = factor(c(rep("levelOne", floor(numDataPoints / 2)), rep("levelTwo", numDataPoints - floor(numDataPoints / 2)))),
    testOffset = runif(numDataPoints, -0.5, 0.5)
  )
  # Set the MCMC control parameters
  mcmcControl <- list(
    numRuns = 10000,
    numChains = 4,
    numBurnIn = 5000,
    thinDensity = 1,
    predictThinDensity = 1
  )
  # Create a set of response values
  meanValues <- applyInverseLink(regCovariates$A * trueRegCoeffs[1] + regCovariates$B * trueRegCoeffs[2] + ifelse(as.integer(regCovariates$C) == 1, 0.0, 1.0) * trueRegCoeffs[3] + trueRegCoeffs[4] + regCovariates$testOffset, "log")
  respValues <- simulateFromErrorFamily(meanValues, scaleVal, "negbinomial")
  # Add the response variable to the data frame
  inputData <- cbind(data.frame(respVar = respValues), regCovariates)
  # Run the regression model
  modelOutput <- glmNIMBLE(respVar ~ A + B + C + offset(testOffset), inputData, list(family = "negbinomial", link = "log"), "none", "negbinomialTest", mcmcControl)
  # Test the residuals of the model
  residualTests <- testResiduals(modelOutput$DHARMaResiduals)
  # Ensure that the residuals are not significantly different from what is expected
  expect_gt(residualTests$uniformity$p.value, 0.05)
  expect_gt(residualTests$outliers$p.value, 0.05)
})
