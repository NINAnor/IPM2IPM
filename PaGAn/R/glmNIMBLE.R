### 1.1. ==== Run Generalised Linear Model in NIMBLE ====
#' @title Run Generalised Linear Model in NIMBLE
#'
#' @description A simple interface to run a generalised linear model using the NIMBLE
#' implementation of MCMC for Bayesian analysis using standard model notation familiar
#' to users of the \code{\link[base::glm]{glm}} function.
#'
#' @param modelFormula A \code{formula} project defining the relationship between the
#' response variable and the predictor variables in a similar manner to that used in
#' the \code{\link[base::glm]{glm}} function
#' @param inputData A \code{data.frame} with dataset containing the variables described in
#' the model formula
#' @param errorFamily A \code{family} object defining the error distribution of the
#' regression model and a possible use of a link function.  This parameter can also be
#' a list with two elements: link and family. Each of these are a \code{character} scalar
#' assigning the link function and error family respectively
#' \code{\link[errorFamilies]{errorFamilies()}} will return a list of error distributions
#' supported by \code{glmNIMBLE} and the potential link functions associated with them.
#' @param regCoeffs A character scalar containing the text \code{"none"},
#' \code{"ridge"}, or \code{"lasso"}.  This determines the type of regularisation to use for the
#' regression coefficients in the model.  \code{"none"} (the default) results in wide normal priors
#' for the regression coefficients.  \code{"ridge"} results in normal priors for each of the
#' regression coefficients regulated by a shared precision parameter which itself has a wide
#' gamma hyperprior.  \code{"lasso"} results in Laplace priors for each of the regression
#' coefficients regulated by a shared rate parameter which itself has a wide gamma hyperprior
#' @param modelSuffix A \code{character} scalar to be used as a suffix to give to nodes in this
#' model.  Defaults to emptry string (\code{""})
#' @param mcmcParams A list containing parameters regulating the Markov chain Monte Carlo
#' algorithm applied in NIMBLE.  This list needs the following named elements: numRuns, numChains,
#' numBurnIn, thinDensity, predictThinDensity
#' @param numCores An \code{integer} scalar denoting the number of cores to use.  MCMC chains
#' will be distributed across the cores.  A value of \code{NA} or \code{0} results in the number
#' of cores being set equal to that returned by \code{\link[parallel::detectCores]{detectCores}}
#'
#' @details \code{mcmcParams} is a list containing the following elements:
#' \itemize{
#'   \item{numSamples}{The number of iterations to sample from the posterior distribution in each chain}
#'   \item{numChains}{The number of chains to simulate in the MCMC algorithm}
#'   \item{numBurnIn}{The number of samples to remove at the beginning of each chain}
#'   \item{thinDensity}{The thinning parameter to use for the model paramters}
#'   \item{predictThinDensity}{The thinning parameter to use for the model predictions}
#' }
#' If any of the elements are missing from the \code{mcmcParams} object then sensible defaults are supplied.
#'
#' @return A \code{list} object containing the following elements:
#' \itemize{
#'   \item{modelDefinition}{The NIMBLE code used in the model definition as would be the output
#'   of \code{\link[nimble]{nimbleCode}}}
#'   \item{compiledModel}{The compiled NIMBLE model object as would be produced by running
#'   \code{\link[nimble]{nimbleModel}} followed by \code{\link[nimble]{compileNimble}}}
#'   \item{compiledMCMC}{The compiled NIMBLE MCMC object as would be produced by running
#'   \code{\link[nimble]{buildMCMC}} followed by \code{\link[nimble]{compileNimble}}}
#'   \item{parameterSamples}{A \code{\link[coda]{mcmc.list}} object containing the samples from
#'   the MCMC analysis}
#'   \item{parameterSummary}{A \code{data.frame} containing summary statistics for each of the
#'   sampled parameters}
#'   \item{predictionSamples}{A \code{\link[coda]{mcmc.list}} object containing the samples of the
#'   mean predictions from the MCMC analysis}
#'   \item{predictionSummary}{A \code{data.frame} containing summary statistics for the mean
#'   predictions}
#'   \item{WAIC}{A scalar containing the Watanabe-Akaike information criterion for the model}
#'   \item{DHARMaResiduals}{A \code{list} of objects as created by \code{\link[DHARMa]{createDHARMa}} that contains
#'   an analysis of residuals of each of the model sub-components.  The first element is the DHARMa analysis for the
#'   overall GPP residuals.  Each element afterwards is a DHARMa analysis for each of the indirect models.}
#'   \item{parameterFigure}{A graphical object (\code{\link[ggplot2::ggplot]{ggplot}}) containing violin plots
#'   for each of the parameters}
#'   \item{rSquared}{A \code{list} object containing the following elements: \code{samples}, a \link[coda]{mcmc.list}
#'   object containing the sampled r squared values and a set of summary statistics for these samples across all chains.
#'   The r squared is calculated according to [Gelman et al. 2019](https://doi.org/10.1080/00031305.2018.1549100)}
#' }
#'
#' @seealso \code{\link[DHARMa::createDHARMa]{createDHARMa}} \code{\link[nimble::nimbleCode]{nimbleCode}}
#' \code{\link[nimble::buildMCMC]{buildMCMC}} \code{\link[glmNIMBLE]{glmNIMBLE}} \code{\link[errorFamilies]{errorFamilies}}
#' \code{\link[base::glm]{glm}}
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#'
glmNIMBLE <- function(modelFormula, inputData, errorFamily = gaussian, regCoeffs = "none", modelSuffix = "", mcmcParams = list(), numCores = 1)  {
  inMCMCParams <- sanityCheckMCMCParameters(mcmcParams)
  # Create the linear model node specifications
  modelNodeDefinitions <- linearModelToNodeDefinition(modelFormula, inputData, errorFamily, regCoeffs, modelSuffix)
  # Create the model code to run in NIMBLE
  modelCode <- nodeDefinitionToNIMBLECode(modelNodeDefinitions$stochasticNodeDef, modelNodeDefinitions$deterministicNodeDef)
  # Retrieve the node names
  nonDataNodeNames <- names(modelNodeDefinitions$stochasticNodeDef)[!(names(modelNodeDefinitions$stochasticNodeDef) %in% names(modelNodeDefinitions$inputData))]
  predictionNodeNames <- names(modelNodeDefinitions$deterministicNodeDef)[grepl("^meanPred", names(modelNodeDefinitions$deterministicNodeDef), perl = TRUE)]
  # Create a set of initial values for the stochastic non-data nodes
  initValues <- setNames(lapply(X = nonDataNodeNames, FUN = function(curName) {
    ifelse(grepl("Coeff$", curName, perl = TRUE), 0.0, 1.0)
  }), nonDataNodeNames)
  # Multiplier for the distributions that have a number of trials
  numTrials <- 1.0
  if(as.character(modelNodeDefinitions$family) == "binomial") {
    numTrials <- modelNodeDefinitions$inputData[[which(grepl("NumTrials$", names(modelNodeDefinitions$inputData), perl = TRUE))[1]]]
  }
  # Set the intercept so that it is the mean value of the response variable
  if(any(paste("intercept", as.character(modelSuffix), "Coeff", sep = "") == nonDataNodeNames)) {
    initValues[[paste("intercept", as.character(modelSuffix), "Coeff", sep = "")]] <- applyLink(
      mean(modelNodeDefinitions$inputData[[modelNodeDefinitions$responseDataNodeNames[1]]] / numTrials, na.rm = TRUE), modelNodeDefinitions$link)
  }
  # Compile and run the model
  modelOutputs <- mcmcNIMBLERun(modelCode, modelNodeDefinitions$inputData, modelNodeDefinitions$inputConstants, nonDataNodeNames, predictionNodeNames, initValues, inMCMCParams, numCores, TRUE)
  uncompiledModel <- modelOutputs$uncompiledModel
  compiledModel <- modelOutputs$compiledModel
  uncompiledMCMC <- modelOutputs$uncompiledMCMC
  compiledMCMC <- modelOutputs$compiledMCMC
  mcmcOutput <- modelOutputs$mcmcOutput
  # Define the model object
  #uncompiledModel <- nimbleModel(modelCode, constants = modelNodeDefinitions$inputConstants, data = modelNodeDefinitions$inputData,
  #  inits = initValues, calculate = TRUE)
  # Compile the model object
  #compiledModel <- compileNimble(uncompiledModel)
  #Create an MCMC object
  #uncompiledMCMC <- buildMCMC(uncompiledModel, enableWAIC = TRUE,
  #  monitors = nonDataNodeNames, monitors2 = predictionNodeNames,
  #  thin = inMCMCParams$thinDensity, thin2 = inMCMCParams$predictThinDensity)
  # Compile the MCMC object
  #compiledMCMC <- compileNimble(uncompiledMCMC, project = uncompiledModel)
  # Run the MCMC
  #mcmcOutput <- runMCMC(compiledMCMC, niter = inMCMCParams$numRuns, nburnin = inMCMCParams$numBurnIn, nchains = inMCMCParams$numChains,
  #  thin = inMCMCParams$thinDensity, thin2 = inMCMCParams$predictThinDensity, samplesAsCodaMCMC = TRUE, WAIC = TRUE, summary = TRUE)
  # Retrieve the offset values (if there are any)
  offsetVals <- 0.0
  if(!is.null(modelNodeDefinitions$inputData[[paste("offset", as.character(modelSuffix), sep = "")]])) {
    offsetVals <- as.double(modelNodeDefinitions$inputData[[paste("offset", as.character(modelSuffix), sep = "")]])
  }
  sampledParams <- as.matrix(do.call(rbind, mcmcOutput$samples))
  covMatrix <- as.matrix(as.data.frame(modelNodeDefinitions$inputData[
    gsub("Coeff$", "", nonDataNodeNames[grepl("Coeff$", nonDataNodeNames, perl = TRUE) & !grepl("^intercept", nonDataNodeNames, perl = TRUE)], perl = TRUE)
  ]))
  if(ncol(covMatrix) <= 0 || nrow(covMatrix) <= 0) {
    covMatrix <- matrix(NA, nrow = max(sapply(X = modelNodeDefinitions$inputData, FUN = length)), ncol = 0)
  }
  # Simulate response values according to the sampled parameter values
  simulatedValues <- apply(X = sampledParams, FUN = function(curRow, covMatrix, inFamily, inLink, inData, inOffset) {
    interceptCoeffVal <- 0.0
    if(any(grepl("^intercept.*Coeff$", names(curRow), perl = TRUE))) {
      interceptCoeffVal <- curRow[grepl("^intercept.*Coeff$", names(curRow), perl = TRUE)]
    }
    if(ncol(covMatrix) <= 0) {
      meanPredVals <- rep(0.0, nrow(covMatrix))
    } else {
      meanPredVals <- covMatrix %*% curRow[paste(colnames(covMatrix), "Coeff", sep = "")]
    }
    # Calculate the mean prediction values
    meanPredVals <- applyInverseLink(
      as.double(meanPredVals + interceptCoeffVal) + inOffset,
      inLink)
    # Retrieve the relevant scale parameter for the error distribution
    scaleValues <- switch(inFamily,
      gaussian = curRow[grepl("Prec$", names(curRow), perl = TRUE)],
      gamma = curRow[grepl("SD$", names(curRow), perl = TRUE)],
      beta = curRow[grepl("Scale$", names(curRow), perl = TRUE)],
      poisson = c(),
      binomial = inData[grepl("NumTrials$", names(inData), perl = TRUE)][[1]],
      negbinomial = curRow[grepl("Scale$", names(curRow), perl = TRUE)])
    # Sample from the relevant distribution
    simulateFromErrorFamily(meanPredVals, scaleValues, inFamily)
  }, covMatrix = covMatrix, inFamily = as.character(modelNodeDefinitions$family), inLink = as.character(modelNodeDefinitions$link),
    inData = modelNodeDefinitions$inputData, inOffset = offsetVals, MARGIN = 1)
  # Calculate the fitted median response
  fittedPred <- mcmcOutput$summary$all.chains[grepl("^meanPred", rownames(mcmcOutput$summary$all.chains), perl = TRUE), "Median"] * numTrials
  # Create a DHARMa object so that model checking can be done
  residAnalysisOb <- createDHARMa(
      simulatedResponse = simulatedValues,
      observedResponse = modelNodeDefinitions$inputData[[modelNodeDefinitions$responseDataNodeNames]],
      fittedPredictedResponse = fittedPred,
      integerResponse = switch(as.character(modelNodeDefinitions$family),
        gaussian = FALSE, gamma = FALSE, beta = FALSE,
        poisson = TRUE, binomial = TRUE, negbinomial = TRUE
      )
  )
  # Create violin plots of the regression coefficients
  regCoeffNames <- nonDataNodeNames[grepl("Coeff$", nonDataNodeNames, perl = TRUE) & !grepl("^intercept", nonDataNodeNames, perl = TRUE)]
  coeffPlot <- ggplot(do.call(rbind, lapply(X = regCoeffNames, FUN = function(curName, inData, inSuffix) {
    data.frame(coeffVal = inData[, curName], covName = rep(gsub(paste(inSuffix, "Coeff$", sep = ""), "", curName, perl = TRUE), nrow(inData)))
  }, inData = do.call(rbind, mcmcOutput$samples), inSuffix = as.character(modelSuffix))), aes(covName, coeffVal)) +
    geom_violin(draw_quantiles = c(0.025, 0.5, 0.975)) + coord_flip() +
    geom_hline(yintercept = 0.0, colour = "grey") + theme_classic() + ylab("Scaled coefficient value") +
    theme(axis.title.y = element_blank(), axis.ticks.y = element_blank())
  # Organise all the outputs into a list
  list(
    modelDefinition = modelCode,
    compiledModel = compiledModel,
    compiledMCMC = compiledMCMC,
    parameterSamples = mcmcOutput$samples,
    parameterSummary = mcmcOutput$summary$all.chains[nonDataNodeNames, ],
    predictionSamples = mcmcOutput$samples2,
    predictionSummary = mcmcOutput$summary$all.chains[!(rownames(mcmcOutput$summary$all.chains) %in% nonDataNodeNames), ],
    WAIC = mcmcOutput$WAIC,
    DHARMaResiduals = residAnalysisOb,
    parameterFigure = coeffPlot,
    stochasticNodeDef = modelNodeDefinitions$stochasticNodeDef,
    deterministicNodeDef = modelNodeDefinitions$deterministicNodeDef,
    modelFormula = as.formula(modelFormula),
    covSummaryStats = modelNodeDefinitions$covSummaryStats,
    modelSuffix = processSuffix(modelSuffix),
    rSquared = bayesRSquared(mcmcOutput$samples2, modelNodeDefinitions$inputData[[modelNodeDefinitions$responseDataNodeNames]])
  )
}

bayesRSquared <- function(predictionChains, observedResponse) {
  # Sample variance calculation function
  vzFunc <- function(zIn) {
    sum((zIn - mean(zIn))^2) / (length(zIn) - 1.0)
  }
  # Calculate the r squared values for each of the sampled predictors across the chains
  sampledRSquared <- do.call(coda::mcmc.list, lapply(X = as.list(predictionChains), FUN = function(curChain, observedResponse) {
    # Retrieve the samples from the current chain as a matrix
    curMatrix <- as.matrix(curChain)
    # Calculate the r squared for each sample of the MCMC chain
    rSquaredVals <- apply(X = curMatrix, FUN = function(curPreds, observedResponse) {
      # Remove any elements that have NA values in the observed response
      valsToUse <- !is.na(observedResponse)
      inCurPreds <- curPreds[valsToUse]
      inObserved <- observedResponse[valsToUse]
      # Calculate the variation in the vector of predicted values
      varPred <- vzFunc(inCurPreds)
      # Calculate the modelled residual variance
      varRes <- vzFunc(inObserved - inCurPreds)
      # Use these to calculate r squared
      varPred / (varPred + varRes)
    }, observedResponse = observedResponse, MARGIN = 1)
    coda::mcmc(data = matrix(rSquaredVals, nrow = length(rSquaredVals), ncol = 1, dimnames = list(NULL, "rSquared")),
      start = attr(curChain, "mcpar")[1], end = attr(curChain, "mcpar")[2], thin = attr(curChain, "mcpar")[3])
  }, observedResponse = observedResponse))
  sampledSummary <- summary(sampledRSquared)
  # Calculate summary statistics of the r squared
  summaryStats <- setNames(
    c(sampledSummary$statistics, sampledSummary$quantiles),
    gsub(" ", ".", tolower(c(names(sampledSummary$statistics), paste("quant", as.double(gsub("%", "", names(sampledSummary$quantiles), fixed = TRUE)) / 100.0, sep = ""))), fixed = TRUE)
  )
  # Arrange the summary statistics and the samples into a list
  c(list(samples = sampledRSquared), as.list(summaryStats))
}

predict.glmNIMBLE <- function(modelOb, newData, type = "link") {
  # Import the model object
  inGLM <- tryCatch(as.list(modelOb), error = function(err) {
    stop("error encountered processing model object: object is not correct type")
  })
  # Retrieve the model formula
  inModForm <- tryCatch(as.formula(inGLM$modelFormula), error = function(err) {
    stop("error encountered processing model object: ", err)
  })
  # Retrieve the centre and scaling information
  cenScaleInfo <- tryCatch(as.matrix(inGLM$covSummaryStats), error = function(err) {
    stop("error encountered processing model object: ", err)
  })
  if(!all(c("mean", "sd") %in% rownames(cenScaleInfo))) {
    stop("error encountered processing model object: centre and scale information is not correctly formatted")
  }
  # Import the new data
  inData <- tryCatch(as.data.frame(newData), error = function(err) {
    stop("error encountered processing prediction data frame: ", err)
  })
  # Centre and scale the new data according to the moment of the training data (if required)
  inData <- as.data.frame(setNames(lapply(X = names(inData), FUN = function(curCol, inData, cenScaleInfo) {
    outVec <- inData[, curCol]
    if(curCol %in% colnames(cenScaleInfo)) {
      if(all(!is.na(cenScaleInfo[, curCol]))) {
        outVec <- (inData[, curCol] - cenScaleInfo["mean", curCol]) / cenScaleInfo["sd", curCol]
      }
    }
    outVec
  }, inData = inData, cenScaleInfo = cenScaleInfo), names(inData)))
  # Retrieve the model suffix (if it has one)
  curSuffix <- ""
  if(!is.null(inGLM$modelSuffix)){
    curSuffix <- inGLM$modelSuffix
  }
  # Retrieve the model matrix
  inModForm <- as.formula(gsub("^.*~", "~", Reduce(paste, deparse(inModForm)), perl = TRUE))
  curModelMatrix <- model.matrix(inModForm, model.frame(inModForm, data = inData, na.action = na.pass))
  # Check whether an intercept term is present in the model matrix
  hasIntercept <- 0 %in% attr(curModelMatrix, "assign")
  if(hasIntercept) {
    # Remove the intercept term from the model matrix if it exists
    curModelMatrix <- curModelMatrix[, attr(curModelMatrix, "assign") != 0, drop = FALSE]
  }
  if(ncol(curModelMatrix) > 0) {
    # Rename the model matrix covariates to NIMBLE-friendly names
    colnames(curModelMatrix) <- paste(makeBUGSFriendlyNames(colnames(curModelMatrix)), curSuffix, sep = "")
  }
  # Retrieve the prediction node of the model
  predNode <- inGLM$deterministicNodeDef
  if(is.null(predNode)) {
    stop("error encountered processing model object: deterministic node specification is missing")
  }
  if(!(paste("meanPred", curSuffix, sep = "") %in% names(predNode))) {
    stop("error encountered processing model object: prediction node is missing")
  }
  predNode <- predNode[[paste("meanPred", curSuffix, sep = "")]]
  # Retrieve the prediction specification text
  modSpecText <- tryCatch(gsub(paste("[1:meanPred", curSuffix, "N]", sep = ""), "", as.character(predNode), fixed = TRUE), error = function(err) {
    stop("error encountered processing model object: ", err)
  })
  if(is.null(inGLM$parameterSamples)) {
    stop("error encountered processing model object: unable to retrieve parameter samples")
  }
  # Retrieve the samples of the parameters as a matrix
  paramSamples <- NULL
  if(class(inGLM$parameterSamples) == "mcmc") {
    paramSamples <- as.data.frame(inGLM$parameterSamples)
  } else {
    paramSamples <- tryCatch(do.call(rbind, lapply(X = inGLM$parameterSamples, FUN = as.data.frame)), error = function(err) {
      stop("error encountered processing model object: parameter samples are in an incorrect format")
    })
  }
  # Retrieve the type of prediction to make
  inType <- tryCatch(as.character(type), error = function(err) {
    stop("error encountered processing prediction: invalid prediction type")
  })
  if(is.null(inType) || length(inType) <= 0) {
    inType <- "link"
  } else if(length(inType) > 1) {
    warning("prediction type specification has length greater than zero: only the first element will be used")
  }
  if(!(inType %in% c("link", "response"))) {
    stop("error encountered processing prediction: invalid prediction type")
  }
  # Retrieve the link specification
  linkFunc <- tryCatch(as.character(attr(predNode, "linkFunction")), error = function(err) {
    stop("error encountered processing model object: ", err)
  })
  if(is.null(linkFunc) || length(linkFunc) <= 0) {
    stop("error encountered processing model object: link function for prediction node not specified")
  } else if(length(linkFunc) > 1) {
    warning("link function specification has length greater than zero: only the first element will be used")
    linkFunc <- linkFunc[1]
  }
  if(!(linkFunc %in% unique(unlist(errorFamilies())))) {
    stop("error encountered processing model object: invalid link function detected")
  }
  # Go through each of the rows of the model frame and calculate the predictions for each
  predOut <- t(apply(X = curModelMatrix, MARGIN = 1, FUN = function(curRow, paramSamples, modSpecText, inType, linkFunc) {
    # Evaluate the model using the current parameter and covarite values
    preds <- eval(parse(text = modSpecText), envir = c(paramSamples, as.list(curRow)))
    if(inType == "response") {
      preds <- applyInverseLink(preds, linkFunc)
    }
    preds
  }, paramSamples = as.list(paramSamples), modSpecText = modSpecText, inType = inType, linkFunc = linkFunc))
  # Return the prediction samples alongside the summary statistics
  list(
    predictionMatrix = predOut,
    summary = t(apply(X = predOut, MARGIN = 1, FUN = function(curRow) {
      setNames(
        c(mean(curRow), median(curRow), sd(curRow), quantile(curRow, c(0.025, 0.975))),
        c("Mean", "Median", "St.Dev.", "95%CI_low", "95%CI_upp")
      )
    }))
  )
}
