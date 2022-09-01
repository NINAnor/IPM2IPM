### 1.1. ==== Run Hierarchical PAR Calibration Model ====
#' @title Run Hierarchical PAR Calibration Model
#'
#' @description Function to fit a calibration model of gross primary productivity (GPP) against
#' photosynthically active radiation (PAR).  The calibration model is hierarchical in nature
#' and the parameters of the calibration curve (the x-assymptote, y-assymptote, and general
#' multiplier) can be linear functions of various covariates.
#'
#' @param inputData A \code{data.frame} containing the data to be used in the calibration model
#' @param lightValues A character scalar of the column name in \code{data} that represents
#' the values of PAR.  Alternatively \code{lightValues} can be a numeric vector with length
#' equal to the number of rows in \code{data} which gives the PAR values directly
#' @param yAsymModel A \code{formula} object giving the relationship between the y-assymptote
#' and a set of covariates.  The left-hand term of the formula should be the term representing
#' GPP.  If the left-hand term is missing then GPP is searched for in the other model
#' formulas \code{xAsymModel} and \code{multiplierModel}.  This parameter can also a be a
#' numeric \code{vector}, in which case the model will be replaced by an 'offset' constant with
#' values given in the vector
#' @param xAsymModel A \code{formula} object giving the relationship between the x-assymptote
#' and a set of covariates.  The left-hand term of the formula should be the term representing
#' GPP.  If the left-hand term is missing then GPP is searched for in the other model
#' formula \code{multiplierModel}.  This parameter can also a be a
#' numeric \code{vector}, in which case the model will be replaced by an 'offset' constant with
#' values given in the vector
#' @param multiplierModel A \code{formula} object giving the relationship between the multiplier
#' and a set of covariates.  The left-hand term of the formula should be the term representing
#' GPP.  This parameter can also a be a numeric \code{vector}, in which case the model will be
#' replaced by an 'offset' constant with values given in the vector.  The default value for this
#' parameter is \code{1} which effectively removes this component from the modelling framework
#' @param indirectComponents A \code{list} object containing a series of model specifications
#' representing each of the indirect model components between the covariates.  Each element in
#' the list should itself be a list containing the following named elements: modelForm, family,
#' priorSpecification, numTrials, and suffix.  These elements are described in more detail
#' below
#' @param regCoeffs A character scalar containing the text \code{"none"},
#' \code{"ridge"}, or \code{"lasso"}.  This determines the type of regularisation to use for the
#' regression coefficients in the model.  \code{"none"} (the default) results in wide normal priors
#' for the regression coefficients.  \code{"ridge"} results in normal priors for each of the
#' regression coefficients regulated by a shared precision parameter which itself has a wide
#' gamma hyperprior.  \code{"lasso"} results in Laplace priors for each of the regression
#' coefficients regulated by a shared rate parameter which itself has a wide gamma hyperprior
#' @param lightStandards A numeric vector containing values of PAR for which predictions are
#' going to be made.  Can be \code{NULL} if no standardisation is required
#' @param mcmcParams A list containing parameters regulating the Markov chain Monte Carlo
#' algorithm applied in NIMBLE.  This list needs the following named elements: numRuns, numChains,
#' numBurnIn, thinDensity, predictThinDensity
#'
#' @details Each element of \code{indirectComponents} is a list containing the following elements:
#' \itemize{
#'   \item{modelFormula}{A \code{formula} object describing the relationships between covariates}
#'   \item{errorFamily}{A \code{family} object descrinig the error family to use to model the indirect
#'   relationships between the covariates and the link function to apply}
#'   \item{regCoeffs}{A character scalar used to determine the prior specification.  This
#'   is applied in the same way as the \code{regCoeffs} input parameter}
#'   \item{modelSuffix}{A suffix to give to nodes in this sub-model.  This can be used to avoid name clashes
#'   in complex models with many indirect components}
#'   \item{mcmcParams}{A list containing the MCMC options}
#' }
#' See \code{\link{glmNIMBLE}} for more information on how these parameters are interpreted.
#'
#' \code{mcmcParams} is a list containing the following elements:
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
#'   of \code{\link[nimble::nimbleCode]{nimbleCode}}}
#'   \item{compiledModel}{The compiled NIMBLE model object as would be produced by running
#'   \code{\link[nimble::nimbleModel]{nimbleModel}} followed by \code{\link[nimble::compileNimble]{compileNimble}}}
#'   \item{compiledMCMC}{The compiled NIMBLE MCMC object as would be produced by running
#'   \code{\link[nimble::buildMCMC]{buildMCMC}} followed by \code{\link[nimble::compileNimble]{compileNimble}}}
#'   \item{parameterSamples}{A \code{\link[coda]{mcmc.list}} object containing the samples from
#'   the MCMC analysis}
#'   \item{parameterSummary}{A \code{data.frame} containing summary statistics for each of the
#'   sampled parameters}
#'   \item{predictionSamples}{A \code{\link[coda]{mcmc.list}} object containing the samples of the
#'   mean predictions from the MCMC analysis}
#'   \item{predictionSummary}{A \code{data.frame} containing summary statistics for the mean
#'   predictions}
#'   \item{WAIC}{A scalar containing the Watanabe-Akaike information criterion for the model}
#'   \item{DHARMaResiduals}{A \copde{list} of objects as created by \code{\link[DHARMa::createDHARMa]{createDHARMa}}
#'   that contains an analysis of residuals of each of the model sub-components.  The first element is the DHARMa
#'   analysis for the overall GPP residuals.  Each element afterwards is a DHARMa analysis for each of the indirect
#'   models.}
#'   \item{parameterFigure}{A graphical object (\code{\link[ggplot2::ggplot]{ggplot}}) containing violin plots
#'   for each of the parameters for the submodels (each in a seperate panel)}
#'   \item{standardSummary}{A \code{list} with a length equal to the number of light standards in the
#'   \code{lightStandards} parameter.  Each element is a \code{data.frame} including the summary information
#'   of the MCMC samples.  This element is absent if \code{lightStandards} is \code{NULL}}
#'   \item{indirectModelOutputs}{A \code{list} containing the outputs from \code{\link{glmNIMBLE}} for each of the
#'   indirect component models specified in \code{indirectComponents}}
#'   \item{rSquared}{A \code{list} object containing the following elements: \code{samples}, a \link[coda]{mcmc.list}
#'   object containing the sampled r squared values and a set of summary statistics for these samples across all chains.
#'   The r squared is calculated according to [Gelman et al. 2019](https://doi.org/10.1080/00031305.2018.1549100)}
#' }
#'
#' @seealso \code{\link[DHARMa::createDHARMa]{createDHARMa}} \code{\link[nimble::nimbleCode]{nimbleCode}}
#' \code{\link[nimble::buildMCMC]{buildMCMC}} \code{\link[glmNIMBLE]{glmNIMBLE}}
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#'
gppLightCurveCorrection <- function(
  inputData,
  lightValues,
  yAsymModel,
  xAsymModel,
  multiplierModel = 1.0,
  indirectComponents = NULL,
  regCoeffs = "none",
  lightStandards = c(),
  mcmcParams = list()
) {
  ## 1.1.1. Sanity test the inputs ----
  inData <- tryCatch(as.data.frame(inputData), error = function(err) {
    stop("error encountered during processing of the input data frame: ", err)
  })
  # Retrieve the light values
  inLightValues <- lightValues
  if(is.character(inLightValues)) {
    # If the input is a character vector then use the first element to lookup the light values from the input data frame
    inLightValues <- tryCatch(as.double(inputData[, inLightValues[1]]), error = function(err) {
      stop("error encountered during processing of the light values: ", err)
    })
  } else {
    # Otherwise use the given value directly
    inLightValues <- tryCatch(as.double(inLightValues), error = function(err) {
      stop("error encountered during processing of the light values: ", err)
    })
  }
  # Recycle the light values so they are the same length as the number of rows in the input data
  inLightValues <- inLightValues[((1:nrow(inputData) - 1) %% length(inLightValues)) + 1]
  if(any(inLightValues <= 0.0) || anyNA(inLightValues)) {
    stop("error encountered during processing of the light values: only positive and non-zero values for the light values are valid")
  }
  # Sanity check the MCMC parameters
  inMCMCParameters <- sanityCheckMCMCParameters(mcmcParams)
  # Retrieve the sub-model formulas as strings
  yAsymModelText <- ""
  yAsymModelFormNew <- NULL
  xAsymModelText <- ""
  xAsymModelFormNew <- NULL
  multiplierModelText <- ""
  multiplierModelFormNew <- NULL
  # Function to tidy long formulas
  tidyLongFormula <- function(inForm) {
    # Remove unneccessary space
    gsub("\\s+", " ", paste(inForm, collapse = ""), perl = TRUE)
  }
  if(!is.numeric(yAsymModel)) {
    yAsymModelText <- tryCatch(tidyLongFormula(deparse(substitute(yAsymModel))), error = function(err) {
      stop("error encountered during processing of y-assymptote sub-model: ", err)
    })
    yAsymModelFormNew <- tryCatch(as.formula(paste("~ ", gsub("^.*~\\s*", "", yAsymModelText, perl = TRUE), sep = "")), error = function(err) {
      stop("error encountered during processing of y-assymptote sub-model: ", err)
    })
  }
  if(!is.numeric(xAsymModel)) {
    xAsymModelText <- tryCatch(tidyLongFormula(deparse(substitute(xAsymModel))), error = function(err) {
      stop("error encountered during processing of x-assymptote sub-model: ", err)
    })
    xAsymModelFormNew <- tryCatch(as.formula(paste("~ ", gsub("^.*~\\s*", "", xAsymModelText, perl = TRUE), sep = "")), error = function(err) {
      stop("error encountered during processing of x-assymptote sub-model: ", err)
    })
  }
  if(!is.numeric(multiplierModel)) {
    multiplierModelText <- tryCatch(tidyLongFormula(deparse(substitute(multiplierModel))), error = function(err) {
      stop("error encountered during processing of the multiplier sub-model: ", err)
    })
    multiplierModelFormNew <- tryCatch(as.formula(paste("~ ", gsub("^.*~\\s*", "", multiplierModelText, perl = TRUE), sep = "")), error = function(err) {
      stop("error encountered during processing of the multiplier sub-model: ", err)
    })
  }
  # Retrieve the GPP variable from the model specifications
  gppVarName <- gsub("^\\s*", "", gsub("\\s*~.*$", "", c(yAsymModelText, xAsymModelText, multiplierModelText), perl = TRUE), perl = TRUE)
  gppVarName <- unique(gppVarName[gppVarName != ""])
  if(length(gppVarName) <= 0) {
    stop("error encountered during processing the model specification for the sub-models: no response variable given in any formula")
  } else if(length(gppVarName) > 1) {
    stop("error encountered during processing the model specification for the sub-models: more than one response variable given in different model formulas")
  }
  gppVarValues <- tryCatch(inputData[, gppVarName], error = function(err) {
    stop("error encountered retrieving response variable values: ", err)
  })
  # Retrieve the light standards to make predictions for
  inLightStandards <- c()
  if(!is.null(lightStandards)) {
    inLightStandards <- tryCatch(as.double(lightStandards), error = function(err) {
      stop("error encountered during processing of the light standards: ", err)
    })
    if(any(inLightStandards <= 0.0) || anyNA(inLightStandards)) {
      stop("error encountered during processing of the light standards: only positive and non-zero values for the light values are valid")
    }
  }
  # Process the indirect modlling elements
  inIndirectComponents <- list()
  if(!is.null(indirectComponents)) {
    inIndirectComponents <- tryCatch(lapply(X = as.list(indirectComponents), FUN = function(curElement, inputData, inMCMCParameters) {
      # Initialise a set of default model specification settings
      retrievedModelSpec <- list(
        modelFormula = NULL,
        inputData = inputData,
        errorFamily = gaussian(),
        regCoeffs = "none",
        modelSuffix = "",
        mcmcParams = inMCMCParameters
      )
      # Import the current list of parameters for the indirect model specification
      inList <- as.list(curElement)
      if(is.null(names(inList))) {
        # If there are no names for the indirect element then fill the parameters based on ordering of parameters
        retrievedModelSpec[1:min(length(retrievedModelSpec), length(inList))] <- inList[1:min(length(retrievedModelSpec), length(inList))]
      } else {
        # Otherwise use the names in the list to fill the elements
        if(!is.null(inList[["modelFormula"]])) {
          retrievedModelSpec$modelFormula <- as.formula(inList[["modelFormula"]])
        }
        if(!is.null(inList[["inputData"]])) {
          retrievedModelSpec$inputData <- as.data.frame(inList[["inputData"]])
        }
        if(!is.null(inList[["errorFamily"]])) {
          retrievedModelSpec$errorFamily <- inList[["errorFamily"]]
        }
        if(!is.null(inList[["regCoeffs"]])) {
          retrievedModelSpec$regCoeffs <- as.character(inList[["regCoeffs"]])
        }
        if(!is.null(inList[["modelSuffix"]])) {
          retrievedModelSpec$modelSuffix <- as.character(inList[["modelSuffix"]])
        }
        if(!is.null(inList[["mcmcParams"]])) {
          retrievedModelSpec$mcmcParams <- as.list(inList[["mcmcParams"]])
        }
      }
      # Retun the model specification in the list
      retrievedModelSpec
    }, inputData = inputData, inMCMCParameters = inMCMCParameters), error = function(err) {
      stop("error encountered during processing of indirect model components: ", err)
    })
  }
  ## 1.1.2. Create model structue ----
  # Function to create a model structure for constant components
  modelConstComponents <- function(numValues, numDataPoint, modelSuffix) {
    list(
      inputData = setNames(list(
        numValues[(1:numDataPoint - 1) %% length(numValues) + 1]
      ), paste("meanPred", modelSuffix, sep = "")),
      inputConstants = list(),
      stochasticNodeDef = list(),
      deterministicNodeDef = list()
    )
  }
  # Create the covariate nodes for each of the sub-components of the light correction model
  yAsymModelNodes <- NULL
  xAsymModelNodes <- NULL
  multiplierModelNodes <- NULL
  if(!is.null(yAsymModelFormNew)) {
    yAsymModelNodes <- linearModelToCovariateNodeDefinition(yAsymModelFormNew, inData, "identity", regCoeffs, "logyAssym")
  } else {
    yAsymModelNodes <- modelConstComponents(log(as.numeric(yAsymModel)), nrow(inData), "logyAssym")
  }
  if(!is.null(xAsymModelFormNew)) {
    xAsymModelNodes <- linearModelToCovariateNodeDefinition(xAsymModelFormNew, inData, "identity", regCoeffs, "logxAssym")
  } else {
    xAsymModelNodes <- modelConstComponents(log(as.numeric(xAsymModel)), nrow(inData), "logxAssym")
  }
  if(!is.null(multiplierModelFormNew)) {
    multiplierModelNodes <- linearModelToCovariateNodeDefinition(multiplierModelFormNew, inData, "identity", regCoeffs, "logmultiplier")
  } else {
    multiplierModelNodes <- modelConstComponents(log(as.numeric(multiplierModel)), nrow(inData), "logmultiplier")
  }
  # Assign the extra nodes to support the model
  extraNodes <- list(
    inputData = setNames(list(
      inLightValues, gppVarValues
    ), c("parValues", "gppValues")),
    inputConstants = setNames(list(
      nrow(inData)
    ), c("NumSamples")),
    stochasticNodeDef = setNames(list(
      structure(
        "dgamma(0.05, 0.005)",
        loopIter = NA, loopMax = NA, vectorSpec = NA),
      structure(
        "dgamma(mean = gppMeanPred[gppIter], sd = gppSD)",
        loopIter = "gppIter", loopMax = "NumSamples", vectorSpec = NA)
    ), c("gppSD", "gppValues")),
    deterministicNodeDef = setNames(list(
      structure(
        "meanPredlogyAssym[1:NumSamples] + log(parValues[1:NumSamples]) + meanPredlogmultiplier[1:NumSamples] - meanPredlogxAssym[1:NumSamples] - log(parValues[1:NumSamples] + exp(meanPredlogxAssym[1:NumSamples] + meanPredlogyAssym[1:NumSamples]))",
        loopIter = NA, loopMax = NA, linkFunction = "log", vectorSpec = "1:NumSamples")
    ), c("gppMeanPred"))
  )
  # Add in extra nodes to account for the GPP standards
  if(length(inLightStandards) > 0) {
    extraNodes$inputData <- c(extraNodes$inputData, setNames(list(
      inLightStandards
    ), c("gppStandards")))
    if(length(inLightStandards) > 1) {
      extraNodes$inputConstants <- c(extraNodes$inputConstants, setNames(list(
        length(inLightStandards)
      ), c("NumStandards")))
      extraNodes$deterministicNodeDef <- c(extraNodes$deterministicNodeDef, setNames(list(
        structure(
          "meanPredlogyAssym[1:NumSamples] + log(gppStandards[standardIter]) + meanPredlogmultiplier[1:NumSamples] - meanPredlogxAssym[1:NumSamples] - log(gppStandards[standardIter] + exp(meanPredlogxAssym[1:NumSamples] + meanPredlogyAssym[1:NumSamples]))",
          loopIter = "standardIter", loopMax = "NumStandards", linkFunction = "log", vectorSpec = "1:NumSamples")
      ), c("standardPred")))
    } else {
      extraNodes$deterministicNodeDef <- c(extraNodes$deterministicNodeDef, setNames(list(
        structure(
          "meanPredlogyAssym[1:NumSamples] + log(gppStandards) + meanPredlogmultiplier[1:NumSamples] - meanPredlogxAssym[1:NumSamples] - log(gppStandards + exp(meanPredlogxAssym[1:NumSamples] + meanPredlogyAssym[1:NumSamples]))",
          loopIter = NA, loopMax = NA, linkFunction = "log", vectorSpec = "1:NumSamples")
      ), c("standardPred")))
    }
  }
  # Aggregate all the nodes to create an entire model
  modelNodeDefinitions <- list(
    inputData = c(yAsymModelNodes$inputData, xAsymModelNodes$inputData, multiplierModelNodes$inputData, extraNodes$inputData),
    inputConstants = c(yAsymModelNodes$inputConstants, xAsymModelNodes$inputConstants, multiplierModelNodes$inputConstants, extraNodes$inputConstants),
    stochasticNodeDef = c(yAsymModelNodes$stochasticNodeDef, xAsymModelNodes$stochasticNodeDef, multiplierModelNodes$stochasticNodeDef, extraNodes$stochasticNodeDef),
    deterministicNodeDef = c(yAsymModelNodes$deterministicNodeDef, xAsymModelNodes$deterministicNodeDef, multiplierModelNodes$deterministicNodeDef, extraNodes$deterministicNodeDef)
  )
  ## 1.1.3. Run the model ----
  # Convert the model node definition to NIMBLE code
  modelCode <- nodeDefinitionToNIMBLECode(modelNodeDefinitions$stochasticNodeDef, modelNodeDefinitions$deterministicNodeDef)
  message("---- Running light curve model ----")
  # Create a set of initial values for the stochastic non-data nodes
  nonDataNodeNames <- names(modelNodeDefinitions$stochasticNodeDef)[
    !(names(modelNodeDefinitions$stochasticNodeDef) %in% names(modelNodeDefinitions$inputData))]
  predictionNodeNames <- "gppMeanPred"
  if(length(inLightStandards) > 0) {
    predictionNodeNames <- c(predictionNodeNames, "standardPred")
  }
  initValues <- setNames(lapply(X = nonDataNodeNames, FUN = function(curName) {
    ifelse(grepl("Coeff$", curName, perl = TRUE), 0.0, 1.0)
  }), nonDataNodeNames)
  # Define the model object
  uncompiledModel <- nimbleModel(modelCode, constants = modelNodeDefinitions$inputConstants, data = modelNodeDefinitions$inputData,
    inits = initValues, calculate = TRUE)
  # Compile the model object
  compiledModel <- compileNimble(uncompiledModel)
  #Create an MCMC object
  uncompiledMCMC <- buildMCMC(uncompiledModel, enableWAIC = TRUE,
    monitors = nonDataNodeNames, monitors2 = c(predictionNodeNames, "gppSD"),
    thin = inMCMCParameters$thinDensity, thin2 = inMCMCParameters$predictThinDensity)
  # Compile the MCMC object
  compiledMCMC <- compileNimble(uncompiledMCMC, project = uncompiledModel)
  # Run the MCMC
  mcmcOutput <- runMCMC(compiledMCMC, niter = inMCMCParameters$numRuns, nburnin = inMCMCParameters$numBurnIn, nchains = inMCMCParameters$numChains,
    thin = inMCMCParameters$thinDensity, thin2 = inMCMCParameters$predictThinDensity, samplesAsCodaMCMC = TRUE, WAIC = TRUE, summary = TRUE)
  ## 1.1.4. Model post-processing ----
  # Simulate responses from the posterior predictive distribution
  simulatedValues <- apply(X = as.matrix(do.call(rbind, mcmcOutput$samples2))[, c(paste("gppMeanPred[", 1:nrow(inputData), "]", sep = ""), "gppSD")], FUN = function(curRow) {
    # Retrieve the mean and standard deviation sampled in the MCMC
    predMeans <- curRow[1:(length(curRow) - 1)]
    predSD <- curRow[length(curRow)]
    rgamma(length(predMeans), shape = (predMeans * predMeans) / (predSD * predSD), scale = (predSD * predSD) / predMeans)
  }, MARGIN = 1)
  # Calculate the fitted median response
  fittedPred <- mcmcOutput$summary$all.chains[paste("gppMeanPred[", 1:nrow(inputData), "]", sep = ""), "Median"]
  # Create a DHARMa object so that model checking can be done
  residAnalysisOb <- createDHARMa(
    simulatedResponse = simulatedValues,
    observedResponse = modelNodeDefinitions$inputData$gppValues,
    fittedPredictedResponse = fittedPred,
    integerResponse = FALSE
  )
  # Create violin plots of the regression coefficients
  regCoeffNames <- nonDataNodeNames[grepl("Coeff$", nonDataNodeNames, perl = TRUE) & !grepl("^intercept", nonDataNodeNames, perl = TRUE)]
  regFrame <- do.call(rbind, lapply(X = regCoeffNames, FUN = function(curName, mcmcSamples) {
    data.frame(
      covName = rep(gsub("(logxAssym|logyAssym|logmultiplier)Coeff$", "", curName, perl = TRUE), nrow(mcmcSamples)),
      submodel = rep(ifelse(grepl("logxAssymCoeff$", curName, perl = TRUE), "X-asymptote model", ifelse(
        grepl("logyAssymCoeff$", curName, perl = TRUE), "Y-asymptote model", "Multiplier model")), nrow(mcmcSamples)),
      coeffVal = mcmcSamples[, curName]
    )
  }, mcmcSamples = do.call(rbind, mcmcOutput$samples)))
  regFrame$covName <- as.factor(regFrame$covName)
  regFrame$submodel <- as.factor(regFrame$submodel)
  coeffPlot <- ggplot(regFrame, aes(covName, coeffVal)) + geom_violin(draw_quantiles = c(0.025, 0.5, 0.975)) + coord_flip() +
    geom_hline(yintercept = 0.0, colour = "grey") + theme_classic() + ylab("Scaled coefficient value") +
    theme(axis.title.y = element_blank(), axis.ticks.y = element_blank()) + facet_grid(. ~ submodel)
  # Organise all the outputs into a list
  outList <- list(
    modelDefinition = modelCode,
    compiledModel = compiledModel,
    compiledMCMC = compiledMCMC,
    parameterSamples = mcmcOutput$samples,
    parameterSummary = mcmcOutput$summary$all.chains[nonDataNodeNames, ],
    predictionSamples = mcmcOutput$samples2,
    predictionSummary = mcmcOutput$summary$all.chains[paste("gppMeanPred[", 1:nrow(inputData), "]", sep = ""), ],
    WAIC = mcmcOutput$WAIC,
    DHARMaResiduals = residAnalysisOb,
    parameterFigure = coeffPlot,
    rSquared = bayesRSquared(mcmcOutput$samples2, modelNodeDefinitions$inputData$gppValues)
  )
  # Rearrange the standardised gpp predictions
  if(length(inLightStandards) > 0) {
    standardSummary <- NULL
    if(length(inLightStandards) > 1) {
      standardSummary <- setNames(lapply(X = 1:length(inLightStandards), FUN = function(curIndex, summaryTable, numRows) {
        summaryTable[paste("standardPred[", 1:numRows, ", ", curIndex, "]", sep = ""), ]
      }, summaryTable = mcmcOutput$summary$all.chains, numRows = nrow(inputData)),
        paste("lightStandard", inLightStandards, sep = ""))
    } else {
      standardSummary <- setNames(
        list(mcmcOutput$summary$all.chains[paste("standardPred[", 1:nrow(inputData), "]", sep = ""), ]),
        paste("lightStandard", inLightStandards, sep = ""))
    }
    outList <- c(outList, list(standardSummary = standardSummary))
  }
  # Run the indirect models
  if(length(inIndirectComponents) > 0) {
    indirectModelOutputs <- lapply(X = inIndirectComponents, FUN = function(curParams) {
      curForm <- as.formula(curParams$modelFormula)
      # Retrieve a name for the current indirect pathway
      modelName <- gsub("\\s*~.*$", "", tidyLongFormula(deparse(substitute(curForm))), perl = TRUE)
      message("---- Running ", modelName, " indirect pathway model ----")
      glmNIMBLE(curParams$modelFormula, curParams$inputData, curParams$errorFamily, curParams$regCoeffs, curParams$modelSuffix, curParams$mcmcParams)
    })
    outList <- c(outList, list(indirectModelOutputs = indirectModelOutputs))
  }
  outList
}
