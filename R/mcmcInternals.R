# Function that returns the error family distributions and the link functions that are available for each of them
errorFamilies <- function() {
  list(
    gaussian = c("identity", "log"),
    gamma = c("log"),
    beta = c("logit", "probit", "cloglog"),
    poisson = c("log"),
    binomial = c("logit", "probit", "cloglog"),
    negbinomial = c("log")
  )
}

mcmcNIMBLERun <- function(modelCode, data, constants, paramNodeNames, predictionNodeNames, inits, mcmcList = list(), numCores = 1, WAIC = TRUE) {
  # Sanity check the MCMC parameters
  inMCMCList <- sanityCheckMCMCParameters(mcmcList)
  # Sanity check the number of cores
  inNumCores <- tryCatch(as.integer(numCores), error = function(err) {
    stop("error encountered during processing of the number of cores to use: ", err)
  })
  if(length(inNumCores) <= 0) {
    stop("error encountered during processing of the number of cores to use: input vector has length 0")
  } else if(length(inNumCores) > 1) {
    warning("length of vector specifying the number of cores to use is greater than one: only the first element will be used")
    inNumCores <- inNumCores[1]
  }
  if(is.na(inNumCores) || inNumCores <= 0) {
    # If the number of cores is NA or equal to less than zero then just set the number
    # of cores equal to the number present in the system
    inNumCores <- parallel::detectCores()
  }
  # Sanity check the WAIC inclusion criterion
  inWAIC <- tryCatch(as.logical(WAIC), error = function(err) {
    stop("error encountered during processing of the WAIC inclusion term: ", err)
  })
  if(length(inWAIC) <= 0) {
    stop("error encountered during processing of the WAIC inclusion term: input vector has length 0")
  } else if(length(inWAIC) > 1) {
    warning("WAIC inclusion terms has a length greater than one: only the first element will be used")
    inWAIC <- inWAIC[1]
  }
  if(is.na(inWAIC)) {
    inWAIC <- FALSE
  }
  # Set the number of cores equal to the number of chains
  inNumCores <- min(inNumCores, inMCMCList$numChains, parallel::detectCores() - 1)
  # Initialise a set of output objects
  mcmcOutput <- NULL
  uncompiledModel <- NULL
  compiledModel <- NULL
  uncompiledMCMC <- NULL
  compiledMCMC <- NULL
  if(inNumCores <= 1) {
    # Define the model object
    uncompiledModel <- nimbleModel(modelCode, constants = constants, data = data, inits = inits, calculate = TRUE)
    # Compile the model object
    compiledModel <- compileNimble(uncompiledModel)
    # Create an MCMC object
    uncompiledMCMC <- buildMCMC(uncompiledModel, enableWAIC = inWAIC, monitors = paramNodeNames, monitors2 = predictionNodeNames, thin = inMCMCList$thinDensity, thin2 = inMCMCList$predictThinDensity)
    # Compile the MCMC object
    compiledMCMC <- compileNimble(uncompiledMCMC, project = uncompiledModel)
    # Run the MCMC
    mcmcOutput <- runMCMC(compiledMCMC, niter = inMCMCList$numRuns + inMCMCList$numBurnIn, nburnin = inMCMCList$numBurnIn, nchains = inMCMCList$numChains, thin = inMCMCList$thinDensity, thin2 = inMCMCList$predictThinDensity, samplesAsCodaMCMC = TRUE, WAIC = inWAIC, summary = TRUE)
  } else {
    # Create a function to run across multiple cores
    parallelRun <- function(procNum, modelCode, data, constants, paramNodeNames, predictionNodeNames, inits, mcmcList, WAIC, tempFiles, chainVec, seedVec) {
      processCompleteTxt <- "Process complete"
      outObject <- NULL
      if(procNum > 0) {
        cat("Initialising process...\n")
        # Initialise the libraries in the local process
        library(coda)
        library(nimble)
        # Define the model object
        cat("Defining model object...\n")
        uncompiledModel <- nimbleModel(modelCode, constants = constants, data = data, inits = inits, calculate = TRUE)
        # Compile the model object
        cat("Compiling model...\n")
        compiledModel <- compileNimble(uncompiledModel)
        # Create an MCMC object
        cat("Creating MCMC object...\n")
        uncompiledMCMC <- buildMCMC(uncompiledModel, enableWAIC = TRUE, monitors = paramNodeNames, monitors2 = predictionNodeNames, thin = mcmcList$thinDensity, thin2 = mcmcList$predictThinDensity)
        # Compile the MCMC object
        cat("Compiling MCMC object...\n")
        compiledMCMC <- compileNimble(uncompiledMCMC, project = uncompiledModel)
        # Run the MCMC
        cat("Running MCMC...\n")
        mcmcOutput <- runMCMC(compiledMCMC, niter = mcmcList$numRuns + mcmcList$numBurnIn, nburnin = mcmcList$numBurnIn, nchains = chainVec[procNum], thin = mcmcList$thinDensity, thin2 = mcmcList$predictThinDensity, samplesAsCodaMCMC = TRUE, WAIC = FALSE, summary = TRUE, setSeed = seedVec[procNum])
        # Print the process complete text
        cat(processCompleteTxt, "\n", sep = "")
        # Create an output object
        outObject <- list(
          uncompiledModel = uncompiledModel,
          compiledModel = compiledModel,
          uncompiledMCMC = uncompiledMCMC,
          compiledMCMC = compiledMCMC,
          mcmcOutput = mcmcOutput
        )
      } else {
        # The first process is a reporting process that reports back on the current status of the other processes
        waitingText <- "Waiting for processes to initialise ...\n"
        cat(waitingText)
        while(any(!sapply(X = tempFiles, FUN = file.exists))) {
          Sys.sleep(10)
        }
        # Retrieve the text in the temporary files
        retrieveProcessText <- function(tempFiles) {
          sapply(X = tempFiles, FUN = function(curFile) {
            paste(readLines(curFile, warn = FALSE), collapse = "\n")
          })
        }
        # Display the output text
        displayOutputText <- function(processText, outMessage) {
          # Clear the buffer of the last status message
          outMsgChar <- nchar(outMessage)
          if(outMsgChar > 0) {
            cat(rep("\r", outMsgChar), sep = "")
          }
          # Create a new status message
          outMessage <- paste("\n--- PROCESS ", 1:length(chainVec), " (", chainVec, " chain", ifelse(chainVec > 1, "s", ""), ") ---\n\n", processText, "\n", sep = "", collapse = "\n")
          cat(outMessage)
          outMessage
        }
        # Intermittently query the status files and report them
        processText <- retrieveProcessText(tempFiles)
        outMessage <- ""
        while(any(!sapply(X = processText, FUN = function(curText, processCompleteTxt) { grepl(paste(processCompleteTxt, "\\s*$", sep = ""), curText, perl = TRUE) }, processCompleteTxt = processCompleteTxt))) {
          # Display the current status
          outMessage <- displayOutputText(processText, outMessage)
          Sys.sleep(20)
          # Retrieve the text in the temporary files
          processText <- retrieveProcessText(tempFiles)
        }
        outMessage <- displayOutputText(processText, outMessage)
      }
      outObject
    }
    # Calculate a balancer to distribute the chains between cores
    chainVec <- rep(ceiling(inMCMCList$numChains / inNumCores), inNumCores)
    overCount <- sum(chainVec) - inMCMCList$numChains
    if(overCount != 0) {
      chainVec[1:overCount] <- chainVec[1:overCount] + ifelse(overCount > 0, -1, 1)
    }
    # Initialise a cluster
    parCluster <- parallel::makeCluster(inNumCores + 1, outfile = "")
    logFiles <- replicate(inNumCores, tempfile())
    clusterApply(cl = parCluster, x = append(list(stdout()), as.list(logFiles)), fun = function(inVal) {
      sink(inVal, type = "output")
      sink(stdout(), type = "message")
    })
    # Call the model with the chains distributed across the cores
    parallelOutputs <- tryCatch(clusterApply(cl = parCluster, x = 0:inNumCores, fun = parallelRun,
      modelCode = modelCode, data = data, constants = constants, paramNodeNames = paramNodeNames, predictionNodeNames = predictionNodeNames,
      inits = inits, mcmcList = inMCMCList, WAIC = inWAIC, tempFiles = logFiles, chainVec = chainVec, seedVec = floor(runif(inNumCores) * .Machine$integer.max)),
      error = function(err) {
        stopCluster(parCluster)
        stop(err)
    })
    # Stop the cluster after the models are run
    stopCluster(parCluster)
    # Stitch together the outputs from the parallel processes
    uncompiledModel <- lapply(X = parallelOutputs[1:inNumCores + 1], FUN = function(curOb) { curOb$uncompiledModel })
    compiledModel <- lapply(X = parallelOutputs[1:inNumCores + 1], FUN = function(curOb) { curOb$compiledModel })
    uncompiledMCMC <- lapply(X = parallelOutputs[1:inNumCores + 1], FUN = function(curOb) { curOb$uncompiledMCMC })
    compiledMCMC <- lapply(X = parallelOutputs[1:inNumCores + 1], FUN = function(curOb) { curOb$compiledMCMC })
    # Create amalgamated coda objects from the runs spread across the processes
    mcmcOutput <- list(
      samples = mcmc.list(setNames(do.call(c, lapply(X = parallelOutputs[1:inNumCores + 1], FUN = function(curOb) { curOb$mcmcOutput$samples })), paste("chain", 1:inMCMCList$numChains, sep = ""))),
      samples2 = mcmc.list(setNames(do.call(c, lapply(X = parallelOutputs[1:inNumCores + 1], FUN = function(curOb) { curOb$mcmcOutput$samples2 })), paste("chain", 1:inMCMCList$numChains, sep = ""))),
      summary = NULL, WAIC = NULL
    )
    # Recreate the summary information for the samples across the chains
    mcmcOutput$summary <- setNames(c(lapply(X = 1:inMCMCList$numChains, FUN = function(inIndex, samplesList, samplesTwoList) {
      rbind(nimble::samplesSummary(as.matrix(samplesList[[inIndex]])), nimble::samplesSummary(as.matrix(samplesTwoList[[inIndex]])))
    }, samplesList = mcmcOutput$samples, samplesTwoList = mcmcOutput$samples2), list(rbind(
      nimble::samplesSummary(do.call(rbind, lapply(X = mcmcOutput$samples, FUN = as.matrix))),
      nimble::samplesSummary(do.call(rbind, lapply(X = mcmcOutput$samples2, FUN = as.matrix)))
    ))), c(paste("chain", 1:inMCMCList$numChains, sep = ""), "all.chains"))
    # Calculate the WAIC using all the samples spread across the processes
    cat("Recompiling model for WAIC calculation...\n")
    tempModel <- nimbleModel(modelCode, constants = constants, data = data, inits = inits, calculate = TRUE)
    tempCModel <- compileNimble(tempModel)
    tempMCMC <- buildMCMC(tempModel, enableWAIC = TRUE, monitors = paramNodeNames, thin = inMCMCList$thinDensity)
    tempCMCMC <- compileNimble(tempMCMC, project = tempModel)
    mcmcOutput$WAIC <- calculateWAIC(do.call(rbind, lapply(X = mcmcOutput$samples, FUN = as.matrix)), model = tempModel)
  }
  list(
    uncompiledModel = uncompiledModel,
    compiledModel = compiledModel,
    uncompiledMCMC = uncompiledMCMC,
    compiledMCMC = compiledMCMC,
    mcmcOutput = mcmcOutput
  )
}

sanityCheckMCMCParameters <- function(inputList) {
  # Default MCMC sampling parameters
  mcmcList <- list(
    numRuns = 10000,
    numChains = 4,
    numBurnIn = 5000,
    thinDensity = 1,
    predictThinDensity = 1
  )
  # Sanity check the MCMC list
  tempList <- tryCatch(as.list(inputList), error = function(err) {
    stop("error encountered during processing of MCMC control parameter list: ", err)
  })
  mcmcList[names(tempList)[names(tempList) %in% names(mcmcList)]] <- tempList
  # Check to ensure the number of numRuns entry is present and correct
  mcmcList$numRuns <- tryCatch(as.integer(mcmcList$numRuns)[1], error = function(err) {
    stop("error encountered during processing of MCMC control parameter (numRuns): ", err)
  })
  if(is.na(mcmcList$numRuns) || mcmcList$numRuns <= 0) {
    stop("error encountered during processing of MCMC control parameter (numRuns): must be an integer value greater than zero")
  }
  # Check to ensure the number of numChains entry is present and correct
  mcmcList$numChains <- tryCatch(as.integer(mcmcList$numChains)[1], error = function(err) {
    stop("error encountered during processing of MCMC control parameter (numChains): ", err)
  })
  if(is.na(mcmcList$numChains) || mcmcList$numChains <= 0) {
    stop("error encountered during processing of MCMC control parameter (numChains): must be an integer value greater than zero")
  }
  # Check to ensure the number of "burn-in" samples entry is present and correct
  mcmcList$numBurnIn <- tryCatch(as.integer(mcmcList$numBurnIn)[1], error = function(err) {
    stop("error encountered during processing of MCMC control parameter (numBurnIn): ", err)
  })
  if(is.na(mcmcList$numBurnIn) || mcmcList$numBurnIn < 0) {
    stop("error encountered during processing of MCMC control parameter (numBurnIn): must be an integer value greater than zero")
  }
  # Check to ensure the thinning parameter is present and correct
  mcmcList$thinDensity <- tryCatch(as.integer(mcmcList$thinDensity)[1], error = function(err) {
    stop("error encountered during processing of MCMC control parameter (thinDensity): ", err)
  })
  if(is.na(mcmcList$thinDensity) || mcmcList$thinDensity <= 0) {
    stop("error encountered during processing of MCMC control parameter (thinDensity): must be an integer value greater than zero")
  }
  # Check to ensure the other thinning parameter is present and correct
  mcmcList$predictThinDensity <- tryCatch(as.integer(mcmcList$predictThinDensity)[1], error = function(err) {
    stop("error encountered during processing of MCMC control parameter (predictThinDensity): ", err)
  })
  if(is.na(mcmcList$predictThinDensity) || mcmcList$predictThinDensity <= 0) {
    stop("error encountered during processing of MCMC control parameter (predictThinDensity): must be an integer value greater than zero")
  }
  mcmcList
}

# Function to ensure that the variable names are valid for BUGS-style model specification
makeBUGSFriendlyNames <- function(inNames) {
  # Remove non-allowed characters from the names
  outNames <- tryCatch(gsub("\\W", "_", as.character(inNames), perl = TRUE), error = function(err) {
    stop("error encountered during processing of input variable names: ", err)
  })
  # Ensure the names does not start with a numeric character
  outNames <- paste(ifelse(
    grepl("^\\d+", outNames, perl = TRUE), "n", ""
  ), outNames, sep = "")
  outNames
}

processSuffix <- function(modelSuffix = "") {
  inSuffix <- ""
  if(!is.null(modelSuffix)) {
    inSuffix <- tryCatch(makeBUGSFriendlyNames(as.character(modelSuffix)), error = function(err) {
      stop("error encountered during processing of model suffix: ", err)
    })
    if(length(inSuffix) == 0) {
      inSuffix <- ""
      warning("input suffix is length zero (default value will be used)")
    } else if(length(inSuffix) > 1) {
      inSuffix <- inSuffix[1]
      warning("input suffix has length greater than one (only the first value will be used)")
    }
    if(is.na(inSuffix)) {
      inSuffix <- ""
      warning("input suffix is NA (defualt value will be used)")
    }
  }
  inSuffix
}

calculateMeanPred <- function(paramVals, modelMatrix, modelOffset = NULL, linkFunction = "identity", modelSuffix = "") {
  # process the parameter values
  inParamVals <- tryCatch(as.double(paramVals), error = function(err) {
    stop("error processing parameter vector: ", err)
  })
  # Process the model matrix
  inMatrix <- tryCatch(as.matrix(modelMatrix), error = function(err) {
    stop("error processing model matrix: ", err)
  })
  # Process the link function input
  inLink <- tryCatch(factor(tolower(as.character(linkFunction)), unique(unlist(errorFamilies()))), error = function(err) {
    stop("error encountered during import of link function: ", err)
  })
  if(is.na(inLink)) {
    stop("error encountered during import of link function: invalid link function selected")
  }
  # Ensure that the processed model matrix has the appropriate column names (if present)
  if(!is.null(names(modelMatrix))) {
    colnames(inMatrix) <- names(modelMatrix)
  }
  if(!is.null(colnames(modelMatrix))) {
    colnames(inMatrix) <- colnames(modelMatrix)
  }
  # Process the input model suffix
  inSuffix <- processSuffix(modelSuffix)
  # Initialise an output values vector
  outValues <- rep(NA, nrow(modelMatrix))
  if(is.null(names(paramVals)) || is.null(colnames(inMatrix))) {
    # If either the parameter vector or the model matrix don't have names then match up coefficients and the covariates
    # based on their position in the matrix/vector
    outValues <- inMatrix %*% inParamVals[(1:ncol(inMatrix) - 1) %% length(inParamVals) + 1]
    if(length(inParamVals) > ncol(inMatrix)) {
      # Add the intercept term (if it is present) - assuming any parameter in the vector with an index larger than the
      # width of the model matrix is an intercept term
      outValues <- outValues + inParamVals[(ncol(inMatrix) + 1):length(inParamVals)]
    }
  } else {
    # If both have names then use those to match up the coefficients and the covariates
    # Find those parameter value names that have entries in the model matrix
    hasEntry <- names(paramVals) %in% colnames(inMatrix)
    outValues <- inMatrix[, names(paramVals)[hasEntry]] %*% inParamVals[hasEntry]
    if(any(!hasEntry)) {
      # If any parameter values do not have an entry in the model matrix then add them as an intercept term
      outValues <- outValues + inParamVals[!hasEntry]
    }
  }
  # Add the offset if it is not NULL then add that to the prediction
  if(!is.null(modelOffset)) {
    outValues <- tryCatch(as.double(modelOffset) + outValues, error = function(err) {
      stop("error encountered during processing of model offset: ", err)
    })
  }
  # Apply the inverse of the link function to scale the predicted values accordingly
  outValues <- tryCatch(switch(as.character(linkFunction),
    indentity = outValues,
    log = exp(outValues),
    logit = exp(outValues) / (1.0 + exp(outValues)),
    probit = pnorm(outValues),
    cloglog = 1.0 - exp(-exp(outValues))
  ), error = function(err) {
    stop("error encountered during application of link function: ", err)
  })
  if(is.null(outValues)) {
    stop("error encountered during application of link function: invalid link function selected")
  }
  outValues
}

applyLink <- function(dataValues, linkFunction = "identity") {
  # Retrieve the data values
  inValue <- tryCatch(as.numeric(dataValues), error = function(err) {
    stop("error encountered applying link function: ", err)
  })
  # Retrieve the link function
  inLink <- tryCatch(as.character(factor(linkFunction, unique(unlist(errorFamilies())))), error = function(err) {
    stop("error encountered importing link function definition: ", err)
  })
  if(any(is.na(inLink))) {
    stop("error encountered applying link function: invalid link function selected")
  }
  # Apply the link function
  switch(inLink,
    identity = inValue,
    log = log(inValue),
    logit = log(inValue / (1.0 - inValue)),
    probit = qnorm(inValue),
    cloglog = log(-log(1.0 - inValue))
  )
}

applyInverseLink <- function(dataValues, linkFunction = "identity") {
  # Retrieve the data values
  inValue <- tryCatch(as.numeric(dataValues), error = function(err) {
    stop("error encountered applying inverse link function: ", err)
  })
  # Retrieve the link function
  inLink <- tryCatch(as.character(factor(linkFunction, unique(unlist(errorFamilies())))), error = function(err) {
    stop("error encountered importing link function definiton: ", err)
  })
  if(any(is.na(inLink))) {
    stop("error encountered applying inverse link function: invalid link function selected")
  }
  # Apply the inverse link function
  switch(inLink,
    identity = inValue,
    log = exp(inValue),
    logit = exp(inValue) / (1.0 + exp(inValue)),
    probit = pnorm(inValue),
    cloglog = 1.0 - exp(-exp(inValue))
  )
}

simulateFromErrorFamily <- function(meanVals, scaleParams, errorFamily = "gaussian") {
  # Import the mean values
  inMeanVals <- tryCatch(as.numeric(meanVals), error = function(err) {
    stop("error encountered during simulation from error distribution: ", err)
  })
  # Import the scale parameters (these have slightly different meaning depending on the error distribution)
  if(is.null(scaleParams)) {
    inScaleParams <- NA
  } else {
    inScaleParams <- tryCatch(as.numeric(scaleParams), error = function(err) {
      stop("error encountered during simulation from error distribution: ", err)
    })
  }
  # Import the error family
  inErrorFamily <- tryCatch(as.character(factor(errorFamily, names(errorFamilies()))), error = function(err) {
    stop("error encountered during simulation from error distribution: ", err)
  })
  if(any(is.na(inErrorFamily))) {
    stop("error encountered during simulation from error dsitribution: invalid error family selected")
  }
  # Sample from the relevant simulation function
  switch(inErrorFamily,
    gaussian = rnorm(length(inMeanVals), inMeanVals, sqrt(1.0 / inScaleParams[1])),
    gamma = rgamma(length(inMeanVals),
      shape = inMeanVals * inMeanVals / (inScaleParams[1] * inScaleParams[1]),
      scale = inScaleParams[1] * inScaleParams[1] / inMeanVals
    ),
    beta = rbeta(length(inMeanVals),
      shape1 = inMeanVals * inScaleParams[1],
      shape2 = (1.0 - inMeanVals) * inScaleParams[1]
    ),
    poisson = rpois(length(inMeanVals), inMeanVals),
    binomial = rbinom(length(inMeanVals), inScaleParams[(1:length(inMeanVals) - 1) %% length(inScaleParams) + 1], inMeanVals),
    negbinomial = rnbinom(length(inMeanVals), 1.0 / inScaleParams[1], 1.0 - 1.0 / (1.0 + inScaleParams[1] * inMeanVals)))
}

# Function to create the node specification for a linear model component
linearModelToCovariateNodeDefinition <- function(modelFormula, inputData, linkFunction = "identity", regCoeffs = "none", modelSuffix = "") {
  # Function to centre and scale the input data frame
  centreAndScale <- function(inCovDataFrame) {
    outFrame <- as.data.frame(lapply(X = as.list(inCovDataFrame), FUN = function(curCol) {
      outVec <- curCol
      # TODO: add in handling routines for binary covariates here
      if(!is.factor(outVec)) {
        if(is.character(outVec)) {
          outVec <- as.factor(outVec)
        } else {
          outVec <- (outVec - mean(outVec, na.rm = TRUE)) / sd(outVec, na.rm = TRUE)
        }
      }
      outVec
    }))
    names(outFrame) <- names(inCovDataFrame)
    outFrame
  }
  # Process the model formula input
  inFormula <- ~ 1
  if(!is.null(modelFormula)) {
    inFormula <- tryCatch(as.formula(modelFormula), error = function(err) {
      stop("error encountered during processing of model specification: ", err)
    })
  }
  # Process the data input and centre and scale it
  inData <- tryCatch(as.data.frame(inputData), error = function(err) {
    stop("error encountered during import of data: ", err)
  })
  # Process the link function input
  inLink <- tryCatch(factor(tolower(as.character(linkFunction)), unique(unlist(errorFamilies()))), error = function(err) {
    stop("error encountered during import of link function: ", err)
  })
  if(is.na(inLink)) {
    stop("error encountered during import of link function: invalid link function selected")
  }
  # Process the regularisation regime for the coefficients
  inReg <- tryCatch(factor(tolower(as.character(regCoeffs)), c("none", "ridge", "lasso")), error = function(err) {
    stop("error encountered during processing of instructions to regularise coefficients: ", err)
  })
  if(is.na(inReg)) {
    stop("error encountered during processing of instructions to regularise coefficients: invalid regularisation technique selected")
  }
  # Process the model suffix
  inSuffix <- processSuffix(modelSuffix)
  # Retrieve the model matrix
  covMeanVals <- sapply(X = as.list(inData), FUN = function(curCol) {
    outVal <- c(NA, NA)
    # TODO: add in handling routines for binary covariates here
    if(!is.factor(curCol) && !is.character(curCol)) {
      outVal <- c(mean(curCol, na.rm = TRUE), sd(curCol, na.rm = TRUE))
    }
    names(outVal) <- c("mean", "sd")
    outVal
  })
  colnames(covMeanVals) <- colnames(inData)
  curModelMatrix <- model.matrix(inFormula, model.frame(inFormula, data = centreAndScale(inData), na.action = na.pass))
  # Check whether an intercept term is present in the model matrix
  hasIntercept <- 0 %in% attr(curModelMatrix, "assign")
  stocNodeDefinitions <- list()
  if(hasIntercept) {
    # Remove the intercept term from the model matrix if it exists
    curModelMatrix <- curModelMatrix[, attr(curModelMatrix, "assign") != 0, drop = FALSE]
    # Add the intercept term to the list of stochastic node definitions
    stocNodeDefinitions <- append(stocNodeDefinitions, setNames(list(
      structure("dnorm(0.0, 0.001)", loopIter = NA, loopMax = NA, vectorSpec = NA)
    ), paste("intercept", inSuffix, "Coeff", sep = "")))
  }
  outData <- list()
  if(ncol(curModelMatrix) > 0) {
    # Ensure that the model matrix covariates have names
    colnames(curModelMatrix) <- paste(makeBUGSFriendlyNames(colnames(curModelMatrix)), inSuffix, sep = "")
    # Assign names for the regularisation parameters
    ridgeParamName <- paste("ridge", inSuffix, "Prec", sep = "")
    lassoParamName <- paste("lasso", inSuffix, "Rate", sep = "")
    # Create an ouput data list with the covariates
    outData <- as.list(as.data.frame(curModelMatrix))
    # Add the regression coefficients to the list of stochastic node definitions
    stocNodeDefinitions <- append(stocNodeDefinitions, setNames(lapply(X = colnames(curModelMatrix), FUN = function(curCovName, ridgeParamName, lassoParamName, inReg) {
      structure(switch(as.character(inReg),
        none = "dnorm(0.0, 0.001)",
        ridge = paste("dnorm(0.0, ", ridgeParamName, ")", sep = ""),
        lasso = paste("ddexp(0.0, ", lassoParamName, ")", sep = "")),
        loopIter = NA, loopMax = NA, vectorSpec = NA)
    }, ridgeParamName = ridgeParamName, lassoParamName = lassoParamName, inReg = inReg), paste(colnames(curModelMatrix), "Coeff", sep = "")))
    # Add the regularisation stochastic node definitions (if needed)
    if(as.character(inReg) == "ridge") {
      stocNodeDefinitions <- append(stocNodeDefinitions, setNames(list(structure("dgamma(0.05, 0.005)", loopIter = NA, loopMax = NA, vectorSpec = NA)), ridgeParamName))
    } else if(as.character(inReg) == "lasso") {
      stocNodeDefinitions <- append(stocNodeDefinitions, setNames(list(structure("dgamma(0.05, 0.005)", loopIter = NA, loopMax = NA, vectorSpec = NA)), lassoParamName))
    }
  }
  # Test to see if there an offset in the model specification
  offsetFrame <- model.frame(inFormula, data = inData, na.action = na.pass)
  offsetText <- ""
  if(!is.null(attr(terms(offsetFrame), which = "offset")) && any(attr(terms(offsetFrame), which = "offset") > 0)) {
    # Retrieve the offset values if they exist
    offsetValues <- as.numeric(model.offset(offsetFrame))
    # Add them to the model specification as data
    outData <- append(outData, setNames(list(offsetValues), paste("offset", inSuffix, sep = "")))
    offsetText <- paste(" + offset", inSuffix, "[1:meanPred", inSuffix, "N]", sep = "")
  }
  if(ncol(curModelMatrix) > 0) {
    # Add the mean predictor to the list of deterministic nodes
    detNodeDefinitions <- setNames(list(
      structure(
        # Set the linear regression relationship
        paste(
          paste(colnames(curModelMatrix), "[1:meanPred", inSuffix, "N] * ", colnames(curModelMatrix), "Coeff", sep = "", collapse = " + "),
          ifelse(hasIntercept, paste(" + intercept", inSuffix, "Coeff", sep = ""), ""),
          offsetText, sep = ""),
        loopIter = NA, loopMax = NA, linkFunction = as.character(inLink), vectorSpec = paste("1:meanPred", inSuffix, "N", sep = "")
      )
    ), paste("meanPred", inSuffix, sep = ""))
  } else {
    detNodeDefinitions <- setNames(list(
      structure(
        # Set the linear relationship (in this case it is intercept-only)
        paste("intercept", inSuffix, "Coeff", sep = ""),
        loopIter = NA, loopMax = NA, linkFunction = as.character(inLink), vectorSpec = paste("1:meanPred", inSuffix, "N", sep = "")
      )
    ), paste("meanPred", inSuffix, sep = ""))
  }
  outConst <- setNames(list(nrow(curModelMatrix)), paste("meanPred", inSuffix, "N", sep = ""))
  list(
    inputData = outData,
    inputConstants = outConst,
    stochasticNodeDef = stocNodeDefinitions,
    deterministicNodeDef = detNodeDefinitions,
    covSummaryStats = covMeanVals
  )
}

linearModelToNodeDefinition <- function(modelFormula, inputData, errorFamily = gaussian, regCoeffs = "none", modelSuffix = "") {
  # Process the model formula input
  inFormula <- ~ 1
  if(!is.null(modelFormula)) {
    inFormula <- tryCatch(as.formula(modelFormula), error = function(err) {
      stop("error encountered during processing of model specification: ", err)
    })
  }
  # Process the data input
  inData <- tryCatch(as.data.frame(inputData), error = function(err) {
    stop("error encountered during import of data: ", err)
  })
  # Process the error family distribution
  testError <- list(link = "identity", family = "gaussian")
  if(is.function(errorFamily)) {
    testError <- tryCatch(as.list(do.call(errorFamily, list())), error = function(err) {
      stop("error encountered during processing of model error distribution: ", err)
    })
  } else if(is.character(errorFamily)) {
    testError <- list(link = NA, family = errorFamily)
  } else {
    testError <- tryCatch(as.list(errorFamily), error = function(err) {
      stop("error encountered during processing of model error distribution: ", err)
    })
  }
  # Retrieve the error distribution and link function from the input
  inFamily <- factor("gaussian", names(errorFamilies()))
  inLink <- factor("identity", unique(unlist(errorFamilies())))
  if(is.null(testError$family)) {
    stop("error encountered during processing of model error distribution: unknown error distribution selected")
  } else {
    inFamily <- tryCatch(factor(tolower(as.character(testError$family)), names(errorFamilies())), error = function(err) {
      stop("error encountered during processing of model error distribution: ", err)
    })
  }
  if(is.na(inFamily)) {
    stop("error encountered during processing of model error distribution: unknown error distribution selected")
  }
  if(is.null(testError$link) || is.na(testError$link)) {
    # Use the default link function for that distribution if it hasn't been set
    inLink <- factor(errorFamilies()[[inFamily]][1], unique(unlist(errorFamilies())))
  } else {
    inLink <- tryCatch(factor(tolower(as.character(testError$link)), unique(unlist(errorFamilies()))), error = function(err) {
      stop("error encountered during processing of model link function: ", err)
    })
    if(is.na(inLink)) {
      stop("error encountered during processing of model link function: unknown link function selected")
    }
  }
  if(!(as.character(inLink) %in% errorFamilies()[[as.character(inFamily)]])) {
    stop("error encountered during processing of model link function: link function is not compatible with error distribution")
  }
  # Process the model suffix
  inSuffix <- processSuffix(modelSuffix)
  # Retrieve the response variable
  modelFrame <- model.frame(inFormula, data = inData, na.action = na.pass)
  if(any(attr(terms(modelFrame), which = "response") == 0)) {
    stop("error encountered during processing of model input: no response variable provided in model formula")
  }
  # Format the response variable information as a list
  responseData <- model.response(modelFrame)
  if(is.matrix(responseData)) {
    responseData <- as.list(as.data.frame(responseData))
  } else {
    responseData <- setNames(list(responseData), all.vars(terms(modelFrame))[attr(terms(modelFrame), which = "response")])
  }
  # Ensure that the response variables have BUGS-friendly names
  names(responseData) <- paste(makeBUGSFriendlyNames(names(responseData)), inSuffix, sep = "")
  # Create the node definitions for the linear model component
  linearDefs <- linearModelToCovariateNodeDefinition(inFormula, inData, inLink, regCoeffs, inSuffix)
  # Initialise a series of lists to store node definitions and constants related to the modelling of the
  # response terms in the model
  stocRespNodeDefinitions <- list()
  detRespNodeDefinitions <- list()
  # Initialise the names of looping variables
  loopIterName <- paste(names(responseData)[1], "Iter", sep = "")
  loopMaxName <- paste(names(responseData)[1], "N", sep = "")
  responseConsts <- setNames(list(length(responseData[[1]])), loopMaxName)
  if(as.character(inFamily) == "gaussian") {
    # Gaussian error family
    # Set relevant parameter names
    precVarName <- paste(names(responseData)[1], "Prec", sep = "")
    # Set a normal distribution node for the data elements and a prior for the precision parameter
    stocRespNodeDefinitions <- setNames(list(
      structure(
        paste("dnorm(meanPred", inSuffix, "[", loopIterName, "], ", precVarName, ")", sep = ""),
        loopIter = loopIterName, loopMax = loopMaxName, vectorSpec = NA
      ),
      structure(
        "dgamma(0.05, 0.005)",
        loopIter = NA, loopMax = NA, vectorSpec = NA
      )
    ), c(names(responseData)[1], precVarName))
  } else if(as.character(inFamily) == "gamma") {
    # Gamma error family
    # Set relevant parameter names
    sdVarName <- paste(names(responseData)[1], "SD", sep = "")
    # Set a gamma distribution node for the data elements and a prior for the standard deviation parameter
    stocRespNodeDefinitions <- setNames(list(
      structure(
        paste("dgamma(mean = meanPred", inSuffix, "[", loopIterName, "], sd = ", sdVarName, ")", sep = ""),
        loopIter = loopIterName, loopMax = loopMaxName, vectorSpec = NA
      ),
      structure(
        "dgamma(0.05, 0.005)",
        loopIter = NA, loopMax = NA, vectorSpec = NA
      )
    ), c(names(responseData)[1], sdVarName))
  } else if(as.character(inFamily) == "beta") {
    # Beta error family (uses the parameterisation described in Douma and Weedon 2019)
    # Set relevant parameter names
    scaleVarName <- paste(names(responseData)[1], "Scale", sep = "")
    # Set a standard deviation deterministic node to link the scale
    detRespNodeDefinitions <- setNames(list(
      structure(
        paste("sqrt(meanPred", inSuffix, "[1:", loopMaxName, "] * (1.0 - meanPred", inSuffix, "[1:", loopMaxName, "]) / (1.0 + ", scaleVarName, "))", sep = ""),
        loopIter = NA, loopMax = NA, linkFunction = NA, vectorSpec = paste("1:", loopMaxName, sep = "")
      )
    ), paste(names(responseData)[1], "SD", sep = ""))
    # Set a beta distribution node for the data elements and set a prior for the scale parameter
    stocRespNodeDefinitions <- setNames(list(
      structure(
        paste("dbeta(mean = meanPred", inSuffix, "[", loopIterName, "], sd = ", names(responseData)[1], "SD[", loopIterName, "])", sep = ""),
        loopIter = loopIterName, loopMax = loopMaxName, vectorSpec = NA
      ),
      structure(
        "dgamma(0.05, 0.005)",
        loopIter = NA, loopMax = NA, vectorSpec = NA
      )
    ), c(names(responseData)[1], scaleVarName))
  } else if(as.character(inFamily) == "poisson") {
    # Poisson error family
    # Set a poisson distribution node for the data elements
    stocRespNodeDefinitions <- setNames(list(
      structure(
        paste("dpois(meanPred", inSuffix, "[", loopIterName, "])", sep = ""),
        loopIter = loopIterName, loopMax = loopMaxName, vectorSpec = NA
      )
    ), names(responseData)[1])
  } else if(as.character(inFamily) == "binomial") {
    # Binomial error family
    if(length(responseData) != 2) {
      stop("error encountered when defining binomial regression model: response variable is not dimensioned sufficiently to calculate the number of trials")
    }
    # Create a number of trials parameter
    responseData[[2]] <- responseData[[1]] + responseData[[2]]
    names(responseData)[2] <- paste(names(responseData)[1], "NumTrials", sep = "")
    # Set a binomial distribution node for the data elements
    stocRespNodeDefinitions <- setNames(list(
      structure(
        paste("dbin(meanPred", inSuffix, "[", loopIterName, "], ", names(responseData)[2], "[", loopIterName, "])", sep = ""),
        loopIter = loopIterName, loopMax = loopMaxName, vectorSpec = NA
      )
    ), names(responseData)[1])
  } else if(as.character(inFamily) == "negbinomial") {
    # Negative binomial error family (uses the parameterisation described in Ver Hoef and Boveng 2007)
    # Set relevant parameter names
    scaleVarName <- paste(names(responseData)[1], "Scale", sep = "")
    # Create a set of deterministic nodes to convert between a mean and scale parameterisation to the canonical parameterisation
    detRespNodeDefinitions <- setNames(list(
      structure(
        paste("1.0 - 1.0 / (1.0 + ", scaleVarName, " * meanPred", inSuffix, "[1:", loopMaxName, "])", sep = ""),
        loopIter = NA, loopMax = NA, linkFunction = NA, vectorSpec = paste("1:", loopMaxName, sep = "")
      ),
      structure(
        paste("1.0 / ", scaleVarName, sep = ""),
        loopIter = NA, loopMax = NA, linkFunction = NA, vectorSpec = NA
      )
    ), paste(names(responseData)[1], c("P", "R"), sep = ""))
    # Set a negative binomial distribution node for the data elements and set a prior for the scale parameter
    stocRespNodeDefinitions <- setNames(list(
      structure(
        paste("dnegbin(", names(responseData)[1], "P[", loopIterName, "], ", names(responseData)[1], "R)", sep = ""),
        loopIter = loopIterName, loopMax = loopMaxName, vectorSpec = NA
      ),
      structure(
        "dgamma(0.05, 0.005)",
        loopIter = NA, loopMax = NA, vectorSpec = NA
      )
    ), c(names(responseData)[1], scaleVarName))
  }
  list(
    inputData = append(linearDefs$inputData, responseData),
    inputConstants = append(linearDefs$inputConstants, responseConsts),
    stochasticNodeDef = append(linearDefs$stochasticNodeDef, stocRespNodeDefinitions),
    deterministicNodeDef = append(linearDefs$deterministicNodeDef, detRespNodeDefinitions),
    link = inLink,
    family = inFamily,
    responseDataNodeNames = names(responseData)[!grepl("NumTrials$", names(responseData), perl = TRUE)],
    covSummaryStats = linearDefs$covSummaryStats
  )
}

nodeDefinitionToNIMBLECode <- function(stochasticNodeDef, deterministicNodeDef) {
  # Process the user
  inStochastic <- tryCatch(as.list(stochasticNodeDef), error = function(err) {
    stop("error encountered when processing stochastic node definitions: ", err)
  })
  inDeterministic <- tryCatch(as.list(deterministicNodeDef), error = function(err) {
    stop("error encountered when processing deterministic node definitions: ", err)
  })
  processNode <- function(curNodeName, isStochastic, nodeList) {
    # Retrieve the current node
    curNode <- tryCatch(as.list(nodeList)[[curNodeName]], error = function(err) {
      stop("error encountered during processing of node ", curNodeName, ": ", err)
    })
    # Initialise attributes for the node
    loopIter <- c()
    loopMax <- c()
    loopStart <- c()
    vectorSpec <- c()
    linkFunction <- c()
    # Import the attributes of the node
    if(!is.null(attr(curNode, which = "loopIter"))) {
      loopIter <- tryCatch(as.character(attr(curNode, which = "loopIter")), error = function(err) {
        stop("error encountered during processing of node ", curNodeName, ": ", err)
      })
      loopIter <- loopIter[!is.na(loopIter)]
    }
    if(!is.null(attr(curNode, which = "loopMax"))) {
      loopMax <- tryCatch(as.character(attr(curNode, which = "loopMax")), error = function(err) {
        stop("error encountered during processing of node ", curNodeName, ": ", err)
      })
      loopMax <- loopMax[!is.na(loopMax)]
    }
    if(!is.null(attr(curNode, which = "loopStart"))) {
      loopStart <- tryCatch(as.character(attr(curNode, which = "loopStart")), error = function(err) {
        stop("error encountered during processing of node ", curNodeName, ": ", err)
      })
      loopStart <- loopStart[!is.na(loopStart)]
    }
    if(!is.null(attr(curNode, which = "vectorSpec"))) {
      vectorSpec <- tryCatch(as.character(attr(curNode, which = "vectorSpec")), error = function(err) {
        stop("error encountered during processing of node ", curNodeName, ": ", err)
      })
      vectorSpec <- vectorSpec[!is.na(vectorSpec)]
    }
    if(!is.null(attr(curNode, which = "linkFunction"))) {
      linkFunction <- tryCatch(as.character(attr(curNode, which = "linkFunction")), error = function(err) {
        stop("error encountered during processing of node ", curNodeName, ": ", err)
      })
      linkFunction <- linkFunction[!is.na(linkFunction)]
    }
    if(length(linkFunction) > 1) {
      stop("error encountered during processing of node ", curNodeName, ": link function specification must have length 1")
    }
    # Convert to node to a character vector
    inNode <- tryCatch(as.character(curNode), error = function(err) {
      stop("error encountered during processing of node ", curNodeName, ": ", err)
    })
    # Add any indexing to the node that might be required
    preText <- curNodeName
    postText <- curNode
    if(length(vectorSpec) > 0 || length(loopIter) > 0) {
      indexText <- c(paste(vectorSpec, collapse = ", ", sep = ""), paste(loopIter, collapse = ", ", sep = ""))
      preText <- paste(preText, "[", paste(indexText[indexText != ""], collapse = ", ", sep = ""), "]", sep = "")
    }
    # Add any link function to the node that might be required (NIMBLE now allows link functions on stochastic nodes as well as deterministic ones)
    if(length(linkFunction) > 0) {
      if(linkFunction == "log") {
         preText <- paste("log(", preText, ")", sep = "")
      } else if(linkFunction == "logit") {
         preText <- paste("logit(", preText, ")", sep = "")
      } else if(linkFunction == "probit") {
         preText <- paste("probit(", preText, ")", sep = "")
      } else if(linkFunction == "cloglog") {
        preText <- paste("cloglog(", preText, ")", sep = "")
      } else if(linkFunction != "identity") {
        stop("error encountered during processing of node ", curNodeName, ": invalid link function selected")
      }
    }
    # Combine the text together to create the text defining the node
    nodeText <- paste("\t", preText, ifelse(isStochastic, " ~ ", " <- "), postText, sep = "")
    # Add in loop information if the node is defined as part of a loop
    if(length(loopIter) > 0) {
      if(length(loopMax) == 0) {
        stop("error encountered during processing of node ", curNodeName, ": loop limit not specified for indexed node")
      }
      if(length(loopStart) == 0) {
        loopStart <- "1"
      }
      nodeText <- paste(paste(rep("\t", length(loopIter)), collapse = ""), nodeText, sep = "")
      loopText <- as.character(t(sapply(X = 1:length(loopIter), FUN = function(curLoopIndex, loopIter, loopMax, numLoops) {
        c(
          paste(paste(rep("\t", curLoopIndex), collapse = ""), "for(", loopIter[curLoopIndex], " in ", loopStart[(curLoopIndex - 1) %% length(loopStart) + 1], ":", loopMax[(curLoopIndex - 1) %% length(loopMax) + 1], ") {", sep = ""),
          paste(paste(rep("\t", numLoops - curLoopIndex + 1), collapse = ""), "}", sep = "")
        )
      }, loopIter = loopIter, loopMax = loopMax, numLoops = length(loopIter))))
      nodeText <- paste(c(loopText[1:(length(loopText) / 2)], nodeText, loopText[(length(loopText) / 2 + 1):length(loopText)]), collapse = "\n")
    }
    nodeText
  }
  # Create the entire model text
  modelText <- paste(
    "nimbleCode({",
    paste(sapply(X = names(inStochastic), FUN = processNode, isStochastic = TRUE, nodeList = inStochastic), collapse = "\n"),
    paste(sapply(X = names(inDeterministic), FUN = processNode, isStochastic = FALSE, nodeList = inDeterministic), collapse = "\n"),
    "})", sep = "\n")
  # Parse the text to create the model code object
  eval(parse(text = modelText))
}
