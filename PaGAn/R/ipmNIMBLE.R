### 1.1. ==== Run Integral Projection Model in NIMBLE ====
#' @title Run Integral Projection Model in NIMBLE
#'
#' @description A interface for the specification of vital rate functions using standard
#' model notation familiar to users of the \code{\link[base::glm]{glm}} function.  This
#' function create a model object that can be analysed using kernel evaluation functions
#' such as those provided in \code{\link[ipmKernelEvaluation]{ipmKernelEvaluation}}
#'
#' @param mcmcParams A list containing parameters regulating the Markov chain Monte Carlo
#' algorithm applied in NIMBLE.  This list needs the following named elements: numRuns, numChains,
#' numBurnIn, thinDensity, predictThinDensity
#' @param inputData A \code{data.frame} with data set containing the variables described in the model formula
#' vital rates.  Can be \code{NULL}, in which case input data must be set for each vital rate separately.  If
#' input data is provided for a vital rate than it supersedes the data provided in this parameter.  This
#' parameter is useful if multiple vital rates share a common data set
#' @param numCores An \code{integer} scalar denoting the number of cores to use.  MCMC chains
#' will be distributed across the cores.  A value of \code{NA} or \code{0} results in the number
#' of cores being set equal to that returned by \code{\link[parallel::detectCores]{detectCores}}
#' @Param ... A set of arguments defining the regression relationships and data sets for each of the individual
#' vital rate functions
#'
#' @seealso \code{\link[DHARMa::createDHARMa]{createDHARMa}} \code{\link[nimble::nimbleCode]{nimbleCode}}
#' \code{\link[nimble::buildMCMC]{buildMCMC}} \code{\link[glmNIMBLE]{glmNIMBLE}} \code{\link[errorFamilies]{errorFamilies}}
#' \code{\link[base::glm]{glm}}
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#'
ipmNIMBLE <- function(mcmcParams = list(), inputData = NULL, numCores = 1, ...) {
  # Sanity check the MCMC parameters
  inMCMCParams <- sanityCheckMCMCParameters(mcmcParams)
  # Sanity check the vital rate specifications
  vitalRateSpecs <- tryCatch(list(...), error = function(err) {
    stop("error encountered during processing of vital rate model specification: ", err)
  })
  if(length(vitalRateSpecs) <= 0) {
    stop("error encountered during processing of vital rate model specification: no vital rates specified")
  } else if(!all(grepl(".", names(vitalRateSpecs), fixed = TRUE))) {
    stop("error encountered during processing of vital rate model specification: arguments do not conform to correct naming convention")
  }
  # Retrieve the names of all the vital rates
  vitalRateNames <- unique(sapply(X = strsplit(names(vitalRateSpecs), ".", fixed = TRUE), FUN = function(curEls) {paste(curEls[1:(length(curEls) - 1)], collapse = ".")}))
  # Run GLMs for each of the separate vital rates
  glmModelOutputs <- setNames(lapply(X = vitalRateNames, FUN = function(curVitalRate, vitalArgs, inputData, inMCMCParams, numCores) {
    message("****** Processing regression model for vital rate ", curVitalRate, " ******")
    # Process the arguments to pass to the vital rate regression model
    if(!(paste(curVitalRate, "modelFormula", sep = ".") %in% names(vitalArgs))) {
      stop("error encountered during processing of vital rate model specification: no model formula given for vital rate ", curVitalRate)
    }
    curInputData <- inputData
    if(paste(curVitalRate, "inputData", sep = ".") %in% names(vitalArgs)) {
      curInputData <- vitalArgs[[paste(curVitalRate, "inputData", sep = ".")]]
    }
    if(is.null(curInputData)) {
      stop("error encountered during processing of vital rate model specification: no input data given for vital rate ", curVitalRate)
    }
    curErrorFamily <- gaussian
    if(paste(curVitalRate, "errorFamily", sep = ".") %in% names(vitalArgs)) {
      curErrorFamily <- vitalArgs[[paste(curVitalRate, "errorFamily", sep = ".")]]
    }
    curRegCoeffs <- "none"
    if(paste(curVitalRate, "regCoeffs", sep = ".") %in% names(vitalArgs)) {
      curRegCoeffs <- vitalArgs[[paste(curVitalRate, "regCoeffs", sep = ".")]]
    }
    # Create a node definition object for the linear model specification
    glmNIMBLE(vitalArgs[[paste(curVitalRate, "modelFormula", sep = ".")]], curInputData, curErrorFamily, curRegCoeffs, curVitalRate, inMCMCParams, numCores)
  }, vitalArgs = vitalRateSpecs, inputData = inputData, inMCMCParams = inMCMCParams, numCores = numCores), vitalRateNames)
  # Return the GLM model outputs
  glmModelOutputs
}

ipmKernelEvaluation <- function(ipmOb, inputData, mappingValues, kernelFunction, ...) {
  # Import the IPM object
  inIPMOb <- tryCatch(as.list(ipmOb), error = function(err) {
    stop("error encountered during processing of IPM object: ", err)
  })
  # Retrieve the names of the vital rate models
  vitalRateNames <- names(inIPMOb)
  if(is.null(vitalRateNames)) {
    stop("error encountered during processing of IPM object: vital rate names not found")
  }
  if(length(inIPMOb) <= 0) {
    stop("error encountered during processing of IPM object: no vital rate models present")
  }
  inKernelFunction <- tryCatch(as.function(kernelFunction), error = function(err) {
    stop("error encountered during processing of kernel function: ", err)
  })
  # Retrieve the arguments of the kernel function
  kernelArgs <- formalArgs(inKernelFunction)
  # Remove those arguments that aren't related to vital rate functions
  kernelArgs <- kernelArgs[sapply(X = kernelArgs, FUN = function(curArg, vitalRateNames) {
    any(sapply(X = vitalRateNames, FUN = function(curVital, curArg) {
      grepl(paste("^", curVital, sep = ""), curArg, perl = TRUE)
    }, curArg = curArg))
  }, vitalRateNames = vitalRateNames)]
  # Retrieve the extra parameters to the function
  extraArgs <- list(...)
  # Retrieve the relevant information from the model to populate the kernel arguments
  inputArgs <- setNames(lapply(X = kernelArgs, FUN = function(kernelArg, inIPMOb, inputData) {
    # Split the argument name
    splitArgName <- strsplit(kernelArg, ".", fixed = TRUE)[[1]]
    splitArgName <- c(splitArgName[1], paste(splitArgName[2:length(splitArgName)], collapse = "."))
    # Retrieve the relevant vital rate
    curVitalRate <- inIPMOb[[splitArgName[1]]]
    if(is.null(curVitalRate)) {
      # Check to ensure that the parameter argument name is valid
      stop("error encountered during processing of kernel function: argument names before the period does not match a vital rate model in the IPM object")
    }
    # Retrieve the parameter samples for the relevant vital rate
    curParamSamples <- do.call(rbind, lapply(X = curVitalRate$parameterSamples, FUN = as.data.frame))
    # Function to retrieve the distributional parameters of the error function
    errorRetrieve <- function(errorType, curParamSamples) {
      outValues <- NULL
      tempValues <- NULL
      # The parameter is implemented in the model as precision so transform it to the appropriate error type
      if(any(grepl("Prec$", names(curParamSamples), perl = TRUE))) {
        tempValues <- curParamSamples[, grepl("Prec$", names(curParamSamples), perl = TRUE)]
        outValues <- switch(errorType,
          prec = tempValues,
          var = 1.0 / tempValues,
          sd = sqrt(1.0 / tempValues),
          NULL)
      # The parameter is implemented in the model as stand deviation so transform it to the appropriate error type
      } else if(any(grepl("SD$", names(curParamSamples), perl = TRUE))) {
        tempValues <- curParamSamples[, grepl("SD$", names(curParamSamples), perl = TRUE)]
        outValues <- switch(errorType,
          prec = 1.0 / (tempValues * tempValues),
          var = tempValues * tempValues,
          sd = tempValues,
          NULL)
      }
      outValues
    }
    # Function to find a parameter in a set of parameter samples
    findParameter <- function(paramName, curParamSamples) {
      outValues <- NULL
      colFind <- grepl(paste(paramName, "$", sep = ""), names(curParamSamples), perl = TRUE)
      if(any(colFind)) {
        outValues <- curParamSamples[, colFind]
      }
      outValues
    }
    # Retrieve the relevant values depending on the nature of the parameter being requested
    outValue <- switch(splitArgName[2],
      response = predict.glmNIMBLE(curVitalRate, inputData, "response")$predictionMatrix,
      linpred = predict.glmNIMBLE(curVitalRate, inputData, "link")$predictionMatrix,
      sd = errorRetrieve("sd", curParamSamples),
      var = errorRetrieve("var", curParamSamples),
      prec = errorRetrieve("prec", curParamSamples),
      scale = findParameter("Scale", curParamSamples),
      NULL
    )
    if(is.null(outValue)) {
      stop("error encountered during processing of kernel function: unable to find parameter corresponding to argument name in the IPM object")
    }
    outValue
  }, inIPMOb = inIPMOb, inputData = inputData), kernelArgs)
  # Iterate over each row of the input data and call the kernel functions
  outMat <- sapply(X = 1:nrow(inputData), FUN = function(curRowIndex, mappingValues, kernelFunction, inputArgs, extraArgs) {
    # Retrieve the arguments relevant to the current row of the input matrix
    curInputArgs <- setNames(lapply(X = inputArgs, FUN = function(curArg, curRowIndex) {
      outVal <- NULL
      if(is.matrix(curArg)) {
        # The argument is a matrix: slice a row
        outVal <- curArg[curRowIndex, ]
      } else {
        # The argument is a vector: return the entire vector
        outVal <- curArg
      }
      outVal
    }, curRowIndex = curRowIndex), names(inputArgs))
    # Create a kernel matrix by calling the kernel function with each mapping value and input arguments
    # retrieved from the IPM object
    kernelMatrix <- t(sapply(X = mappingValues, FUN = function(curMapVal, curInputArgs, kernelFunction, extraArgs) {
      do.call(kernelFunction, c(list(x = curMapVal), curInputArgs, extraArgs))
    }, curInputArgs = curInputArgs, kernelFunction = kernelFunction, extraArgs = extraArgs))
    dimnames(kernelMatrix) <- list(paste("mappedValue", 1:length(mappingValues), sep = ""), paste("sample", 1:ncol(kernelMatrix), sep = ""))
    kernelMatrix
  }, mappingValues = mappingValues, kernelFunction = kernelFunction, inputArgs = inputArgs, extraArgs = extraArgs, simplify = "array")
  dimnames(outMat)[[3]] <- paste("dataRow", 1:(dim(outMat)[3]), sep = "")
  # Reorder the matrix so that samples are on the third dimension
  outMat <- aperm(outMat, c(1, 3, 2))
  # Return the matrices and calculated the expected value
  list(
    samples = outMat,
    expected = apply(X = outMat, FUN = mean, MARGIN = c(1, 2))
  )
}

ipmKernelAnalysis <- function(fullKernel) {
  inKernel <- fullKernel
  if(is.list(inKernel)) {
    inKernel <- inKernel$samples
    if(is.null(inKernel)) {
      stop("input kernel is a list with no \"samples\" element")
    }
  }
  if(is.null(dim(inKernel)) || length(dim(inKernel)) != 3) {
    stop("input kernel does not have a correct dimension structure")
  }
  # Calculate the assymptotic population growth from the real part of the dominant eigenvalue
  assymPopGrowth <- apply(X = inKernel, FUN = function(curMat) {
    Re(eigen(curMat)$values[1])
  }, MARGIN = 3)
  # Estimate the assymptotic stable distribution
  stableDist <- apply(X = inKernel, FUN = function(curMat) {
    outEigen <- Re(eigen(curMat)$vectors[, 1])
    outEigen / sum(outEigen)
  }, MARGIN = 3)
  # Estimate the assymptotic repoductive value
  reproVal <- apply(X = inKernel, FUN = function(curMat) {
    outEigen <- Re(eigen(t(curMat))$vectors[, 1])
    outEigen / outEigen[1]
  }, MARGIN = 3)
  # Create a set of summary statistics
  createSummaryStats <- function(inVec) {
    makeStats <- function(inRow) {
      setNames(
        c(mean(inRow), sd(inRow), median(inRow), quantile(inRow, c(0.025, 0.975))),
        c("mean", "sd", "median", "lower95cred", "upper95cred")
      )
    }
    outVal <- NA
    if(is.null(dim(inVec)) || length(dim(inVec)) != 2) {
      outVal <- makeStats(inVec)
    } else {
      outVal <- t(apply(X = inVec, FUN = makeStats, MARGIN = 1))
    }
    outVal
  }
  list(
    assymPopGrowth = assymPopGrowth,
    stableDist = stableDist,
    reproVal = reproVal,
    assymPopGrowthSummary = createSummaryStats(assymPopGrowth),
    stableDistSummary = createSummaryStats(stableDist),
    reproValSummary = createSummaryStats(reproVal)
  )
}
