## 1. ------ DEFINE INTERNAL FUNCTIONS ------

### 1.1. ==== Change Variable Names to BUGS-Friendly Versions ====
#' @title Change Variable Names to BUGS-Friendly Versions
#'
#' @description Not all potential variable names used in R can be used in BUGS code without
#' producing a syntax error.  This function changes a vector of input names and ensures that
#' they are valid BUGS variable names.
#'
#' @param inName A character vector of variable names
#'
#' @return A character vector of BUGS-compliant variable names
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @keywords internal
#'
setBUGSVariableName <- function(inName) {
  outName <- tryCatch(as.character(inName), error = function(err) {
    stop(paste("invalid parameter name:", err, sep = " "))
  })
  # Remove any non-word characters in the name
  outName <- gsub("\\W+", "_", outName, perl = TRUE)
  # Replace any digit values with a text equivalent (this is to ensure that numbers in variable names aren't parsed incorrectly)
  outName <- gsub("0", "Zero", outName, fixed = TRUE)
  outName <- gsub("1", "One", outName, fixed = TRUE)
  outName <- gsub("2", "Two", outName, fixed = TRUE)
  outName <- gsub("3", "Three", outName, fixed = TRUE)
  outName <- gsub("4", "Four", outName, fixed = TRUE)
  outName <- gsub("5", "Five", outName, fixed = TRUE)
  outName <- gsub("6", "Six", outName, fixed = TRUE)
  outName <- gsub("7", "Seven", outName, fixed = TRUE)
  outName <- gsub("8", "Eight", outName, fixed = TRUE)
  outName <- gsub("9", "Nine", outName, fixed = TRUE)
  outName
}

### 1.2. ==== Create NIMBLE Model Structures for the Multinomial Ecosystem State Model ====
#' @title Create NIMBLE Model Structures for the Multinomial Ecosystem State Model
#'
#' @description This function takes a set of model specifications for the three sub-model
#' components of the multinomial ecosystem state model and generates a set of structures
#' for initialisation of a NIMBLE model
#'
#' @param stateValModels A formula describing the regression relationship between the mean
#' ecosystem state value and the covariates for all ecosystem states, or a list of formulas
#' with each element giving the regression relationship between the mean ecosystem state
#' value and the covariates for each ecosystem state (ordered according to intercept of the
#' ecosystem state value on the y-axis).  The ecosystem state variable must be given on the
#' left-hand side of the formula.
#' @param stateProbModels A formula describing the regression relationship between the probability
#' of the existence of ecosystem state and the covariates for all ecosystem states, or a list
#' of formulas with each element giving the regression relationship between the probability
#' of the existence of ecosystem state and the covariates for each ecosystem state (ordered
#' according to intercept of the ecosystem state value on the y-axis).  The ecosystem state
#' variable must be given on the left-hand side of the formula.
#' @param statePrecModels A formula describing the regression relationship between the variance
#' in the ecosystem state value and the covariates for all ecosystem states, or a list of
#' formulas with each element giving the regression relationship between the variance in the
#' ecosystem state value and the covariates for each ecosystem state (ordered according to
#' intercept of the ecosystem state value on the y-axis).  The ecosystem state variable must
#' be given on the left-hand side of the formula.
#' @param inputData A data frame containing the covariate information and ecosystem state
#' variable.
#' @param numStates An integer scalar containing the number of stable states in the ecosystem.
#' If any of the \code{stateValModels}, \code{stateProbModels}, or \code{statePrecModels} parameters
#' is a list then \code{numStates} can be omitted and is therefore set to the length of the
#' list.
#' @param stateValError A description of the error distribution and link function to be used
#' in the model describing the ecosystem state value.  This can be from the \link[stats]{family}
#' specification or \code{character} scalar with the following possible values: \code{"gaussian"},
#' \code{"gamma"}, \code{"beta"}, \code{"negbinomial"}, or \code{"betabinomial"}.
#' @param setPriors A named list of prior distributions. Distribution are specified using character
#' strings. If sublists are not provided, values from the list are distributed to all sublist items
#' allowing to specify several priors at once. Sublist items are \code{"int"}, for specification of
#' priors on intercept parameters, and \code{"pred"}, from specification of priors on predictor
#' parameters. \code{"int"} is followed by \code{"1"} or \code{"2"} marking priors for the first
#' intercept and all the other intercepts respectively. For full structure of the list see default
#' values. Prior \code{"stateVal$Int2"} should allow only positive values to ensure distinctness of
#' states.
#'
#' @return A list containing the following components:
#' \itemize{
#' \item{\code{modelText}}{A character scalar containing the text of the model specification in
#' the NIMBLE BUGS dialect}
#' \item{\code{modelCode}}{A \code{nimbleCode} object containing the model specification}
#' \item{\code{constants}}{A list containing the constants to be provided to NIMBLE}
#' \item{\code{data}}{A list containing the data to be provided to NIMBLE}
#' \item{\code{errorModel}}{A factor containing the error model used for the specification of
#' the error distribution for the ecoystem state variable}
#' \item{\code{linkFunction}}{A factor containing the link function used in the specification
#' of the ecosystem state variable models}
#' \item{\code{initialValues}}{A list containing potential initial values for each of the stochastic
#' nodes in the NIMBLE model specification}
#' }
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @keywords internal
#'
modelSpecificationMultinomialEcosystemState <- function(
  stateValModels,
  stateProbModels,
  statePrecModels,
  inputData,
  numStates = NULL,
  stateValError = gaussian,
  setPriors = list(
    stateVal = list(
      int1 = "dnorm(0.0, 0.001)",
      int2 = "dgamma(0.001, 0.001)",
      pred = "dnorm(0.0, 0.001)"),
    stateProb = list(
      int2 = "dnorm(0.0, 0.001)",
      pred = "dnorm(0.0, 0.001)"),
    statePrec = list(
      int = "dnorm(0.0, 0.001)",
      pred = "dnorm(0.0, 0.001)"))
) {
  # Small helper function to test whether a variable is a formula
  is.formula <- function(inVal){
    inherits(inVal, "formula")
  }
  # The supported error distributions
  supportedError <- c("gaussian", "gamma", "beta", "negbinomial", "betabinomial")
  # The suported link functions
  supportedLink <- c("identity", "log", "logit", "probit", "cloglog")
  # Sanity test the error distribution for the state error
  inStateValError <- stateValError
  inLinkFunction <- NA
  if(is.function(inStateValError)) {
    # If it is a function then call it to retrieve the family object
    inStateValError <- stateValError()
  }
  if(is.factor(inStateValError)) {
    # If it is a factor then convert it to a string
    inStateValError <- as.character(inStateValError)
  }
  if(is.character(inStateValError)) {
    # If it is a string then use the default link function associated with the specified family
    if(tolower(inStateValError[1]) == "gaussian") {
      inStateValError <- factor("gaussian", levels = supportedError)
      inLinkFunction <- factor("identity", levels = supportedLink)
    } else if(tolower(inStateValError[1]) == "gamma") {
      inStateValError <- factor("gamma", levels = supportedError)
      inLinkFunction <- factor("log", levels = supportedLink)
    } else if(tolower(inStateValError[1]) == "beta") {
      inStateValError <- factor("beta", levels = supportedError)
      inLinkFunction <- factor("logit", levels = supportedLink)
    } else if(tolower(inStateValError[1]) == "negbinomial") {
      inStateValError <- factor("negbinomial", levels = supportedError)
      inLinkFunction <- factor("log", levels = supportedLink)
    } else if(tolower(inStateValError[1]) == "betabinomial") {
      inStateValError <- factor("betabinomial", levels = supportedError)
      inLinkFunction <- factor("logit", levels = supportedLink)
    } else {
      stop("selected error family is not supported")
    }
  } else if(class(inStateValError) == "family") {
    # If it is a "family" object then use the set error distribution and link function
    inLinkFunction <- factor(tolower(inStateValError$link), levels = supportedLink)
    inStateValError <- factor(tolower(inStateValError$family), levels = supportedError)
    if(is.na(inStateValError)) {
      stop("selected error family is not supported")
    }
    if(is.na(inLinkFunction)) {
      stop("selected link function is not supported")
    }
  } else {
    stop("error family specification is of invalid type")
  }
  inStateValModels <- stateValModels
  inStateProbModels <- stateProbModels
  inStatePrecModels <- statePrecModels
  inNumStates <- 1
  if(is.null(numStates)) {
    # If the maximum number of states is not set then find the largest model component to set the number
    # of states from
    inNumStates <- as.integer(max(c(
      ifelse(is.list(inStateValModels), length(inStateValModels), 1),
      ifelse(is.list(inStateProbModels), length(inStateProbModels), 1),
      ifelse(is.list(inStatePrecModels), length(inStatePrecModels), 1)
    )))
  } else {
    inNumStates <- tryCatch(as.integer(numStates), error = function(err) {
      stop("invalid entry for the number of states: ", err)
    })
  }
  if(inNumStates <= 0) {
    # Ensure that the number of states is greater than zero
    stop("invalid entry for the number of states: values equal or less than zero given")
  }
  # Define a function to facilitate the recycling and type specification of the model specification lists
  recycleModelFormulae <- function(modelInput, numStates, allowNull = FALSE) {
    inModelInput <- modelInput
    if((allowNull && is.null(modelInput)) || is.formula(modelInput)) {
      # If the model input is a formula or NULL then recycle that value to form a list
      inModelInput <- replicate(numStates, modelInput, simplify = FALSE)
    } else {
      # Otherwise force the input to be a list
      inModelInput <- tryCatch(as.list(modelInput), error = function(err) {
        stop("invalid entry for the model specification list: ", err)
      })
    }
    # Ensure that the list is of at least length one
    if(length(inModelInput)) {
      inModelInput <- lapply(X = inModelInput[((1:numStates - 1) %% length(inModelInput)) + 1], FUN = function(curEntry, allowNull) {
        outVal <- NULL
        if(is.null(curEntry)) {
          if(!allowNull) {
            stop("error encountered during the processing of the model specification list: NULL entries encountered")
          }
        } else {
          outVal <- tryCatch(as.formula(curEntry), error = function(err) {
            stop("error encountered during the processing of the model specification list: ", err)
          })
        }
        outVal
      }, allowNull = allowNull)
    } else {
      stop("invalid entry for the model specification list: list is empty")
    }
    inModelInput
  }
  # Recycle the model formulae so that they are of the correct length and type
  inStateValModels <- recycleModelFormulae(inStateValModels, inNumStates, FALSE)
  inStateProbModels <- recycleModelFormulae(inStateProbModels, inNumStates, FALSE)
  inStatePrecModels <- recycleModelFormulae(inStatePrecModels, inNumStates, TRUE)
  # Retrieve the list of formula strings for each of the models
  formulaStrings <- matrix(unlist(lapply(X = c(inStateValModels, inStateProbModels, inStatePrecModels), FUN = function(curForm) {
    outText <- NA
    if(!is.null(curForm)) {
      outText <- Reduce(paste, deparse(curForm))
    }
    outText
  })), ncol = 3, dimnames = list(NULL, c("stateVal", "stateProb", "statePrec")))
  # Retrieve the names of any response variables mentioned in any of the models
  respVariables <- unique(as.vector(gsub("\\s*~.*$", "", formulaStrings, perl = TRUE)))
  respVariables <- respVariables[!is.na(respVariables) & respVariables != ""]
  if(length(respVariables) != 1) {
    stop("invalid entry for the response variable: only one variable name must be present on the left-hand side of the formulae")
  }
  # Retrieve the names of the predictor variables mentioned in any of the models
  predVariables <- unique(gsub("^.*~\\s*", "", formulaStrings, perl = TRUE))
  predVariables <- predVariables[!is.na(predVariables) & predVariables != ""]
  predVariables <- unlist(lapply(X = predVariables, FUN = function(curVars) {
    strsplit(curVars, "\\s*\\+\\s*", perl = TRUE)[[1]]
  }))
  # Create a formula with the entrie set of predictor variables
  fullFormula <- as.formula(paste(respVariables, "~", paste(predVariables, collapse = " + "), sep = " "))
  # Retrieve the model matrix
  modelMatrix <- tryCatch(model.matrix(fullFormula, model.frame(fullFormula, inputData, na.action = NULL)), error = function(err) {
    stop("error thrown during construction of the model matrix: ", err)
  })
  # Remove the intercept term in the model matrix
  modelMatrix <- modelMatrix[, colnames(modelMatrix) != "(Intercept)", drop = FALSE]
  # Retrieve the model response variable
  respValues <- model.response(model.frame(fullFormula, inputData, na.action = NULL))
  numTrials <- NULL
  if(!is.null(dim(respValues))) {
    if(length(dim(respValues)) == 2) {
      if(ncol(respValues == 1)) {
        # Treat a single column in a same way as a dimensionless
        respValues <- tryCatch(as.numeric(respValues), error = function(err) {
          stop("error thrown during processing of response variable: ", err)
        })
        numTrials <- rep(NA, length(respValues))
      } else if(ncol(respValues) == 2) {
        # Check to ensure that the values in these columns take appropriate values
        if(any(respValues < 0.0)) {
          stop("a two-column response variable must have only positive values")
        }
        numTrials <- tryCatch(as.numeric(apply(X = as.matrix(respValues), FUN = sum, na.rm = TRUE, MARGIN = 1)), error = function(err) {
          stop("error thrown during processing of the number of trials: ", err)
        })
        respValues <- tryCatch(as.numeric(respValues[, 1]), error = function(err) {
          stop("error thrown during processing of response variable: ", err)
        })
      } else {
        stop("response variable can only have one or two columns")
      }
    } else {
      stop("response variable has invalid dimension structure")
    }
  } else {
    # Process the response variable
    respValues <- tryCatch(as.numeric(respValues), error = function(err) {
      stop("error thrown during processing of response variable: ", err)
    })
    numTrials <- rep(NA, length(respValues))
  }
  if(length(respValues) != nrow(modelMatrix)) {
    stop("response variable does not have the same length as the predictor variables")
  }
  # Convert the name of the response variable to something BUGS-friendly
  respVariablesBUGS <- setBUGSVariableName(respVariables)
  # Rename the covariates from the model matrix to something BUGS-friendly
  covariatesBUGS <- setBUGSVariableName(colnames(modelMatrix))
  # Very occasionally the conversion of the covariate names results in duplicates: this line protects from the possibility
  covariatesBUGS <- setNames(
    ifelse(duplicated(covariatesBUGS), paste(covariatesBUGS, setBUGSVariableName(as.character(1:length(covariatesBUGS))), sep = "_"), covariatesBUGS),
    colnames(modelMatrix))
  colnames(modelMatrix) <- covariatesBUGS
  # Set the link prefix and suffix for the state value model
  linkPrefix <- ""
  linkSuffix <- ""
  if(inLinkFunction != "identity") {
    linkPrefix <- paste(inLinkFunction, "(", sep = "")
    linkSuffix <- ")"
  }
  # Function to retrieve the covariate names from each sub-model formula
  getCovNames <- function(curFormula, inputData, covariatesBUGS) {
    outNames <- NA
    if(!is.na(curFormula)) {
      # Retrieve the model covariates
      modelCovars <- tryCatch(colnames(model.matrix(as.formula(curFormula), model.frame(as.formula(curFormula), inputData, na.action = NULL))), error = function(err) {
        stop("error thrown during construction of sub-model matrix: ", err)
      })
      modelCovars <- modelCovars[modelCovars != "(Intercept)"]
      outNames <- c("intercept", covariatesBUGS[modelCovars])
      if(anyNA(outNames)) {
        stop("error thrown during construction of sub-model matrix: covariates present in a sub-model that are not found in the full model")
      }
    }
    outNames
  }
  # Check priors
  if (!all(grepl("^d.+\\(.*)$", unlist(setPriors))))
    stop("unexpected prior specification")
  # Helper function which takes call of the list and modify/amends its items on all levels
  callModify <- function(oldCall, mods){
    callObj <- eval(oldCall)
    # Recursive function which modifies items of a list (and add new ones)
    recMod <- function(target, mod){
      for (i in names(mod)){
        target[i] <- if (is.list(mod[[i]]))
          list(recMod(target[[i]], mod[[i]])) else mod[[i]]
      }
      target
    }
    # Recursive function which copies structure of the list propagating items to lower levels
    recStr <- function(structure, content){
      for (i in names(structure)){
        contentItem <- if (is.null(names(content))) content else content[[i]]
        structure[i] <- if (is.list(structure[[i]]))
          list(recStr(structure[[i]], contentItem)) else contentItem
      }
      structure
    }
    recStr(callObj, recMod(callObj, mods))
  }
  # Prepare full object specifying them priors
  inSetPriors <- callModify(formals(modelSpecificationMultinomialEcosystemState)$setPriors, setPriors)
  # Initialise a vector to store potential initial values for the model
  initialValues <- as.character(c())
  # If the model has more than one state (99% of the times this function will be called) then create the relevant model strings
  if(inNumStates > 1) {
    # Create a matrix of strings containing the model specification
    modelStrings <- t(sapply(X = 1:inNumStates, FUN = function(curState, formulaStrings, inputData, covariatesBUGS, stateValError, linkFunction) {
      # Create a string version of the current state number
      stateString <- tolower(setBUGSVariableName(curState))
      # Retrieve the covariate names for each of the sub-models
      stateValCovs <- getCovNames(formulaStrings[curState, 1], inputData, covariatesBUGS)
      stateValCovs_nonIntercept <- stateValCovs[stateValCovs != "intercept"]
      stateProbCovs <- getCovNames(formulaStrings[curState, 2], inputData, covariatesBUGS)
      statePrecCovs <- getCovNames(formulaStrings[curState, 3], inputData, covariatesBUGS)
      # Prepare auxiliary vectors of priors
      auxPriorsPrec <- c(inSetPriors$statePrec$int, rep(inSetPriors$statePrec$pred, length(statePrecCovs) - 1))
      auxPriorsProb <- c(inSetPriors$stateProb$int2, rep(inSetPriors$stateProb$pred, length(stateProbCovs) - 1))
      # Set the prior text for the state value model
      priorStateVal <- paste(
        paste("\t# Set priors for the state variable value model for state ", stateString, sep = ""),
        if (length(stateValCovs_nonIntercept) > 0) paste("\t", stateValCovs_nonIntercept, "_stateVal[", curState, "] ~ ", inSetPriors$stateVal$pred, sep = "", collapse = "\n"),
        # The intercept of the first state has a normal prior.  All other states are forced to have positive priors in order to ensure that the
        # state labels are ordered and that MCMC doesn't just do state relabelling.
        paste("\tintercept_stateVal[", curState, "] ~ ", ifelse(curState == 1, inSetPriors$stateVal$int1, inSetPriors$stateVal$int2), sep = ""),
        sep = "\n")
      # Set the prior text for the state probability model
      priorStateProb <- "\t# The first state probability model is a baseline model so has no parameters"
      if(curState > 1) {
        priorStateProb <- paste(
          paste("\t# Set priors for the state probability model for state ", stateString, sep = ""),
          paste("\t", stateProbCovs, "_stateProb[", curState, "] ~", auxPriorsProb, sep = "", collapse = "\n"),
          sep = "\n")
      }
      # Set the prior text for the state precision model
      # Intially assume a simple multiplier model
      priorStatePrec <- paste(
        paste("\t# Set priors for the state precision model for state ", stateString, " (simple multiplier model)", sep = ""),
        paste("\tlinStateProb_statePrec[", curState, "] ~", inSetPriors$statePrec$pred, sep = ""),
        paste("\tintercept_statePrec[", curState, "] ~", inSetPriors$statePrec$int, sep = ""),
        sep = "\n")
      if(!is.na(formulaStrings[curState, 3])) {
        # If a formula has been specified for the precision model then use a linear sub-model instead
        priorStatePrec <- paste(
          paste("\t# Set priors for the state precision model for state ", stateString, sep = ""),
          paste("\t", statePrecCovs, "_statePrec[", curState, "] ~", auxPriorsPrec, sep = "", collapse = "\n"),
          sep = "\n")
      }
      # Set the model specification text for the state value model
      likelihoodStateVal <- paste(
        paste("\t\t# Set the model specification for the state value for state ", stateString, sep = ""),
        paste("\t\t", linkPrefix, "linStateVal[dataIter, ", curState, "]", linkSuffix, " <- ",
          ifelse(curState > 1, paste("sum(intercept_stateVal[1:", curState, "])", sep = ""), "intercept_stateVal[1]"), " * intercept[dataIter]", if (length(stateValCovs_nonIntercept > 0)) "+",
          if (length(stateValCovs_nonIntercept > 0)) paste(stateValCovs_nonIntercept, "_stateVal[", curState, "] * ", stateValCovs_nonIntercept, "[dataIter]", sep = "", collapse = " + "), sep = ""),
        sep = "\n")
      # Set the model specification text for the state probability model
      likelihoodStateProb <- paste("\t\t# Set the model specification for the state probability model for state ", stateString, sep = "")
      if(curState > 1) {
        likelihoodStateProb <- paste(
          likelihoodStateProb,
          paste("\t\tlog(linStateProb[dataIter, ", curState, "]) <- ", paste(stateProbCovs, "_stateProb[", curState, "] * ", stateProbCovs, "[dataIter]", sep = "", collapse = " + "), sep = ""),
          sep = "\n")
      } else {
        likelihoodStateProb <- paste(
          likelihoodStateProb,
          paste("\t\tlinStateProb[dataIter, ", curState, "] <- 1.0", sep = ""),
          sep = "\n")
      }
      # Set the model specification text for the state precision model
      # Initially assume a simple multiplier model
      likelihoodStatePrec <- paste(
        paste("\t\t# Set the model specification for the state precision model for state ", stateString, sep = ""),
        paste("\t\tlog(linStatePrec[dataIter, ", curState, "]) <- intercept_statePrec[", curState, "] + linStateProb_statePrec[", curState, "] * linStateProb[dataIter, ", curState, "] / sum(linStateProb[dataIter, 1:numStates])", sep = ""),
        sep = "\n")
      if(!is.na(formulaStrings[curState, 3])) {
        # If a formula has been specified for the precision model then use a linear sub-model instead
        likelihoodStatePrec <- paste(
          paste("\t\t# Set the model specification for the state precision model for state ", stateString, sep = ""),
          paste("\t\tlog(linStatePrec[dataIter, ", curState, "]) <- ", paste(statePrecCovs, "_statePrec[", curState, "] * ", statePrecCovs, "[dataIter]", sep = "", collapse = " + "), sep = ""),
          sep = "\n")
      }
      # Set an output vector with the appropriate names
      setNames(
        c(priorStateVal, priorStateProb, priorStatePrec, likelihoodStateVal, likelihoodStateProb, likelihoodStatePrec),
        c("priorValModel", "priorProbModel", "priorPrecModel", "likelihoodValModel", "likelihoodProbModel", "likelihoodPrecModel"))
    }, formulaStrings = formulaStrings, inputData = inputData, covariatesBUGS = covariatesBUGS, stateValError = as.character(inStateValError), linkFunction = as.character(inLinkFunction)))
    # Assign the error distribution for the ecosystem state value model
    errorStrings <- paste("\t\t# Set the error specification model for the state value sub-model", switch(as.character(inStateValError),
      "gaussian" = paste(
        "\t\t", respVariablesBUGS, "[dataIter] ~ dnormStateValueMembership(linStateVal[dataIter, 1:numStates], linStatePrec[dataIter, 1:numStates], linStateProb[dataIter, 1:numStates])",
      sep = ""),
      "gamma" = paste(
        "\t\t", respVariablesBUGS, "[dataIter] ~ dgammaStateValueMembership(linStateVal[dataIter, 1:numStates], linStatePrec[dataIter, 1:numStates], linStateProb[dataIter, 1:numStates])",
        sep = ""),
      "beta" = paste(
        "\t\t", respVariablesBUGS, "[dataIter] ~ dbetaStateValueMembership(linStateVal[dataIter, 1:numStates], linStatePrec[dataIter, 1:numStates], linStateProb[dataIter, 1:numStates])",
        sep = ""),
      "negbinomial" = paste(
        "\t\t", respVariablesBUGS, "[dataIter] ~ dnegbinStateValueMembership(linStateVal[dataIter, 1:numStates], linStatePrec[dataIter, 1:numStates], linStateProb[dataIter, 1:numStates])",
        sep = ""),
      "betabinomial" = paste(
        "\t\t", respVariablesBUGS, "[dataIter] ~ dbetabinStateValueMembership(linStateVal[dataIter, 1:numStates], linStatePrec[dataIter, 1:numStates], linStateProb[dataIter, 1:numStates], numTrials[dataIter])",
        sep = "")
    ), sep = "\n")
    # Create a vector of potential initial values for the model parameters
    initialValues <- unlist(lapply(X = 1:inNumStates, FUN = function(curState, formulaStrings, inputData, covariatesBUGS) {
      # Retrieve the covariate names for each of the sub-models
      stateValCovs <- getCovNames(formulaStrings[curState, 1], inputData, covariatesBUGS)
      stateValCovs_nonIntercept <- stateValCovs[stateValCovs != "intercept"]
      stateProbCovs <- getCovNames(formulaStrings[curState, 2], inputData, covariatesBUGS)
      statePrecCovs <- getCovNames(formulaStrings[curState, 3], inputData, covariatesBUGS)
      # Initialise a vecotr of output values
      outValuesNames <- c("intercept", stateValCovs_nonIntercept)
      outValues <- setNames(c(
        rnorm(length(stateValCovs_nonIntercept), 0.0, 4.0),
        ifelse(curState > 1, abs(rnorm(1, 0.0, 4.0)), rnorm(1, 0.0, 4.0))
      ), paste(outValuesNames, "_stateVal[", curState, "]", sep = ""))
      if(curState > 1) {
        # Add the the probability sub-model parameters if the current state is greater than 1
        outValues <- c(outValues, setNames(
          rnorm(length(stateProbCovs), 0.0, 4.0),
          paste(stateProbCovs, "_stateProb[", curState, "]", sep = "")
        ))
      }
      # Add the precision sub-model parameters
      if(is.na(formulaStrings[curState, 3])) {
        outValues <- c(outValues, setNames(
          rnorm(2, 0.0, 4.0),
          paste(c("intercept_statePrec", "linStateProb_statePrec"), "[", curState, "]", sep = "")
        ))
      } else {
        outValues <- c(outValues, setNames(
          rnorm(length(statePrecCovs), 0.0, 4.0),
          paste(statePrecCovs, "_statePrec[", curState, "]", sep = "")
        ))
      }
    }, formulaStrings = formulaStrings, inputData = inputData, covariatesBUGS = covariatesBUGS))
  } else {
    # If the model has only one state then completely remove the multinomial state latent state components
    # Retrieve the covariate names for each of the sub-models
    stateValCovs <- getCovNames(formulaStrings[1, 1], inputData, covariatesBUGS)
    statePrecCovs <- getCovNames(formulaStrings[1, 3], inputData, covariatesBUGS)
    # Create a matrix of model text
    auxPriorsVal <- c(inSetPriors$stateVal$int1, rep(inSetPriors$stateVal$pred, length(stateValCovs) - 1))
    auxPriorsPrec <- c(inSetPriors$statePrec$int, rep(inSetPriors$statePrec$pred, length(statePrecCovs) - 1))
    modelStrings <- matrix(nrow = 1, dimnames = list(NULL, c("priorValModel", "priorProbModel", "priorPrecModel", "likelihoodValModel", "likelihoodProbModel", "likelihoodPrecModel")), data = c(
      # Set the prior for the state value model
      paste(
        "\t# Set priors for the state variable value model",
        paste("\t", stateValCovs, "_stateVal ~", auxPriorsVal, sep = "", collapse = "\n"),
      sep = "\n"),
      # Set the prior for the state probability model: there is no state probability model because there is only one state
      "\t# There is no state probability model because there is only one state",
      # Set the prior for the state precision model
      ifelse(is.na(formulaStrings[1, 3]),
        paste("\t# Set priors for the state precision model\n\tintercept_statePrec ~", inSetPriors$statePrec$int),
        paste(
          "\t# Set priors for the state precision model",
          paste("\t", statePrecCovs, "_statePrec ~", auxPriorsPrec, sep = "", collapse = "\n"),
        sep = "\n")
      ),
      # Set the model specification text for the state value model
      paste(
        "\t\t# Set the model specification for the state value",
        paste("\t\t", linkPrefix, "linStateVal[dataIter]", linkSuffix, " <- ", paste(stateValCovs, "_stateVal * ", stateValCovs, "[dataIter]", sep = "", collapse = " + "), sep = ""),
        sep = "\n"),
      # Set the model specification text for the state probability model: there is no state probability model because there is only one state
      "\t\t# There is no state probability model because there is only one state",
      # Set teh model specification text for the state precision model
      ifelse(is.na(formulaStrings[1, 3]),
        "\t\t# Set the model specification for the state precision model\n\t\tlog(linStatePrec[dataIter]) <- intercept_statePrec",
        paste(
          "\t\t# Set the model specification for the state precision model",
          paste("\t\tlog(linStatePrec[dataIter]) <- ", paste(statePrecCovs, "_statePrec * ", statePrecCovs, "[dataIter]", sep = "", collapse = " + "), sep = ""),
          sep = "\n")
      )
    ))
    # Assign the error distribution for the ecosystem state precision model
    errorStrings <- paste("\t\t# Set the error specification model for the state value sub-model", switch(as.character(inStateValError),
      "gaussian" = paste("\t\t", respVariablesBUGS, "[dataIter] ~ dnorm(linStateVal[dataIter], linStatePrec[dataIter])", sep = ""),
      "gamma" = paste("\t\t", respVariablesBUGS, "[dataIter] ~ dgamma(mean = linStateVal[dataIter], sd = pow(linStatePrec[dataIter], -0.5))", sep = ""),
      "beta" = paste("\t\t", respVariablesBUGS, "[dataIter] ~ dbeta(mean = linStateVal[dataIter], sd = pow(linStatePrec[dataIter], -0.5))", sep = ""),
      "negbinomial" = paste("\t\t", respVariablesBUGS, "[dataIter] ~ dnegbin(\n\t\t\t1.0 - linStateVal[dataIter] * linStatePrec[dataIter], \n\t\t\tlinStateVal[dataIter] * linStateVal[dataIter] * linStatePrec[dataIter] / (1.0 - linStateVal[dataIter] * linStatePrec[dataIter]))", sep = ""),
      "betabinomial" = paste("\t\t", respVariablesBUGS, "[dataIter] ~ dbetabin(mean = linStateVal[dataIter], prec = linStatePrec[dataIter], size = numTrials[dataIter])")
    ), sep = "\n")
    # Create a vector of potential initial values for the model parameters
    initialValues <- setNames(rnorm(length(stateValCovs), 0.0, 4.0), paste(stateValCovs, "_stateVal", sep = ""))
    if(is.na(formulaStrings[1, 3])) {
      initialValues <- c(initialValues, setNames(rnorm(1, 0.0, 4.0), "intercept_statePrec"))
    } else {
      initialValues <- c(initialValues, setNames(rnorm(length(statePrecCovs), 0.0, 4.0), paste(statePrecCovs, "_statePrec", sep = "")))
    }
  }
  # Create NIMBLE model code
  modelText <- paste(
    "nimbleCode({",
    "\n\t#### PRIOR SPECIFICATION ####",
    "\n\t# Set priors for the ecosystem state value sub-model ----",
    paste(modelStrings[, "priorValModel"], collapse = "\n"),
    "\n\t# Set priors for the ecosystem state probability sub-model ----",
    paste(modelStrings[, "priorProbModel"], collapse = "\n"),
    "\n\t# Set priors for the ecosystem state precision sub-model ----",
    paste(modelStrings[, "priorPrecModel"], collapse = "\n"),
    "\n\t#### MODEL STRUCTURE ####",
    "\n\t# Iterate over each data point",
    "\tfor(dataIter in 1:numData) {",
    "\n\t\t# Set model structure for the ecosystem state value sub-model ----",
    paste(modelStrings[, "likelihoodValModel"], collapse = "\n"),
    "\n\t\t# Set model structure for the ecosystem state probability sub-model ----",
    paste(modelStrings[, "likelihoodProbModel"], collapse = "\n"),
    "\n\t\t# Set model structure for the ecosystem state precision sub-model ----",
    paste(modelStrings[, "likelihoodPrecModel"], collapse = "\n"),
    "\n\t\t# Set the error model for the ecosystem state value sub-model ----",
    errorStrings,
    "\t}",
    "})",
    sep = "\n")
  # Parse the NIMBLE model code to create a code object
  modelCode <- eval(parse(text = modelText))
  # Create a set of constants to use in the model
  modelConstants <- append(list(
    numData = nrow(modelMatrix),
    numStates = inNumStates,
    intercept = rep(1.0, nrow(modelMatrix))
  ), as.list(as.data.frame(modelMatrix)))
  if(as.character(inStateValError) == "betabinomial") {
    # Add the number of trials if the betabinomial error distribution is being used
    modelConstants <- append(modelConstants, list(numTrials = numTrials))
  }
  # Restructrue the initial values as a list
  vectorNames <- unique(gsub("\\[.*$", "", names(initialValues), perl = TRUE))
  initialValuesList <- setNames(lapply(X = vectorNames, FUN = function(curVecName, initialValues, inNumStates) {
    # Initialise an output vector
    outVec <- rep(0.0, inNumStates)
    if(inNumStates > 1) {
      # Retrieve the covairates with the current vector
      curCovs <- initialValues[grepl(paste("^", curVecName, "\\[\\d+\\]$", sep = ""), names(initialValues), perl = TRUE)]
      # Fill in the covariate values in the respective places
      outVec[as.integer(gsub("\\]$", "", gsub("^.*\\[", "", names(curCovs), perl = TRUE), perl = TRUE))] <- as.double(curCovs)
    } else {
      outVec <- initialValues[curVecName]
    }
    outVec
  }, initialValues = initialValues, inNumStates = inNumStates), vectorNames)
  # Create a set of data to use in the model
  modelData <- list(response = respValues)
  names(modelData) <- respVariablesBUGS
  list(
    modelText = modelText,
    modelCode = modelCode,
    constants = modelConstants,
    data = modelData,
    errorModel = inStateValError,
    linkFunction = inLinkFunction,
    initialValues = initialValuesList
  )
}

## 2. ------ DEFINE MODELLING FUNCTIONS ------

### 2.1. ==== Simulate Instances of the Multinomial Ecosystem State Model ====
#' @title Simulate Instances for the Multinomial Ecosystem State Model
#'
#' @description This function takes a set of model specifications for the three sub-model
#' components of the multinomial ecosystem state model and generates a simulation from this
#' model specification
#'
#' @param numSims The number of simulations to draw from the multinomial ecosystem state model.
#' @param stateValModels A formula describing the regression relationship between the mean
#' ecosystem state value and the covariates for all ecosystem states, or a list of formulas
#' with each element giving the regression relationship between the mean ecosystem state
#' value and the covariates for each ecosystem state (ordered according to intercept of the
#' ecosystem state value on the y-axis).  The ecosystem state variable must be given on the
#' left-hand side of the formula.
#' @param stateProbModels A formula describing the regression relationship between the probability
#' of the existence of ecosystem state and the covariates for all ecosystem states, or a list
#' of formulas with each element giving the regression relationship between the probability
#' of the existence of ecosystem state and the covariates for each ecosystem state (ordered
#' according to intercept of the ecosystem state value on the y-axis).  The ecosystem state
#' variable must be given on the left-hand side of the formula.
#' @param statePrecModels A formula describing the regression relationship between the variance
#' in the ecosystem state value and the covariates for all ecosystem states, or a list of
#' formulas with each element giving the regression relationship between the variance in the
#' ecosystem state value and the covariates for each ecosystem state (ordered according to
#' intercept of the ecosystem state value on the y-axis).  The ecosystem state variable must
#' be given on the left-hand side of the formula.
#' @param inputData A data frame containing the covariate information and ecosystem state
#' variable.
#' @param numStates An integer scalar containing the number of stable states in the ecosystem.
#' If any of the \code{stateValModels}, \code{stateProbModels}, or \code{statePrecModels} parameters
#' is a list then \code{numStates} can be omitted and is therefore set to the length of the
#' list.
#' @param stateValError A description of the error distribution and link function to be used
#' in the model describing the ecosystem state value.  This can be from the \link[stats]{family}
#' specification or \code{character} scalar with the following possible values: \code{"gaussian"},
#' \code{"gamma"}, \code{"beta"}, \code{"negbinomial"}, or \code{"betabinomial"}.
#' @param coefficientValues A list containing the values of the coefficients to use in the
#' simulation.
#'
#' @return A list containg a vector of simulated values for each stochastic node.  In addition the
#' following elements are appended to the list:
#' \itemize{
#' \item{\code{linStateVal}}{A matrix containing the predicted ecosystem state values at each row
#' of the input covariate data.frame.  Each column represents the predicted ecosystem state value
#' for each stable state}
#' \item{\code{linStatePrec}}{A matrix containing the predicted precision of the ecosystem state
#' value at each row of the input covariate data.frame.  Each column represents the predicted
#' precision of the ecosystem state value for each stable state}
#' \item{\code{linStateProb}}{A matrix containing the probability of each ecosystem state existing
#' at each row of the input covariate data.frame.  Each column represents the probability of each
#' ecosystem state existing for each stable state}
#' }
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
#'
simulateMultinomialEcosystemState <- function(
  numSims,
  stateValModels,
  stateProbModels,
  statePrecModels,
  inputData,
  numStates = NULL,
  stateValError = gaussian,
  coefficientValues = NULL
) {
  # Ensure that the coefficient values are correctly specified
  inCoefficientValues <- list()
  if(!is.null(coefficientValues)) {
    inCoefficientValues <- tryCatch(as.list(coefficientValues), error = function(err) {
      stop("error encountered during processing of the coefficient value list: ", err)
    })
  }
  # Ensure that the number of simulations input is correctly specified
  inNumSims <- tryCatch(as.integer(numSims)[1], error = function(err) {
    stop("error encountered during processing of the number of simulations: ", err)
  })
  if(inNumSims <= 0) {
    stop("error encountered during processing of the number of simulations: number of simulations requested is less than 1")
  }
  # Create a NIMBLE model specification
  modelSpecification <- modelSpecificationMultinomialEcosystemState(stateValModels, stateProbModels, statePrecModels, inputData, numStates, stateValError)
  # Initialise the data with some dummy (but plausable) values for the relevant error model.  This is so that the NIMBLE
  # model is fully initialised (avoids some error messages and ensures NIMBLE starts from a valid state)
  modelSpecification$data[[1]] <- rep(switch(as.character(modelSpecification$errorModel),
    "gaussian" = 0.0,
    "gamma" = 0.1,
    "beta" = 0.5,
    "negbinomial" = 0,
    "betabinomial" = 0
    ), length(modelSpecification$data[[1]]))
  # If coefficient values have been provided then use those in the model initialisation
  if(length(inCoefficientValues) > 0 && !is.null(names(inCoefficientValues))) {
    # Overwrite the list of initial values with the input coefficients (where appropriate)
    modelSpecification$initialValues <- setNames(lapply(X = names(modelSpecification$initialValues), FUN = function(curName, curInits, inputInits) {
      # Initialise the current output
      outValues <- curInits[[curName]]
      if(curName %in% names(inputInits) && length(outValues) > 0) {
        outValues <- ifelse(is.na(inputInits[[curName]][(1:length(outValues) - 1) %% length(inputInits) + 1]), outValues, inputInits[[curName]][(1:length(outValues) - 1) %% length(inputInits) + 1])
      }
      outValues
    }, curInits = modelSpecification$initialValues, inputInits = inCoefficientValues), names(modelSpecification$initialValues))
    # Check the variable names to ensure that they are present in the model
    isInModel <- names(inCoefficientValues) %in% names(modelSpecification$initialValues)
    if(any(!isInModel)) {
      warning("some coefficient names are not present in the model: ", print(names(inCoefficientValues)[!isInModel], collapse = ", "))
    }
  }
  # Initialise the NIMBLE model
  modelObject <- nimbleModel(modelSpecification$modelCode, constants = modelSpecification$constants, data = modelSpecification$data, inits = modelSpecification$initialValues)
  # Specify a function to simulate data (of a particular data node)
  singleSimulation <- function(curDataName, modelObject, modelSpecification) {
    # Simulate the data for the model (overwrites the previously set dummy data)
    simulate(modelObject, names(modelSpecification$data), includeData = TRUE)
    # Retrieve the simulated data
    modelObject[[curDataName]]
  }
  # Initialise an output object with the simulations of the data
  outObject <- setNames(lapply(X = names(modelSpecification$data), FUN = function(curDataName, modelObject, modelSpecification, numSims) {
    # Call simulation function for each requested simulation
    outVec <- do.call(cbind, replicate(numSims, singleSimulation(curDataName, modelObject, modelSpecification), simplify = FALSE))
    # Ensure that the output object has two dimensions
    if(is.null(dim(outVec)) || length(dim(outVec)) <= 1) {
      dim(outVec) <- c(length(outVec), 1)
    }
    outVec
  }, modelObject = modelObject, modelSpecification = modelSpecification, numSims = inNumSims), names(modelSpecification$data))
  # Retrieve the linear predictors for the value and precision sub-models also
  outObject <- append(outObject, list(
    linStateVal = modelObject[["linStateVal"]],
    linStatePrec = modelObject[["linStatePrec"]]
  ))
  # Ensure consistent dimensions of the linear predictor objects
  if(is.null(dim(outObject$linStateVal)) || length(dim(outObject$linStateVal)) == 1) {
    dim(outObject$linStateVal) <- c(length(outObject$linStateVal), 1)
  }
  if(is.null(dim(outObject$linStatePrec)) || length(dim(outObject$linStatePrec)) == 1) {
    dim(outObject$linStatePrec) <- c(length(outObject$linStatePrec), 1)
  }
  # Add the probability model outputs if they are missing
  if(is.null(modelObject$linStateProb)) {
    outObject <- append(outObject, list(linStateProb = rep(1.0, length(outObject$linStateVal))))
    dim(outObject$linStateProb) <- c(length(outObject$linStateVal), 1)
  # Otherwise normalise the probability values and add it to the output list
  } else {
    outObject <- append(outObject, list(linStateProb = t(apply(X = modelObject$linStateProb, FUN = function(curRow) {
      curRow / sum(curRow, na.rm = TRUE)
    }, MARGIN = 1))))
  }
  outObject
}

### 2.2. ==== Fit the Multinomial Ecosystem State Model ====
#' @title Fit the Multinomial Ecosystem State Model
#'
#' @description This function takes a set of model specifications for the three sub-model
#' components of the multinomial ecosystem state model and fits them to a data set.
#'
#' @param stateValModels A formula describing the regression relationship between the mean
#' ecosystem state value and the covariates for all ecosystem states, or a list of formulas
#' with each element giving the regression relationship between the mean ecosystem state
#' value and the covariates for each ecosystem state (ordered according to intercept of the
#' ecosystem state value on the y-axis).  The ecosystem state variable must be given on the
#' left-hand side of the formula.
#' @param stateProbModels A formula describing the regression relationship between the probability
#' of the existence of ecosystem state and the covariates for all ecosystem states, or a list
#' of formulas with each element giving the regression relationship between the probability
#' of the existence of ecosystem state and the covariates for each ecosystem state (ordered
#' according to intercept of the ecosystem state value on the y-axis).  The ecosystem state
#' variable must be given on the left-hand side of the formula.
#' @param statePrecModels A formula describing the regression relationship between the variance
#' in the ecosystem state value and the covariates for all ecosystem states, or a list of
#' formulas with each element giving the regression relationship between the variance in the
#' ecosystem state value and the covariates for each ecosystem state (ordered according to
#' intercept of the ecosystem state value on the y-axis).  The ecosystem state variable must
#' be given on the left-hand side of the formula.
#' @param inputData A data frame containing the covariate information and ecosystem state
#' variable.
#' @param numStates An integer scalar containing the number of stable states in the ecosystem.
#' If any of the \code{stateValModels}, \code{stateProbModels}, or \code{statePrecModels} parameters
#' is a list then \code{numStates} can be omitted and is therefore set to the length of the
#' list.
#' @param stateValError A description of the error distribution and link function to be used
#' in the model describing the ecosystem state value.  This can be from the \link[stats]{family}
#' specification or \code{character} scalar with the following possible values: \code{"gaussian"},
#' \code{"gamma"}, \code{"beta"}, \code{"negbinomial"}, or \code{"betabinomial"}.
#' @param mcmcIter An integer scalar providing the number of post-burn-in samples to draw from the
#' MCMC sampler (per chain).
#' @param mcmcBurnIn An integer scalar providing the number of MCMC samples to use for the
#' adaption or burn-in portion of the process (per chain).
#' @param mcmcChains An integer scalar giving the number of MCMC chains to use.
#' @param mcmcThin An integer scalar giving the thinning frequency in the MCMC chains.  For example,
#' a value of \code{4} results in every fourth values being retained.
#' @param setPriors A named list of prior distributions. Distribution are specified using character
#' strings. If sublists are not provided, values from the list are distributed to all sublist items
#' allowing to specify several priors at once. Sublist items are \code{"int"}, for specification of
#' priors on intercept parameters, and \code{"pred"}, from specification of priors on predictor
#' parameters. \code{"int"} is followed by \code{"1"} or \code{"2"} marking priors for the first
#' intercept and all the other intercepts respectively. For full structure of the list see default
#' values. Prior \code{"stateVal$Int2"} should allow only positive values to ensure distinctness of
#' states.
#' @param setInit list of initial values which overwrites generated ones.
#'
#' @return A list containing the following components:
#' \itemize{
#' \item{\code{mcmcSamples}}{An \link[coda]{mcmc} object if \code{mcmcChains == 1} or \link[coda]{mcmc.list}
#' object if \code{mcmcChains > 1} containing the sampled values}
#' \item{\code{compiledModel}}{An R interface object containing an interface for the compiled NIMBLE model.
#' This is the same output as the \link[nimble]{compileNimble} function applied to the model object}
#' \item{\code{modelText}}{A character scalar containing the text of the model specification in
#' the NIMBLE BUGS dialect}
#' \item{\code{modelCode}}{A \code{nimbleCode} object containing the model specification}
#' \item{\code{constants}}{A list containing the constants provided to NIMBLE}
#' \item{\code{data}}{A list containing the data provided to NIMBLE}
#' \item{\code{errorModel}}{A factor containing the error model used for the specification of
#' the error distribution for the ecoystem state variable}
#' \item{\code{linkFunction}}{A factor containing the link function used in the specification
#' of the ecosystem state variable models}
#' \item{\code{initialValues}}{A list containing potential initial values used for each of the stochastic
#' nodes in the NIMBLE model specification}
#' }
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @keywords internal
#'
fitMultinomialEcosystemState <- function(
  stateValModels,
  stateProbModels,
  statePrecModels,
  inputData,
  numStates = NULL,
  stateValError = gaussian,
  mcmcIters = 10000,
  mcmcBurnin = 5000,
  mcmcChains = 4,
  mcmcThin = 1,
  setPriors = list(
    stateVal = list(
      int1 = "dnorm(0.0, 0.001)",
      int2 = "dgamma(0.001, 0.001)",
      pred = "dnorm(0.0, 0.001)"),
    stateProb = list(
      int2 = "dnorm(0.0, 0.001)",
      pred = "dnorm(0.0, 0.001)"),
    statePrec = list(
      int = "dnorm(0.0, 0.001)",
      pred = "dnorm(0.0, 0.001)")),
  setInit = NULL
) {
  # Create a NIMBLE model specification
  modelSpecification <- modelSpecificationMultinomialEcosystemState(stateValModels, stateProbModels, statePrecModels, inputData, numStates, stateValError, setPriors)
  # Change initial values if provided
  if (!is.null(setInit)) modelSpecification$initialValues <- setInit
  modelObject <- nimbleModel(modelSpecification$modelCode, constants = modelSpecification$constants, data = modelSpecification$data, inits = modelSpecification$initialValues)
  # Build the MCMC object and compile it
  varsToMonitor <- c(modelObject$getVarNames(), "linStateVal", "linStatePrec")
  if (grepl("linStateProb", modelSpecification$modelText)) varsToMonitor <- c(varsToMonitor, "linStateProb")
  mcmcObject <- buildMCMC(modelObject, enableWAIC = TRUE, monitors = varsToMonitor)
  mcmcObjectCompiled <- compileNimble(mcmcObject, modelObject)
  # Run the MCMC
  mcmcOutput <- runMCMC(mcmcObjectCompiled$mcmcObject, niter = mcmcIters, nburnin = mcmcBurnin, thin = mcmcThin, nchains = mcmcChains, WAIC = TRUE, samplesAsCodaMCMC = TRUE)
  # Structure the compiled model, the MCMC samples, and the model specification into a list
  out <- append(list(mcmcSamples = mcmcOutput, compiledModel = mcmcObjectCompiled), modelSpecification)
  class(out) <- "PaGAnmesm"
  out
}

## 3. ------ DEFINE GENERAL MODEL METHODS ------

### 3.1. ==== Plot results of Multinomial Ecosystem State Model ====
#' @title Plot results of Multinomial Ecosystem State Model
#'
#' @description This function plots results of multinomial ecosystem state model on the
#' current graphical device.
#'
#' @param form formula, such as y ~ pred, specifying variables to be plotted
#' @param mod an object of class "mesm"
#' @param yaxis vector of values to be marked on y-axis
#' @param transCol logical value indicating usage of transparent colours
#' @param addWAIC logical value indication display of WAIC in upper right corner of the plot
#' @param setCol vector of colours to be used for states
#' @param drawXaxis logical value indicating whether values should be marked on x-axis
#' @param SDmult scalar multiplying visualized standard deviation (to make lines for small standard deviation visible)
#' @param byChain logical value indicating whether to plot states for each chain
#' @param ... additional arguments passed to plot
#'
#' @return Returns invisibly a list containing posterior means of state value
#' coefficients for each chain used in plotting.
#'
#' @author Adam Klimes
#' @export
#'
plot.PaGAnmesm <- function(form, mod, yaxis, transCol = TRUE, addWAIC = FALSE,
                      setCol = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"),
                      drawXaxis = TRUE, SDmult = 1, byChains = TRUE, ...) {
  resp <- mod$data[[1]]
  dat <- data.frame(mod$data, mod$constants[sapply(mod$constants, length) ==
                                              length(resp)])
  svar <- labels(terms(form))
  svar <- svar[svar %in% names(dat)]
  auxRange <- max(resp) - min(resp)
  invlink <- switch(as.character(mod$linkFunction), identity = function(x) x, log = exp,
                    logit = function(x) exp(x)/(1+exp(x)))
  par(mai = c(0.8,0.8,0.1,0.1))
  plot(form, data = dat, ylim = c(min(resp) - 0.05 * auxRange, max(resp) + 0.3 * auxRange),
       yaxs = "i", axes = FALSE, ...)
  box(bty = "l")
  if (drawXaxis) axis(1)
  axis(2, labels = yaxis, at = yaxis)
  axis(2, labels = 0:1, at = c(max(resp) + 0.1 * auxRange, max(resp) + 0.25 * auxRange))
  abline(h = max(resp) + 0.05 * auxRange, lwd = 3)
  abline(h = max(resp) + 0.1 * auxRange, lty = 2)
  abline(h = max(resp) + 0.25 * auxRange, lty = 2)
  if (addWAIC) text(par("usr")[2] - (par("usr")[2] - par("usr")[1]) * 0.2, max(resp) + 0.175 * auxRange, paste("WAIC:", round(mod$mcmcSamples$WAIC$WAIC, 1)))
  parsTab <- summary.mesm(mod, byChains = byChains, absInt = TRUE, digit = NULL)
  auxLines <- function(parsChain, dat, mod){
    nstates <- mod$constants$numStates
    xx <- seq(min(dat[, svar]), max(dat[, svar]), length.out = 100)
    ind <- NULL
    cNames <- rownames(parsChain)
    if (nstates > 1) {
      ind <- paste0("[", 1:nstates, "]")
      probInt <- parsChain[paste0("intercept_stateProb", ind), "mean"]
    }
    valInt <- parsChain[paste0("intercept_stateVal", ind), "mean"]
    precInt <- parsChain[paste0("intercept_statePrec", ind), "mean"]
    valCov <- if (paste0(svar, "_stateVal", ind[1]) %in% cNames) parsChain[paste0(svar, "_stateVal", ind), "mean"] else rep(0, nstates)
    precCov <- if (paste0(svar, "_statePrec", ind[1]) %in% cNames) parsChain[paste0(svar, "_statePrec", ind), "mean"] else rep(0, nstates)
    probCov <- if (paste0(svar, "_stateProb", ind[1]) %in% cNames) parsChain[paste0(svar, "_stateProb", ind), "mean"] else rep(0, nstates)
    if (nstates > 1) {
      probVals <- as.matrix(data.frame(Map(function(int, cov) exp(int + cov * xx), probInt, probCov)))
      probVals[is.na(probVals)] <- 1
      probVals <- probVals / rowSums(probVals)
      probVals[is.nan(probVals)] <- 1
      }
    for (i in 1:nstates){
      cols <- setCol[i]
      if (nstates > 1) {
        lines(xx, max(resp) + 0.1 * auxRange + probVals[, i] * 0.15 * auxRange, col = setCol[i], lwd = 3)
        if (transCol) {
          rgbVec <- col2rgb(cols)[, 1]
          cols <- rgb(rgbVec[1], rgbVec[2], rgbVec[3], alpha = 40 + probVals[, i] * 215, maxColorValue = 255)
        }
      }
      sdVals <- 1 / sqrt(exp(precInt[i] + precCov[i] * xx))
      yEst <- do.call(invlink, list(valInt[i] + valCov[i] * xx))
      segments(head(xx, -1), head(yEst, -1), x1 = tail(xx, -1), y1 = tail(yEst, -1), col = cols, lwd = 3)
      lines(xx, do.call(invlink, list(valInt[i] + valCov[i] * xx + sdVals * SDmult)), col = setCol[i], lty = 2, lwd = 1)
      lines(xx, do.call(invlink, list(valInt[i] + valCov[i] * xx - sdVals * SDmult)), col = setCol[i], lty = 2, lwd = 1)
    }
  }
  invisible(lapply(parsTab, auxLines, dat, mod))
}

### 3.2. ==== Summary of Multinomial Ecosystem State Model ====
#' @title Summarize Multinomial Ecosystem State Model
#'
#' @description This function calculates posterior quantiles of parameters of
#' Multinomial Ecosystem State Model across all chains or for each chain separately
#'
#' @param object an object of class "mesm"
#' @param byChains logical value indicating if the summary should be calculated for each chain separately
#' @param digit integer specifying the number of decimal places to be used. Use \code{"NULL"} for no rounding.
#' @param absInt logical value indicating if intercepts for state values should be absolute (by default, they represent differences)
#'
#' @return Returns data.frame of quantiles of posterior of parameters
#'
#' @author Adam Klimes
#' @export
#'
summary.PaGAnmesm <- function(object, byChains = FALSE, digit = 4, absInt = FALSE){
  varsSamples <- lapply(object$mcmcSamples$samples,
    function(x) x[, !grepl(paste0("^lifted|^linState|^", names(object$data)), colnames(x))])
  if (!byChains) varsSamples <- list(do.call(rbind, varsSamples))
  sepInt <- function(samp){
    scol <- grepl("intercept_stateVal", colnames(samp))
    samp[, scol] <- t(apply(samp, 1, function(x, scol) cumsum(x[scol]), scol))
    samp
  }
  if (absInt) varsSamples <- lapply(varsSamples, sepInt)
  auxSummary <- function(x)
    c(mean = mean(x), sd = sd(x), quantile(x, c(0.025,0.25,0.75,0.975), na.rm = TRUE))
  out <- lapply(varsSamples, function(x) t(apply(x, 2, auxSummary)))
  # if (length(out) == 1) out <- out[[1]]
  if (!is.null(digit)) out <- lapply(out, round, digit)
  out
}

### 3.3. ==== Predict ecosystem characteristics based on Multinomial Ecosystem State Model ====
#' @title Predict from Multinomial Ecosystem State Model
#'
#' @description This function calculates probability curves for ecosystems based on Multinomial Ecosystem State Model
#'
#' @param mod an object of class "mesm"
#' @param newdata dataframe of predictor values of ecosystems to be predicted.
#'   If not provided, prediction is done for modelled data.
#' @param samples number of samples to take along the respons variable
#'
#' @return A list containing the following components:
#' \itemize{
#' \item{\code{sampledResp}}{A numeric vector of samples along response variable}
#' \item{\code{probCurves}}{A data frame of probability curves for each observation}
#' \item{\code{tipPoints}}{A list of tipping points for each observation. Border values are never included}
#' \item{\code{stableStates}}{A list of stable states for each observation}
#' \item{\code{obsDat}}{A data frame containing values of response variable,
#' distance to closest tipping point and stable state for each observation}
#' }
#'
#' @author Adam Klimes
#' @export
#'
predict.PaGAnmesm <- function(mod, newdata = NULL, samples = 1000){
  if (is.null(newdata)) newdata <- as.data.frame(mod$constants[-(1:3)])
  form <- formula(paste("~", colnames(newdata)[1]))
  slices <- slice.mesm(form, mod, value = newdata, byChains = FALSE, doPlot = FALSE, samples = samples)
  probCurve <- as.data.frame(slices[[1]])
  names(probCurve) <- paste0("obs", seq_along(probCurve))
  resp <- slices$resp
  respVal <- mod$data[[1]]
  getMin <- function(x, resp, extremes = TRUE) {
    id <- findMin(x)
    if (!extremes) id <- id[!id %in% c(1, length(x))]
    out <- resp[id] * (1 - id %% 1) + resp[min(id + 1, length(resp))] * id %% 1
    if (length(out) == 0) out <- NA
    out
  }
  tipPoints <- lapply(probCurve, getMin, resp, extremes = FALSE)
  stableStates <- lapply(-probCurve, getMin, resp)
  auxDist <- function(x, target) min(abs(x - target))
  obsDat <- data.frame(
    respVal = respVal,
    distToTip = unlist(Map(auxDist, respVal, tipPoints)),
    distToState = unlist(Map(auxDist, respVal, stableStates)))
  list(sampledResp = resp,
              probCurves = probCurve,
              tipPoints = tipPoints,
              stableStates = stableStates,
              obsDat = obsDat)
}

## 4. ------ DEFINE HELPER FUNCTIONS ------

### 4.1. ==== Find positions of local minima in a vector ====
#' @title Find positions of local minima in a vector
#'
#' @description Finds position of all local minima in a vector including start
#'   and end point. For flat minima (identical subsequent values), it denotes middle point
#'
#' @param x numeric vector
#'
#' @return Positions of minima in x
#'
#' @author Adam Klimes
#' @keywords internal
#'
findMin <- function(x){
  dfXin <- diff(x)
  seqCount <- diff(c(0, which(dfXin != 0), length(x)))
  Nflat <- rep(seqCount, seqCount) - 1
  xClear <- x[c(TRUE,  dfXin != 0)]
  dfX <- diff(xClear)
  loc <- which(diff(sign(dfX)) == 2) + 1
  if (dfX[1] > 0) loc <- c(1, loc)
  if (tail(dfX, 1) < 0) loc <- c(loc, length(xClear))
  inLoc <- seq_along(x)[c(TRUE, dfXin != 0)][loc]
  inLoc[inLoc %in% which(dfXin == 0)] <-
    0.5 * Nflat[inLoc[inLoc %in% which(dfXin == 0)]] + inLoc[inLoc %in% which(dfXin == 0)]
  inLoc
}

### 4.2. ==== Plot slice from Multinomial Ecosystem State Model ====
#' @title Plot slice from Multinomial Ecosystem State Model
#'
#' @description This function plots probability density for given predictor value
#'
#' @param form formula with one predictor specifying which variables to plot
#' @param mod an object of class "mesm"
#' @param value numeric vector of values of the preditor specified by
#'   \code{"form"} where the slice is done or data.frame with values of predictors
#'   in named columns
#' @param byChains logical value indicating if slice should be done for each chain separately
#' @param xlab string used as label for x-axis
#' @param doPlot logical value indicating if plotting should be done
#' @param setCol vector of colours to be used for visualization of estimated states
#' @param plotEst logical value indicating if estimated states should be visualized
#' @param xaxis logical value indicating if values should be marked on x-axis
#' @param addEcos logical value indicating if ecosystems within \code{"ecosTol"} from \code{"value"} should be visualized on the line
#' @param ecosTol scalar specifying range of predictor from the \code{"value"} to select ecosystems to be visualized
#' @param samples scalar specifying number of samples to be taken along predictor
#'
#' @return Returns list of plotted values
#'
#' @author Adam Klimes
#' @export
#'
sliceMESM <- function(form, mod, value = 0, byChains = TRUE, xlab = "", doPlot = TRUE,
                       setCol = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"),
                       plotEst = TRUE, xaxis = TRUE, addEcos = FALSE, ecosTol = 0.1, samples = 1000){
  resp <- mod$data[[1]]
  parsTab <- summary.mesm(mod, byChains = byChains, absInt = TRUE, digit = NULL)
  svar <- labels(terms(form))
  Nstates <- mod$constants$numStates
  invlink <- switch(as.character(mod$linkFunction), identity = function(x) x, log = exp,
                    logit = function(x) exp(x)/(1+exp(x)))
  xx <- seq(min(resp), max(resp), length.out = samples)
  if (addEcos) {
    pred <- mod$constants[[svar]]
    xx <- c(xx, resp[abs(pred - value) < ecosTol])
  }
  if (is.null(dim(value))) value <- matrix(value, length(value), dimnames = list(NULL, svar))
  value <- as.matrix(value)
  plotSlice <- function(pars, value, mod){
    getPars <- function(curState, pars, value){
      auxExtract <- function(toGet, curState, pars, value){
        pos <- match(paste0(c("intercept",colnames(value)), "_", toGet, "[", curState, "]"), rownames(pars))
        pars[pos[1], "mean"] + as.vector(value %*% pars[pos[-1], "mean"])
      }
      est <- auxExtract("stateVal", curState, pars, value)
      prec <- auxExtract("statePrec", curState, pars, value)
      prob <- auxExtract("stateProb", curState, pars, value)
      cbind(est = do.call(invlink, list(est)), sd = 1 / sqrt(exp(prec)), prob = prob)
    }
    parsVal <- vapply(1:Nstates, getPars, FUN.VALUE = array(0, dim = c(nrow(value), 3)), pars, value)
    parsVal[, "prob", 1] <- rep(0, nrow(value))
    parsVal[, "prob", ] <- exp(parsVal[, "prob", ]) / rowSums(exp(parsVal[, "prob", , drop = FALSE]))
    rownames(parsVal) <- paste0("value", 1:nrow(value))
    aux <- parsVal[, "est", ]
    auxDim <- c(nrow(value), Nstates, 2)
    parsD <- switch(as.character(mod$errorModel),
                    gaussian = array(c(parsVal[, "est", ], parsVal[, "sd", ]), dim = auxDim),
                    gamma = array(c(aux^2/parsVal[, "sd", ]^2, aux/parsVal[, "sd", ]^2), dim = auxDim),
                    beta = array(c(aux*(aux*(1-aux)/parsVal[, "sd", ]^2-1), (aux*(1-aux)/parsVal[, "sd", ]^2-1)*(1-aux)), dim = auxDim))
    dfun <- switch(as.character(mod$errorModel),
                   gaussian = dnorm,
                   gamma = dgamma,
                   beta = dbeta)
    dens <- apply(parsD, 1, apply, 1, function(pars) do.call(dfun, list(xx, pars[1], pars[2])), simplify = FALSE)
    dens <- Map(function(den, prob) den * rep(prob, each = nrow(den)), dens, data.frame(t(matrix(parsVal[, "prob", ], nrow = nrow(value)))))
    dens <- lapply(dens, rowSums)
    densSt <- lapply(dens, function(x) x / max(x))
    if (doPlot){
      lines(xx[1:samples], densSt[[1]][1:samples])
      if (plotEst) {
        rgbVec <- col2rgb(setCol)
        cols <- rgb(rgbVec[1, ], rgbVec[2, ], rgbVec[3, ], alpha = 40 + parsVal[1, "prob", ] * 215, maxColorValue = 255)
        abline(v = parsVal[1, "est", ], lty = 2, lwd = 3, col = cols)
      }
      if (addEcos) points(tail(xx, -samples), tail(densSt[[1]], -samples), pch = 16)
    }
    densSt
  }
  if (doPlot){
    plot(range(resp), c(1, 0), type = "n", ylab = "Potential energy", xlab = xlab, ylim = c(1, 0), axes = FALSE, yaxs = "i")
    if (xaxis) axis(1)
    axis(2, labels = 0:5/5, at = 5:0/5, las = 2)
    box(bty = "l")
  }
  out <- lapply(parsTab, plotSlice, value, mod)
  invisible(c(out, list(resp = xx)))
}

### 4.3. ==== Probability landscape from Multinomial Ecosystem State Model ====
#' @title Plot probability landscape from Multinomial Ecosystem State Model
#'
#' @description This function plots probability landscape for given predictor
#'
#' @param form formula with one predictor specifying which variables to plot
#' @param mod an object of class "mesm"
#' @param addPoints logical value indicating if ecosystems should be visualized
#' @param addMinMax logical value indicating if stable states and tipping points should be visualized
#' @param ... parameters passed to image()
#'
#' @return Returns Probability density (scaled to [0,1]) matrix.
#'
#' @author Adam Klimes
#' @export
#'
landscapeMESM <- function(form, mod, addPoints = TRUE, addMinMax = TRUE, ...){
  svar <- labels(terms(form))
  resp <- mod$data[[1]]
  pred <- mod$constants[[svar]]
  grad <- seq(min(pred), max(pred), length.out = 500)
  slices <- slice.mesm(form, mod, value = grad, byChains = FALSE, doPlot = FALSE)
  mat <- do.call(cbind, slices[[1]])
  image(t(mat), ...)
  plotMinMax <- function(matCol, xCoors) {
    yCoors <- seq(0, 1, length.out = nrow(mat))
    mins <- findMin(matCol)
    maxs <- findMin(-matCol)
    points(rep(xCoors, length(maxs)), yCoors[maxs], pch = 16, col = "red", cex = 0.5)
    points(rep(xCoors, length(mins)), yCoors[mins], pch = 16, col = "blue", cex = 0.5) #-yCoors[mins]+1
  }
  if (addMinMax) Map(plotMinMax, data.frame(-mat+1), seq(0, 1, length.out = ncol(mat)))
  stRange <- function(x) (x - min(x)) / max(x - min(x))
  if (addPoints) points(stRange(pred), stRange(resp), cex = 0.4, pch = 16)
  invisible(mat)
}

### 4.4. ==== Fit Multinomial Ecosystem State Model using rasters ====
#' @title Fit Multinomial Ecosystem State Model using rasters
#'
#' @description Wrapper function to fit Multinomial Ecosystem State Model using
#'raster layers and export results as rasters
#'
#' @param resp A raster of response variable
#' @param preds Named list of rasters used as predictors
#' @param subsample A scalar denoting number of randomly sampled cells used for modelling. NULL for no subsampling
#' @param numStates A scalar denoting number of distributions to fit in the mixture
#' @param stateValError A description of the error distribution and link function to be used
#' in the model describing the ecosystem state value.  This can be from the \link[stats]{family}
#' specification or \code{character} scalar with the following possible values: \code{"gaussian"},
#' \code{"gamma"}, \code{"beta"}, \code{"negbinomial"}, or \code{"betabinomial"}.
#' @param transResp A function to be used for transformation of response variable
#' @param mcmcChains An integer scalar giving the number of MCMC chains to use
#'
#' @return A list containing the following components:
#' \itemize{
#' \item{\code{}}{}
#' \item{\code{}}{}
#' }
#'
#' @author Adam Klimes
#' @export
#'
fitRasterMESM <- function(resp, preds, subsample = NULL, numStates = 4, stateValError = gaussian,
                          transResp = function(x) x, mcmcChains = 2){
  library(raster) # add to package libraries?
  # rasters checking - resolution, extent, projection
  checkFun <- function(resp, preds, fun) {
    all(vapply(preds, function(x, ref) identical(fun(x), ref),
               fun(resp), FUN.VALUE = FALSE))
  }
  if (!checkFun(resp, preds, res)) stop("Resolution of rasters has to be identical")
  if (!checkFun(resp, preds, extent)) stop("Extent of rasters has to be identical")
  if (!checkFun(resp, preds, projection)) stop("Projection of rasters has to be identical")
  # data preparation
  dat <- data.frame(resp = transResp(getValues(resp)), lapply(preds, getValues))
  selID <- which(!apply(is.na(dat), 1, any))
  if (!is.null(subsample)) selID <- sample(selID, subsample)
  datSel <- dat[selID, ]
  st <- function(x, y = x) (x - mean(y)) / sd(y)
  stinv <- function(x, y) x * sd(y) + mean(y)
  datSelSt <- data.frame(resp = datSel$resp, lapply(datSel[, -1, drop = FALSE], st))
  # model preparation
  form <- formula(paste("resp ~", paste(colnames(dat)[-1], collapse = " + ")))
  predInit <- lapply(dat[, -1, drop = FALSE], function(x) rep(0.01, numStates))
  setInit <- c(list(intercept_stateVal = c(0, rep(0.01, numStates - 1)),
    intercept_statePrec = rep(2, numStates),
    intercept_stateProb = c(0, rep(0.01, numStates - 1))),
    predInit, predInit, predInit)
  names(setInit)[-(1:3)] <- paste0(rep(colnames(dat)[-1], each = 3), c("_stateVal", "_statePrec", "_stateProb"))
  # model fit
  mod <- fitMultinomialEcosystemState(
    stateValModels = form,
    stateProbModels = form,
    statePrecModels = form,
    stateValError = stateValError,
    inputData = datSelSt,
    numStates = numStates,
    mcmcChains = mcmcChains,
    setInit = setInit,
    setPriors = list(stateVal = list(int1 = "dnorm(0, 1)", pred = "dnorm(0, 10)"),
                     statePrec = list(int = "dnorm(0, 1)", pred = "dnorm(0, 10)"))
  )
  # results inverse transformation
  inverseSt.mesm <- function(mod, datOrig){
    mod$constants[-(1:3)] <- datOrig[-1]
    getInt <- function(cf, xvars, st){
      cf[1] + sum(cf[-1] * vapply(xvars, function(x) st(0, x), FUN.VALUE = 1.1))
    }
    getSlope <- function(newSlope, xvar){
      (- newSlope * st(0, xvar)) / stinv(0, xvar)
    }
    selState <- function(x, state) x[, grepl(paste0("\\[",state,"\\]$"), colnames(x)), drop = FALSE]
    auxInvSt <- function(mod, datOrig, stateType, chain) {
      Nstates <- mod$constants$numStates
      samp <- mod$mcmcSamples$samples[[chain]]
      ints <- samp[, paste0("intercept_state", stateType, "[", 1:Nstates, "]"), drop = FALSE]
      if (stateType == "Val") ints <- t(apply(ints, 1, cumsum))
      slopes <- samp[, paste0(rep(names(mod$constants)[-(1:3)], each = Nstates), "_state", stateType, "[", 1:Nstates, "]"), drop = FALSE]
      origInt <- data.frame(lapply(1:Nstates, function(curState, ints, slopes) apply(cbind(selState(ints, curState), selState(slopes, curState)), 1, getInt, xvars = datOrig[, -1, drop = FALSE], st = st), ints, slopes))
      if (stateType == "Val") origInt <- t(apply(origInt, 1, function(x) x - c(0, head(cumsum(x), -1))))
      origSlopes <- Map(getSlope, data.frame(samp[, paste0(rep(names(mod$constants)[-(1:3)], each = Nstates), "_state", stateType, "[", 1:Nstates, "]"), drop = FALSE]),
                        xvar = datOrig[, rep(2:ncol(datOrig), each = Nstates), drop = FALSE])
      mod$mcmcSamples$samples[[chain]][, paste0("intercept_state", stateType, "[", 1:Nstates, "]")] <- as.matrix(origInt)
      mod$mcmcSamples$samples[[chain]][, paste0(rep(names(mod$constants)[-(1:3)], each = Nstates), "_state", stateType, "[", 1:Nstates, "]")] <- as.matrix(data.frame(origSlopes))
      mod
    }
    for (chain in seq_along(mod$mcmcSamples$samples)){
      for (stateType in c("Val", "Prec", "Prob")){
        mod <- auxInvSt(mod, datOrig, stateType, chain)
      }
    }
    mod
  }
  modISt <- inverseSt.mesm(mod, datSel)

  # raster reconstruction - to be used for model output
  distToState <- rep(NA, nrow(dat))
  distToState[selID] <- predict.mesm(mod)$obsDat$distToState
  precar <- rep(NA, nrow(dat))
  precar[selID] <- predict.mesm(mod)$obsDat$distToTip
  dToStateR <- raster(matrix(distToState, nrow = dim(resp)[1], ncol = dim(resp)[2], byrow = TRUE), template = resp)
  precarR <- raster(matrix(precar, nrow = dim(resp)[1], ncol = dim(resp)[2], byrow = TRUE), template = resp)
  out <- list(mod = mod, modISt = modISt, dToStateR = dToStateR, precarR = precarR)
}
