## 1. ------ DEFINE A LIST OF DISTRIBUTIONS AND CONSTANTS ------
### 1.1. ==== Initialise the distribution list ====
# Initialise a list of distributions that are defined in this file.  Note that the distributions are not
# registered just yet (that is done at the end of the source file).
distributionList <- list(
  ## 1.1.1. Define the dbetabin distribution ----
  dbetabin = list(
    # Define the BUGS code to call the distribution
    BUGSdist = "dbetabin(shape1, shape2, size, mean, prec)",
    # Define the way it will be called
    Rdist = "dbetabin(shape1 = mean * (prec * mean * (size - mean) - 1.0) / (size - prec * mean * (size - mean)), shape2 = (size - mean) * (prec * mean * (size - mean) - 1.0) / (size - prec * mean * (size - mean)), size = size)",
    # Define how the alternative parameterisations relate to the default parameterisation
    altParams = c(
      "mean = (size * shape1) / (shape1 + shape2)",
      "prec = pow(shape1 + shape2, 2.0) * (shape1 + shape2 + 1.0) / (size * shape1 * shape2 * (shape1 + shape2 + size))",
      "size = size"
    ),
    # Set the input and output types and dimension structure
    types = c(
      "value = double(0)", "shape1 = double(0)", "shape2 = double(0)",
      "size = double(0)"),
    # It is a discrete-valued distribution
    discrete = TRUE,
    # Define the cumulative probability and quantile function availability
    pqAvail = FALSE,
    # Specify the range of values
    range = c(0, Inf)
  )
  ## 1.1.2. Define the dnormStateValueMembership distribution ----
  , dnormStateValueMembership = list(
    # Define the BUGS code to call the distribution
    BUGSdist = "dnormStateValueMembership(stateVal, statePrec, stateProb)",
    # Set the input and output types and dimension structure
    types = c(
      "value = double(0)", "stateVal = double(1)", "statePrec = double(1)",
      "stateProb = double(1)"),
    # Define the cumulative probability and quantile function availability
    pqAvail = FALSE,
    mixedSizes = TRUE   # Turn off warnings about possible dimension mismatch
  )
  ## 1.1.3. Define the dgammaStateValueMembership distribution ----
  , dgammaStateValueMembership = list(
  # Define the BUGS code to call the distribution
  BUGSdist = "dgammaStateValueMembership(stateVal, statePrec, stateProb)",
    # Set the input and output types and dimension structure
    types = c(
      "value = double(0)", "stateVal = double(1)", "statePrec = double(1)",
      "stateProb = double(1)"),
    # Define the cumulative probability and quantile function availability
    pqAvail = FALSE,
    # Specify the range of values
    range = c(0, Inf),
    mixedSizes = TRUE   # Turn off warnings about possible dimension mismatch
  )
  ## 1.1.4. Define the dbetaStateValueMembership distribution ----
  , dbetaStateValueMembership = list(
    # Define the BUGS code to call the distribution
    BUGSdist = "dbetaStateValueMembership(stateVal, statePrec, stateProb)",
    # Set the input and output types and dimension structure
    types = c(
      "value = double(0)", "stateVal = double(1)", "statePrec = double(1)",
      "stateProb = double(1)"),
    # Define the cumulative probability and quantile function availability
    pqAvail = FALSE,
    # Specify the range of values
    range = c(0, 1),
    mixedSizes = TRUE   # Turn off warnings about possible dimension mismatch
  )
  ## 1.1.5. Define the dnegbinStateValueMembership distribution ----
  , dnegbinStateValueMembership = list(
    # Define the BUGS code to call the distribution
    BUGSdist = "dnegbinStateValueMembership(stateVal, statePrec, stateProb)",
    # Set the input and output types and dimension structure
    types = c(
      "value = double(0)", "stateVal = double(1)", "statePrec = double(1)",
      "stateProb = double(1)"),
    # Define the cumulative probability and quantile function availability
    pqAvail = FALSE,
    # It is a discrete-valued distribution
    discrete = TRUE,
    # Specify the range of values
    range = c(0, Inf),
    mixedSizes = TRUE   # Turn off warnings about possible dimension mismatch
  )
  ## 1.1.6. Define the dbetabinStateValueMembership distribution ----
  , dbetabinStateValueMembership = list(
    # Define the BUGS code to call the distribution
    BUGSdist = "dbetabinStateValueMembership(stateVal, statePrec, stateProb, numTrials)",
    # Set the input and output types and dimension structure
    types = c(
      "value = double(0)", "stateVal = double(1)", "statePrec = double(1)",
      "stateProb = double(1)", "numTrials = double(0)"),
    # Define the cumulative probability and quantile function availability
    pqAvail = FALSE,
    # It is a discrete-valued distribution
    discrete = TRUE,
    # Specify the range of values
    range = c(0, Inf),
    mixedSizes = TRUE   # Turn off warnings about possible dimension mismatch
  )
)
### 1.2. ==== Deregister any previously registered distributions ====
# In some versions of NIMBLE the distribution registeration proceedures can return an error if the distributions
# have already been previously registered.  This little section of code deregisters the functions to ensure no
# errors are thrown if this file is sourced more than once.
# This seems to be required from NIMBLE 0.6-12 onwards
if(exists("distributions", nimbleUserNamespace)) {
  # Find those distributions that are defined in this source file and see if they are already registered
  distributionNames <- names(distributionList)
  isDefinedDist <- distributionNames %in% nimbleUserNamespace$distributions$namesVector
  if(any(isDefinedDist)) {
    # If the distributions are registered then deregister them
    deregisterDistributions(distributionNames[isDefinedDist])
  }
}
### 1.3. ==== NIMBLE constant helper function ====
#' @title Create Scalar Constant for NIMBLE
#'
#' @description Function to create global constants that can be used in NIMBLE functions.
#' As of NIMBLE 0.6-12 there doesn't appear to be official support for global constants
#' but this instead creates functions that can be used in place of constants.
#'
#' @param inValue A scalar value to set as a constant
#' @param isDouble A logical scalar.  If \code{TRUE} then the constant is a numeric scalar
#' and if \code{FALSE} then the constant is an integer.
#'
#' @return An object of the type returned by \code{nimbleFunction} that can be called inside
#' other NIMBLE code to retrieve the value of the global constant.
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @keywords internal
#'
nimbleScalarConstant <- function(inValue, isDouble = TRUE) {
  # Define whether the constant is an integer or a double
  retType <- ifelse(as.logical(isDouble)[1], "double", "integer")
  # Create a string containing the call to create a nimble function that returns the constant
  # This will also resolve 'inValue'
  functionText <- paste(
    "nimbleFunction(run = function() {",
    paste("\treturnType(", retType, "(0))", sep = ""),
    paste("\treturn(", inValue, ")", sep = ""),
    "})", sep = "\n")
  # Parse the string and create the nimble function
  eval(parse(text = functionText))
}
### 1.4. ==== Define constants used in NIMBLE functions ====
maxInt <- nimbleScalarConstant(.Machine$integer.max, FALSE)           # Maximum integer size
maxDouble <- nimbleScalarConstant(.Machine$double.xmax, TRUE)         # Maximum double size
minNonZeroDouble <- nimbleScalarConstant(.Machine$double.xmin, TRUE)  # Minimum positive double that is non-zero

## 2. ------ DEFINE NIMBLE CUSTOM UTILITY FUNCTIONS USED IN THE NIMBLE DISTRIBUTIONS ------
# This section is for general utility functions that are shared by multiple distributions
### 2.1. ==== Parameter value testing function (scalar double) ====
#' @title Test a Scalar Input Parameter (Double)
#'
#' @description Internal function to test whether an input scalar parameter
#' has a valid value.
#'
#' @param inValue A scalar numeric value to test the value of.
#' @param valRange A numeric vector giving the range of values that the input
#' parameter can take.  If it is of length one then the value is tested whether
#' it is larger than \code{valRange}.  If it is of length two then the values of
#' \code{inValue} is tested between the values of \code{valRange}.
#'
#' @return The input value.
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @keywords internal
#'
parameterTest_scalarDouble <- nimbleFunction(
  run = function(
    inValue = double(0),
    valRange = double(1)
  ) {
    # Define the return type for the function
    returnType(double(0))
    # Sanity test the input parameters
    flagValue <- 0
    if(dim(valRange)[1] == 0) {
      # If the value range is zero-length then all values are valid
      flagValue <- 1
    } else if(dim(valRange)[1] == 1) {
      # If the value range is length one then treat the input as a minimum value
      flagValue <- inValue >= valRange[1]
    } else if(dim(valRange)[1] == 2) {
      # If the value range is length two then treat the input as an acceptable range
      flagValue <- inValue >= min(valRange) & inValue <= max(valRange)
    } else {
      stop("invalid dimensionality for the vector defining the range of values")
    }
    # Throw an error if the value is invalid
    if(flagValue == 0) {
      stop("values in the input parameter are outside of the valid range")
    }
    return(inValue)
  }
)
### 2.2. ==== Parameter value testing function (vector double) ====
#' @title Test a Vector Input Parameter (Double)
#'
#' @description Internal function to test whether an input vector parameter
#' has valid values and dimension.
#'
#' @param inValue A numeric vector to test the values of.
#' @param valRange A numeric vector giving the range of values that the input
#' parameter can take.  If it is of length one then the value is tested whether
#' it is larger than \code{valRange}.  If it is of length two then the values of
#' \code{inValue} is tested between the values of \code{valRange}.
#' @param lengthRange An integer vector giving the range of sizes the input
#' vector can take.  If it is of length one, then the vector must be of the
#' length given.  If it is of length two, then the vector length must be within
#' the given values.
#' @param lengthHandle A scalar integer.  If the value is 0 an error is given
#' if the length requirements are not met.  If the value is 1 then the values
#' are recycled to the correct length.  If the value is 2 then the vector is
#' padded with zeros to the correct length.
#'
#' @return The input value (recycled and padded if appropriate).
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @keywords internal
#'
parameterTest_vectorDouble <- nimbleFunction(
  run = function(
    inValue = double(1),
    valRange = double(1),
    lengthRange = integer(1),
    lengthHandle = integer(0, default = 0)
  ) {
    outValue <- inValue
    # Define the return type for the function
    returnType(double(1))
    # Sanity test the input parameters
    # Retrieve the length of the input parameter
    valLength <- dim(inValue)[1]
    if(dim(lengthRange)[1] == 1) {
      # Length range is only a single value so test if the input length is that value
      if(lengthHandle == 0 & valLength != lengthRange[1]) {
        stop("input parameter vector has invalid length")
      } else if(lengthHandle == 1) {
        # Recycle the elements of the vector if requested
        outValue <- numeric(length = lengthRange[1], value = inValue, recycle = TRUE)
      } else if(lengthHandle == 2) {
        # Pad the elements of the vector if requested
        outValue <- numeric(length = lengthRange[1], value = inValue, fillZeros = TRUE, recycle = FALSE)
      } else if(lengthHandle >= 3) {
        stop("invalid handler value for the parameter vector length")
      }
    } else if(dim(lengthRange)[1] == 2) {
      # Length range is a vector so test if the input length falls within those values
      if(lengthHandle == 0 & valLength < min(lengthRange) | valLength > max(lengthRange)) {
        stop("input parameter vector has invalid length")
      } else if(lengthHandle == 1) {
        # Recycle the elements of the vector if requested
        if(valLength < min(lengthRange)) {
          outValue <- numeric(length = min(lengthRange), value = inValue, recycle = TRUE)
        } else {
          outValue <- numeric(length = min(valLength, max(lengthRange)), value = inValue, recycle = TRUE)
        }
      } else if(lengthHandle == 2) {
        # Pad the elements of the vector if requested
        if(valLength < min(lengthRange)) {
          outValue <- numeric(length = min(lengthRange), value = inValue, fillZeros = TRUE, recycle = FALSE)
        } else {
          outValue <- numeric(length = min(valLength, max(lengthRange)), value = inValue, fillZeros = TRUE, recycle = FALSE)
        }
      } else if(lengthHandle >= 3) {
        stop("invalid handler value for the parameter vector length")
      }
    } else {
      # The length specification vector is itself of invalid length
      stop("length specification vector has invalid length")
    }
    # Test to ensure that the values in the vector fall within set values
    flagValue <- 0
    if(dim(valRange)[1] == 0) {
      # If the value range is zero-length then all values are valid
      flagValue <- 1
    } else if(dim(valRange)[1] == 1) {
      # If the value range is length one then treat the input as a minimum value
      flagValue <- prod(outValue >= rep(valRange[1], valLength))
    } else if(dim(valRange)[1] == 2) {
      # If the value range is length two then treat the input as an acceptable range
      flagValue <- prod(outValue >= rep(min(valRange), valLength) & outValue <= rep(max(valRange), valLength))
    } else {
      stop("invalid dimensionality for the vector defining the range of values")
    }
    # Throw an error if there are values in the input vector that are invalid
    if(flagValue == 0) {
      stop("values in the input parameter vector are outside of the valid range")
    }
    # Return the transformed input vector
    return(outValue)
  }
)
### 2.3. ==== Parameter value testing function (matrix double) ====
#' @title Test a Matrix Input Parameter (Double)
#'
#' @description Internal function to test whether an input matrix parameter
#' has valid values and dimension.
#'
#' @param inValue A numeric matrix to test the values of.
#' @param valRange A numeric vector giving the range of values that the input
#' parameter can take.  If it is of length one then the value is tested whether
#' it is larger than \code{valRange}.  If it is of length two then the values of
#' \code{inValue} is tested between the values of \code{valRange}.
#' @param lengthRow An integer vector giving the range of sizes the input
#' rows can take.  If it is of length one, then the vector must be of the
#' length given.  If it is of length two, then the vector rows must be within
#' the given values.
#' @param lengthCol An integer vector giving the range of sizes the input
#' columns can take.  If it is of length one, then the vector must be of the
#' length given.  If it is of length two, then the vector columns must be within
#' the given values.
#' @param lengthHandle A scalar integer.  If the value is 0 an error is given
#' if the length requirements are not met.  If the value is 1 then the values
#' are recycled to the correct length.  If the value is 2 then the vector is
#' padded with zeros to the correct length.
#'
#' @return The input value (recycled and padded if appropriate).
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @keywords internal
#'
parameterTest_matrixDouble <- nimbleFunction(
  run = function(
    inValue = double(2),
    valRange = double(1),
    lengthRow = integer(1),
    lengthCol = integer(1),
    lengthHandle = integer(0, default = 0)
  ) {
    # Define the return type for the function
    returnType(double(2))
    # Sanity test the input parameters
    # Retrieve the number of rows and columns of the input parameter
    valRows <- dim(inValue)[1]
    valCols <- dim(inValue)[2]
    outRows <- valRows
    outCols <- valCols
    # Ensure that the input matrix has a correct number of rows
    if(dim(lengthRow)[1] == 1) {
      outRows <- lengthRow
      if(lengthHandle == 0 & valRows != outRows) {
        stop("input parameter matrix has an invalid number of rows")
      }
    } else if(dim(lengthRow)[1] == 2) {
      outRows <- min(max(valRows, min(lengthRow)), max(lengthRow))
      if(lengthHandle == 0 & (valRows < min(lengthRow) | valRows > max(lengthRow))) {
        stop("input parameter matrix has an invalid number of rows")
      }
    } else {
      # The length specification vector is itself of invalid length
      stop("row length specification vector has invalid length")
    }
    # Ensure that the input matrix has a correct number of columns
    if(dim(lengthCol)[1] == 1) {
      outCols <- lengthCol
      if(lengthHandle == 0 & valCols != outCols) {
        stop("input parameter matrix has an invalid number of columns")
      }
    } else if(dim(lengthCol)[1] == 2) {
      outCols <- min(max(valCols, min(lengthCol)), max(lengthCol))
      if(lengthHandle == 0 & (valCols < min(lengthCol) | valRows > max(lengthCol))) {
        stop("input parameter matrix has an invalid number of columns")
      }
    } else {
      # The length specification vector is itself of invalid length
      stop("column length specification vector has invalid length")
    }
    # Create an output matrix
    # Initialise a vector of values
    matVals <- rep(0.0, outRows * outCols)
    # Iterate over each of the columns and rows in the output matrix
    for(curCol in 1:outCols) {
      for(curRow in 1:outRows) {
        # Find the corresponding value in the input matrix (recycling if neccessary)
        curVal <- inValue[(curRow - 1) %% valRows + 1, (curCol - 1) %% valCols + 1]
        # Set the value to zero if the cell coordinates fall outside the cell coordinates
        # of the input matrix and the matrix is to be padded with zeros
        if(lengthHandle == 2 & (curRow > valRows | curCol > valCols)) {
          curVal <- 0.0
        }
        # Set the relevant value in the output matrix
        matVals[curRow + (curCol - 1) * outRows] <- curVal
      }
    }
    # Initialise the output matrix (whilst testing the matrix values)
    outValue <- matrix(value = parameterTest_vectorDouble(matVals, valRange, outRows * outCols, 0), nrow = outRows, ncol = outCols, type = "double")
    # Return the transformed input matrix
    return(outValue)
  }
)

## 3. ------ DEFINE CUSTOM DISTRIBUTION FUNCTIONS ------

### 3.1. ==== Create NIMBLE functions for the beta-binomial distribution ====
## 3.1.1. dbetabin density function ----
#' @title Probability Density of the Beta-Binomial Distribution ----
#'
#' @description Calcualte the probability density of the beta-binomial distribution.
#'
#' @param x A scalar integer containing the value to calculate a density for
#' @param shape1 A scalar value containing the first shape parameter of the composite
#' beta distribution
#' @param shape2 A scalar value containing the second shape parameter of the composite
#' beta distribution
#' @param size A scalar integer containing the number of trials from the composite
#' binomial distribution
#' @param log If \code{TRUE} then return the log density
#'
#' @return A scalar containing the probability density
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
#'
dbetabin <- nimbleFunction(
  run = function(
    x = double(0),
    shape1 = double(0),
    shape2 = double(0),
    size = double(0),
    log = integer(0, default = 0)
  ) {
    # Specify the return type dimensionality
    returnType(double(0))
    # Initialise an output value
    outProb <- 0.0
    if(log) {
      outProb <- -Inf
    }
    # Ensure that the parameterisation is valid
    if(shape1 > 0.0 & shape2 > 0.0 & x >= 0 & size >= 0 & x <= size) {
      if(size == 0 & x == 0) {
        # Catch the edge-case where both n and the number of samples is 0
        outProb <- 0.0
      } else {
        # Calculate the density using the log-gamma representation of the density function
        outProb <- lgamma(size + 1.0) + lgamma(x + shape1) + lgamma(size - x + shape2) + lgamma(shape1 + shape2) -
          lgamma(x + 1.0) - lgamma(size - x + 1.0) - lgamma(size + shape1 + shape2) - lgamma(shape1) - lgamma(shape2)
      }
      if(!log) {
        outProb <- exp(outProb)
      }
    }
    return(outProb)
  }
)
## 3.1.2. dbetabin simulation function ----
#' @title Probability Density of the Beta-Binomial Distribution ----
#'
#' @description Generate a varibale from the beta-binomial distribution.
#'
#' @param n A scalar integer containing the number of variables to simulate (NIMBLE
#' currently only allows n = 1)
#' @param shape1 A scalar value containing the first shape parameter of the composite
#' beta distribution
#' @param shape2 A scalar value containing the second shape parameter of the composite
#' beta distribution
#' @param size A scalar integer containing the number of trials from the composite
#' binomial distribution
#'
#' @return A scalar containing the simulated variable
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
#'
rbetabin <- nimbleFunction(
  run = function(
    n = integer(0),
    shape1 = double(0),
    shape2 = double(0),
    size = double(0)
  ) {
    # Specify the return type of the function
    returnType(double(0))
    # Trap the situations
    if(shape1 <= 0.0 | shape2 <= 0.0 | size < 0) {
      stop("invalid parameterisation of the shape and/or size parameters")
    }
    # Ensure that only one sample is requested
    if(n <= 0) {
      stop("the number of requested samples must be above zero")
    } else if(n > 1) {
      print("this function only allows n = 1; using n = 1")
    }
    # Initialise an output value
    outVal <- 0
    if(size > 0) {
      outVal <- rbinom(1, size, rbeta(1, shape1, shape2))
    }
    return(outVal)
  }
)

### 3.2. ==== Create NIMBLE functions for the ecosystem-state value distribution (normal error) ====
## 3.2.1. dnormStateValueMembership density function ----
#' @title Probability Density of an Ecosystem-State Value Distribution (Normal Error)
#'
#' @description Calculate the probability density of an ecosystem-state value distribution with a
#' normal error distribution for the ecosystem state variable.
#'
#' @param x A scalar value containing the value of the ecosystem state variable
#' @param stateVal A numeric vector containing the mean values of the ecosystem state variable
#' associated with each alternative stable state
#' @param statePrec A numeric vector containing the precision (reciprocal of the variance) of the
#' error of the ecosystem state variable associated with each alternative stable state
#' @param stateProb A numeric vector containing the probabilities of being in each ecosystem stable
#' state (internally normalised to one)
#' @param log If \code{TRUE} then return the log density
#'
#' @return A scalar containing the probability density
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
#'
dnormStateValueMembership <- nimbleFunction(
  run = function(
    x = double(0),
    stateVal = double(1),
    statePrec = double(1),
    stateProb = double(1),
    log = integer(0, default = 0)
  ) {
    # Specify the return type dimensionality
    returnType(double(0))
    # Retrieve the number of states in the model
    numStates <- max(dim(stateVal)[1], max(dim(statePrec)[1], dim(statePrec)[1]))
    if(numStates < 2) {
      stop("invalid number of states (there must be at least two stable states)")
    }
    # Recycle the input vector if neccessary and check to ensure the values are valid
    inStateVal <- parameterTest_vectorDouble(stateVal, numeric(length = 2, value = c(-maxDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStatePrec <- parameterTest_vectorDouble(statePrec, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStateProb <- parameterTest_vectorDouble(stateProb, numeric(length = 2, value = c(0.0, maxDouble())), integer(length = 1, value = numStates), 1)
    # Calculate the total of the state probabilities
    fullProb <- sum(inStateProb[1:numStates])
    # Intialise a vector of conditional probabilities for each state
    condStateProb <- numeric(length = numStates)
    for(stateIter in 1:numStates) {
      condStateProb[stateIter] <- inStateProb[stateIter] * dnorm(x, inStateVal[stateIter], pow(inStatePrec[stateIter], -0.5), FALSE)
    }
    # Intialise the output probability
    outProb <- sum(condStateProb[1:numStates]) / fullProb
    if(log) {
      outProb <- log(outProb)
    }
    return(outProb)
  }
)
## 3.2.2. dnormStateValueMembership simulation function ----
#' @title Simulate a Variable from an Ecosystem-State Value Distribution (Normal Error)
#'
#' @description Simulate a variable according to an ecosystem-state value distribution with a
#' normal error distribution for the ecosystem state variable.
#'
#' @param n An integer scalar containing the number of variables to simulate (NIMBLE currently
#' only allows n = 1)
#' @param stateVal A numeric vector containing the mean values of the ecosystem state variable
#' associated with each alternative stable state
#' @param statePrec A numeric vector containing the precision (reciprocal of the variance) of the
#' error of the ecosystem state variable associated with each alternative stable state
#' @param stateProb A numeric vector containing the probabilities of being in each ecosystem stable
#' state (internally normalised to one)
#'
#' @return A scalar containing the simulated ecosystem state variable
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
#'
rnormStateValueMembership <- nimbleFunction(
  run = function(
    n = integer(0),
    stateVal = double(1),
    statePrec = double(1),
    stateProb = double(1)
  ) {
    # Specify the return type dimensionality
    returnType(double(0))
    # Retrieve the number of states in the model
    numStates <- max(dim(stateVal)[1], max(dim(statePrec)[1], dim(statePrec)[1]))
    if(numStates < 2) {
      stop("invalid number of states (there must be at least two stable states)")
    }
    # Ensure that only one sample is requested
    if(n <= 0) {
      stop("the number of requested samples must be above zero")
    } else if(n > 1) {
      print("this function only allows n = 1; using n = 1")
    }
    # Recycle the input vector if neccessary and check to ensure the values are valid
    inStateVal <- parameterTest_vectorDouble(stateVal, numeric(length = 2, value = c(-maxDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStatePrec <- parameterTest_vectorDouble(statePrec, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStateProb <- parameterTest_vectorDouble(stateProb, numeric(length = 2, value = c(0.0, maxDouble())), integer(length = 1, value = numStates), 1)
    # Draw a state to be in
    outState <- rcat(1, inStateProb)
    # Using this state: draw the state value
    return(rnorm(1, inStateVal[outState], pow(inStatePrec[outState], -0.5)))
  }
)

### 3.3. ==== Create NIMBLE functions for the ecosystem-state value distribution (gamma error) ====
## 3.3.1. dgammaStateValueMembership density function ----
#' @title Probability Density of an Ecosystem-State Value Distribution (Gamma Error)
#'
#' @description Calculate the probability density of an ecosystem-state value distribution with a
#' gamma error distribution for the ecosystem state variable.
#'
#' @param x A scalar value containing the value of the ecosystem state variable
#' @param stateVal A numeric vector containing the mean values of the ecosystem state variable
#' associated with each alternative stable state
#' @param statePrec A numeric vector containing the precision (reciprocal of the variance) of the
#' error of the ecosystem state variable associated with each alternative stable state
#' @param stateProb A numeric vector containing the probabilities of being in each ecosystem stable
#' state (internally normalised to one)
#' @param log If \code{TRUE} then return the log density
#'
#' @return A scalar containing the probability density
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
#'
dgammaStateValueMembership <- nimbleFunction(
  run = function(
    x = double(0),
    stateVal = double(1),
    statePrec = double(1),
    stateProb = double(1),
    log = integer(0, default = 0)
  ) {
    # Specify the return type dimensionality
    returnType(double(0))
    # Retrieve the number of states in the model
    numStates <- max(dim(stateVal)[1], max(dim(statePrec)[1], dim(statePrec)[1]))
    if(numStates < 2) {
      stop("invalid number of states (there must be at least two stable states)")
    }
    # Recycle the input vector if neccessary and check to ensure the values are valid
    inStateVal <- parameterTest_vectorDouble(stateVal, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStatePrec <- parameterTest_vectorDouble(statePrec, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStateProb <- parameterTest_vectorDouble(stateProb, numeric(length = 2, value = c(0.0, maxDouble())), integer(length = 1, value = numStates), 1)
    # Calculate the total of the state probabilities
    fullProb <- sum(inStateProb[1:numStates])
    # Reparameterise in terms of the canoconical NIMBLE parameterisation
    inShape <- inStatePrec[1:numStates] * inStateVal[1:numStates] * inStateVal[1:numStates]
    inScale <- 1.0 / (inStatePrec[1:numStates] * inStateVal[1:numStates])
    # Intialise a vector of conditional probabilities for each state
    condStateProb <- numeric(length = numStates)
    for(stateIter in 1:numStates) {
      condStateProb[stateIter] <- inStateProb[stateIter] * dgamma(x, inShape[stateIter], inScale[stateIter], log = FALSE)
    }
    # Intialise the output probability
    outProb <- sum(condStateProb[1:numStates]) / fullProb
    if(log) {
      outProb <- log(outProb)
    }
    return(outProb)
  }
)
## 3.3.2. dgammaStateValueMembership simulation function ----
#' @title Simulate a Variable from an Ecosystem-State Value Distribution (Gamma Error)
#'
#' @description Simulate a variable according to an ecosystem-state value distribution with a
#' gamma error distribution for the ecosystem state variable.
#'
#' @param n An integer scalar containing the number of variables to simulate (NIMBLE currently
#' only allows n = 1)
#' @param stateVal A numeric vector containing the mean values of the ecosystem state variable
#' associated with each alternative stable state
#' @param statePrec A numeric vector containing the precision (reciprocal of the variance) of the
#' error of the ecosystem state variable associated with each alternative stable state
#' @param stateProb A numeric vector containing the probabilities of being in each ecosystem stable
#' state (internally normalised to one)
#'
#' @return A scalar containing the simulated ecosystem state variable
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
#'
rgammaStateValueMembership <- nimbleFunction(
  run = function(
    n = integer(0),
    stateVal = double(1),
    statePrec = double(1),
    stateProb = double(1)
  ) {
    # Specify the return type dimensionality
    returnType(double(0))
    # Retrieve the number of states in the model
    numStates <- max(dim(stateVal)[1], max(dim(statePrec)[1], dim(statePrec)[1]))
    if(numStates < 2) {
      stop("invalid number of states (there must be at least two stable states)")
    }
    # Ensure that only one sample is requested
    if(n <= 0) {
      stop("the number of requested samples must be above zero")
    } else if(n > 1) {
      print("this function only allows n = 1; using n = 1")
    }
    # Recycle the input vector if neccessary and check to ensure the values are valid
    inStateVal <- parameterTest_vectorDouble(stateVal, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStatePrec <- parameterTest_vectorDouble(statePrec, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStateProb <- parameterTest_vectorDouble(stateProb, numeric(length = 2, value = c(0.0, maxDouble())), integer(length = 1, value = numStates), 1)
    # Reparameterise in terms of the canoconical NIMBLE parameterisation
    inShape <- inStatePrec[1:numStates] * inStateVal[1:numStates] * inStateVal[1:numStates]
    inScale <- 1.0 / (inStatePrec[1:numStates] * inStateVal[1:numStates])
    # Draw a state to be in
    outState <- rcat(1, inStateProb)
    # Using this state: draw the state value
    return(rgamma(1, inShape[outState], inScale[outState]))
  }
)

### 3.4. ==== Create NIMBLE functions for the ecosystem-state value distribution (beta error) ====
## 3.4.1. dbetaStateValueMembership density function ----
#' @title Probability Density of an Ecosystem-State Value Distribution (Beta Error)
#'
#' @description Calculate the probability density of an ecosystem-state value distribution with a
#' beta error distribution for the ecosystem state variable.
#'
#' @param x A scalar value containing the value of the ecosystem state variable
#' @param stateVal A numeric vector containing the mean values of the ecosystem state variable
#' associated with each alternative stable state
#' @param statePrec A numeric vector containing the precision (reciprocal of the variance) of the
#' error of the ecosystem state variable associated with each alternative stable state
#' @param stateProb A numeric vector containing the probabilities of being in each ecosystem stable
#' state (internally normalised to one)
#' @param log If \code{TRUE} then return the log density
#'
#' @return A scalar containing the probability density
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
#'
dbetaStateValueMembership <- nimbleFunction(
  run = function(
    x = double(0),
    stateVal = double(1),
    statePrec = double(1),
    stateProb = double(1),
    log = integer(0, default = 0)
  ) {
    # Specify the return type dimensionality
    returnType(double(0))
    # Retrieve the number of states in the model
    numStates <- max(dim(stateVal)[1], max(dim(statePrec)[1], dim(statePrec)[1]))
    if(numStates < 2) {
      stop("invalid number of states (there must be at least two stable states)")
    }
    # Recycle the input vector if neccessary and check to ensure the values are valid
    inStateVal <- parameterTest_vectorDouble(stateVal, numeric(length = 2, value = c(minNonZeroDouble(), 1.0 - minNonZeroDouble())), integer(length = 1, value = numStates), 1)
    inStatePrec <- parameterTest_vectorDouble(statePrec, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStateProb <- parameterTest_vectorDouble(stateProb, numeric(length = 2, value = c(0.0, maxDouble())), integer(length = 1, value = numStates), 1)
    # Calculate the total of the state probabilities
    fullProb <- sum(inStateProb[1:numStates])
    # Reparameterise in terms of the canoconical NIMBLE parameterisation
    inAlpha <- inStateVal[1:numStates] * inStateVal[1:numStates] * (1.0 - inStateVal[1:numStates]) * inStatePrec[1:numStates] - inStateVal[1:numStates]
    inBeta <- inStateVal[1:numStates] * pow(1.0 - inStateVal[1:numStates], 2.0) * inStatePrec[1:numStates] + inStateVal[1:numStates] - 1.0
    # Intialise a vector of conditional probabilities for each state
    condStateProb <- numeric(length = numStates)
    for(stateIter in 1:numStates) {
      condStateProb[stateIter] <- inStateProb[stateIter] * dbeta(x, inAlpha[stateIter], inBeta[stateIter], FALSE)
    }
    # Intialise the output probability
    outProb <- sum(condStateProb[1:numStates]) / fullProb
    if(log) {
      outProb <- log(outProb)
    }
    return(outProb)
  }
)
## 3.4.2. dbetaStateValueMembership simulation function ----
#' @title Simulate a Variable from an Ecosystem-State Value Distribution (Beta Error)
#'
#' @description Simulate a variable according to an ecosystem-state value distribution with a
#' beta error distribution for the ecosystem state variable.
#'
#' @param n An integer scalar containing the number of variables to simulate (NIMBLE currently
#' only allows n = 1)
#' @param stateVal A numeric vector containing the mean values of the ecosystem state variable
#' associated with each alternative stable state
#' @param statePrec A numeric vector containing the precision (reciprocal of the variance) of the
#' error of the ecosystem state variable associated with each alternative stable state
#' @param stateProb A numeric vector containing the probabilities of being in each ecosystem stable
#' state (internally normalised to one)
#'
#' @return A scalar containing the simulated ecosystem state variable
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
#'
rbetaStateValueMembership <- nimbleFunction(
  run = function(
    n = integer(0),
    stateVal = double(1),
    statePrec = double(1),
    stateProb = double(1)
  ) {
    # Specify the return type dimensionality
    returnType(double(0))
    # Retrieve the number of states in the model
    numStates <- max(dim(stateVal)[1], max(dim(statePrec)[1], dim(statePrec)[1]))
    if(numStates < 2) {
      stop("invalid number of states (there must be at least two stable states)")
    }
    # Ensure that only one sample is requested
    if(n <= 0) {
      stop("the number of requested samples must be above zero")
    } else if(n > 1) {
      print("this function only allows n = 1; using n = 1")
    }
    # Recycle the input vector if neccessary and check to ensure the values are valid
    inStateVal <- parameterTest_vectorDouble(stateVal, numeric(length = 2, value = c(minNonZeroDouble(), 1.0 - minNonZeroDouble())), integer(length = 1, value = numStates), 1)
    inStatePrec <- parameterTest_vectorDouble(statePrec, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStateProb <- parameterTest_vectorDouble(stateProb, numeric(length = 2, value = c(0.0, maxDouble())), integer(length = 1, value = numStates), 1)
    # Reparameterise in terms of the canoconical NIMBLE parameterisation
    inAlpha <- inStateVal[1:numStates] * inStateVal[1:numStates] * (1.0 - inStateVal[1:numStates]) * inStatePrec[1:numStates] - inStateVal[1:numStates]
    inBeta <- inStateVal[1:numStates] * pow(1.0 - inStateVal[1:numStates], 2.0) * inStatePrec[1:numStates] + inStateVal[1:numStates] - 1.0
    # Draw a state to be in
    outState <- rcat(1, inStateProb)
    # Using this state: draw the state value
    return(rbeta(1, inAlpha[outState], inBeta[outState]))
  }
)

### 3.5. ==== Create NIMBLE functions for the ecosystem-state value distribution (negative binomial error) ====
## 3.5.1. dnegbinStateValueMembership density function ----
#' @title Probability Density of an Ecosystem-State Value Distribution (Negative Binomial Error)
#'
#' @description Calculate the probability density of an ecosystem-state value distribution with a
#' negative binomial error distribution for the ecosystem state variable.
#'
#' @param x A scalar value containing the value of the ecosystem state variable
#' @param stateVal A numeric vector containing the mean values of the ecosystem state variable
#' associated with each alternative stable state
#' @param statePrec A numeric vector containing the precision (reciprocal of the variance) of the
#' error of the ecosystem state variable associated with each alternative stable state
#' @param stateProb A numeric vector containing the probabilities of being in each ecosystem stable
#' state (internally normalised to one)
#' @param log If \code{TRUE} then return the log density
#'
#' @return A scalar containing the probability density
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
#'
dnegbinStateValueMembership <- nimbleFunction(
  run = function(
    x = double(0),
    stateVal = double(1),
    statePrec = double(1),
    stateProb = double(1),
    log = integer(0, default = 0)
  ) {
    # Specify the return type dimensionality
    returnType(double(0))
    # Retrieve the number of states in the model
    numStates <- max(dim(stateVal)[1], max(dim(statePrec)[1], dim(statePrec)[1]))
    if(numStates < 2) {
      stop("invalid number of states (there must be at least two stable states)")
    }
    # Recycle the input vector if neccessary and check to ensure the values are valid
    inStateVal <- parameterTest_vectorDouble(stateVal, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStatePrec <- parameterTest_vectorDouble(statePrec, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStateProb <- parameterTest_vectorDouble(stateProb, numeric(length = 2, value = c(0.0, maxDouble())), integer(length = 1, value = numStates), 1)
    # Calculate the total of the state probabilities
    fullProb <- sum(inStateProb[1:numStates])
    # Reparameterise in terms of the canoconical NIMBLE parameterisation
    inProb <- 1.0 - inStateVal[1:numStates] * inStatePrec[1:numStates]
    inSize <- inStateVal[1:numStates] * inStateVal[1:numStates] * inStatePrec[1:numStates] / inProb[1:numStates]
    # Intialise a vector of conditional probabilities for each state
    condStateProb <- numeric(length = numStates)
    for(stateIter in 1:numStates) {
      # Test to ensure that the parameterisation is valid
      if(inProb[stateIter] < 0.0 | inProb[stateIter] > 1.0 | inSize[stateIter] < 0.0) {
        if(log) {
          return(-Inf)
        } else {
          return(0.0)
        }
      }
      # Calculate the conditional probability
      condStateProb[stateIter] <- inStateProb[stateIter] * dnbinom(x, inSize[stateIter], inProb[stateIter], FALSE)
    }
    # Intialise the output probability
    outProb <- sum(condStateProb[1:numStates]) / fullProb
    if(log) {
      outProb <- log(outProb)
    }
    return(outProb)
  }
)
## 3.5.2. dnegbinStateValueMembership simulation function ----
#' @title Simulate a Variable from an Ecosystem-State Value Distribution (Negative Binomial Error)
#'
#' @description Simulate a variable according to an ecosystem-state value distribution with a
#' negative binomial error distribution for the ecosystem state variable.
#'
#' @param n An integer scalar containing the number of variables to simulate (NIMBLE currently
#' only allows n = 1)
#' @param stateVal A numeric vector containing the mean values of the ecosystem state variable
#' associated with each alternative stable state
#' @param statePrec A numeric vector containing the precision (reciprocal of the variance) of the
#' error of the ecosystem state variable associated with each alternative stable state
#' @param stateProb A numeric vector containing the probabilities of being in each ecosystem stable
#' state (internally normalised to one)
#'
#' @return A scalar containing the simulated ecosystem state variable
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
#'
rnegbinStateValueMembership <- nimbleFunction(
  run = function(
    n = integer(0),
    stateVal = double(1),
    statePrec = double(1),
    stateProb = double(1)
  ) {
    # Specify the return type dimensionality
    returnType(double(0))
    # Retrieve the number of states in the model
    numStates <- max(dim(stateVal)[1], max(dim(statePrec)[1], dim(statePrec)[1]))
    if(numStates < 2) {
      stop("invalid number of states (there must be at least two stable states)")
    }
    # Ensure that only one sample is requested
    if(n <= 0) {
      stop("the number of requested samples must be above zero")
    } else if(n > 1) {
      print("this function only allows n = 1; using n = 1")
    }
    # Recycle the input vector if neccessary and check to ensure the values are valid
    inStateVal <- parameterTest_vectorDouble(stateVal, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStatePrec <- parameterTest_vectorDouble(statePrec, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStateProb <- parameterTest_vectorDouble(stateProb, numeric(length = 2, value = c(0.0, maxDouble())), integer(length = 1, value = numStates), 1)
    # Reparameterise in terms of the canoconical NIMBLE parameterisation
    inProb <- 1.0 - inStateVal[1:numStates] * inStatePrec[1:numStates]
    inSize <- inStateVal[1:numStates] * inStateVal[1:numStates] * inStatePrec[1:numStates] / inProb[1:numStates]
    # Draw a state to be in
    outState <- rcat(1, inStateProb)
    # Test to ensure that the parameterisation is valid
    if(inProb[outState] < 0.0 | inProb[outState] > 1.0 | inSize[outState] < 0.0) {
      stop("invalid parameter values given for the negative binomial distribution")
    }
    # Using this state: draw the state value
    return(rnbinom(1, inSize[outState], inProb[outState]))
  }
)

### 3.6. ==== Create NIMBLE functions for the ecosystem-state value distribution (beta-binomial error) ====
## 3.6.1. dbetabinStateValueMembership density function ----
#' @title Probability Density of an Ecosystem-State Value Distribution (Beta-Binomial Error)
#'
#' @description Calculate the probability density of an ecosystem-state value distribution with a
#' beta-binomial error distribution for the ecosystem state variable.
#'
#' @param x A scalar value containing the value of the ecosystem state variable
#' @param stateVal A numeric vector containing the mean values of the ecosystem state variable
#' associated with each alternative stable state
#' @param statePrec A numeric vector containing the precision (reciprocal of the variance) of the
#' error of the ecosystem state variable associated with each alternative stable state
#' @param stateProb A numeric vector containing the probabilities of being in each ecosystem stable
#' state (internally normalised to one)
#' @param numTrials A scalar containing the number of trials
#' @param log If \code{TRUE} then return the log density
#'
#' @return A scalar containing the probability density
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
#'
dbetabinStateValueMembership <- nimbleFunction(
  run = function(
    x = double(0),
    stateVal = double(1),
    statePrec = double(1),
    stateProb = double(1),
    numTrials = double(0),
    log = integer(0, default = 0)
  ) {
    # Specify the return type dimensionality
    returnType(double(0))
    # Retrieve the number of states in the model
    numStates <- max(dim(stateVal)[1], max(dim(statePrec)[1], dim(statePrec)[1]))
    if(numStates < 2) {
      stop("invalid number of states (there must be at least two stable states)")
    }
    # Recycle the input vector if neccessary and check to ensure the values are valid
    inStateVal <- parameterTest_vectorDouble(stateVal, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStatePrec <- parameterTest_vectorDouble(statePrec, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStateProb <- parameterTest_vectorDouble(stateProb, numeric(length = 2, value = c(0.0, maxDouble())), integer(length = 1, value = numStates), 1)
    # Test the number of trials input
    inNumTrials <- numTrials
    if(inNumTrials < 0) {
      if(log) {
        return(-Inf)
      } else {
        return(0.0)
      }
    }
    # Calculate the total of the state probabilities
    fullProb <- sum(inStateProb[1:numStates])
    # Reparameterise in terms of the canoconical NIMBLE parameterisation
    inAlpha <- inStateVal[1:numStates] * (inStatePrec[1:numStates] * inStateVal[1:numStates] * (inNumTrials - inStateVal[1:numStates]) - 1.0) / (inNumTrials - inStatePrec[1:numStates] * inStateVal[1:numStates] * (inNumTrials - inStateVal[1:numStates]))
    inBeta <- (inNumTrials - inStateVal[1:numStates]) * (inStatePrec[1:numStates] * inStateVal[1:numStates] * (inNumTrials - inStateVal[1:numStates]) - 1.0) / (inNumTrials - inStatePrec[1:numStates] * inStateVal[1:numStates] * (inNumTrials - inStateVal[1:numStates]))
    # Intialise a vector of conditional probabilities for each state
    condStateProb <- numeric(length = numStates)
    for(stateIter in 1:numStates) {
      condStateProb[stateIter] <- inStateProb[stateIter] * dbetabin(x, inAlpha[stateIter], inBeta[stateIter], inNumTrials, FALSE)
    }
    # Intialise the output probability
    outProb <- sum(condStateProb[1:numStates]) / fullProb
    if(log) {
      outProb <- log(outProb)
    }
    return(outProb)
  }
)
## 3.6.2. dbetabinStateValueMembership simulation function ----
#' @title Simulate a Variable from an Ecosystem-State Value Distribution (Beta-Binomial Error)
#'
#' @description Simulate a variable according to an ecosystem-state value distribution with a
#' beta-binomial error distribution for the ecosystem state variable.
#'
#' @param n An integer scalar containing the number of variables to simulate (NIMBLE currently
#' only allows n = 1)
#' @param stateVal A numeric vector containing the mean values of the ecosystem state variable
#' associated with each alternative stable state
#' @param statePrec A numeric vector containing the precision (reciprocal of the variance) of the
#' error of the ecosystem state variable associated with each alternative stable state
#' @param stateProb A numeric vector containing the probabilities of being in each ecosystem stable
#' state (internally normalised to one)
#' @param numTrials A scalar containing the number of trials
#'
#' @return A scalar containing the simulated ecosystem state variable
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
#'
rbetabinStateValueMembership <- nimbleFunction(
  run = function(
    n = integer(0),
    stateVal = double(1),
    statePrec = double(1),
    stateProb = double(1),
    numTrials = double(0)
  ) {
    # Specify the return type dimensionality
    returnType(double(0))
    # Retrieve the number of states in the model
    numStates <- max(dim(stateVal)[1], max(dim(statePrec)[1], dim(statePrec)[1]))
    if(numStates < 2) {
      stop("invalid number of states (there must be at least two stable states)")
    }
    # Ensure that only one sample is requested
    if(n <= 0) {
      stop("the number of requested samples must be above zero")
    } else if(n > 1) {
      print("this function only allows n = 1; using n = 1")
    }
    # Recycle the input vector if neccessary and check to ensure the values are valid
    inStateVal <- parameterTest_vectorDouble(stateVal, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStatePrec <- parameterTest_vectorDouble(statePrec, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStateProb <- parameterTest_vectorDouble(stateProb, numeric(length = 2, value = c(0.0, maxDouble())), integer(length = 1, value = numStates), 1)
    inNumTrials <- parameterTest_scalarDouble(numTrials, numeric(length = 2, value = c(0, maxDouble())))
    # Reparameterise in terms of the canoconical NIMBLE parameterisation
    inAlpha <- inStateVal[1:numStates] * (inStatePrec[1:numStates] * inStateVal[1:numStates] * (inNumTrials - inStateVal[1:numStates]) - 1.0) / (inNumTrials - inStatePrec[1:numStates] * inStateVal[1:numStates] * (inNumTrials - inStateVal[1:numStates]))
    inBeta <- (inNumTrials - inStateVal[1:numStates]) * (inStatePrec[1:numStates] * inStateVal[1:numStates] * (inNumTrials - inStateVal[1:numStates]) - 1.0) / (inNumTrials - inStatePrec[1:numStates] * inStateVal[1:numStates] * (inNumTrials - inStateVal[1:numStates]))
    # Draw a state to be in
    outState <- rcat(1, inStateProb)
    # Using this state: draw the state value
    return(rbetabin(1, inAlpha[outState], inBeta[outState], inNumTrials))
  }
)

## 4. ------ REGISTER CUSTOM DISTRIBUTION FUNCTIONS ------
# Register the distributions with NIMBLE (the distribution list is defined at the top of the source file)
registerDistributions(distributionList)

## 5. ------ ADD IN DEBUG INFORMATION ------
# If the 'debugNimbleFunctions' flag is set then add debug information to each of the functions.  This is to aid
# development of the functions and stop processing during calculation steps in NIMBLE
if(get0("debugNimbleFunctions", ifnotfound = FALSE)) {
  debug(parameterTest_scalarDouble)
  debug(parameterTest_vectorDouble)
  debug(parameterTest_matrixDouble)
  debug(dnegbin)
  debug(rnegbin)
  debug(dnormStateValueMembership)
  debug(rnormStateValueMembership)
  debug(dgammaStateValueMembership)
  debug(rgammaStateValueMembership)
  debug(dbetaStateValueMembership)
  debug(rbetaStateValueMembership)
  debug(dnegbinStateValueMembership)
  debug(rnegbinStateValueMembership)
  debug(dbetabinStateValueMembership)
  debug(rbetabinStateValueMembership)
}
