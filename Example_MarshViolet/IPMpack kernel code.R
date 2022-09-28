library(IPMpack)

makeIPMPmatrix





function (nEnvClass = 1, nBigMatrix = 50, minSize = -1, maxSize = 50,
          chosenCov = data.frame(covariate = 1), growObj, survObj,
          discreteTrans = 1, integrateType = "midpoint", correction = "none")
{
  if (class(growObj) == "growthObjPois" | class(growObj) ==
      "growthObjNegBin")
    print("warning: IPMs not appropriate with discrete growth processes")
  b <- minSize + c(0:nBigMatrix) * (maxSize - minSize)/nBigMatrix
  y <- 0.5 * (b[1:nBigMatrix] + b[2:(nBigMatrix + 1)])
  h <- y[2] - y[1]
  if (integrateType == "midpoint") {
    get.matrix <- t(outer(y, y, growSurv, cov = chosenCov,
                          growthObj = growObj, survObj = survObj)) * h
  }
  if (integrateType == "cumul") {
    get.matrix.cum <- t(outer(y, b, growthCum, cov = chosenCov,
                              growthObj = growObj))
    get.matrix <- get.matrix.cum[2:(nBigMatrix + 1), ] -
      get.matrix.cum[1:nBigMatrix, ]
    get.matrix <- t(t(get.matrix) * surv(size = y, cov = chosenCov,
                                         survObj = survObj))
  }
  if (correction == "constant") {
    nvals <- colSums(get.matrix, na.rm = TRUE)
    loc0 <- which(nvals == 0, arr.ind = TRUE)
    if (length(loc0) > 0) {
      print("warnings - columns that sum to 0 or that have NAs - assuming survival is along the diagonal; plot your Pmatrix to check it")
      get.matrix[, loc0] <- 0
      get.matrix[cbind(loc0, loc0)] <- surv(size = y[loc0],
                                            cov = chosenCov, survObj = survObj)
    }
    nvals <- colSums(get.matrix, na.rm = TRUE)
    get.matrix <- t((t(get.matrix)/nvals) * surv(size = y,
                                                 cov = chosenCov, survObj = survObj))
  }
  if (correction == "discretizeExtremes") {
    tooLow <- growthCum(y, b[1], cov = chosenCov, growthObj = growObj)
    tooHigh <- 1 - growthCum(y, b[length(b)], cov = chosenCov,
                             growthObj = growObj)
    get.matrix[1, ] <- get.matrix[1, ] + tooLow * surv(size = y,
                                                       cov = chosenCov, survObj = survObj)
    get.matrix[nBigMatrix, ] <- get.matrix[nBigMatrix, ] +
      tooHigh * surv(size = y, cov = chosenCov, survObj = survObj)
  }
  rc <- new("IPMmatrix", nDiscrete = 0, nEnvClass = 1,
            nBigMatrix = nBigMatrix, nrow = 1 * nBigMatrix, ncol = 1 *
              nBigMatrix, meshpoints = y, env.index = rep(1:nEnvClass,
                                                          each = nBigMatrix), names.discrete = "")
  rc[, ] <- get.matrix
  if (class(discreteTrans) == "discreteTrans") {
    nDisc <- ncol(discreteTrans@meanToCont)
    moveToDiscrete <- predict(discreteTrans@moveToDiscrete,
                              data.frame(size = y, size2 = (y * y)), type = "response")
    cont.to.cont <- get.matrix * matrix(1 - moveToDiscrete,
                                        nrow = nBigMatrix, ncol = nBigMatrix, byrow = TRUE)
    disc.to.disc <- discreteTrans@discreteTrans[1:nDisc,
                                                1:nDisc]
    disc.to.cont <- matrix(0, ncol = nDisc, nrow = nBigMatrix)
    cont.to.disc <- matrix(0, nrow = nDisc, ncol = nBigMatrix)
    for (j in 1:nDisc) {
      tmp <- dnorm(y, discreteTrans@meanToCont[j], discreteTrans@sdToCont[j]) *
        h
      if (correction == "constant")
        tmp <- tmp/sum(tmp)
      tmp[which(is.na(tmp))] <- 0
      disc.to.cont[, j] <- discreteTrans@discreteTrans["continuous",
                                                       j] * tmp
      if (discreteTrans@discreteTrans[j, "continuous"] >
          0) {
        cont.to.disc[j, ] <- surv(y, chosenCov, survObj) *
          moveToDiscrete * discreteTrans@discreteTrans[j,
                                                       "continuous"]/sum(discreteTrans@discreteTrans[1:nDisc,
                                                                                                     "continuous"])
      }
    }
    get.disc.matrix <- rbind(cbind(disc.to.disc, cont.to.disc),
                             cbind(disc.to.cont, cont.to.cont))
    rc <- new("IPMmatrix", nDiscrete = nDisc, nEnvClass = 1,
              nBigMatrix = nBigMatrix, nrow = 1 * nBigMatrix +
                nDisc, ncol = 1 * nBigMatrix + nDisc, meshpoints = y,
              env.index = rep(1:nEnvClass, each = nBigMatrix),
              names.discrete = rownames(discreteTrans@discreteTrans)[1:nDisc])
    rc[, ] <- get.disc.matrix
  }
  return(rc)
}
