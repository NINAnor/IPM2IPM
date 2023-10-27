# Constants needed in the model
#	numParticles - The number of samples from the posterior
#	numTime	- The number of time periods that you have population-level data for
#	numVitalParams - The number of vital rate parameters that are sampled per particle
#	numTraitBins - The number of trait bins used in the intregral projection approximation
#	vitalRateParams[1:numParticles, 1:numVitalParams] - A matrix of vital rate parameters sampled from
#		the posterior distribution
#
# Data needed in the model
# 	numObsIndividuals[1:numTime] - A vector of the number of individuals observed in each time period
#   obsTraitDensity[1:numTraitBins, 1:numTime] - A matrix of counts of trait values that fall within
#		  each of the trait bins in each of the time periods
#	
nimbleCode{
	# Set a prior for the weight of the different particles (samples from the posterior)
	for(partIter in 1:numParticles) {
		particlePriorWeight[partIter] <- 1.0 / numParticles 
	}
	# Sample the index for the particle
	particleIndex ~ dcat(particlePriorWeight[1:numParticles])
	
	# Calculate the projection matrix using the current particle index and the associated vital rate parameter values
	for(traitIterOne in 1:numTraitBins) {
		for(traitIterTwo in 1:numTraitBins) {
			traitProj[traitIterOne, traitIterTwo] <- projCalc(traitIterOne, traitIterTwo, vitalRateParams[particleIndex, 1:numVitalParam])
		}
		# Set a prior for the relative density of trait values in time period one
		traitDensity[traitIterOne, 1] ~ dgamma(0.001, 0.001)
	}
	# Use the projection matrix and the relative density of trait values in time period one to estimate the density values across all
	# other time periods
	for(timeIter in 2:numTime) {
		# Calculate the trait density based the projection matrix
		traitDensity[1:numTraitBins, timeIter] <- traitDensity[1:numTraitBins, timeIter - 1] %*% traitProj[1:numTraitBins, 1:numTraitBins]
		obsTraitDensity[1:numTraitBins, timeIter] ~ dmulti(traitDensity[1:numTraitBins, timeIter], numObsIndividuals[timeIter])
		# Link the total density to the number of individuals caught in that time period (provides scaling for the density function)
		totDensity[timeIter] <- sum(traitDensity[1:numTraitBins, timeIter])
		numObsIndividuals[timeIter] ~ dpois(totDensity[timeIter])
	}
	obsTraitDensity[1:numTraitBins, 1] ~ dmulti(traitDensity[1:numTraitBins, 1])
	# Link the total density to the number of individuals caught for the first time period (provides scaling for the density function)
	totDensity[1] <- sum(traitDensity[1:numTraitBins, timeIter])
	numObsIndividuals[1] ~ dpois(totDensity[timeIter])
}