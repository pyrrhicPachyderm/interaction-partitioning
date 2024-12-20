#include <math.h>
#include <numeric> //Gives std::accumulate.
#include "utils/array.hpp"
#include "rjsolver.hpp"

#define MAX_TRANS_MODEL_JUMP_PROBABILITY 0.9
#define INITIAL_APPROXIMATE_POSTERIOR_VARIANCE_MULTIPLIER 0.01
#define DESIRED_ACCEPTANCE_RATE 0.23 //For within-model jumps.
#define MAX_JUMP_VARIANCE_MULTIPLIER_CHANGE 2.0 //Controls how quickly dialIn changes the jump variance multiplier.
#define MAX_POSTERIOR_VARIANCE_CHANGE_FACTOR 100.0 //Controls how quickly dialIn changes the estimates of posterior variance.

//The probability of any particular jump.
template<typename ErrDistT> double ReversibleJumpSolver<ErrDistT>::getTransModelJumpProbability(GroupingType groupingType, MoveType moveType) const {
	return getTransModelJumpProbability(groupingType, moveType, false);
}

template<typename ErrDistT> double ReversibleJumpSolver<ErrDistT>::getTransModelJumpProbability(GroupingType groupingType, MoveType moveType, bool reverse) const {
	GroupingSizeSet sourceGroupingSizes = getGroupingSizeSet(currentGroupings);
	return transModelJumpProbabilityMultiplier * getUnscaledTransModelJumpProbability(sourceGroupingSizes, groupingType, moveType, reverse);
}

//The number of possible moves of a specified type.
template<typename ErrDistT> size_t ReversibleJumpSolver<ErrDistT>::getNumTransModelJumps(GroupingType groupingType, MoveType moveType) const {
	return moveType == MERGE ? currentGroupings[groupingType].getNumMerges() : currentGroupings[groupingType].getNumSplits();
}

template<typename ErrDistT> double ReversibleJumpSolver<ErrDistT>::getUnscaledTransModelJumpProbability(GroupingSizeSet sourceGroupingSizes, GroupingSizeSet destGroupingSizes) const {
	double sourceHyperprior = hyperprior.getDensity(sourceGroupingSizes);
	double destHyperprior = hyperprior.getDensity(destGroupingSizes);
	return std::min(1.0, destHyperprior/sourceHyperprior);
}

template<typename ErrDistT> double ReversibleJumpSolver<ErrDistT>::getUnscaledTransModelJumpProbability(GroupingSizeSet sourceGroupingSizes, GroupingType groupingType, MoveType moveType, bool reverse) const {
	GroupingSizeSet destGroupingSizes = sourceGroupingSizes;
	destGroupingSizes[groupingType] += moveType == MERGE ? -1 : 1;
	return reverse ?
		getUnscaledTransModelJumpProbability(destGroupingSizes, sourceGroupingSizes) :
		getUnscaledTransModelJumpProbability(sourceGroupingSizes, destGroupingSizes);
}

template<typename ErrDistT> double ReversibleJumpSolver<ErrDistT>::getUnscaledTransModelJumpProbability(GroupingSizeSet sourceGroupingSizes, GroupingType groupingType, MoveType moveType) const {
	return getUnscaledTransModelJumpProbability(sourceGroupingSizes, groupingType, moveType, false);
}

template<typename ErrDistT> double ReversibleJumpSolver<ErrDistT>::findUnscaledMaxTransModelJumpProbability(GroupingSizeSet groupingSizes) const {
	double result = 0.0;
	for(size_t groupingType = 0; groupingType < NUM_GROUPING_TYPES; groupingType++) {
		if(!isChangingGroupings[groupingType]) continue;
		for(size_t moveType = 0; moveType < NUM_MOVE_TYPES; moveType++) {
			size_t maxNumMoves = moveType == MERGE ?
				Grouping::getNumMerges(groupingSizes[groupingType]) :
				Grouping::getMaxNumSplits(this->data.getNumSpecies((GroupingType)groupingType), groupingSizes[groupingType]);
			result += maxNumMoves * getUnscaledTransModelJumpProbability(groupingSizes, (GroupingType)groupingType, (MoveType)moveType);
		}
	}
	return result;
}

template<typename ErrDistT> double ReversibleJumpSolver<ErrDistT>::findUnscaledMaxTransModelJumpProbability(GroupingSizeSet groupingSizes, size_t recursionLevel) const {
	//We need to loop over each type of model that's changing; this is a NUM_GROUPING_TYPES-times nested loop.
	//However, whether each loop exists is conditional on isChangingGroupings.
	//So the neatest way to do this is recursion.
	//recursionLevel is an index into a GroupingSet; the index to loop over in this recursion.
	//groupingSizes is only determined for indices below recursionLevel.
	
	if(recursionLevel == NUM_GROUPING_TYPES) {
		return findUnscaledMaxTransModelJumpProbability(groupingSizes);
	}
	
	if(!isChangingGroupings[recursionLevel]) {
		//We're not looping at this level; just set the actual grouping we're using at this level.
		groupingSizes[recursionLevel] = getGrouping((GroupingType)recursionLevel).getNumGroups();
		return findUnscaledMaxTransModelJumpProbability(groupingSizes, recursionLevel+1);
	}
	
	double result = 0.0;
	for(size_t i = 1; i <= this->data.getNumSpecies((GroupingType)recursionLevel); i++) {
		groupingSizes[recursionLevel] = i;
		result = std::max(result, findUnscaledMaxTransModelJumpProbability(groupingSizes, recursionLevel+1));
	}
	return result;
}

template<typename ErrDistT> double ReversibleJumpSolver<ErrDistT>::findTransModelJumpProbabilityMultiplier() const {
	//This is the value of c used in trans-model jump probabilities, as used in section 4.3 of Green 1995.
	
	//The GroupingIndexSet to pass the recursive function doesn't matter.
	//It's only determined for indices below recursionLevel, which is zero.
	double unscaledMaxTransModelJumpProbability = findUnscaledMaxTransModelJumpProbability(GroupingSizeSet(), 0);
	
	return MAX_TRANS_MODEL_JUMP_PROBABILITY / unscaledMaxTransModelJumpProbability;
}

template<typename ErrDistT> double ReversibleJumpSolver<ErrDistT>::guessInitialGrowthRateApproximatePosteriorVariance() const {
	return pow(this->data.guessGrowthRateMagnitude(), 2) * INITIAL_APPROXIMATE_POSTERIOR_VARIANCE_MULTIPLIER;
}

template<typename ErrDistT> double ReversibleJumpSolver<ErrDistT>::guessInitialCompetitionCoefficientApproximatePosteriorVariance() const {
	return pow(this->data.guessCompetitionCoefficientMagnitude(), 2) * INITIAL_APPROXIMATE_POSTERIOR_VARIANCE_MULTIPLIER;
}

template<typename ErrDistT> ReversibleJumpSolver<ErrDistT>::AdditionalParametersVector ReversibleJumpSolver<ErrDistT>::guessInitialAdditionalParametersApproximatePosteriorVariance() const {
	//Multiply guessInitialAdditionalParameters by INITIAL_APPROXIMATE_POSTERIOR_VARIANCE_MULTIPLIER;
	return array_map(
		[] (double i) -> double {
			return i * INITIAL_APPROXIMATE_POSTERIOR_VARIANCE_MULTIPLIER;
		},
		this->guessInitialAdditionalParameters()
	);
}

template<typename ErrDistT> double ReversibleJumpSolver<ErrDistT>::getJumpVarianceMultiplier(JumpType jumpType) const {
	if(jumpType == WITHIN_JUMP) {
		return withinModelJumpVarianceMultiplier;
	} else if(jumpType == MERGE_JUMP || jumpType == SPLIT_JUMP) {
		return transModelJumpVarianceMultiplier;
	} else __builtin_unreachable();
}

template<typename ErrDistT> Distribution<double> ReversibleJumpSolver<ErrDistT>::getGrowthRateJumpDistribution(JumpType jumpType) const {
	return Distribution<double>(new Distributions::Normal(0.0, growthRateApproximatePosteriorVariance * getJumpVarianceMultiplier(jumpType)));
}

template<typename ErrDistT> Distribution<double> ReversibleJumpSolver<ErrDistT>::getCompetitionCoefficientJumpDistribution(JumpType jumpType) const {
	return Distribution<double>(new Distributions::Normal(0.0, competitionCoefficientApproximatePosteriorVariance * getJumpVarianceMultiplier(jumpType)));
}

template<typename ErrDistT> std::array<Distribution<double>, ReversibleJumpSolver<ErrDistT>::NUM_ADDITIONAL_PARAMETERS> ReversibleJumpSolver<ErrDistT>::getAdditionalParametersJumpDistribution(JumpType jumpType) const {
	return array_map(
		[this, jumpType] (size_t index) -> Distribution<double> {
			return Distribution<double>(new Distributions::Normal(0.0, additionalParametersApproximatePosteriorVariance[index] * getJumpVarianceMultiplier(jumpType)));
		},
		make_index_array<NUM_ADDITIONAL_PARAMETERS>()
	);
}

template<typename ErrDistT> Distribution<double> ReversibleJumpSolver<ErrDistT>::getTransModelJumpDistribution(GroupingType groupingType, JumpType jumpType) const {
	if(groupingType == GROWTH) {
		return getGrowthRateJumpDistribution(jumpType);
	} else if(groupingType == ROW || groupingType == COL) {
		return getCompetitionCoefficientJumpDistribution(jumpType);
	} else __builtin_unreachable();
}

template<typename ErrDistT> double ReversibleJumpSolver<ErrDistT>::getRandomProbability() {
	//A random double in [0,1).
	return Distributions::Uniform(0,1).getRandom(randomGenerator);
}

template<typename ErrDistT> double ReversibleJumpSolver<ErrDistT>::proposeTransModelJump(GroupingType groupingType, MoveType moveType, size_t index) {
	if(moveType == MERGE) proposedJumpType = MERGE_JUMP;
	else if(moveType == SPLIT) proposedJumpType = SPLIT_JUMP;
	else __builtin_unreachable();
	
	Grouping newGrouping = moveType == MERGE ?
		currentGroupings[groupingType].getMerge(index) :
		currentGroupings[groupingType].getSplit(index);
	GroupingMove groupingMove = GroupingMove(currentGroupings[groupingType], newGrouping);
	
	proposedGroupings = currentGroupings;
	proposedGroupings[groupingType] = newGrouping;
	
	Distribution<double> randomVariableDistribution = getTransModelJumpDistribution(groupingType, proposedJumpType);
	
	proposedParameters = currentParameters;
	double logAcceptanceRatio = proposedParameters.moveModel(groupingType, moveType, groupingMove, randomVariableDistribution, randomGenerator);
	
	logAcceptanceRatio += log(getTransModelJumpProbability(groupingType, moveType, false));
	logAcceptanceRatio -= log(getTransModelJumpProbability(groupingType, moveType, true));
	
	return logAcceptanceRatio;
}

template<typename ErrDistT> double ReversibleJumpSolver<ErrDistT>::proposeWithinModelJump() {
	proposedGroupings = currentGroupings;
	proposedParameters = currentParameters;
	
	proposedJumpType = WITHIN_JUMP;
	
	proposedParameters.moveParameters(getGrowthRateJumpDistribution(proposedJumpType), getCompetitionCoefficientJumpDistribution(proposedJumpType), getAdditionalParametersJumpDistribution(proposedJumpType), randomGenerator);
	
	//TODO: If using a non-symmetric jumping density, the jumping density component of the acceptance ratio may not be 1.
	return 1.0;
}

template<typename ErrDistT> double ReversibleJumpSolver<ErrDistT>::proposeJump() {
	double selector = getRandomProbability();
	
	//Iterate over each type of trans-model jump, and see if we do it.
	//Instead of accumulating some sum of previous probabilities, we'll just decrement selector.
	for(size_t groupingType = 0; groupingType < NUM_GROUPING_TYPES; groupingType++) {
		if(!isChangingGroupings[groupingType]) continue;
		for(size_t moveType = 0; moveType < NUM_MOVE_TYPES; moveType++) {
			//Get the number of jumps, and the probability of each jump, and multiply to find the total.
			double probability = getTransModelJumpProbability((GroupingType)groupingType, (MoveType)moveType);
			size_t numJumps = getNumTransModelJumps((GroupingType)groupingType, (MoveType)moveType);
			double totalProbability = probability * numJumps;
			if(selector < totalProbability) {
				size_t index = selector / probability;
				return proposeTransModelJump((GroupingType)groupingType, (MoveType)moveType, index);
			} else {
				selector -= totalProbability;
			}
		}
	}
	
	return proposeWithinModelJump();
}

template<typename ErrDistT> void ReversibleJumpSolver<ErrDistT>::acceptJump() {
	currentGroupings = proposedGroupings;
	currentParameters = proposedParameters;
}

template<typename ErrDistT> void ReversibleJumpSolver<ErrDistT>::rejectJump() {
	return;
}

template<typename ErrDistT> double ReversibleJumpSolver<ErrDistT>::getLogLikelihoodRatio() const {
	return this->getLogLikelihood(proposedParameters, proposedGroupings) - this->getLogLikelihood(currentParameters, currentGroupings);
}

template<typename ErrDistT> double ReversibleJumpSolver<ErrDistT>::getLogPriorRatio() const {
	double logHyperpriorRatio = hyperprior.getLogDensity(proposedGroupings) - hyperprior.getLogDensity(currentGroupings);
	double logParametersPriorRatio = parametersPrior.getLogDensity(proposedParameters) - parametersPrior.getLogDensity(currentParameters);
	return logHyperpriorRatio + logParametersPriorRatio;
}

template<typename ErrDistT> bool ReversibleJumpSolver<ErrDistT>::makeJump(bool canTransModelJump) {
	double logJumpingDensityRatio = canTransModelJump ? proposeJump() : proposeWithinModelJump();
	
	double logPriorRatio = getLogPriorRatio();
	double logLikelihoodRatio = getLogLikelihoodRatio();
	double logAcceptanceRatio = logPriorRatio + logLikelihoodRatio + logJumpingDensityRatio;
	double acceptanceRatio = exp(logAcceptanceRatio);
	
	double selector = getRandomProbability();
	if(selector < acceptanceRatio) {
		acceptJump();
		return true;
	} else {
		rejectJump();
		return false;
	}
}

template<typename ErrDistT> void ReversibleJumpSolver<ErrDistT>::burnIn(size_t numJumps, bool canTransModelJump) {
	for(size_t i = 0; i < numJumps; i++) {
		makeJump(canTransModelJump);
	}
}

//A pair of functions to raise or lower a jump variance multiplier.
//They enact a maximum proportional change of MAX_JUMP_VARIANCE_MULTIPLIER_CHANGE.
//They enact no change with discrepantProportion = 0, up to full change with discrepantProportion = 1.
//These edit the jump variance in-place.
static void raiseJumpVarianceMultiplier(double &multiplier, double discrepantProportion) {
	multiplier *= 1.0 + discrepantProportion * (MAX_JUMP_VARIANCE_MULTIPLIER_CHANGE - 1.0);
}
static void lowerJumpVarianceMultiplier(double &multiplier, double discrepantProportion) {
	multiplier /= 1.0 + discrepantProportion * (MAX_JUMP_VARIANCE_MULTIPLIER_CHANGE - 1.0);
}

//A function to adjust an approximatePosteriorVariance.
//It sets currentVar equal to newVar,
//unless doing so would be a proportional change greater than MAX_POSTERIOR_VARIANCE_CHANGE_FACTOR.
//In which case it simply adjusts by MAX_POSTERIOR_VARIANCE_CHANGE_FACTOR.
//This damps oscillation in approximate posterior variance, and prevents it reaching zero.
//It edits the variable in-place.
static void adjustApproximatePosteriorVariance(double &currVar, double newVar) {
	if(newVar > currVar) {
		currVar = std::min(newVar, currVar * MAX_POSTERIOR_VARIANCE_CHANGE_FACTOR);
	} else {
		currVar = std::max(newVar, currVar / MAX_POSTERIOR_VARIANCE_CHANGE_FACTOR);
	}
}

//A function to calculate approximate posterior variances for dialing in.
//It takes a matrix where each row is a draw from the posterior, and each column is one parameter (after ungrouping).
static double getAverageColumnwiseVariance(Eigen::MatrixXd mat) {
	double sum = 0.0;
	for(auto col : mat.colwise()) {
		Eigen::VectorXd residuals = col - Eigen::VectorXd::Constant(col.size(), col.mean());
		double variance = residuals.dot(residuals) / (residuals.size() - 1);
		sum += variance;
	}
	return sum / mat.cols();
}

template<typename ErrDistT> void ReversibleJumpSolver<ErrDistT>::dialIn(size_t jumpsPerDial, size_t numDials) {
	//Dials in the jumping variances.
	//Combining two different dialling in methods here.
	//The first to balance the sizes of jumps in different variables against one another.
	//The second to balance overall sizes of jumps.
	//The first tries to estimate the variance of each posterior.
	//The second aims for a given within-model acceptance ratio.
	
	for(size_t dialIndex = 0; dialIndex < numDials; dialIndex++) {
		//For the first method.
		//We want to track the variance of each growth rate/competition coefficient (after ungrouping) separately.
		//So we build matrices for them.
		Eigen::MatrixXd growthRates = Eigen::MatrixXd(jumpsPerDial, this->data.getNumRowSpecies());
		Eigen::MatrixXd competitionCoefficients = Eigen::MatrixXd(jumpsPerDial, this->data.getNumRowSpecies() * this->data.getNumColSpecies());
		std::array<Eigen::VectorXd, NUM_ADDITIONAL_PARAMETERS> additionalParameters;
		additionalParameters.fill(Eigen::VectorXd(jumpsPerDial));
		
		//For the second method.
		Eigen::Array<double, NUM_JUMP_TYPES, 1> numProposals = Eigen::Array<double, NUM_JUMP_TYPES, 1>::Zero();
		Eigen::Array<double, NUM_JUMP_TYPES, 1> numAccepts = Eigen::Array<double, NUM_JUMP_TYPES, 1>::Zero();
		
		for(size_t j = 0; j < jumpsPerDial; j++) {
			bool accepted = makeJump();
			numProposals[proposedJumpType] += 1.0;
			if(accepted) numAccepts[proposedJumpType] += 1.0;
			
			//To get consistency, we need to ungroup the parameters.
			AugmentedParameters<NUM_ADDITIONAL_PARAMETERS> ungroupedParameters(currentParameters, currentGroupings);
			
			growthRates.row(j) = ungroupedParameters.getGrowthRates();
			competitionCoefficients.row(j) = ungroupedParameters.getCompetitionCoefficients().reshaped();
			for(size_t i = 0; i < NUM_ADDITIONAL_PARAMETERS; i++) {
				additionalParameters[i][j] = ungroupedParameters.getAdditionalParameter(i);
			}
		}
		
		adjustApproximatePosteriorVariance(growthRateApproximatePosteriorVariance, getAverageColumnwiseVariance(growthRates));
		adjustApproximatePosteriorVariance(competitionCoefficientApproximatePosteriorVariance, getAverageColumnwiseVariance(competitionCoefficients));
		for(size_t i = 0; i < NUM_ADDITIONAL_PARAMETERS; i++) {
			adjustApproximatePosteriorVariance(additionalParametersApproximatePosteriorVariance[i], getAverageColumnwiseVariance(additionalParameters[i]));
		}
		
		Eigen::Array<double, NUM_JUMP_TYPES, 1> acceptanceRates = numAccepts / numProposals;
		
		//If the acceptance rate is too high, we want to raise the withinModelJumpVarianceMultiplier, and vice versa.
		if(acceptanceRates[WITHIN_JUMP] > DESIRED_ACCEPTANCE_RATE) {
			double discrepantProportion = (acceptanceRates[WITHIN_JUMP] - DESIRED_ACCEPTANCE_RATE) / (1.0 - DESIRED_ACCEPTANCE_RATE);
			raiseJumpVarianceMultiplier(withinModelJumpVarianceMultiplier, discrepantProportion);
		} else {
			double discrepantProportion = (DESIRED_ACCEPTANCE_RATE - acceptanceRates[WITHIN_JUMP]) / DESIRED_ACCEPTANCE_RATE;
			lowerJumpVarianceMultiplier(withinModelJumpVarianceMultiplier, discrepantProportion);
		}
	}
}

template<typename ErrDistT> void ReversibleJumpSolver<ErrDistT>::resetChain() {
	currentGroupings = initialGroupings;
	currentParameters = initialParameters;
}

//Explicitly instantiate.
template class ReversibleJumpSolver<Distributions::Normal>;
template class ReversibleJumpSolver<Distributions::Gamma2>;
template class ReversibleJumpSolver<Distributions::DiscreteWrapper<Distributions::NegativeBinomial2>>;
