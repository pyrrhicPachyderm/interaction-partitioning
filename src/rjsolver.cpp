#include <numeric> //Gives std::accumulate.
#include "utils/array.hpp"
#include "rjsolver.hpp"

#define MAX_TRANS_MODEL_JUMP_PROBABILITY 0.9
#define DESIRED_ACCEPTANCE_RATE 0.23 //For within-model jumps.
#define MAX_JUMP_VARIANCE_MULTIPLIER_CHANGE 2.0 //Controls how quickly dialIn changes the jump variance multiplier.

//The probability of any particular jump.
double ReversibleJumpSolver::getTransModelJumpProbability(GroupingType groupingType, MoveType moveType) const {
	return getTransModelJumpProbability(groupingType, moveType, false);
}

double ReversibleJumpSolver::getTransModelJumpProbability(GroupingType groupingType, MoveType moveType, bool reverse) const {
	GroupingSizeSet sourceGroupingSizes = getGroupingSizeSet(currentGroupings);
	return transModelJumpProbabilityMultiplier * getUnscaledTransModelJumpProbability(sourceGroupingSizes, groupingType, moveType, reverse);
}

//The number of possible moves of a specified type.
size_t ReversibleJumpSolver::getNumTransModelJumps(GroupingType groupingType, MoveType moveType) const {
	return moveType == MERGE ? currentGroupings[groupingType].getNumMerges() : currentGroupings[groupingType].getNumSplits();
}

double ReversibleJumpSolver::getUnscaledTransModelJumpProbability(GroupingSizeSet sourceGroupingSizes, GroupingSizeSet destGroupingSizes) const {
	double sourceHyperprior = hyperprior.getDensity(sourceGroupingSizes);
	double destHyperprior = hyperprior.getDensity(destGroupingSizes);
	return std::min(1.0, destHyperprior/sourceHyperprior);
}

double ReversibleJumpSolver::getUnscaledTransModelJumpProbability(GroupingSizeSet sourceGroupingSizes, GroupingType groupingType, MoveType moveType, bool reverse) const {
	GroupingSizeSet destGroupingSizes = sourceGroupingSizes;
	destGroupingSizes[groupingType] += moveType == MERGE ? -1 : 1;
	return reverse ?
		getUnscaledTransModelJumpProbability(destGroupingSizes, sourceGroupingSizes) :
		getUnscaledTransModelJumpProbability(sourceGroupingSizes, destGroupingSizes);
}

double ReversibleJumpSolver::getUnscaledTransModelJumpProbability(GroupingSizeSet sourceGroupingSizes, GroupingType groupingType, MoveType moveType) const {
	return getUnscaledTransModelJumpProbability(sourceGroupingSizes, groupingType, moveType, false);
}

double ReversibleJumpSolver::getUnscaledMaxTransModelJumpProbability(GroupingSizeSet groupingSizes) const {
	double result = 0.0;
	for(size_t groupingType = 0; groupingType < NUM_GROUPING_TYPES; groupingType++) {
		if(!isChangingGroupings[groupingType]) continue;
		for(size_t moveType = 0; moveType < NUM_MOVE_TYPES; moveType++) {
			size_t maxNumMoves = moveType == MERGE ?
				Grouping::getNumMerges(groupingSizes[groupingType]) :
				Grouping::getMaxNumSplits(data.getNumSpecies((GroupingType)groupingType), groupingSizes[groupingType]);
			result += maxNumMoves * getUnscaledTransModelJumpProbability(groupingSizes, (GroupingType)groupingType, (MoveType)moveType);
		}
	}
	return result;
}

double ReversibleJumpSolver::getUnscaledMaxTransModelJumpProbability(GroupingSizeSet groupingSizes, size_t recursionLevel) const {
	//We need to loop over each type of model that's changing; this is a NUM_GROUPING_TYPES-times nested loop.
	//However, whether each loop exists is conditional on isChangingGroupings.
	//So the neatest way to do this is recursion.
	//recursionLevel is an index into a GroupingSet; the index to loop over in this recursion.
	//groupingSizes is only determined for indices below recursionLevel.
	
	if(recursionLevel == NUM_GROUPING_TYPES) {
		return getUnscaledMaxTransModelJumpProbability(groupingSizes);
	}
	
	if(!isChangingGroupings[recursionLevel]) {
		//We're not looping at this level; just set the actual grouping we're using at this level.
		groupingSizes[recursionLevel] = getGrouping((GroupingType)recursionLevel).getNumGroups();
		return getUnscaledMaxTransModelJumpProbability(groupingSizes, recursionLevel+1);
	}
	
	double result = 0.0;
	for(size_t i = 1; i <= data.getNumSpecies((GroupingType)recursionLevel); i++) {
		groupingSizes[recursionLevel] = i;
		result = std::max(result, getUnscaledMaxTransModelJumpProbability(groupingSizes, recursionLevel+1));
	}
	return result;
}

double ReversibleJumpSolver::getTransModelJumpProbabilityMultiplier() const {
	//This is the value of c used in trans-model jump probabilities, as used in section 4.3 of Green 1995.
	
	//The GroupingIndexSet to pass the recursive function doesn't matter.
	//It's only determined for indices below recursionLevel, which is zero.
	double unscaledMaxTransModelJumpProbability = getUnscaledMaxTransModelJumpProbability(GroupingSizeSet(), 0);
	
	return MAX_TRANS_MODEL_JUMP_PROBABILITY / unscaledMaxTransModelJumpProbability;
}

Distribution<double> ReversibleJumpSolver::getGrowthRateJumpDistribution() const {
	return Distribution<double>(new Distributions::Normal(0.0, growthRateApproximatePosteriorVariance * jumpVarianceMultiplier));
}

Distribution<double> ReversibleJumpSolver::getCompetitionCoefficientJumpDistribution() const {
	return Distribution<double>(new Distributions::Normal(0.0, competitionCoefficientApproximatePosteriorVariance * jumpVarianceMultiplier));
}

std::array<Distribution<double>, ReversibleJumpSolver::NUM_ADDITIONAL_PARAMETERS> ReversibleJumpSolver::getAdditionalParametersJumpDistribution() const {
	return array_map(
		[this] (size_t index) -> Distribution<double> {
			return Distribution<double>(new Distributions::Normal(0.0, additionalParametersApproximatePosteriorVariance[index] * jumpVarianceMultiplier));
		},
		make_index_array<NUM_ADDITIONAL_PARAMETERS>()
	);
}

Distribution<double> ReversibleJumpSolver::getTransModelJumpDistribution(GroupingType groupingType) const {
	if(groupingType == GROWTH) {
		return getGrowthRateJumpDistribution();
	} else if(groupingType == ROW || groupingType == COL) {
		return getCompetitionCoefficientJumpDistribution();
	} else __builtin_unreachable();
}

static double getRandomProbability() {
	//A random double in [0,1).
	return Distributions::Uniform(0,1).getRandom();
}

double ReversibleJumpSolver::proposeTransModelJump(GroupingType groupingType, MoveType moveType, size_t index) {
	Grouping newGrouping = moveType == MERGE ?
		currentGroupings[groupingType].getMerge(index) :
		currentGroupings[groupingType].getSplit(index);
	GroupingMove groupingMove = GroupingMove(currentGroupings[groupingType], newGrouping);
	
	proposedGroupings = currentGroupings;
	proposedGroupings[groupingType] = newGrouping;
	
	Distribution<double> randomVariableDistribution = getTransModelJumpDistribution(groupingType);
	
	proposedParameters = currentParameters;
	double acceptanceRatio = proposedParameters.moveModel(groupingType, moveType, groupingMove, randomVariableDistribution);
	
	acceptanceRatio *= getTransModelJumpProbability(groupingType, moveType, false)
		/ getTransModelJumpProbability(groupingType, moveType, true);
	
	if(moveType == MERGE) proposedJumpType = MERGE_JUMP;
	else if(moveType == SPLIT) proposedJumpType = SPLIT_JUMP;
	else __builtin_unreachable();
	
	return acceptanceRatio;
}

double ReversibleJumpSolver::proposeWithinModelJump() {
	proposedGroupings = currentGroupings;
	proposedParameters = currentParameters;
	
	proposedParameters.moveParameters(getGrowthRateJumpDistribution(), getCompetitionCoefficientJumpDistribution(), getAdditionalParametersJumpDistribution());
	
	proposedJumpType = WITHIN_JUMP;
	
	//TODO: If using a non-symmetric jumping density, the jumping density component of the acceptance ratio may not be 1.
	return 1.0;
}

double ReversibleJumpSolver::proposeJump() {
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

void ReversibleJumpSolver::acceptJump() {
	currentGroupings = proposedGroupings;
	currentParameters = proposedParameters;
}

void ReversibleJumpSolver::rejectJump() {
	return;
}

Distribution<double> ReversibleJumpSolver::getErrorDistribution(const AugmentedParameters<NUM_ADDITIONAL_PARAMETERS> &parameters) const {
	//This uses the additional parameter, and is used in calculating the likelihood.
	double errorVariance = parameters.getAdditionalParameter(0);
	return Distribution(new Distributions::Normal(0, errorVariance));
}

double ReversibleJumpSolver::getLikelihoodRatio() {
	//Having a getLikelihood() function and calling it for source and destination
	//just results in it returning 0 twice, as it multiplies several hundred small numbers together, and gets too small.
	//So we must calculate the likelihood ratio, getting the ratio as we go.
	Eigen::VectorXd sourceResiduals = getResiduals(currentParameters, currentGroupings);
	Eigen::VectorXd destResiduals = getResiduals(proposedParameters, proposedGroupings);
	Distribution<double> sourceErrorDistribution = getErrorDistribution(currentParameters);
	Distribution<double> destErrorDistribution = getErrorDistribution(proposedParameters);
	
	double likelihoodRatio = 1.0;
	for(size_t i = 0; i < (size_t)sourceResiduals.size(); i++) {
		double sourceLikelihood = sourceErrorDistribution.getDensity(sourceResiduals[i]);
		double destLikelihood = destErrorDistribution.getDensity(destResiduals[i]);
		likelihoodRatio *= destLikelihood / sourceLikelihood;
	}
	return likelihoodRatio;
}

double ReversibleJumpSolver::getPriorRatio() const {
	double hyperpriorRatio = hyperprior.getDensity(proposedGroupings) / hyperprior.getDensity(currentGroupings);
	double parametersPriorRatio = parametersPrior.getDensity(proposedParameters) / parametersPrior.getDensity(currentParameters);
	return hyperpriorRatio * parametersPriorRatio;
}

bool ReversibleJumpSolver::makeJump(bool canTransModelJump) {
	double jumpingDensityRatio = canTransModelJump ? proposeJump() : proposeWithinModelJump();
	
	double priorRatio = getPriorRatio();
	double likelihoodRatio = getLikelihoodRatio();
	double acceptanceRatio = priorRatio * likelihoodRatio * jumpingDensityRatio;
	
	
	double selector = getRandomProbability();
	if(selector < acceptanceRatio) {
		acceptJump();
		return true;
	} else {
		rejectJump();
		return false;
	}
}

void ReversibleJumpSolver::burnIn(size_t numJumps, bool canTransModelJump) {
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

void ReversibleJumpSolver::dialIn(size_t jumpsPerDial, size_t numDials) {
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
		Eigen::MatrixXd growthRates = Eigen::MatrixXd(jumpsPerDial, data.getNumRowSpecies());
		Eigen::MatrixXd competitionCoefficients = Eigen::MatrixXd(jumpsPerDial, data.getNumRowSpecies() * data.getNumColSpecies());
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
		growthRateApproximatePosteriorVariance = getAverageColumnwiseVariance(growthRates);
		competitionCoefficientApproximatePosteriorVariance = getAverageColumnwiseVariance(competitionCoefficients);
		for(size_t i = 0; i < NUM_ADDITIONAL_PARAMETERS; i++) {
			additionalParametersApproximatePosteriorVariance[i] = getAverageColumnwiseVariance(additionalParameters[i]);
		}
		
		Eigen::Array<double, NUM_JUMP_TYPES, 1> acceptanceRates = numAccepts / numProposals;
		
		//If the acceptance rate is too high, we want to raise the withinModelJumpVarianceMultiplier, and vice versa.
		if(acceptanceRates[WITHIN_JUMP] > DESIRED_ACCEPTANCE_RATE) {
			double discrepantProportion = (acceptanceRates[WITHIN_JUMP] - DESIRED_ACCEPTANCE_RATE) / (1.0 - DESIRED_ACCEPTANCE_RATE);
			raiseJumpVarianceMultiplier(jumpVarianceMultiplier, discrepantProportion);
		} else {
			double discrepantProportion = (DESIRED_ACCEPTANCE_RATE - acceptanceRates[WITHIN_JUMP]) / DESIRED_ACCEPTANCE_RATE;
			lowerJumpVarianceMultiplier(jumpVarianceMultiplier, discrepantProportion);
		}
	}
}

void ReversibleJumpSolver::resetChain() {
	currentGroupings = initialGroupings;
	currentParameters = initialParameters;
}
