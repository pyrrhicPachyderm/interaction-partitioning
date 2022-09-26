#include <numeric> //Gives std::accumulate.
#include <random>
#include "utils/array.hpp"
#include "rjsolver.hpp"

#define RANDOM_SEED 42
#define MAX_TRANS_MODEL_JUMP_PROBABILITY 0.9
#define DESIRED_ACCEPTANCE_RATE 0.23 //For within-model jumps.
#define MAX_JUMP_VARIANCE_MULTIPLIER_CHANGE 2.0 //Controls how quickly dialIn changes the jump variance multiplier.

static std::default_random_engine randomNumberGenerator(RANDOM_SEED);

double ReversibleJumpSolver::getTransModelJumpProbability(GroupingIndexSet sourceGroupingIndices, GroupingIndexSet destGroupingIndices) const {
	return getTransModelJumpProbability(sourceGroupingIndices, destGroupingIndices, transModelJumpProbabilityMultiplier);
}

double ReversibleJumpSolver::getTransModelJumpProbability(GroupingIndexSet sourceGroupingIndices, GroupingIndexSet destGroupingIndices, double multiplier) const {
	double sourceHyperprior = hyperprior.getDensity(getGroupings(sourceGroupingIndices));
	double destHyperprior = hyperprior.getDensity(getGroupings(destGroupingIndices));
	return multiplier * std::min(1.0, destHyperprior/sourceHyperprior);
}

std::vector<double> ReversibleJumpSolver::getTransModelJumpProbabilities(GroupingType groupingType, MoveType moveType) const {
	return getTransModelJumpProbabilities(groupingType, moveType, currentGroupings, transModelJumpProbabilityMultiplier);
}

std::vector<double> ReversibleJumpSolver::getTransModelJumpProbabilities(GroupingType groupingType, MoveType moveType, GroupingIndexSet groupingIndices, double multiplier) const {
	const std::vector<size_t> &destIndices = groupingLattice.getMoveDests(moveType, groupingIndices[groupingType]);
	std::vector<double> probabilities;
	
	GroupingIndexSet newGroupingIndices = groupingIndices;
	for(size_t i = 0; i < destIndices.size(); i++) {
		newGroupingIndices[groupingType] = destIndices[i];
		probabilities.push_back(getTransModelJumpProbability(groupingIndices, newGroupingIndices, multiplier));
	}
	
	return probabilities;
}

double ReversibleJumpSolver::getUnscaledTotalTransModelJumpProbability(GroupingIndexSet groupingIndices) const {
	double result = 0.0;
	for(size_t groupingType = 0; groupingType < NUM_GROUPING_TYPES; groupingType++) {
		if(!isChangingGroupings[groupingType]) continue;
		for(size_t moveType = 0; moveType < NUM_MOVE_TYPES; moveType++) {
			std::vector<double> probabilities = getTransModelJumpProbabilities((GroupingType)groupingType, (MoveType)moveType, groupingIndices, 1.0);
			result += std::accumulate(probabilities.begin(), probabilities.end(), 0.0);
		}
	}
	return result;
}

double ReversibleJumpSolver::getUnscaledMaxTransModelJumpProbability(size_t recursionLevel, GroupingIndexSet groupingIndices) const {
	//We need to loop over each type of model that's changing; this is a NUM_GROUPING_TYPES-times nested loop.
	//However, whether each loop exists is conditional on isChangingGroupings.
	//So the neatest way to do this is recursion.
	//recursionLevel is an index into a GroupingSet; the index to loop over in this recursion.
	//groupingIndices is only determined for indices below recursionLevel.
	
	if(recursionLevel == NUM_GROUPING_TYPES) {
		return getUnscaledTotalTransModelJumpProbability(groupingIndices);
	}
	
	if(!isChangingGroupings[recursionLevel]) {
		//We're not looping at this level; just set the actual grouping we're using at this level.
		groupingIndices[recursionLevel] = currentGroupings[recursionLevel];
		return getUnscaledMaxTransModelJumpProbability(recursionLevel+1, groupingIndices);
	}
	
	double result = 0.0;
	for(size_t i = 0; i < groupingLattice.getNumGroupings(); i++) {
		groupingIndices[recursionLevel] = i;
		result = std::max(result, getUnscaledMaxTransModelJumpProbability(recursionLevel+1, groupingIndices));
	}
	return result;
}

double ReversibleJumpSolver::getTransModelJumpProbabilityMultiplier() const {
	//This is the value of c used in trans-model jump probabilities, as given defined in Green 1995 and the report.
	
	//The GroupingIndexSet to pass the recursive function doesn't matter.
	//It's only determined for indices below recursionLevel, which is zero.
	double unscaledMaxTransModelJumpProbability = getUnscaledMaxTransModelJumpProbability(0, GroupingIndexSet());
	
	return MAX_TRANS_MODEL_JUMP_PROBABILITY / unscaledMaxTransModelJumpProbability;
}

GroupingSet ReversibleJumpSolver::getGroupings(GroupingIndexSet groupingIndices) const {
	return array_map(
		[this, groupingIndices] (size_t index) -> Grouping {
			return groupingLattice.getGrouping(groupingIndices[index]);
		},
		make_index_array<NUM_GROUPING_TYPES>()
	);
}

ReversibleJumpSolver::GroupingIndexSet ReversibleJumpSolver::getGroupingIndices(GroupingSet groupings) const {
	return array_map(
		[this, groupings] (size_t index) -> size_t {
			return groupingLattice.getIndex(groupings[index]);
		},
		make_index_array<NUM_GROUPING_TYPES>()
	);
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

static double getVariance(std::vector<double> nums) {
	Eigen::Map<Eigen::VectorXd> vec(&nums[0], nums.size());
	Eigen::VectorXd residuals = vec - Eigen::VectorXd::Constant(vec.size(), vec.mean());
	return residuals.dot(residuals) / residuals.size();
}

static double getRandomProbability() {
	//A random double in [0,1).
	return Distributions::Uniform(0,1).getRandom();
}

//We need to use this, rather than direct setting, to ensure it also dirties the groupings.
//Further, this requires that this function is called after any changes to current or proposed groupings are performed.
void ReversibleJumpSolver::setIsProposing(bool b) {
	if(isProposing == b) return;
	isProposing = b;
	for(size_t i = 0; i < NUM_GROUPING_TYPES; i++) {
		if(currentGroupings[i] != proposedGroupings[i]) dirtyGrouping((GroupingType)i);
	}
}

double ReversibleJumpSolver::proposeTransModelJump(GroupingType groupingType, MoveType moveType, size_t adjIndex) {
	//adjIndex is the index into the adjacency list of the relevant part of groupingLattice.
	size_t newGroupingIndex = groupingLattice.getMoveDest(moveType, currentGroupings[groupingType], adjIndex);
	GroupingMove groupingMove = groupingLattice.getMove(moveType, currentGroupings[groupingType], adjIndex);
	
	proposedGroupings = currentGroupings;
	proposedGroupings[groupingType] = newGroupingIndex;
	
	Distribution<double> randomVariableDistribution = getTransModelJumpDistribution(groupingType);;
	
	proposedParameters = currentParameters;
	double acceptanceRatio = proposedParameters.moveModel(groupingType, moveType, groupingMove, randomVariableDistribution);
	
	acceptanceRatio *= getTransModelJumpProbability(currentGroupings, proposedGroupings) / getTransModelJumpProbability(proposedGroupings, currentGroupings);
	
	if(moveType == MERGE) proposedJumpType = MERGE_JUMP;
	else if(moveType == SPLIT) proposedJumpType = SPLIT_JUMP;
	else __builtin_unreachable();
	setIsProposing(true);
	
	return acceptanceRatio;
}

double ReversibleJumpSolver::proposeWithinModelJump() {
	proposedGroupings = currentGroupings;
	proposedParameters = currentParameters;
	
	proposedParameters.moveParameters(getGrowthRateJumpDistribution(), getCompetitionCoefficientJumpDistribution(), getAdditionalParametersJumpDistribution());
	
	proposedJumpType = WITHIN_JUMP;
	setIsProposing(true);
	
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
			std::vector<double> probabilities = getTransModelJumpProbabilities((GroupingType)groupingType, (MoveType)moveType);
			for(size_t adjIndex = 0; adjIndex < probabilities.size(); adjIndex++) {
				selector -= probabilities[adjIndex];
				if(selector < 0) {
					return proposeTransModelJump((GroupingType)groupingType, (MoveType)moveType, adjIndex);
				}
			}
		}
	}
	
	return proposeWithinModelJump();
}

void ReversibleJumpSolver::acceptJump() {
	currentGroupings = proposedGroupings;
	currentParameters = proposedParameters;
	setIsProposing(false);
}

void ReversibleJumpSolver::rejectJump() {
	setIsProposing(false);
}

Distribution<double> ReversibleJumpSolver::getErrorDistribution() const {
	//This uses the additional parameter, and is used in calculating the likelihood.
	double errorVariance = getParameters().getAdditionalParameter(0);
	return Distribution(new Distributions::Normal(0, errorVariance));
}

Eigen::VectorXd ReversibleJumpSolver::getResiduals() {
	return Solver::getResiduals(isProposing ? (Parameters)proposedParameters : (Parameters)currentParameters);
}

double ReversibleJumpSolver::getLikelihoodRatio(Eigen::VectorXd sourceResiduals, Eigen::VectorXd destResiduals, Distribution<double> sourceErrorDistribution,  Distribution<double> destErrorDistribution) {
	//Having a getLikelihood() function and calling it for source and destination
	//just results in it returning 0 twice, as it multiplies several hundred small numbers together, and gets too small.
	//So we must calculate the likelihood ratio, getting the ratio as we go.
	double likelihoodRatio = 1.0;
	for(size_t i = 0; i < (size_t)sourceResiduals.size(); i++) {
		double sourceLikelihood = sourceErrorDistribution.getDensity(sourceResiduals[i]);
		double destLikelihood = destErrorDistribution.getDensity(destResiduals[i]);
		likelihoodRatio *= destLikelihood / sourceLikelihood;
	}
	return likelihoodRatio;
}

double ReversibleJumpSolver::getPriorDensity() const {
	double hyperpriorDensity = hyperprior.getDensity(getGroupings());
	
	//TODO: There should also be the prior of all the parameters here.
	//Currently, the competition coefficients have a prior.
	//But I'm trying some sort of improper flat prior for the growth rates; this is only really valid as long as their are no trans-model jumps changing the number of growth rates.
	//Lastly, some prior on the variance might be suitable.
	double parameterPriorDensity = 1.0;
	
	Eigen::MatrixXd competitionCoefficients = getParameters().getCompetitionCoefficients();
	for(size_t i = 0; i < (size_t)competitionCoefficients.rows(); i++) {
		for(size_t j = 0; j < (size_t)competitionCoefficients.cols(); j++) {
			parameterPriorDensity *= Distributions::Normal(0.0, competitionCoefficientPriorVariance).getDensity(competitionCoefficients(i,j));
		}
	}
	
	return hyperpriorDensity * parameterPriorDensity;
}

bool ReversibleJumpSolver::makeJump(bool canTransModelJump) {
	double sourcePrior = getPriorDensity();
	Eigen::VectorXd sourceResiduals = getResiduals();
	Distribution<double> sourceErrorDistribution = getErrorDistribution();
	
	double jumpingDensityRatio = canTransModelJump ? proposeJump() : proposeWithinModelJump();
	
	double destPrior = getPriorDensity();
	Eigen::VectorXd destResiduals = getResiduals();
	Distribution<double> destErrorDistribution = getErrorDistribution();
	
	double priorRatio = destPrior / sourcePrior;
	double likelihoodRatio = getLikelihoodRatio(sourceResiduals, destResiduals, sourceErrorDistribution, destErrorDistribution);
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


void ReversibleJumpSolver::dialIn(size_t jumpsPerDial, size_t numDials) {
	//Dials in the jumping variances.
	//Combining two different dialling in methods here.
	//The first to balance the sizes of jumps in different variables against one another.
	//The second to balance overall sizes of jumps.
	//The first tries to estimate the variance of each posterior.
	//The second aims for a given within-model acceptance ratio, and tries to equalise merge and split acceptance ratios.
	
	for(size_t i = 0; i < numDials; i++) {
		//For the first method.
		std::vector<double> growthRates;
		std::vector<double> competitionCoefficients;
		std::array<std::vector<double>, NUM_ADDITIONAL_PARAMETERS> additionalParameters;
		
		//For the second method.
		Eigen::Array<double, NUM_JUMP_TYPES, 1> numProposals = Eigen::Array<double, NUM_JUMP_TYPES, 1>::Zero();
		Eigen::Array<double, NUM_JUMP_TYPES, 1> numAccepts = Eigen::Array<double, NUM_JUMP_TYPES, 1>::Zero();
		
		for(size_t j = 0; j < jumpsPerDial; j++) {
			bool accepted = makeJump();
			numProposals[proposedJumpType] += 1.0;
			if(accepted) numAccepts[proposedJumpType] += 1.0;
			
			//We can't rely on there being a certain number of growth rates or competition coefficients.
			//And it's altogether too much work to tally multiple sets of growth rates or competition coefficients, and take the variance of each set.
			//So we just take the first of each.
			growthRates.push_back(getParameters().getGrowthRate(0));
			competitionCoefficients.push_back(getParameters().getCompetitionCoefficient(0,0));
			for(size_t i = 0; i < NUM_ADDITIONAL_PARAMETERS; i++) {
				additionalParameters[i].push_back(getParameters().getAdditionalParameter(i));
			}
		}
		growthRateApproximatePosteriorVariance = getVariance(growthRates);
		competitionCoefficientApproximatePosteriorVariance = getVariance(competitionCoefficients);
		for(size_t i = 0; i < NUM_ADDITIONAL_PARAMETERS; i++) {
			additionalParametersApproximatePosteriorVariance[i] = getVariance(additionalParameters[i]);
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
