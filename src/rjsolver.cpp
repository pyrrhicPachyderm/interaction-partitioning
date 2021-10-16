#include <numeric> //Gives std::accumulate.
#include <random>
#include "rjsolver.hpp"

#define RANDOM_SEED 42
#define MAX_TRANS_MODEL_JUMP_PROBABILITY 0.9

static std::default_random_engine randomNumberGenerator(RANDOM_SEED);

double ReversibleJumpSolver::getTransModelJumpProbability(GroupingIndexSet sourceGroupingIndices, GroupingIndexSet destGroupingIndices) const {
	return getTransModelJumpProbability(sourceGroupingIndices, destGroupingIndices, transModelJumpProbabilityMultiplier);
}

double ReversibleJumpSolver::getTransModelJumpProbability(GroupingIndexSet sourceGroupingIndices, GroupingIndexSet destGroupingIndices, double multiplier) const {
	double sourceHyperprior = hyperpriorFunc(getGroupings(sourceGroupingIndices));
	double destHyperprior = hyperpriorFunc(getGroupings(destGroupingIndices));
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

double ReversibleJumpSolver::aicHyperprior(GroupingSet groupings) {
	size_t numParameters = groupings[GROWTH].getNumGroups() + groupings[ROW].getNumGroups() * groupings[COL].getNumGroups();
	return exp(-numParameters);
}

GroupingSet ReversibleJumpSolver::getGroupings(GroupingIndexSet groupingIndices) const {
	//TODO: This should be done based on NUM_GROUPING_TYPES, but constructing a std::array is hard.
	//Perhaps look at https://stackoverflow.com/a/32175958
	return GroupingSet({groupingLattice.getGrouping(groupingIndices[GROWTH]), groupingLattice.getGrouping(groupingIndices[ROW]), groupingLattice.getGrouping(groupingIndices[COL])});
}

ReversibleJumpSolver::GroupingIndexSet ReversibleJumpSolver::getGroupingIndices(GroupingSet groupings) const {
	//TODO: Should also be done based on NUM_GROUPING_TYPES.
	return GroupingIndexSet({groupingLattice.getIndex(groupings[GROWTH]), groupingLattice.getIndex(groupings[ROW]), groupingLattice.getIndex(groupings[COL])});
}

static double getRandomProbability() {
	//A random double in [0,1).
	return std::uniform_real_distribution(0.0, 1.0)(randomNumberGenerator);
}

static double getRandomNormal(double variance) {
	return std::normal_distribution(0.0, sqrt(variance))(randomNumberGenerator);
}

static double getNormalDensity(double variance, double residual) {
	return 1.0 / sqrt(2 * M_PI * variance) * exp(-0.5 * residual*residual / variance);
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
	
	double jumpVariance = groupingType == GROWTH ? growthRateJumpVariance : competitionCoefficientJumpVariance;
	RandomVariableFunc getRandomVariable = std::bind(getRandomNormal, jumpVariance);
	RandomVariableDensityFunc getRandomVariableDensity = std::bind(getNormalDensity, jumpVariance, std::placeholders::_1);
	
	proposedParameters = currentParameters;
	double acceptanceRatio = proposedParameters.moveModel(groupingType, moveType, groupingMove, getRandomVariable, getRandomVariableDensity);
	
	acceptanceRatio *= getTransModelJumpProbability(currentGroupings, proposedGroupings) / getTransModelJumpProbability(proposedGroupings, currentGroupings);
	
	setIsProposing(true);
	return acceptanceRatio;
}

double ReversibleJumpSolver::proposeWithinModelJump() {
	proposedGroupings = currentGroupings;
	proposedParameters = currentParameters;
	
	RandomVariableFunc getGrowthRateJump = std::bind(getRandomNormal, growthRateJumpVariance);
	RandomVariableFunc getCompetitionCoefficientJump = std::bind(getRandomNormal, competitionCoefficientJumpVariance);
	std::array<RandomVariableFunc, 1> getAdditionalParameterJumps = {std::bind(getRandomNormal, varianceJumpVariance)};
	
	proposedParameters.moveParameters(getGrowthRateJump, getCompetitionCoefficientJump, getAdditionalParameterJumps);
	
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
