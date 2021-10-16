#include <numeric> //Gives std::accumulate.
#include <random>
#include "rjsolver.hpp"

#define RANDOM_SEED 42
#define MAX_TRANS_MODEL_JUMP_PROBABILITY 0.9

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
	std::vector<double> probabilities(destIndices.size());
	
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

static double generateJump(double variance) {
	static std::default_random_engine rng(RANDOM_SEED);
	return std::normal_distribution(0.0, sqrt(variance))(rng);
}

static double getJumpDensity(double variance, double jumpSize) {
	return 1.0 / sqrt(2 * M_PI * variance) * exp(-0.5 * jumpSize*jumpSize / variance);
}

double ReversibleJumpSolver::proposeTransModelJump(GroupingType groupingType, MoveType moveType, size_t adjIndex) {
	//adjIndex is the index into the adjacency list of the relevant part of groupingLattice.
	isProposing = true;
	
	size_t newGroupingIndex = groupingLattice.getMoveDest(moveType, currentGroupings[groupingType], adjIndex);
	GroupingMove groupingMove = groupingLattice.getMove(moveType, currentGroupings[groupingType], adjIndex);
	
	proposedGroupings = currentGroupings;
	proposedGroupings[groupingType] = newGroupingIndex;
	
	double jumpVariance = groupingType == GROWTH ? growthRateJumpVariance : competitionCoefficientJumpVariance;
	RandomVariableFunc getRandomVariable = std::bind(generateJump, jumpVariance);
	RandomVariableDensityFunc getRandomVariableDensity = std::bind(getJumpDensity, jumpVariance, std::placeholders::_1);
	
	proposedParameters = currentParameters;
	double acceptanceRatio = proposedParameters.moveModel(groupingType, moveType, groupingMove, getRandomVariable, getRandomVariableDensity);
	
	return acceptanceRatio;
}
