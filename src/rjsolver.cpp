#include "rjsolver.hpp"

double aicHyperprior(GroupingSet groupings) {
	size_t numParameters = groupings[GROWTH].getNumGroups() + groupings[ROW].getNumGroups() * groupings[COL].getNumGroups();
	return exp(-numParameters);
}

GroupingSet ReversibleJumpSolver::getGroupings(GroupingIndexSet groupingIndices) const {
	//TODO: This should be done based on NUM_GROUPING_TYPES, but constructing a std::array is hard.
	//Perhaps look at https://stackoverflow.com/a/32175958
	return GroupingSet({groupingLattice.getGrouping(groupingIndices[GROWTH]), groupingLattice.getGrouping(groupingIndices[ROW]), groupingLattice.getGrouping(groupingIndices[COL])});
}
