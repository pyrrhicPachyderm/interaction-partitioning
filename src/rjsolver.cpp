#include "rjsolver.hpp"

double aicHyperprior(GroupingSet groupings) {
	size_t numParameters = groupings[GROWTH].getNumGroups() + groupings[ROW].getNumGroups() * groupings[COL].getNumGroups();
	return exp(-numParameters);
}
