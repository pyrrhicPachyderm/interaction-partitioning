#include "rjsolver.hpp"

double aicHyperprior(Grouping growthGrouping, Grouping rowGrouping, Grouping colGrouping) {
	size_t numParameters = growthGrouping.getNumGroups() + rowGrouping.getNumGroups() * colGrouping.getNumGroups();
	return exp(-numParameters);
}
