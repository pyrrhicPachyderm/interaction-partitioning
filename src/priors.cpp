#include "priors.hpp"

double Hyperprior::flatFunc(const GroupingSet &groupings) {
	return 1.0;
}

double Hyperprior::aicFunc(const GroupingSet &groupings) {
	size_t numParameters = groupings[GROWTH].getNumGroups() + groupings[ROW].getNumGroups() * groupings[COL].getNumGroups();
	return exp(-1.0 * (double)numParameters);
}
