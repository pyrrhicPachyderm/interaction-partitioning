#include "priors.hpp"

double Hyperprior::flatFunc(const GroupingSizeSet &groupingSizes) {
	return 1.0;
}

double Hyperprior::aicFunc(const GroupingSizeSet &groupingSizes) {
	size_t numParameters = groupingSizes[GROWTH] + groupingSizes[ROW] * groupingSizes[COL];
	return exp(-1.0 * (double)numParameters);
}
