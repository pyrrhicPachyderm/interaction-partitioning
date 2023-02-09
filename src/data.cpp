#include "data.hpp"

size_t Data::findNumFocals(std::vector<size_t> focals) {
	//We want the number of unique focals.
	//We sort and remove duplicates, then take the length.
	std::sort(focals.begin(), focals.end());
	auto last = std::unique(focals.begin(), focals.end());
	return last - focals.begin();
}

bool Data::areFocalsFirst() const {
	//Check whether the maximum numbered focal is lower than the number of focals.
	return *std::max_element(focal.begin(), focal.end()) < numRowSpecies;
}

double Data::getResponseMean() const {
	return response.mean();
}

double Data::getResponseVariance() const {
	Eigen::VectorXd residuals = response - Eigen::VectorXd::Constant(response.size(), response.mean());
	return residuals.dot(residuals) / residuals.size();
}

double Data::guessGrowthRate() const {
	//We might assume that all the species are in one group, and that all competition coefficients are zero.
	//This gives us the average observed response.
	return response.mean();
}

double Data::guessCompetitionCoefficientMagnitude() const {
	//We will assume that with all species present at average density, growth halts.
	//This gives us 1, divided by the average density, divided by the number of species.
	return 1.0 / design.mean() / numColSpecies;
}

double Data::guessErrorVariance() const {
	//Simply return the variance of the response variable.
	return getResponseVariance();
}
