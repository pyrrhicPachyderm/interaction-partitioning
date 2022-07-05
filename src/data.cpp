#include "data.hpp"

double Data::getResponseMean() const {
	return response.mean();
}

double Data::getResponseVariance() const {
	Eigen::VectorXd residuals = response - Eigen::VectorXd::Constant(response.size(), response.mean());
	return residuals.dot(residuals) / residuals.size();
}

double Data::guessGrowthRate() const {
	//we might assume that all the species are in one group, and that all competition coefficients are zero.
	//This gives us the average observed response.
	//If this is total, rather than per capita, we must divide by the average species density in the design.
	double guess = response.mean();
	if(!isPerCapita) guess /= design.mean();
	return guess;
}

double Data::guessCompetitionCoefficientMagnitude() const {
	//We will assume that with all species present at average density, growth halts.
	//This gives us 1, divided by the average density, divided by the number of species.
	return 1.0 / design.mean() / numSpecies;
}
