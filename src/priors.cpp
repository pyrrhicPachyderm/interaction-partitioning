#include <math.h>
#include "priors.hpp"

double Hyperpriors::Base::getDensity(const GroupingSizeSet &groupingSizes) const {
	//A base function that may be overwritten if the subclass can do it better.
	return exp(getLogDensity(groupingSizes));
}

double Hyperpriors::Flat::getDensity(const GroupingSizeSet &groupingSizes) const {
	return 1.0;
}
double Hyperpriors::Flat::getLogDensity(const GroupingSizeSet &groupingSizes) const {
	return 0.0;
}

double Hyperpriors::AIC::getLogDensity(const GroupingSizeSet &groupingSizes) const {
	size_t numParameters = groupingSizes[GROWTH] + groupingSizes[ROW] * groupingSizes[COL];
	return -1.0 * (double)numParameters;
}

double ParametersPrior::getLogDensity(const Parameters &parameters) const {
	double logDensity = 0.0;
	
	Eigen::VectorXd growthRates = parameters.getGrowthRates();
	for(size_t i = 0; i < (size_t)growthRates.size(); i++) {
		logDensity += growthRatePrior.getLogDensity(growthRates(i));
	}
	
	Eigen::MatrixXdRowMajor competitionCoefficients = parameters.getCompetitionCoefficients();
	for(size_t i = 0; i < (size_t)competitionCoefficients.rows(); i++) {
		for(size_t j = 0; j < (size_t)competitionCoefficients.cols(); j++) {
			logDensity += competitionCoefficientPrior.getLogDensity(competitionCoefficients(i,j));
		}
	}
	
	return logDensity;
}

template<size_t nAug> double AugmentedParametersPrior<nAug>::getLogDensity(const AugmentedParameters<nAug> &parameters) const {
	double logDensity = ParametersPrior::getLogDensity((Parameters)parameters); //Call the parent class function.
	
	for(size_t i = 0; i < nAug; i++) {
		logDensity += additionalParameterPriors[i].getLogDensity(parameters.getAdditionalParameter(i));
	}
	
	return logDensity;
}

//Explicitly instantiate.
template class AugmentedParametersPrior<1>;
