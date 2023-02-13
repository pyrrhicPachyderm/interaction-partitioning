#include "priors.hpp"

double Hyperprior::flatFunc(const GroupingSizeSet &groupingSizes) {
	return 1.0;
}

double Hyperprior::aicFunc(const GroupingSizeSet &groupingSizes) {
	size_t numParameters = groupingSizes[GROWTH] + groupingSizes[ROW] * groupingSizes[COL];
	return exp(-1.0 * (double)numParameters);
}

double ParametersPrior::getDensity(Parameters parameters) const {
	double density = 1.0;
	
	Eigen::VectorXd growthRates = parameters.getGrowthRates();
	for(size_t i = 0; i < (size_t)growthRates.size(); i++) {
		density *= growthRatePrior.getDensity(growthRates(i));
	}
	
	Eigen::MatrixXdRowMajor competitionCoefficients = parameters.getCompetitionCoefficients();
	for(size_t i = 0; i < (size_t)competitionCoefficients.rows(); i++) {
		for(size_t j = 0; j < (size_t)competitionCoefficients.cols(); j++) {
			density *= competitionCoefficientPrior.getDensity(competitionCoefficients(i,j));
		}
	}
	
	return density;
}

template<size_t nAug> double AugmentedParametersPrior<nAug>::getDensity(AugmentedParameters<nAug> parameters) const {
	double density = ParametersPrior::getDensity((Parameters)parameters); //Call the parent class function.
	
	for(size_t i = 0; i < nAug; i++) {
		density *= additionalParameterPriors[i].getDensity(parameters.getAdditionalParameter(i));
	}
	
	return density;
}

//Explicitly instantiate.
template class AugmentedParametersPrior<1>;
