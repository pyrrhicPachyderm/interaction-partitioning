#include "priors.hpp"

double Hyperprior::flatFunc(const GroupingSizeSet &groupingSizes) {
	return 1.0;
}

double Hyperprior::aicFunc(const GroupingSizeSet &groupingSizes) {
	size_t numParameters = groupingSizes[GROWTH] + groupingSizes[ROW] * groupingSizes[COL];
	return exp(-1.0 * (double)numParameters);
}

std::vector<double> ParametersPrior::getDensities(const Parameters &parameters) const {
	std::vector<double> densities;
	
	Eigen::VectorXd growthRates = parameters.getGrowthRates();
	for(size_t i = 0; i < (size_t)growthRates.size(); i++) {
		densities.push_back(growthRatePrior.getDensity(growthRates(i)));
	}
	
	Eigen::MatrixXdRowMajor competitionCoefficients = parameters.getCompetitionCoefficients();
	for(size_t i = 0; i < (size_t)competitionCoefficients.rows(); i++) {
		for(size_t j = 0; j < (size_t)competitionCoefficients.cols(); j++) {
			densities.push_back(competitionCoefficientPrior.getDensity(competitionCoefficients(i,j)));
		}
	}
	
	return densities;
}

template<size_t nAug> std::vector<double> AugmentedParametersPrior<nAug>::getDensities(const AugmentedParameters<nAug> &parameters) const {
	std::vector<double> densities = ParametersPrior::getDensities((Parameters)parameters); //Call the parent class function.
	
	for(size_t i = 0; i < nAug; i++) {
		densities.push_back(additionalParameterPriors[i].getDensity(parameters.getAdditionalParameter(i)));
	}
	
	return densities;
}

//Explicitly instantiate.
template class AugmentedParametersPrior<1>;
