#include "parameters.hpp"
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

Eigen::MatrixXd Data::getColGroupedDesign(const Grouping &grouping) const {
	Eigen::MatrixXd colGroupedDesign = Eigen::MatrixXd::Zero(design.rows(), grouping.getNumGroups());
	
	for(size_t obs = 0; obs < (size_t)design.rows(); obs++) {
		for(size_t sp = 0; sp < (size_t)design.cols(); sp++) {
			colGroupedDesign(obs, grouping.getGroup(sp)) += design(obs, sp);
		}
	}
	
	return colGroupedDesign;
}

Eigen::VectorXd Data::getPredictions(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const {
	Eigen::VectorXd predictions = Eigen::VectorXd::Zero(response.size());
	
	Eigen::MatrixXd colGroupedDesign = getColGroupedDesign(groupings[COL]);
	
	for(size_t obs = 0; obs < numObservations; obs++) {
		size_t focalGrowthGroup = groupings[GROWTH].getGroup(focal[obs]);
		size_t focalRowGroup = groupings[ROW].getGroup(focal[obs]);
		
		double focalDensity = 1.0;
		double focalGrowthRate = parameters.getGrowthRate(focalGrowthGroup);
		Eigen::VectorXd densities = colGroupedDesign.row(obs);
		Eigen::VectorXd competitionCoefficients = parameters.getCompetitionCoefficientsRow(focalRowGroup);
		
		predictions[obs] = model.getDerivative(focalDensity, focalGrowthRate, densities, competitionCoefficients);
	}
	
	return predictions;
}

Eigen::VectorXd Data::getResiduals(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const {
	return getResponse() - getPredictions(model, parameters, groupings);
}

Data::Jacobian Data::getPredictionsJacobian(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const {
	Data::Jacobian jacobian = Eigen::MatrixXd::Zero(numObservations, parameters.getNumParameters());
	
	Eigen::MatrixXd colGroupedDesign = getColGroupedDesign(groupings[COL]);
	
	//Note that this is the Jacobian of the residuals, not of the predicted values.
	//As such, it is negated, compared to the predicted values.
	for(size_t obs = 0; obs < numObservations; obs++) {
		size_t focalGrowthGroup = groupings[GROWTH].getGroup(focal[obs]);
		size_t focalRowGroup = groupings[ROW].getGroup(focal[obs]);
		
		double focalDensity = 1.0;
		double focalGrowthRate = parameters.getGrowthRate(focalGrowthGroup);
		Eigen::VectorXd densities = colGroupedDesign.row(obs);
		Eigen::VectorXd competitionCoefficients = parameters.getCompetitionCoefficientsRow(focalRowGroup);
		
		//First, the derivatives with respect to the growth rates.
		//If it's not the growth rate of the observation's focal species, this is zero.
		//So there will be only one per observation.
		jacobian(obs, parameters.getAsVectorGrowthRateIndex(focalGrowthGroup)) = model.getGrowthRateJacobian(focalDensity, focalGrowthRate, densities, competitionCoefficients);
		
		//Second, the derivatives with respect to the competition coefficients.
		//If it's not a competition coefficient *on* the focal species, this is zero.
		//So there will be a number per row equal to the number of column groups.
		for(size_t colGroup = 0; colGroup < groupings[COL].getNumGroups(); colGroup++) {
			jacobian(obs, parameters.getAsVectorCompetitionCoefficientIndex(focalRowGroup, colGroup)) = model.getCompetitionCoefficientJacobian(focalDensity, focalGrowthRate, densities, competitionCoefficients, colGroup);
		}
	}
	
	return jacobian;
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
	Eigen::VectorXd residuals = response - Eigen::VectorXd::Constant(response.size(), response.mean());
	return residuals.dot(residuals) / residuals.size();
}
