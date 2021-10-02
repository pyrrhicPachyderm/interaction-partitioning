#include "solver.hpp"

size_t Solver::getGrowthRateIndex(size_t growthGroup) const {
	return growthGroup;
}

double Solver::getGrowthRate(Solver::ParameterVector parameters, size_t growthGroup) const {
	return parameters(getGrowthRateIndex(growthGroup));
}

Eigen::VectorXd Solver::getGrowthRates(Solver::ParameterVector parameters) const {
	return parameters.block(0, getGrowthRateIndex(0), 1, growthGrouping.getNumGroups());
}

size_t Solver::getCompetitionCoefficientIndex(size_t rowGroup, size_t colGroup) const {
	return growthGrouping.getNumGroups() + colGrouping.getNumGroups() * rowGroup + colGroup;
}

double Solver::getCompetitionCoefficient(Solver::ParameterVector parameters, size_t rowGroup, size_t colGroup) const {
	return parameters(getCompetitionCoefficientIndex(rowGroup, colGroup));
}

Eigen::VectorXd Solver::getCompetitionCoefficientsRow(Solver::ParameterVector parameters, size_t rowGroup) const {
	return parameters.block(0, getCompetitionCoefficientIndex(rowGroup, 0), 1, colGrouping.getNumGroups());
}

Eigen::MatrixXd Solver::getColGroupedDesign() const {
	//TODO: Memoise this.
	Eigen::MatrixXd colGroupedDesign = Eigen::MatrixXd::Zero(data.numObservations, colGrouping.getNumGroups());
	
	for(size_t obs = 0; obs < data.numObservations; obs++) {
		for(size_t sp = 0; sp < data.numSpecies; sp++) {
			colGroupedDesign(obs, colGrouping.getGroup(sp)) += data.getDesign()(obs, sp);
		}
	}
	
	return colGroupedDesign;
}

Solver::ParameterVector Solver::getInitialParameterValues() const {
	//We need somewhat reasonable guesses for the growth rates and the competition coefficients.
	//It is reasonable to guess that all the competition coefficients are zero.
	double competitionCoefficient = 0.0;
	//As for the growth rates, we might assume that all the species are in one group, and that all competition coefficients are zero.
	//This gives us, approximately, the average observed response divided by the average species density in the design.
	double growthRate = data.getResponse().mean() / data.getDesign().mean();
	
	size_t numGrowthRates = growthGrouping.getNumGroups();
	size_t numCompetitionCoefficients = rowGrouping.getNumGroups() * colGrouping.getNumGroups();
	
	Solver::ParameterVector parameters = Eigen::VectorXd(numGrowthRates + numCompetitionCoefficients);
	for(size_t i = 0; i < numGrowthRates; i++) {
		parameters[i] = growthRate;
	}
	for(size_t i = numGrowthRates; i < numGrowthRates + numCompetitionCoefficients; i++) {
		parameters[i] = competitionCoefficient;
	}
	
	return parameters;
}

Eigen::VectorXd Solver::getPredictions(ParameterVector parameters) const {
	Eigen::VectorXd predictions = Eigen::VectorXd::Zero(data.numObservations);
	
	Eigen::MatrixXd colGroupedDesign = getColGroupedDesign();
	
	for(size_t obs = 0; obs < data.numObservations; obs++) {
		size_t focal = data.getFocal()[obs];
		size_t focalGrowthGroup = growthGrouping.getGroup(focal);
		size_t focalRowGroup = rowGrouping.getGroup(focal);
		
		double focalGrowthRate = getGrowthRate(parameters, focalGrowthGroup);
		double focalDensity = data.getDesign()(obs, focal);
		
		double intrinsicGrowth = focalGrowthRate * focalDensity;
		double totalCompetition = getCompetitionCoefficientsRow(parameters, focalRowGroup).dot(colGroupedDesign.row(obs));
		double prediction = intrinsicGrowth * (1.0 - totalCompetition);
		predictions[obs] = prediction;
	}
	
	return predictions;
}

Eigen::VectorXd Solver::getResiduals(ParameterVector parameters) const {
	return data.getResponse() - getPredictions(parameters);
}
