#include "solver.hpp"

double Solver::getGrowthRate(Solver::ParameterVector parameters, size_t growthGroup) const {
	return parameters(growthGroup);
}

Eigen::VectorXd Solver::getGrowthRates(Solver::ParameterVector parameters) const {
	return parameters.block(0, 0, 1, growthGrouping.getNumGroups());
}

double Solver::getCompetitionCoefficient(Solver::ParameterVector parameters, size_t rowGroup, size_t colGroup) const {
	return parameters(growthGrouping.getNumGroups() + colGrouping.getNumGroups() * rowGroup + colGroup);
}

Eigen::VectorXd Solver::getCompetitionCoefficientsRow(Solver::ParameterVector parameters, size_t rowGroup) const {
	return parameters.block(0, growthGrouping.getNumGroups() + colGrouping.getNumGroups() * rowGroup, 1, colGrouping.getNumGroups());
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
