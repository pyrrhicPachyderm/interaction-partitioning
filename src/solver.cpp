#include "solver.hpp"

void Solver::calculateColGroupedDesign() {
	colGroupedDesign = Eigen::MatrixXd::Zero(data.numObservations, getColGrouping().getNumGroups());
	
	for(size_t obs = 0; obs < data.numObservations; obs++) {
		for(size_t sp = 0; sp < data.numSpecies; sp++) {
			colGroupedDesign(obs, getColGrouping().getGroup(sp)) += data.getDesign()(obs, sp);
		}
	}
	
	isDirtyColGroupedDesign = false;
}

Eigen::MatrixXd Solver::getColGroupedDesign() {
	if(isDirtyColGroupedDesign) calculateColGroupedDesign();
	return colGroupedDesign;
}

Eigen::VectorXd Solver::getPredictions(const Parameters &parameters) {
	Eigen::VectorXd predictions = Eigen::VectorXd::Zero(data.numObservations);
	
	Eigen::MatrixXd colGroupedDesign = getColGroupedDesign();
	
	for(size_t obs = 0; obs < data.numObservations; obs++) {
		size_t focal = data.getFocal()[obs];
		size_t focalGrowthGroup = getGrowthGrouping().getGroup(focal);
		size_t focalRowGroup = getRowGrouping().getGroup(focal);
		
		double focalGrowthRate = parameters.getGrowthRate(focalGrowthGroup);
		double focalDensity = data.getDesign()(obs, focal);
		
		double intrinsicGrowth = focalGrowthRate;
		if(!data.isPerCapita) intrinsicGrowth *= focalDensity;
		double totalCompetition = parameters.getCompetitionCoefficientsRow(focalRowGroup).dot(colGroupedDesign.row(obs));
		double prediction = intrinsicGrowth * (1.0 - totalCompetition);
		predictions[obs] = prediction;
	}
	
	return predictions;
}

Eigen::VectorXd Solver::getResiduals(const Parameters &parameters) {
	return data.getResponse() - getPredictions(parameters);
}
