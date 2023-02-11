#include "solver.hpp"

Eigen::MatrixXd Solver::getColGroupedDesign() const {
	Eigen::MatrixXd colGroupedDesign = Eigen::MatrixXd::Zero(data.getNumObservations(), getGrouping(COL).getNumGroups());
	
	for(size_t obs = 0; obs < data.getNumObservations(); obs++) {
		for(size_t sp = 0; sp < data.getNumColSpecies(); sp++) {
			colGroupedDesign(obs, getGrouping(COL).getGroup(sp)) += data.getDesign()(obs, sp);
		}
	}
	
	return colGroupedDesign;
}

Eigen::VectorXd Solver::getPredictions(const Parameters &parameters) const {
	Eigen::VectorXd predictions = Eigen::VectorXd::Zero(data.getNumObservations());
	
	Eigen::MatrixXd colGroupedDesign = getColGroupedDesign();
	
	for(size_t obs = 0; obs < data.getNumObservations(); obs++) {
		size_t focal = data.getFocal()[obs];
		size_t focalGrowthGroup = getGrouping(GROWTH).getGroup(focal);
		size_t focalRowGroup = getGrouping(ROW).getGroup(focal);
		
		double focalGrowthRate = parameters.getGrowthRate(focalGrowthGroup);
		
		double intrinsicGrowth = focalGrowthRate;
		double totalCompetition = parameters.getCompetitionCoefficientsRow(focalRowGroup).dot(colGroupedDesign.row(obs));
		double prediction = intrinsicGrowth * (1.0 - totalCompetition);
		predictions[obs] = prediction;
	}
	
	return predictions;
}

Eigen::VectorXd Solver::getResiduals(const Parameters &parameters) const {
	return data.getResponse() - getPredictions(parameters);
}
