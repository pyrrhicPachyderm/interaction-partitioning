#include "nls.hpp"
#include "solver.hpp"

#define RELATIVE_TOLERANCE 1e-6

size_t Solver::getGrowthRateIndex(size_t growthGroup) const {
	return growthGroup;
}

double Solver::getGrowthRate(Solver::ParameterVector parameters, size_t growthGroup) const {
	return parameters(getGrowthRateIndex(growthGroup));
}

Eigen::VectorXd Solver::getGrowthRates(Solver::ParameterVector parameters) const {
	return parameters.segment(getGrowthRateIndex(0), growthGrouping.getNumGroups());
}

size_t Solver::getCompetitionCoefficientIndex(size_t rowGroup, size_t colGroup) const {
	return growthGrouping.getNumGroups() + colGrouping.getNumGroups() * rowGroup + colGroup;
}

double Solver::getCompetitionCoefficient(Solver::ParameterVector parameters, size_t rowGroup, size_t colGroup) const {
	return parameters(getCompetitionCoefficientIndex(rowGroup, colGroup));
}

Eigen::VectorXd Solver::getCompetitionCoefficientsRow(Solver::ParameterVector parameters, size_t rowGroup) const {
	return parameters.segment(getCompetitionCoefficientIndex(rowGroup, 0), colGrouping.getNumGroups());
}

void Solver::calculateColGroupedDesign() {
	colGroupedDesign = Eigen::MatrixXd::Zero(data.numObservations, colGrouping.getNumGroups());
	
	for(size_t obs = 0; obs < data.numObservations; obs++) {
		for(size_t sp = 0; sp < data.numSpecies; sp++) {
			colGroupedDesign(obs, colGrouping.getGroup(sp)) += data.getDesign()(obs, sp);
		}
	}
	
	isDirtyColGroupedDesign = false;
}

Eigen::MatrixXd Solver::getColGroupedDesign() {
	if(isDirtyColGroupedDesign) calculateColGroupedDesign();
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

Solver::ParameterVector Solver::getParameterTolerances() const {
	//We need somewhat reasonable guesses for the magnitudes of the growth rates and the competition coefficients.
	//We will then multiply these by the RELATIVE_TOLERANCE.
	//For the growth rates, we will use the same guess as for the initial values.
	double growthRateTolerance = data.getResponse().mean() / data.getDesign().mean() * RELATIVE_TOLERANCE;
	//For the competition coefficients, we will assume that with all species present at average density, growth halts.
	//This gives us 1, divided by the square of average density, divided by the number of species.
	double competitionCoefficientTolerance = 1.0 / pow(data.getResponse().mean(), 2.0) / data.numSpecies * RELATIVE_TOLERANCE;
	
	size_t numGrowthRates = growthGrouping.getNumGroups();
	size_t numCompetitionCoefficients = rowGrouping.getNumGroups() * colGrouping.getNumGroups();
	
	Solver::ParameterVector tolerances = Eigen::VectorXd(numGrowthRates + numCompetitionCoefficients);
	for(size_t i = 0; i < numGrowthRates; i++) {
		tolerances[i] = growthRateTolerance;
	}
	for(size_t i = numGrowthRates; i < numGrowthRates + numCompetitionCoefficients; i++) {
		tolerances[i] = competitionCoefficientTolerance;
	}
	
	return tolerances;
}

Eigen::VectorXd Solver::getPredictions(ParameterVector parameters) {
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

Eigen::VectorXd Solver::getResiduals(ParameterVector parameters) {
	return data.getResponse() - getPredictions(parameters);
}

Solver::Jacobian Solver::getJacobian(ParameterVector parameters) {
	Solver::Jacobian jacobian = Eigen::MatrixXd::Zero(data.numObservations, parameters.size());
	
	Eigen::MatrixXd colGroupedDesign = getColGroupedDesign();
	
	//Note that this is the Jacobian of the residuals, not of the predicted values.
	//As such, it is negated, compared to the predicted values.
	for(size_t obs = 0; obs < data.numObservations; obs++) {
		size_t focal = data.getFocal()[obs];
		size_t focalGrowthGroup = growthGrouping.getGroup(focal);
		size_t focalRowGroup = rowGrouping.getGroup(focal);
		
		double focalGrowthRate = getGrowthRate(parameters, focalGrowthGroup);
		double focalDensity = data.getDesign()(obs, focal);
		
		//First, the derivatives with respect to the growth rates.
		//If it's not the growth rate of the observation's focal species, this is zero.
		//So there will be only one per observation.
		//This one will be equal to the prediction, divided by the growth rate itself.
		double totalCompetition = getCompetitionCoefficientsRow(parameters, focalRowGroup).dot(colGroupedDesign.row(obs));
		double derivative = focalDensity * (1.0 - totalCompetition);
		jacobian(obs, getGrowthRateIndex(focalGrowthGroup)) = -derivative;
		
		//Second, the derivatives with respect to the competition coefficients.
		//If it's not a competition coefficient *on* the focal species, this is zero.
		//So there will be a number per row equal to the number of column groups.
		//This will be the negative of the focal growth rate, times the focal density, times the column group density.
		for(size_t colGroup = 0; colGroup < colGrouping.getNumGroups(); colGroup++) {
			double derivative = - focalGrowthRate * focalDensity * colGroupedDesign(obs, colGroup);
			jacobian(obs, getCompetitionCoefficientIndex(focalRowGroup, colGroup)) = -derivative;
		}
	}
	
	return jacobian;
}

Solver::ParameterVector Solver::solve() {
	ResidualsFunc residualsFunc = std::bind(&Solver::getResiduals, this, std::placeholders::_1);
	JacobianFunc jacobianFunc = std::bind(&Solver::getJacobian, this, std::placeholders::_1);
	return gaussNewtonNLS(residualsFunc, jacobianFunc, getInitialParameterValues(), getParameterTolerances());
}
