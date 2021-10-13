#include "nls.hpp"
#include "solver.hpp"

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

Eigen::VectorXd Solver::getPredictions(const Parameters &parameters) {
	Eigen::VectorXd predictions = Eigen::VectorXd::Zero(data.numObservations);
	
	Eigen::MatrixXd colGroupedDesign = getColGroupedDesign();
	
	for(size_t obs = 0; obs < data.numObservations; obs++) {
		size_t focal = data.getFocal()[obs];
		size_t focalGrowthGroup = growthGrouping.getGroup(focal);
		size_t focalRowGroup = rowGrouping.getGroup(focal);
		
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

Eigen::VectorXd Solver::getResidualsFromVector(const Eigen::VectorXd &parameterVector) {
	return getResiduals(Parameters(parameterVector, growthGrouping, rowGrouping, colGrouping));
}

Solver::Jacobian Solver::getJacobian(const Parameters &parameters) {
	Solver::Jacobian jacobian = Eigen::MatrixXd::Zero(data.numObservations, parameters.getNumParameters());
	
	Eigen::MatrixXd colGroupedDesign = getColGroupedDesign();
	
	//Note that this is the Jacobian of the residuals, not of the predicted values.
	//As such, it is negated, compared to the predicted values.
	for(size_t obs = 0; obs < data.numObservations; obs++) {
		size_t focal = data.getFocal()[obs];
		size_t focalGrowthGroup = growthGrouping.getGroup(focal);
		size_t focalRowGroup = rowGrouping.getGroup(focal);
		
		double focalGrowthRate = parameters.getGrowthRate(focalGrowthGroup);
		double focalDensity = data.getDesign()(obs, focal);
		
		//First, the derivatives with respect to the growth rates.
		//If it's not the growth rate of the observation's focal species, this is zero.
		//So there will be only one per observation.
		//This one will be equal to the prediction, divided by the growth rate itself.
		double totalCompetition = parameters.getCompetitionCoefficientsRow(focalRowGroup).dot(colGroupedDesign.row(obs));
		double derivative = 1.0 - totalCompetition;
		if(!data.isPerCapita) derivative *= focalDensity;
		jacobian(obs, parameters.getAsVectorGrowthRateIndex(focalGrowthGroup)) = -derivative;
		
		//Second, the derivatives with respect to the competition coefficients.
		//If it's not a competition coefficient *on* the focal species, this is zero.
		//So there will be a number per row equal to the number of column groups.
		//This will be the negative of the focal growth rate, times the focal density (if not per capita), times the column group density.
		for(size_t colGroup = 0; colGroup < colGrouping.getNumGroups(); colGroup++) {
			double derivative = - focalGrowthRate * colGroupedDesign(obs, colGroup);
			if(!data.isPerCapita) derivative *= focalDensity;
			jacobian(obs, parameters.getAsVectorCompetitionCoefficientIndex(focalRowGroup, colGroup)) = -derivative;
		}
	}
	
	return jacobian;
}

Solver::Jacobian Solver::getJacobianFromVector(const Eigen::VectorXd &parameterVector) {
	return getJacobian(Parameters(parameterVector, growthGrouping, rowGrouping, colGrouping));
}

void Solver::calculateSolution() {
	ResidualsFunc residualsFunc = std::bind(&Solver::getResidualsFromVector, this, std::placeholders::_1);
	JacobianFunc jacobianFunc = std::bind(&Solver::getJacobianFromVector, this, std::placeholders::_1);
	
	Eigen::VectorXd initialParameterVector = Parameters(data, growthGrouping, rowGrouping, colGrouping).getAsVector();
	Eigen::VectorXd parameterTolerances = Parameters::getTolerances(data, growthGrouping, rowGrouping, colGrouping);
	
	Eigen::VectorXd solutionVector = gaussNewtonNLS(residualsFunc, jacobianFunc, initialParameterVector, parameterTolerances);
	solution = Parameters(solutionVector, growthGrouping, rowGrouping, colGrouping);
	isDirtySolution = false;
}

Parameters Solver::getSolution() {
	if(isDirtySolution) calculateSolution();
	return solution;
}

Eigen::VectorXd Solver::getSolutionPredictions() {
	//TODO: Memoise.
	return getPredictions(getSolution());
}

Eigen::VectorXd Solver::getSolutionResiduals() {
	//TODO: Memoise.
	return getResiduals(getSolution());
}

double Solver::getDeviance() {
	//Returns the deviance.
	//That is, the negative of twice the log likelihood.
	Parameters parameters = getSolution();
	Eigen::VectorXd residuals = getSolutionResiduals();
	double sumOfSquares = residuals.dot(residuals);
	
	//First, we must estimate the variance of the residuals.
	//They should already be centred around zero, so we shouldn't need to subtract the mean.
	//TODO: Double-check that this is the correct formula and degrees of freedom.
	double variance = sumOfSquares / (residuals.size() - parameters.getNumParameters());
	
	//The likelihood is the product of the likelihoods for each residual.
	//So the log likelihood is the sum of log likelihoods.
	//The likelihood for a given data point (noting that these are reisudals, so the mean is zero) is:
	//(2 * pi * variance) ^ (-1/2) * exp(-1/2 * residual^2 / variance)
	//So the log likelihood is:
	//-1/2 * log(2 * pi * variance) - 1/2 * residual^2 / variance
	//And the total log likehood is just the sum of these.
	//Note that the deviance is twice the negative log likelihood.
	//This cancels the factor of -1/2.
	double deviance = residuals.size() * log(2 * M_PI * variance) + sumOfSquares / variance;
	
	return deviance;
}

double Solver::getAIC() {
	Parameters parameters = getSolution();
	double deviance = getDeviance();
	
	double aic = 2 * parameters.getNumParameters() + deviance;
	
	return aic;
}

double Solver::getR2() {
	//Returns R^2, the coefficient of determination.
	
	Eigen::VectorXd response = data.getResponse();
	Eigen::VectorXd normalisedResponse = response - Eigen::VectorXd::Constant(response.size(), response.mean());
	double totalSS = normalisedResponse.dot(normalisedResponse);
	
	Eigen::VectorXd residuals = getSolutionResiduals();
	double residualSS = residuals.dot(residuals);
	
	return 1.0 - (residualSS / totalSS);
}
