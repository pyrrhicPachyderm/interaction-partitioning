#include "nls.hpp"
#include "mlsolver.hpp"

Eigen::VectorXd MaximumLikelihoodSolver::getResidualsFromVector(const Eigen::VectorXd &parameterVector) const {
	return getResiduals(Parameters(parameterVector, groupings));
}

MaximumLikelihoodSolver::Jacobian MaximumLikelihoodSolver::getJacobian(const Parameters &parameters) const {
	MaximumLikelihoodSolver::Jacobian jacobian = Eigen::MatrixXd::Zero(data.getNumObservations(), parameters.getNumParameters());
	
	Eigen::MatrixXd colGroupedDesign = getColGroupedDesign();
	
	//Note that this is the Jacobian of the residuals, not of the predicted values.
	//As such, it is negated, compared to the predicted values.
	for(size_t obs = 0; obs < data.getNumObservations(); obs++) {
		size_t focal = data.getFocal()[obs];
		size_t focalGrowthGroup = groupings[GROWTH].getGroup(focal);
		size_t focalRowGroup = groupings[ROW].getGroup(focal);
		
		double focalGrowthRate = parameters.getGrowthRate(focalGrowthGroup);
		
		//First, the derivatives with respect to the growth rates.
		//If it's not the growth rate of the observation's focal species, this is zero.
		//So there will be only one per observation.
		//This one will be equal to the prediction, divided by the growth rate itself.
		double totalCompetition = parameters.getCompetitionCoefficientsRow(focalRowGroup).dot(colGroupedDesign.row(obs));
		double derivative = 1.0 - totalCompetition;
		jacobian(obs, parameters.getAsVectorGrowthRateIndex(focalGrowthGroup)) = -derivative;
		
		//Second, the derivatives with respect to the competition coefficients.
		//If it's not a competition coefficient *on* the focal species, this is zero.
		//So there will be a number per row equal to the number of column groups.
		//This will be the negative of the focal growth rate times the column group density.
		for(size_t colGroup = 0; colGroup < groupings[COL].getNumGroups(); colGroup++) {
			double derivative = - focalGrowthRate * colGroupedDesign(obs, colGroup);
			jacobian(obs, parameters.getAsVectorCompetitionCoefficientIndex(focalRowGroup, colGroup)) = -derivative;
		}
	}
	
	return jacobian;
}

MaximumLikelihoodSolver::Jacobian MaximumLikelihoodSolver::getJacobianFromVector(const Eigen::VectorXd &parameterVector) const {
	return getJacobian(Parameters(parameterVector, groupings));
}

void MaximumLikelihoodSolver::calculateSolution() {
	ResidualsFunc residualsFunc = std::bind(&MaximumLikelihoodSolver::getResidualsFromVector, this, std::placeholders::_1);
	JacobianFunc jacobianFunc = std::bind(&MaximumLikelihoodSolver::getJacobianFromVector, this, std::placeholders::_1);
	
	Eigen::VectorXd initialParameterVector = Parameters(data, groupings).getAsVector();
	Eigen::VectorXd parameterTolerances = Parameters::getTolerances(data, groupings);
	
	Eigen::VectorXd solutionVector = gaussNewtonNLS(residualsFunc, jacobianFunc, initialParameterVector, parameterTolerances);
	solution = Parameters(solutionVector, groupings);
	isDirtySolution = false;
}

Parameters MaximumLikelihoodSolver::getSolution() {
	if(isDirtySolution) calculateSolution();
	return solution;
}

Eigen::VectorXd MaximumLikelihoodSolver::getSolutionPredictions() {
	//TODO: Memoise.
	return getPredictions(getSolution());
}

Eigen::VectorXd MaximumLikelihoodSolver::getSolutionResiduals() {
	//TODO: Memoise.
	return getResiduals(getSolution());
}

double MaximumLikelihoodSolver::getDeviance() {
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
	//The likelihood for a given data point (noting that these are residuals, so the mean is zero) is:
	//(2 * pi * variance) ^ (-1/2) * exp(-1/2 * residual^2 / variance)
	//So the log likelihood is:
	//-1/2 * log(2 * pi * variance) - 1/2 * residual^2 / variance
	//And the total log likehood is just the sum of these.
	//Note that the deviance is twice the negative log likelihood.
	//This cancels the factor of -1/2.
	double deviance = residuals.size() * log(2 * M_PI * variance) + sumOfSquares / variance;
	
	return deviance;
}

double MaximumLikelihoodSolver::getAIC() {
	Parameters parameters = getSolution();
	double deviance = getDeviance();
	
	double aic = 2 * parameters.getNumParameters() + deviance;
	
	return aic;
}

double MaximumLikelihoodSolver::getAICc() {
	Parameters parameters = getSolution();
	double n = data.getNumObservations();
	double p = parameters.getNumParameters();
	
	double aic = getAIC();
	double aicc = aic + (double)(2 * p * p + 2 * p) / (n - p - 1);
	
	return aicc;
}

double MaximumLikelihoodSolver::getR2() {
	//Returns R^2, the coefficient of determination.
	
	Eigen::VectorXd response = data.getResponse();
	Eigen::VectorXd normalisedResponse = response - Eigen::VectorXd::Constant(response.size(), response.mean());
	double totalSS = normalisedResponse.dot(normalisedResponse);
	
	Eigen::VectorXd residuals = getSolutionResiduals();
	double residualSS = residuals.dot(residuals);
	
	return 1.0 - (residualSS / totalSS);
}
