#include <nlopt.hpp>
#include "nls.hpp"
#include "mlsolver.hpp"

#define NLOPT_RELATIVE_TOLERANCE 2e-3

template<typename SolverT> SolverT::ParametersT MaximumLikelihoodSolver<SolverT>::getSolution() {
	//Returns Parameters or AugmentedParameters as appropriate.
	if(isDirtySolution) calculateSolution();
	return solution;
}

template<typename SolverT> Parameters MaximumLikelihoodSolver<SolverT>::getSolutionParameters() {
	//Coerces to Parameters (not AugmentedParameters).
	return (Parameters)getSolution();
}

template<> double MaximumLikelihoodSolver<Solver>::getSolutionAdditionalParameter(size_t i) {
	assert(false);
}

template<typename SolverT> double MaximumLikelihoodSolver<SolverT>::getSolutionAdditionalParameter(size_t i) {
	//This will only work if SolverT is a type of GeneralisedSolver (the standard Solver case is handled above).
	//But alas, you can't partially specialise methods.
	return getSolution().getAdditionalParameter(i);
}

template<typename SolverT> Eigen::VectorXd MaximumLikelihoodSolver<SolverT>::parametersToVector(const SolverT::ParametersT &parameters) const {
	return parameters.getAsVector();
}

template<typename SolverT> SolverT::ParametersT MaximumLikelihoodSolver<SolverT>::vectorToParameters(const Eigen::VectorXd &vector) const {
	return ParametersT(vector, this->groupings);
}

Eigen::VectorXd GaussNewtonSolver::getResidualsFromVector(const Eigen::VectorXd &parameterVector) const {
	return getResiduals(vectorToParameters(parameterVector), groupings);
}

Jacobian GaussNewtonSolver::getResidualsJacobianFromVector(const Eigen::VectorXd &parameterVector) const {
	return getResidualsJacobian(vectorToParameters(parameterVector), groupings);
}

void GaussNewtonSolver::calculateSolution() {
	ResidualsFunc residualsFunc = std::bind(&GaussNewtonSolver::getResidualsFromVector, this, std::placeholders::_1);
	JacobianFunc jacobianFunc = std::bind(&GaussNewtonSolver::getResidualsJacobianFromVector, this, std::placeholders::_1);
	
	Eigen::VectorXd initialParameterVector = parametersToVector(Parameters(data, groupings));
	Eigen::VectorXd parameterTolerances = Parameters::getTolerances(data, groupings);
	
	Eigen::VectorXd solutionVector = gaussNewtonNLS(residualsFunc, jacobianFunc, initialParameterVector, parameterTolerances);
	solution = vectorToParameters(solutionVector);
	isDirtySolution = false;
}

template<typename ErrDistT> double NLoptSolver<ErrDistT>::getLogLikelihoodFromVector(const std::vector<double> &parametersVector) {
	return this->getLogLikelihood(this->vectorToParameters(Eigen::Map<const Eigen::VectorXd>(parametersVector.data(), parametersVector.size())), this->groupings);
}

template<typename ErrDistT> double NLoptSolver<ErrDistT>::optimisationFunc(const std::vector<double>& parametersVector, std::vector<double>& grad, void* solver) {
	//We ignore grad. We'll need to fix this if we ever use non-derivative-free optimisation.
	return ((NLoptSolver<ErrDistT>*)solver)->getLogLikelihoodFromVector(parametersVector);
}

template<typename ErrDistT> void NLoptSolver<ErrDistT>::calculateSolution() {
	Eigen::VectorXd parametersVectorEigen = this->parametersToVector(ParametersT(this->data, this->groupings, this->guessInitialAdditionalParameters()));
	std::vector<double> parametersVector = std::vector<double>(parametersVectorEigen.data(), parametersVectorEigen.data() + parametersVectorEigen.size()); //Take a copy as a std::vector.
	
	nlopt::opt optimiser = nlopt::opt(nlopt::LN_SBPLX, parametersVector.size());
	optimiser.set_max_objective(NLoptSolver<ErrDistT>::optimisationFunc, (void*)this);
	optimiser.set_xtol_rel(NLOPT_RELATIVE_TOLERANCE);
	
	double finalLogLikelihood;
	optimiser.optimize(parametersVector, finalLogLikelihood);
	
	this->solution = this->vectorToParameters(Eigen::Map<const Eigen::VectorXd>(parametersVector.data(), parametersVector.size()));
	this->isDirtySolution = false;
}

template<typename SolverT> Eigen::VectorXd MaximumLikelihoodSolver<SolverT>::getSolutionPredictions() {
	return this->getPredictions(getSolutionParameters(), groupings);
}

template<typename SolverT> Eigen::VectorXd MaximumLikelihoodSolver<SolverT>::getSolutionResiduals() {
	return this->getResiduals(getSolutionParameters(), groupings);
}

template<typename SolverT> size_t MaximumLikelihoodSolver<SolverT>::getNumParameters() {
	ParametersT parameters = getSolution();
	return parameters.getNumParameters();
}

template<typename SolverT> double MaximumLikelihoodSolver<SolverT>::getAIC() {
	double p = getNumParameters();
	double deviance = getDeviance();
	
	double aic = 2 * p + deviance;
	
	return aic;
}

template<typename SolverT> double MaximumLikelihoodSolver<SolverT>::getAICc() {
	double n = this->data.getNumObservations();
	double p = getNumParameters();
	
	double aic = getAIC();
	double aicc = aic + (double)(2 * p * p + 2 * p) / (n - p - 1);
	
	return aicc;
}

template<typename SolverT> double MaximumLikelihoodSolver<SolverT>::getR2() {
	//Returns R^2, the coefficient of determination.
	
	Eigen::VectorXd response = this->getObservations();
	Eigen::VectorXd normalisedResponse = response - Eigen::VectorXd::Constant(response.size(), response.mean());
	double totalSS = normalisedResponse.dot(normalisedResponse);
	
	Eigen::VectorXd residuals = getSolutionResiduals();
	double residualSS = residuals.dot(residuals);
	
	return 1.0 - (residualSS / totalSS);
}

double GaussNewtonSolver::getDeviance() {
	//Returns the deviance.
	//That is, the negative of twice the log likelihood.
	Parameters parameters = getSolutionParameters();
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

template<typename ErrDistT> double NLoptSolver<ErrDistT>::getDeviance() {
	//Returns the deviance.
	//That is, the negative of twice the log likelihood.
	return -2 * this->getLogLikelihood(this->getSolution(), this->groupings);
}

//Explicitly instantiate.
template class NLoptSolver<Distributions::Gamma2>;
template class NLoptSolver<Distributions::DiscreteWrapper<Distributions::NegativeBinomial2>>;
