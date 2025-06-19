#include <nlopt.hpp>
#include "nls.hpp"
#include "mlsolver.hpp"

#define NLOPT_RELATIVE_TOLERANCE 5e-5

template<typename SolverT> SolverT::ParametersT MaximumLikelihoodSolver<SolverT>::getSolution(bool isNull) {
	//Returns Parameters or AugmentedParameters as appropriate.
	if(isNull) {
		if(isDirtyNullSolution) calculateSolution(isNull);
		return nullSolution;
	} else {
		if(isDirtySolution) calculateSolution(isNull);
		return solution;
	}
}

template<typename SolverT> Parameters MaximumLikelihoodSolver<SolverT>::getSolutionParameters(bool isNull) {
	//Coerces to Parameters (not AugmentedParameters).
	return (Parameters)getSolution(isNull);
}

template<> double MaximumLikelihoodSolver<Solver>::getSolutionAdditionalParameter(bool isNull, size_t i) {
	assert(false);
}

template<typename SolverT> double MaximumLikelihoodSolver<SolverT>::getSolutionAdditionalParameter(bool isNull, size_t i) {
	//This will only work if SolverT is a type of GeneralisedSolver (the standard Solver case is handled above).
	//But alas, you can't partially specialise methods.
	return getSolution(isNull).getAdditionalParameter(i);
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

void GaussNewtonSolver::calculateSolution(bool isNull) {
	assert(!isNull); //TODO: Make this work if isNull.
	
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

template<typename ErrDistT> void NLoptSolver<ErrDistT>::calculateSolution(bool isNull) {
	ParametersT parameters = ParametersT(this->data, this->groupings, this->guessInitialAdditionalParameters());
	Eigen::VectorXd parametersVectorEigen = this->parametersToVector(parameters);
	std::vector<double> parametersVector = std::vector<double>(parametersVectorEigen.data(), parametersVectorEigen.data() + parametersVectorEigen.size()); //Take a copy as a std::vector.
	
	nlopt::opt optimiser = nlopt::opt(nlopt::LN_SBPLX, parametersVector.size());
	optimiser.set_max_objective(NLoptSolver<ErrDistT>::optimisationFunc, (void*)this);
	optimiser.set_ftol_rel(NLOPT_RELATIVE_TOLERANCE);
	
	if(isNull) {
		//Force all competition coefficients to be zero.
		std::vector<double> lowerBounds = optimiser.get_lower_bounds();
		std::vector<double> upperBounds = optimiser.get_upper_bounds();
		for(size_t i = 0; i < this->groupings[ROW].getNumGroups(); i++) {
			for(size_t j = 0; j < this->groupings[COL].getNumGroups(); j++) {
				size_t index = parameters.getAsVectorCompetitionCoefficientIndex(i, j);
				lowerBounds[index] = 0.0;
				upperBounds[index] = 0.0;
			}
		}
		optimiser.set_lower_bounds(lowerBounds);
		optimiser.set_upper_bounds(upperBounds);
	}
	
	double finalLogLikelihood;
	optimiser.optimize(parametersVector, finalLogLikelihood);
	
	ParametersT result = this->vectorToParameters(Eigen::Map<const Eigen::VectorXd>(parametersVector.data(), parametersVector.size()));
	
	if(isNull) {
		this->nullSolution = result;
		this->isDirtyNullSolution = false;
	} else {
		this->solution = result;
		this->isDirtySolution = false;
	}
}

template<typename SolverT> Eigen::VectorXd MaximumLikelihoodSolver<SolverT>::getSolutionPredictions(bool isNull) {
	return this->getPredictions(getSolutionParameters(isNull), groupings);
}

template<typename SolverT> Eigen::VectorXd MaximumLikelihoodSolver<SolverT>::getSolutionResiduals(bool isNull) {
	return this->getResiduals(getSolutionParameters(isNull), groupings);
}

template<typename SolverT> size_t MaximumLikelihoodSolver<SolverT>::getNumParameters(bool isNull) {
	ParametersT parameters = getSolution(isNull);
	if(isNull) {
		return parameters.getNumParameters() - parameters.getNumCompetitionCoefficients();
	} else {
		return parameters.getNumParameters();
	}
}

template<typename SolverT> double MaximumLikelihoodSolver<SolverT>::getAIC(bool isNull) {
	double p = getNumParameters(isNull);
	double deviance = getDeviance(isNull);
	
	double aic = 2 * p + deviance;
	
	return aic;
}

template<typename SolverT> double MaximumLikelihoodSolver<SolverT>::getAICc(bool isNull) {
	double n = this->data.getNumObservations();
	double p = getNumParameters(isNull);
	
	double aic = getAIC(isNull);
	double aicc = aic + (double)(2 * p * p + 2 * p) / (n - p - 1);
	
	return aicc;
}

template<typename SolverT> double MaximumLikelihoodSolver<SolverT>::getR2(bool isNull) {
	//Returns R^2, the coefficient of determination.
	
	Eigen::VectorXd response = this->getObservations();
	Eigen::VectorXd normalisedResponse = response - Eigen::VectorXd::Constant(response.size(), response.mean());
	double totalSS = normalisedResponse.dot(normalisedResponse);
	
	Eigen::VectorXd residuals = this->getSolutionResiduals(isNull);
	double residualSS = residuals.dot(residuals);
	
	return 1.0 - (residualSS / totalSS);
}

template<typename SolverT> double MaximumLikelihoodSolver<SolverT>::getMcFaddenR2() {
	//Returns McFadden's pseudo-R^2.
	return 1.0 - (this->getDeviance(false) / this->getDeviance(true));
}

double GaussNewtonSolver::getDeviance(bool isNull) {
	//Returns the deviance.
	//That is, the negative of twice the log likelihood.
	Parameters parameters = getSolutionParameters(isNull);
	Eigen::VectorXd residuals = getSolutionResiduals(isNull);
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

template<typename ErrDistT> double NLoptSolver<ErrDistT>::getDeviance(bool isNull) {
	//Returns the deviance.
	//That is, the negative of twice the log likelihood.
	return -2 * this->getLogLikelihood(this->getSolution(isNull), this->groupings);
}

//Explicitly instantiate.
template class NLoptSolver<Distributions::Gamma2>;
template class NLoptSolver<Distributions::DiscreteWrapper<Distributions::NegativeBinomial2>>;
