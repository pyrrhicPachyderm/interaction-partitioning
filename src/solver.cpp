#include "solver.hpp"

template<> GeneralisedSolver<Distributions::Normal>::AdditionalParametersVector GeneralisedSolver<Distributions::Normal>::guessInitialAdditionalParameters() const {
	//Normal error variance.
	return {{data.guessErrorVariance()}};
}
template<> GeneralisedSolver<Distributions::Gamma2>::AdditionalParametersVector GeneralisedSolver<Distributions::Gamma2>::guessInitialAdditionalParameters() const {
	//Gamma2 dispersion parameter.
	return {{1.0}};
}
template<> GeneralisedSolver<Distributions::DiscreteWrapper<Distributions::NegativeBinomial2>>::AdditionalParametersVector GeneralisedSolver<Distributions::DiscreteWrapper<Distributions::NegativeBinomial2>>::guessInitialAdditionalParameters() const {
	//NegativeBinomial2 dispersion parameter.
	return {{1.0}};
}
