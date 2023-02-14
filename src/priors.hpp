#ifndef PRIORS_HPP
#define PRIORS_HPP

#include "parameters.hpp"

//It is a limitation of the fact the way RJSolver::findTransModelJumpProbabilityMultiplier() works
//that hyperpriors may only depend on the number of groups in a grouping,
//not the full structure of the grouping.
//This is necessary to find the relevant constant in polynomial time.
//This limitation enforced by only passing the relevant information to a HyperpriorFunc.

class Hyperprior {
	protected:
		typedef std::function<double(const GroupingSizeSet &groupingSizes)> HyperpriorFunc;
		HyperpriorFunc hyperpriorFunc;
		
		Hyperprior(HyperpriorFunc hyperpriorFunc):
			hyperpriorFunc(hyperpriorFunc) {};
		
		static double flatFunc(const GroupingSizeSet &groupingSizes);
		static double aicFunc(const GroupingSizeSet &groupingSizes);
	public:
		static Hyperprior flat() {return Hyperprior(flatFunc);}
		static Hyperprior aic() {return Hyperprior(aicFunc);}
		
		double getDensity(const GroupingSizeSet &groupingSizes) const {return hyperpriorFunc(groupingSizes);}
		double getDensity(const GroupingSet &groupings) const {return getDensity(getGroupingSizeSet(groupings));}
};

class ParametersPrior {
	protected:
		Distribution<double> growthRatePrior;
		Distribution<double> competitionCoefficientPrior;
	public:
		ParametersPrior(Distribution<double> growthRatePrior, Distribution<double> competitionCoefficientPrior):
			growthRatePrior(growthRatePrior), competitionCoefficientPrior(competitionCoefficientPrior) {};
		
		std::vector<double> getDensities(Parameters parameters) const;
};

template<size_t nAug> class AugmentedParametersPrior : public ParametersPrior {
	protected:
		std::array<Distribution<double>, nAug> additionalParameterPriors;
	public:
		AugmentedParametersPrior(Distribution<double> growthRatePrior, Distribution<double> competitionCoefficientPrior, std::array<Distribution<double>, nAug> additionalParameterPriors):
			ParametersPrior(growthRatePrior, competitionCoefficientPrior), additionalParameterPriors(additionalParameterPriors) {};
		
		std::vector<double> getDensities(AugmentedParameters<nAug> parameters) const;
};

#endif
