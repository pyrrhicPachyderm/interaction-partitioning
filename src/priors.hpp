#ifndef PRIORS_HPP
#define PRIORS_HPP

#include "parameters.hpp"

//It is a limitation of the fact the way RJSolver::findTransModelJumpProbabilityMultiplier() works
//that hyperpriors may only depend on the number of groups in a grouping,
//not the full structure of the grouping.
//This is necessary to find the relevant constant in polynomial time.
//This limitation enforced by only passing the relevant information to a Hyperprior subclass.

namespace Hyperpriors {
	class Base {
		public:
			virtual double getLogDensity(const GroupingSizeSet &groupingSizes) const = 0;
			virtual double getDensity(const GroupingSizeSet &groupingSizes) const;
			double getLogDensity(const GroupingSet &groupings) const {return getLogDensity(getGroupingSizeSet(groupings));}
			double getDensity(const GroupingSet &groupings) const {return getDensity(getGroupingSizeSet(groupings));}
			
			//Virtual destructor, as this is an abstract class.
			virtual ~Base() {};
	};
	
	class Flat : public Base {
		public:
			double getDensity(const GroupingSizeSet &groupingSizes) const override;
			double getLogDensity(const GroupingSizeSet &groupingSizes) const override;
	};
	
	class AIC : public Base {
		public:
			double getLogDensity(const GroupingSizeSet &groupingSizes) const override;
	};
}

class Hyperprior {
	protected:
		std::shared_ptr<const Hyperpriors::Base> p;
	public:
		Hyperprior() = default;
		Hyperprior(const Hyperpriors::Base *p):
			p(p) {}
		
		double getDensity(const GroupingSizeSet &groupingSizes) const {return p->getDensity(groupingSizes);}
		double getDensity(const GroupingSet &groupings) const {return p->getDensity(groupings);}
		double getLogDensity(const GroupingSizeSet &groupingSizes) const {return p->getLogDensity(groupingSizes);}
		double getLogDensity(const GroupingSet &groupings) const {return p->getLogDensity(groupings);}
};

class ParametersPrior {
	protected:
		Distribution<double> growthRatePrior;
		Distribution<double> competitionCoefficientPrior;
	public:
		ParametersPrior(Distribution<double> growthRatePrior, Distribution<double> competitionCoefficientPrior):
			growthRatePrior(growthRatePrior), competitionCoefficientPrior(competitionCoefficientPrior) {};
		
		double getLogDensity(const Parameters &parameters) const;
};

template<size_t nAug> class AugmentedParametersPrior : public ParametersPrior {
	protected:
		std::array<Distribution<double>, nAug> additionalParameterPriors;
	public:
		AugmentedParametersPrior(Distribution<double> growthRatePrior, Distribution<double> competitionCoefficientPrior, std::array<Distribution<double>, nAug> additionalParameterPriors):
			ParametersPrior(growthRatePrior, competitionCoefficientPrior), additionalParameterPriors(additionalParameterPriors) {};
		
		double getLogDensity(const AugmentedParameters<nAug> &parameters) const;
};

#endif
