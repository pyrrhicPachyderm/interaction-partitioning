#ifndef PRIORS_HPP
#define PRIORS_HPP

#include "parameters.hpp"

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
		
		double getDensity(GroupingSizeSet groupingSizes) const {return hyperpriorFunc(groupingSizes);}
		double getDensity(GroupingSet groupings) const {return getDensity(getGroupingSizeSet(groupings));}
};

#endif
