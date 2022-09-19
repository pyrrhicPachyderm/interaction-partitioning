#ifndef PRIORS_HPP
#define PRIORS_HPP

#include "parameters.hpp"

class Hyperprior {
	protected:
		typedef std::function<double(const GroupingSet &groupings)> HyperpriorFunc;
		HyperpriorFunc hyperpriorFunc;
		
		Hyperprior(HyperpriorFunc hyperpriorFunc):
			hyperpriorFunc(hyperpriorFunc) {};
		
		static double flatFunc(const GroupingSet &groupings);
		static double aicFunc(const GroupingSet &groupings);
	public:
		static Hyperprior flat() {return Hyperprior(flatFunc);}
		static Hyperprior aic() {return Hyperprior(aicFunc);}
		
		double getDensity(GroupingSet groupings) const {return hyperpriorFunc(groupings);}
};

#endif
