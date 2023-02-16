#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP

#include <math.h>
#include <random>
#include <vector>
#include <memory>

//This sets the type of random engine that will be used throughout the program.
typedef std::default_random_engine RandomGenerator;

namespace Distributions {
	template<typename DomainValue> class Base {
		public:
			virtual double getDensity(DomainValue x) const = 0;
			virtual double getLogDensity(DomainValue x) const {return log(getDensity(x));}; //A base function that may be overwritten if there's a better way for a given distribution.
			virtual DomainValue getRandom(RandomGenerator &generator) const = 0;
			
			virtual ~Base() {};
	};
	
	class Uniform : public Base<double> {
		protected:
			double min;
			double max;
		public:
			Uniform(double min, double max):
				min(min), max(max) {};
			
			double getDensity(double x) const override;
			double getRandom(RandomGenerator &generator) const override;
	};
	
	class Normal : public Base<double> {
		protected:
			double mean;
			double variance;
		public:
			Normal(double mean, double variance):
				mean(mean), variance(variance) {};
			
			double getDensity(double x) const override;
			double getLogDensity(double x) const override;
			double getRandom(RandomGenerator &generator) const override;
	};
	
	class InverseGamma : public Base<double> {
		protected:
			double shape;
			double scale;
		public:
			InverseGamma(double shape, double scale):
				shape(shape), scale(scale) {};
			
			double getDensity(double x) const override;
			double getRandom(RandomGenerator &generator) const override;
	};
}

//A wrapper class to deal with maintaining the pointer required for runtime polymorphism.
template<typename DomainValue> class Distribution {
	protected:
		std::shared_ptr<const Distributions::Base<DomainValue>> d;
	public:
		Distribution(const Distributions::Base<DomainValue> *d):
			d(d) {}
		
		double getDensity(DomainValue x) const {return d->getDensity(x);};
		double getLogDensity(DomainValue x) const {return d->getLogDensity(x);};
		DomainValue getRandom(RandomGenerator &generator) const {return d->getRandom(generator);};
};

#endif
