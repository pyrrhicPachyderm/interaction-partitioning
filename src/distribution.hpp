#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP

#include <vector>
#include <memory>

namespace Distributions {
	template<typename DomainValue> class Base {
		public:
			virtual double getDensity(DomainValue x) const = 0;
			virtual DomainValue getRandom() const = 0;
			
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
			double getRandom() const override;
	};
	
	class Normal : public Base<double> {
		protected:
			double mean;
			double variance;
		public:
			Normal(double mean, double variance):
				mean(mean), variance(variance) {};
			
			double getDensity(double x) const override;
			double getRandom() const override;
	};
	
	class InverseGamma : public Base<double> {
		protected:
			double shape;
			double scale;
		public:
			InverseGamma(double shape, double scale):
				shape(shape), scale(scale) {};
			
			double getDensity(double x) const override;
			double getRandom() const override;
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
		DomainValue getRandom() const {return d->getRandom();};
};

#endif
