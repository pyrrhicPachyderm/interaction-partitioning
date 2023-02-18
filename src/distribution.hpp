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
		protected:
			bool hasValidParameters = true;
			virtual bool isInDomain(DomainValue x) const {return true;};
			
			virtual double calculateDensity(DomainValue x) const = 0;
			virtual double calculateLogDensity(DomainValue x) const {return log(calculateDensity(x));}; //A base function that may be overwritten if there's a better way for a given distribution.
		public:
			Base() = default;
			Base(bool hasValidParameters):
				hasValidParameters(hasValidParameters) {};
			
			double getDensity(DomainValue x) const {
				if(!hasValidParameters || !isInDomain(x)) return 0.0;
				return calculateDensity(x);
			}
			double getLogDensity(DomainValue x) const {
				if(!hasValidParameters || !isInDomain(x)) return -INFINITY;
				return calculateLogDensity(x);
			}
			//TODO: What to do about getRandom if !hasValidParameters?
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
			
			double getRandom(RandomGenerator &generator) const override;
		protected:
			double calculateDensity(double x) const override;
	};
	
	class Normal : public Base<double> {
		protected:
			double mean;
			double variance;
		public:
			Normal(double mean, double variance):
				Base(variance > 0), mean(mean), variance(variance) {};
			
			double getRandom(RandomGenerator &generator) const override;
		protected:
			double calculateDensity(double x) const override;
			double calculateLogDensity(double x) const override;
	};
	
	class InverseGamma : public Base<double> {
		protected:
			double shape;
			double scale;
			
			bool isInDomain(double x) const override {return x > 0.0;};
		public:
			InverseGamma(double shape, double scale):
				Base(shape > 0 && scale > 0), shape(shape), scale(scale) {};
			
			double getRandom(RandomGenerator &generator) const override;
		protected:
			double calculateDensity(double x) const override;
			double calculateLogDensity(double x) const override;
	};
	
	class Gamma : public Base<double> {
		protected:
			double shape;
			double scale;
			
			bool isInDomain(double x) const override {return x > 0.0;};
		public:
			Gamma(double shape, double scale):
				Base(shape > 0 && scale > 0), shape(shape), scale(scale) {};
			
			double getRandom(RandomGenerator &generator) const override;
		protected:
			double calculateDensity(double x) const override;
			double calculateLogDensity(double x) const override;
	};
	
	class Gamma2 : public Gamma {
		//Parameterised by a mean and dispersion parameter.
		public:
			Gamma2(double mean, double dispersion):
				Gamma(1 / dispersion, mean * dispersion) {};
	};
	
	class NegativeBinomial : public Base<int> {
		protected:
			double r;
			double p;
			
			bool isInDomain(int x) const override {return x >= 0;};
		public:
			NegativeBinomial(double r, double p):
				Base(r > 0 && p >= 0 && p <= 1), r(r), p(p) {};
			
			int getRandom(RandomGenerator &generator) const override;
		protected:
			double calculateDensity(int x) const override;
			double calculateLogDensity(int x) const override;
	};
	
	class NegativeBinomial2 : public NegativeBinomial {
		//Parameterised by a mean and dispersion parameter.
		public:
			NegativeBinomial2(double mean, double dispersion):
				NegativeBinomial(1 / dispersion, (1 / dispersion) / (mean + 1 / dispersion)) {};
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
