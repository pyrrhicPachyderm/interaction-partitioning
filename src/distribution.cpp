#include <assert.h>
#include "distribution.hpp"

//NB: lgamma() is the log of the gamma function, while tgamma() is the true gamma function.
//gamma() should be avoided; it refers to lgamma() on my machine, but apparently not all machines.

double Distributions::Uniform::getDensity(double x) const {
	if(x >= min && x < max) return 1.0 / (max - min);
	else return 0;
}

double Distributions::Uniform::getRandom(RandomGenerator &generator) const {
	return std::uniform_real_distribution(min, max)(generator);
}

double Distributions::Normal::getDensity(double x) const {
	double residual = x - mean;
	return 1.0 / sqrt(2 * M_PI * variance) * exp(-0.5 * residual*residual / variance);
}

double Distributions::Normal::getLogDensity(double x) const {
	double residual = x - mean;
	return -0.5 * (log(variance) + log(2 * M_PI) + residual*residual / variance);
}

double Distributions::Normal::getRandom(RandomGenerator &generator) const {
	return std::normal_distribution(mean, sqrt(variance))(generator);
}

double Distributions::InverseGamma::getDensity(double x) const {
	return pow(scale, shape) / tgamma(shape) * pow(1 / x, shape + 1) * exp(-scale / x);
}

double Distributions::InverseGamma::getLogDensity(double x) const {
	return log(scale) * shape - lgamma(shape) - log(x) * (shape + 1) - scale / x;
}

double Distributions::InverseGamma::getRandom(RandomGenerator &generator) const {
	//TODO: Implement.
	assert(false);
}

double Distributions::Gamma::getDensity(double x) const {
	return 1 / (tgamma(shape) * pow(scale, shape)) * pow(x, shape - 1) * exp(-x / scale);
}

double Distributions::Gamma::getLogDensity(double x) const {
	return - lgamma(shape) - log(scale) * shape + log(x) * (shape - 1) - x / scale;
}

double Distributions::Gamma::getRandom(RandomGenerator &generator) const {
	return std::normal_distribution(shape, scale)(generator);
}
