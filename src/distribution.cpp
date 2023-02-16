#include <assert.h>
#include "distribution.hpp"

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

double Distributions::InverseGamma::getRandom(RandomGenerator &generator) const {
	//TODO: Implement.
	assert(false);
}
