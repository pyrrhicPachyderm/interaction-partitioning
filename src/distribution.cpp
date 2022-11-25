#include <math.h>
#include <random>
#include "distribution.hpp"

#define RANDOM_SEED 42

static std::default_random_engine randomNumberGenerator(RANDOM_SEED);

double Distributions::Uniform::getDensity(double x) const {
	if(x >= min && x < max) return 1.0 / (max - min);
	else return 0;
}

double Distributions::Uniform::getRandom() const {
	return std::uniform_real_distribution(min, max)(randomNumberGenerator);
}

double Distributions::Normal::getDensity(double x) const {
	double residual = x - mean;
	return 1.0 / sqrt(2 * M_PI * variance) * exp(-0.5 * residual*residual / variance);
}

double Distributions::Normal::getRandom() const {
	return std::normal_distribution(mean, sqrt(variance))(randomNumberGenerator);
}
