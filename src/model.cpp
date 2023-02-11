#include "model.hpp"

double Models::LotkaVolterra::getDerivative(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients) const {
	return focalDensity * focalGrowthRate * (1.0 - densities.dot(competitionCoefficients));
}
