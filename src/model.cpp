#include "model.hpp"

double Models::LotkaVolterra::getDerivative(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients) const {
	return focalDensity * focalGrowthRate * (1.0 - densities.dot(competitionCoefficients));
}

double Models::LotkaVolterra::getGrowthRateJacobian(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients) const {
	return focalDensity * (1.0 - densities.dot(competitionCoefficients));
}

double Models::LotkaVolterra::getCompetitionCoefficientJacobian(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients, size_t index) const {
	return focalDensity * focalGrowthRate * (- densities[index]);
}
