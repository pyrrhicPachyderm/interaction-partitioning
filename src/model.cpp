#include "model.hpp"

Eigen::VectorXd Models::Base::getDerivatives(const Eigen::VectorXd &densities, const Eigen::VectorXd &growthRates, const Eigen::MatrixXd &competitionCoefficients) const {
	Eigen::VectorXd derivatives(densities.size());
	for(size_t i = 0; i < (size_t)derivatives.size(); i++) {
		derivatives[i] = getDerivative(densities[i], growthRates[i], densities, competitionCoefficients.row(i));
	}
	return derivatives;
}

double Models::LotkaVolterra::getDerivative(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients) const {
	return focalDensity * focalGrowthRate * (1.0 - densities.dot(competitionCoefficients));
}

double Models::LotkaVolterra::getGrowthRateJacobian(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients) const {
	return focalDensity * (1.0 - densities.dot(competitionCoefficients));
}

double Models::LotkaVolterra::getCompetitionCoefficientJacobian(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients, size_t index) const {
	return focalDensity * focalGrowthRate * (- densities[index]);
}
