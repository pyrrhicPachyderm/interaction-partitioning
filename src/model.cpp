#include "model.hpp"

Eigen::VectorXd Models::Base::getResponses(const Eigen::VectorXd &densities, const Eigen::VectorXd &growthRates, const Eigen::MatrixXdRowMajor &competitionCoefficients) const {
	Eigen::VectorXd derivatives(densities.size());
	for(size_t i = 0; i < (size_t)derivatives.size(); i++) {
		derivatives[i] = getResponse(densities[i], growthRates[i], densities, competitionCoefficients.row(i));
	}
	return derivatives;
}

double Models::LotkaVolterra::getResponse(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients) const {
	return focalDensity * focalGrowthRate * (1.0 - densities.dot(competitionCoefficients));
}

double Models::LotkaVolterra::getGrowthRateJacobian(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients) const {
	return focalDensity * (1.0 - densities.dot(competitionCoefficients));
}

double Models::LotkaVolterra::getCompetitionCoefficientJacobian(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients, size_t index) const {
	return focalDensity * focalGrowthRate * (- densities[index]);
}

double Models::BevertonHolt::getResponse(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients) const {
	return focalDensity * focalGrowthRate / (1.0 + densities.dot(competitionCoefficients));
}

double Models::BevertonHolt::getGrowthRateJacobian(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients) const {
	return focalDensity / (1.0 + densities.dot(competitionCoefficients));
}

double Models::BevertonHolt::getCompetitionCoefficientJacobian(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients, size_t index) const {
	return - focalDensity * focalGrowthRate * densities[index] / pow(1.0 + densities.dot(competitionCoefficients), 2);
}
