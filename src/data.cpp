#include "data.hpp"

double Data::getResponseMean() const {
	return response.mean();
}

double Data::getResponseVariance() const {
	Eigen::VectorXd residuals = response - Eigen::VectorXd::Constant(response.size(), response.mean());
	return residuals.dot(residuals) / residuals.size();
}
