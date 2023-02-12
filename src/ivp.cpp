#include "ivp.hpp"

Eigen::VectorXd IVPStepFuncs::forwardEuler(IVPDerivativeFunc f, Eigen::VectorXd y0, double t0, double h) {
	return y0 + h * f(t0, y0);
}

Eigen::VectorXd IVPStepFuncs::rungeKuttaFour(IVPDerivativeFunc f, Eigen::VectorXd y0, double t0, double h) {
	Eigen::VectorXd k1 = f(t0, y0);
	Eigen::VectorXd k2 = f(t0 + h/2, y0 + h/2 * k1);
	Eigen::VectorXd k3 = f(t0 + h/2, y0 + h/2 * k2);
	Eigen::VectorXd k4 = f(t0 + h, y0 + h * k3);
	return y0 + h * (k1 + 2*k2 + 2*k3 + k4)/6;
}

Eigen::VectorXd solveIVP(IVPDerivativeFunc f, Eigen::VectorXd y0, double t0, double tMax, size_t numSteps, IVPStepFunc stepFunc) {
	double h = (tMax - t0) / numSteps;
	for(size_t i = 0; i < numSteps; i++) {
		double t = t0 + i * h;
		y0 = stepFunc(f, y0, t, h);
	}
	return y0;
}
