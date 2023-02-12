#ifndef IVP_HPP
#define IVP_HPP

#include <Eigen/Core>

//A module to solve initial value problems of the form.
//dy/dt = f(t, y(t))

//This is dy/dt = f(t, y(t))
typedef std::function<Eigen::VectorXd(double t, const Eigen::VectorXd &y)> IVPDerivativeFunc;

//A function for stepping the solution to an IVP, from t=t0 to t=t0+h.
typedef std::function<Eigen::VectorXd(IVPDerivativeFunc f, Eigen::VectorXd y0, double t0, double h)> IVPStepFunc;

namespace IVPStepFuncs {
	extern Eigen::VectorXd forwardEuler(IVPDerivativeFunc f, Eigen::VectorXd y0, double t0, double h);
	extern Eigen::VectorXd rungeKuttaFour(IVPDerivativeFunc f, Eigen::VectorXd y0, double t0, double h);
}

//The function to apply an IVPStepFunc the appropriate number of times to solve an IVP.
extern Eigen::VectorXd solveIVP(IVPDerivativeFunc f, Eigen::VectorXd y0, double t0, double tMax, size_t numSteps, IVPStepFunc stepFunc);

#endif
