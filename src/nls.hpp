#ifndef NLS_HPP
#define NLS_HPP

#include "types.hpp"

//These functions take a vector of parameters, and return:
//a vector of residuals, and
//the Jacobian of the residuals with respect to the parameters.
typedef std::function<Eigen::VectorXd(const Eigen::VectorXd &parameters)> ResidualsFunc;
typedef std::function<Jacobian(const Eigen::VectorXd &parameters)> JacobianFunc;

//Methods require a function to calculate the residuals, a function to calculate the Jacobian, an initial guess at the parameters, and tolerances for each parameter.
extern Eigen::VectorXd gaussNewtonNLS(ResidualsFunc residualsFunc, JacobianFunc jacobianFunc, Eigen::VectorXd parameters, const Eigen::VectorXd &tolerances);

#endif
