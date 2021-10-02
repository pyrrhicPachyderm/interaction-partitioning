#ifndef NLS_HPP
#define NLS_HPP

#include <Eigen/Core>

//These functions take a vector of parameters, and return:
//a vector of residuals, and
//the Jacobian of the residuals with respect to the parameters.
typedef std::function<Eigen::VectorXd(Eigen::VectorXd parameters)> ResidualsFunc;
typedef std::function<Eigen::MatrixXd(Eigen::VectorXd parameters)> JacobianFunc;

//Methods require a function to calculate the residuals, a function to calculate the Jacobian, an initial guess at the parameters, and tolerances for each parameter.
extern Eigen::VectorXd gaussNewtonNLS(ResidualsFunc residualsFunc, JacobianFunc jacobianFunc, Eigen::VectorXd parameters, Eigen::VectorXd tolerances);

#endif
