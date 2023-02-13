#include <Eigen/SVD>
#include "nls.hpp"

static bool checkTolerances(const Eigen::VectorXd &oldParameters, const Eigen::VectorXd &newParameters, const Eigen::VectorXd &tolerances) {
	for(size_t i = 0; i < (size_t)tolerances.size(); i++) {
		if(abs(newParameters[i] - oldParameters[i]) > tolerances[i]) return false;
	}
	return true;
}

static Eigen::VectorXd gaussNewtonShift(ResidualsFunc residualsFunc, JacobianFunc jacobianFunc, const Eigen::VectorXd &parameters) {
	//See https://eigen.tuxfamily.org/dox/group__LeastSquares.html for a comparison of methods.
	Eigen::VectorXd residuals = residualsFunc(parameters);
	Jacobian jacobian = jacobianFunc(parameters);
	return jacobian.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(residuals);
}

Eigen::VectorXd gaussNewtonNLS(ResidualsFunc residualsFunc, JacobianFunc jacobianFunc, Eigen::VectorXd parameters, const Eigen::VectorXd &tolerances) {
	while(true) {
		Eigen::VectorXd shift = gaussNewtonShift(residualsFunc, jacobianFunc, parameters);
		Eigen::VectorXd newParameters = parameters - shift;
		//TODO: Check for divergence, as described on the Wikipedia page.
		if(checkTolerances(parameters, newParameters, tolerances)) return newParameters;
		parameters = newParameters;
	}
}
