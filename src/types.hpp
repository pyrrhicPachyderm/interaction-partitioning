#ifndef TYPES_HPP
#define TYPES_HPP

#include <Eigen/Core>

//A Jacobian matrix.
//It has a number of rows equal to the number of response variables, and
//a number of columns equal to the number of parameters.
typedef Eigen::MatrixXd Jacobian;

namespace Eigen {
	//As MatrixXd, but row-major.
	typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXdRowMajor;
}

#endif
