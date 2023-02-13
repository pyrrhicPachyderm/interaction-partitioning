#ifndef TYPES_HPP
#define TYPES_HPP

#include <Eigen/Core>

namespace Eigen {
	//As MatrixXd, but row-major.
	typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXdRowMajor;
}

#endif
