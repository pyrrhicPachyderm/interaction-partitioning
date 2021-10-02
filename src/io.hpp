#ifndef IO_HPP
#define IO_HPP

#include <Eigen/Core>

extern std::vector<size_t> readIndexVector(const char *filename);
extern Eigen::VectorXd readDoubleVector(const char *filename);
extern Eigen::MatrixXd readDoubleMatrix(const char *filename);

#endif
