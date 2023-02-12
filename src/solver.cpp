#include "solver.hpp"

Eigen::MatrixXd Solver::getColGroupedDesign() const {
	return data.getColGroupedDesign(getGrouping(COL));
}
