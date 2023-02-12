#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "grouping.hpp"
#include "data.hpp"
#include "parameters.hpp"

class Solver {
	protected:
		Data data;
	public:
		Solver(Data data): data(data) {};
		
		const Grouping getGrouping(GroupingType groupingType) const {
			return getGroupings()[groupingType];
		}
		virtual const GroupingSet &getGroupings() const = 0;
	protected:
		Eigen::MatrixXd getColGroupedDesign() const;
		
		Eigen::VectorXd getPredictions(const Parameters &parameters) const;
		Eigen::VectorXd getResiduals(const Parameters &parameters) const;
};

#endif
