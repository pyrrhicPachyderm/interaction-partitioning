#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "grouping.hpp"
#include "data.hpp"
#include "parameters.hpp"

class Solver {
	protected:
		Model model;
		Data data;
	public:
		Solver(Model model, Data data): model(model), data(data) {};
		
		const Grouping getGrouping(GroupingType groupingType) const {
			return getGroupings()[groupingType];
		}
		virtual const GroupingSet &getGroupings() const = 0;
	protected:
		Eigen::MatrixXd getColGroupedDesign() const;
		
		Eigen::VectorXd getPredictions(const Parameters &parameters) const {
			return data.getPredictions(model, parameters, getGroupings());
		}
		Eigen::VectorXd getResiduals(const Parameters &parameters) const {
			return data.getResiduals(model, parameters, getGroupings());
		}
		Data::Jacobian getPredictionsJacobian(const Parameters &parameters) const {
			return data.getPredictionsJacobian(model, parameters, getGroupings());
		}
		Data::Jacobian getResidualsJacobian(const Parameters &parameters) const {
			return data.getResidualsJacobian(model, parameters, getGroupings());
		}
};

#endif
