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
	protected:
		const Eigen::VectorXd &getObservations() {
			return data.getObservations();
		}
		Eigen::VectorXd getPredictions(const Parameters &parameters, const GroupingSet &groupings) const {
			return data.getPredictions(model, parameters, groupings);
		}
		Eigen::VectorXd getResiduals(const Parameters &parameters, const GroupingSet &groupings) const {
			return data.getResiduals(model, parameters, groupings);
		}
		Jacobian getPredictionsJacobian(const Parameters &parameters, const GroupingSet &groupings) const {
			return data.getPredictionsJacobian(model, parameters, groupings);
		}
		Jacobian getResidualsJacobian(const Parameters &parameters, const GroupingSet &groupings) const {
			return data.getResidualsJacobian(model, parameters, groupings);
		}
};

#endif
