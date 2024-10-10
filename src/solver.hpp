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
		constexpr static size_t NUM_ADDITIONAL_PARAMETERS = 0;
		typedef Parameters ParametersT;
		
		Solver(Model model, Data data): model(model), data(data) {};
	protected:
		const Eigen::VectorXd &getObservations() const {
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

template<typename ErrDistT> class GeneralisedSolver : public Solver {
	//GeneralisedSolver is templated on the type of error distribution, which must inherit from Distributions::Base<double>.
	//The first parameter in the constructor of the given distribution MUST be the mean (e.g. use Gamma2 instead of Gamma).
	//Every other parameter will be estimated, and is stored in the additionalParameters of an AugmentedParameters object.
	static_assert(std::is_base_of_v<Distributions::Base<double>, ErrDistT>);
	public:
		//TODO: Use std::tuple_size<refl::ctor_as_tuple<ErrDistT>>{} - 1 for NUM_ADDITIONAL_PARAMETERS.
		//This works for Normal and Gamma2, but not DiscreteWrapper<NegativeBinomial>.
		//I'm not sure why; refl is terribly complicated, after all.
		constexpr static size_t NUM_ADDITIONAL_PARAMETERS = 1;
		typedef AugmentedParameters<NUM_ADDITIONAL_PARAMETERS> ParametersT;
		typedef ParametersT::AdditionalParametersVector AdditionalParametersVector;
		
		using Solver::Solver;
	protected:
		double getLogLikelihood(const AugmentedParameters<NUM_ADDITIONAL_PARAMETERS> &parameters, const GroupingSet &groupings) const {
			const Eigen::VectorXd &observations = getObservations();
			Eigen::VectorXd predictions = getPredictions(parameters, groupings);
			
			double result = 0.0;
			for(size_t i = 0; i < (size_t)observations.size(); i++) {
				result += std::make_from_tuple<ErrDistT>(std::tuple_cat(std::tuple(predictions[i]), parameters.getAdditionalParameters())).getLogDensity(observations[i]);
			}
			return result;
		};
};

#endif
