#ifndef MAXIMUM_LIKELIHOOD_SOLVER_HPP
#define MAXIMUM_LIKELIHOOD_SOLVER_HPP

#include "solver.hpp"

class MaximumLikelihoodSolverInterface {
	public:
		virtual Parameters getSolutionParameters() = 0;
		virtual double getSolutionAdditionalParameter(size_t i) = 0;
		
		virtual Eigen::VectorXd getSolutionPredictions() = 0;
		virtual Eigen::VectorXd getSolutionResiduals() = 0;
		
		virtual double getDeviance() = 0;
		virtual size_t getNumParameters() = 0;
		virtual double getAIC() = 0;
		virtual double getAICc() = 0;
		virtual double getR2() = 0;
};

template<typename SolverT> class MaximumLikelihoodSolver : public SolverT, public MaximumLikelihoodSolverInterface {
	static_assert(std::is_base_of_v<Solver, SolverT>);
	protected:
		GroupingSet groupings;
	public:
		MaximumLikelihoodSolver(Model model, Data data):
			SolverT(model, data), groupings({Grouping(data.getNumRowSpecies()), Grouping(data.getNumRowSpecies()), Grouping(data.getNumColSpecies())}) {};
		MaximumLikelihoodSolver(Model model, Data data, GroupingSet groupings):
			SolverT(model, data),
			groupings(groupings)
		{
			for(size_t i = 0; i < NUM_GROUPING_TYPES; i++) {
				assert(groupings[i].numSpecies == data.getNumSpecies((GroupingType)i));
			}
		};
	protected:
		SolverT::ParametersT solution;
		bool isDirtySolution = true;
	public:
		//Functions to update groupings.
		//Can use reset(), separate(), or advance().
		template<typename T> T updateGrouping(GroupingType groupingType, T (Grouping::*updateFunc)()) {
			isDirtySolution = true;
			return (groupings[groupingType].*updateFunc)();
		}
		
		const GroupingSet &getGroupings() const {return groupings;}
		const Grouping &getGrouping(GroupingType groupingType) const {return groupings[groupingType];}
	protected:
		virtual void calculateSolution() = 0;
	public:
		Parameters getSolutionParameters() override;
		double getSolutionAdditionalParameter(size_t i) override;
		
		Eigen::VectorXd getSolutionPredictions() override;
		Eigen::VectorXd getSolutionResiduals() override;
		
		virtual double getDeviance() = 0;
		size_t getNumParameters() override;
		double getAIC() override;
		double getAICc() override;
		double getR2() override; //TODO: Check whether this needs to be different with a non-normal error distribution.
};

class GaussNewtonSolver : public MaximumLikelihoodSolver<Solver> {
	public:
		using MaximumLikelihoodSolver<Solver>::MaximumLikelihoodSolver;
	protected:
		Eigen::VectorXd getResidualsFromVector(const Eigen::VectorXd &parameterVector) const;
		Jacobian getResidualsJacobianFromVector(const Eigen::VectorXd &parameterVector) const;
		
		void calculateSolution() override;
	public:
		double getDeviance() override;
};

#endif
