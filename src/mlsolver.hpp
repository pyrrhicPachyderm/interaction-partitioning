#ifndef MAXIMUM_LIKELIHOOD_SOLVER_HPP
#define MAXIMUM_LIKELIHOOD_SOLVER_HPP

#include "solver.hpp"

class MaximumLikelihoodSolver : public Solver {
	protected:
		GroupingSet groupings;
	public:
		MaximumLikelihoodSolver(Model model, Data data):
			Solver(model, data), groupings({Grouping(data.getNumRowSpecies()), Grouping(data.getNumRowSpecies()), Grouping(data.getNumColSpecies())}) {};
		MaximumLikelihoodSolver(Model model, Data data, Grouping growthGrouping, Grouping rowGrouping, Grouping colGrouping):
			Solver(model, data), groupings({growthGrouping, rowGrouping, colGrouping})
		{
			for(size_t i = 0; i < NUM_GROUPING_TYPES; i++) {
				assert(groupings[i].numSpecies == data.getNumSpecies((GroupingType)i));
			}
		}
	protected:
		Parameters solution;
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
		Eigen::VectorXd getResidualsFromVector(const Eigen::VectorXd &parameterVector) const;
		Datasets::Jacobian getResidualsJacobianFromVector(const Eigen::VectorXd &parameterVector) const;
		
		void calculateSolution();
	public:
		Parameters getSolution();
		Eigen::VectorXd getSolutionPredictions();
		Eigen::VectorXd getSolutionResiduals();
		double getDeviance();
		double getAIC();
		double getAICc();
		double getR2();
};

#endif
