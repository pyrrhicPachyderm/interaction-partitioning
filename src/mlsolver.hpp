#ifndef MAXIMUM_LIKELIHOOD_SOLVER_HPP
#define MAXIMUM_LIKELIHOOD_SOLVER_HPP

#include "solver.hpp"

class MaximumLikelihoodSolver : public Solver {
	protected:
		GroupingSet groupings;
	public:
		MaximumLikelihoodSolver(Data data):
			Solver(data), groupings({Grouping(data.getNumRowSpecies()), Grouping(data.getNumRowSpecies()), Grouping(data.getNumColSpecies())}) {};
		MaximumLikelihoodSolver(Data data, Grouping growthGrouping, Grouping rowGrouping, Grouping colGrouping):
			Solver(data), groupings({growthGrouping, rowGrouping, colGrouping})
		{
			for(size_t i = 0; i < NUM_GROUPING_TYPES; i++) {
				assert(groupings[i].numSpecies == data.getNumSpecies((GroupingType)i));
			}
		}
		
		//The Jacobian will be represented simply as a matrix.
		//It has a number of rows equal to the number of observations, and
		//a number of columns equal to the number of parameters.
		//It treats parameters in the order given by Parameters::getAsVector().
		//TODO: Do this more nicely.
		typedef Eigen::MatrixXd Jacobian;
	protected:
		Parameters solution;
		bool isDirtySolution = true;
		
		//Functions to say that particular elements have changed and mark appropriate things as dirty.
		void dirtyDataSubclass() override {
			isDirtySolution = true;
		}
		void dirtyGroupingSubclass(GroupingType groupingType) override {
			isDirtySolution = true;
		}
	public:
		//Functions to update groupings.
		//Can use reset(), separate(), or advance().
		template<typename T> T updateGrouping(GroupingType groupingType, T (Grouping::*updateFunc)()) {
			dirtyGrouping(groupingType);
			return (groupings[groupingType].*updateFunc)();
		}
		
		//Functions to retrieve groupings.
		const Grouping &getGrouping(GroupingType groupingType) const override {
			return groupings[groupingType];
		}
		GroupingSet getGroupings() const {
			return groupings;
		}
	protected:
		Eigen::VectorXd getResidualsFromVector(const Eigen::VectorXd &parameterVector);
		
		Jacobian getJacobian(const Parameters &parameters);
		Jacobian getJacobianFromVector(const Eigen::VectorXd &parameterVector);
		
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
