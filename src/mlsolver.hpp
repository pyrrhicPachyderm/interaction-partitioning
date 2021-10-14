#ifndef MAXIMUM_LIKELIHOOD_SOLVER_HPP
#define MAXIMUM_LIKELIHOOD_SOLVER_HPP

#include "solver.hpp"

class MaximumLikelihoodSolver : public Solver {
	protected:
		Grouping growthGrouping;
		Grouping rowGrouping;
		Grouping colGrouping;
	public:
		MaximumLikelihoodSolver(Data data):
			Solver(data), growthGrouping(Grouping(data.numSpecies)), rowGrouping(Grouping(data.numSpecies)), colGrouping(Grouping(data.numSpecies)) {};
		MaximumLikelihoodSolver(Data data, Grouping growthGrouping, Grouping rowGrouping, Grouping colGrouping):
			Solver(data), growthGrouping(growthGrouping), rowGrouping(rowGrouping), colGrouping(colGrouping)
		{
			assert(growthGrouping.numSpecies == data.numSpecies);
			assert(rowGrouping.numSpecies == data.numSpecies);
			assert(colGrouping.numSpecies == data.numSpecies);
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
		void dirtyDataSubclass() {
			isDirtySolution = true;
		}
		void dirtyGrowthGroupingSubclass() {
			isDirtySolution = true;
		}
		void dirtyRowGroupingSubclass() {
			isDirtySolution = true;
		}
		void dirtyColGroupingSubclass() {
			isDirtySolution = true;
		}
	public:
		//Functions to update groupings.
		//Can use reset(), separate(), or advance().
		template<typename T> T updateGrowthGrouping(T (Grouping::*updateFunc)()) {
			dirtyGrowthGrouping();
			return (growthGrouping.*updateFunc)();
		}
		template<typename T> T updateRowGrouping(T (Grouping::*updateFunc)()) {
			dirtyRowGrouping();
			return (rowGrouping.*updateFunc)();
		}
		template<typename T> T updateColGrouping(T (Grouping::*updateFunc)()) {
			dirtyColGrouping();
			return (colGrouping.*updateFunc)();
		}
		
		//Functions to retrieve groupings.
		const Grouping &getGrowthGrouping() const {
			return growthGrouping;
		}
		const Grouping &getRowGrouping() const {
			return rowGrouping;
		}
		const Grouping &getColGrouping() const {
			return colGrouping;
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
		double getR2();
};

#endif
