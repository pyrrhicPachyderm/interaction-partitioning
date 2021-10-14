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
	protected:
		Eigen::MatrixXd colGroupedDesign;
		bool isDirtyColGroupedDesign = true;
		
		//We will have functions to say that particular elements have changed and mark appropriate things as dirty.
		//There needs to be a set for the superclass, Solver, as well as ones for subclasses.
		//The subclass ones will be define pure virtual here, and will need overriding.
		//The superclass ones will dirty elements of the superclass, then call the subclass ones.
		virtual void dirtyDataSubclass() = 0;
		virtual void dirtyGrowthGroupingSubclass() = 0;
		virtual void dirtyRowGroupingSubclass() = 0;
		virtual void dirtyColGroupingSubclass() = 0;
		
		void dirtyData() {
			isDirtyColGroupedDesign = true;
			dirtyDataSubclass();
		}
		void dirtyGrowthGrouping() {
			dirtyGrowthGroupingSubclass();
		}
		void dirtyRowGrouping() {
			dirtyRowGroupingSubclass();
		}
		void dirtyColGrouping() {
			isDirtyColGroupedDesign = true;
			dirtyColGroupingSubclass();
		}
	public:
		//Functions to retrieve groupings.
		//Subclasses each store groupings in their own way, so these are pure virtual.
		virtual const Grouping &getGrowthGrouping() const = 0;
		virtual const Grouping &getRowGrouping() const = 0;
		virtual const Grouping &getColGrouping() const = 0;
	protected:
		void calculateColGroupedDesign();
		Eigen::MatrixXd getColGroupedDesign();
		
		Eigen::VectorXd getPredictions(const Parameters &parameters);
		Eigen::VectorXd getResiduals(const Parameters &parameters);
};

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
