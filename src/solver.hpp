#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "grouping.hpp"
#include "data.hpp"
#include "parameters.hpp"

class Solver {
	protected:
		Data data;
		Grouping growthGrouping;
		Grouping rowGrouping;
		Grouping colGrouping;
	public:
		Solver(Data data):
			data(data), growthGrouping(Grouping(data.numSpecies)), rowGrouping(Grouping(data.numSpecies)), colGrouping(Grouping(data.numSpecies)) {};
		Solver(Data data, Grouping growthGrouping, Grouping rowGrouping, Grouping colGrouping):
			data(data), growthGrouping(growthGrouping), rowGrouping(rowGrouping), colGrouping(colGrouping)
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
		Eigen::MatrixXd colGroupedDesign;
		bool isDirtyColGroupedDesign = true;
		
		Parameters solution;
		bool isDirtySolution = true;
		
		//Functions to say that particular elements have changed and mark appropriate things as dirty.
		void dirtyData() {
			isDirtySolution = true;
			isDirtyColGroupedDesign = true;
		}
		void dirtyGrowthGrouping() {
			isDirtySolution = true;
		}
		void dirtyRowGrouping() {
			isDirtySolution = true;
		}
		void dirtyColGrouping() {
			isDirtySolution = true;
			isDirtyColGroupedDesign = true;
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
		Grouping getGrowthGrouping() const {
			return growthGrouping;
		}
		Grouping getRowGrouping() const {
			return rowGrouping;
		}
		Grouping getColGrouping() const {
			return colGrouping;
		}
	protected:
		void calculateColGroupedDesign();
		Eigen::MatrixXd getColGroupedDesign();
		
		Eigen::VectorXd getPredictions(const Parameters &parameters);
		Eigen::VectorXd getResiduals(const Parameters &parameters);
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
