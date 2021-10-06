#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "grouping.hpp"
#include "data.hpp"

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
		
		//To ease matrix manipulation, the set of parameters will just be represented as a vector.
		//In reality, it will consist of:
		//one growth rate for each growth rate group, then
		//one competition coefficient for each pair of a row group and a column group, in row-major order.
		typedef Eigen::VectorXd ParameterVector;
		
		//Similarly, the Jacobian will be represented simply as a matrix.
		//It has a number of rows equal to the number of observations, and
		//a number of columns equal to the number of parameters, as given above.
		typedef Eigen::MatrixXd Jacobian;
	protected:
		Eigen::MatrixXd colGroupedDesign;
		bool isDirtyColGroupedDesign = true;
		
		ParameterVector solution;
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
		size_t getGrowthRateIndex(size_t growthGroup) const;
		double getGrowthRate(const ParameterVector &parameters, size_t growthGroup) const;
		Eigen::VectorXd getGrowthRates(const ParameterVector &parameters) const;
		size_t getCompetitionCoefficientIndex(size_t rowGroup, size_t colGroup) const;
		double getCompetitionCoefficient(const ParameterVector &parameters, size_t rowGroup, size_t colGroup) const;
		Eigen::VectorXd getCompetitionCoefficientsRow(const ParameterVector &parameters, size_t rowGroup) const;
		
		void calculateColGroupedDesign();
		Eigen::MatrixXd getColGroupedDesign();
		
		ParameterVector getInitialParameterValues() const;
		ParameterVector getParameterTolerances() const;
		
		Eigen::VectorXd getPredictions(const ParameterVector &parameters);
		Eigen::VectorXd getResiduals(const ParameterVector &parameters);
		
		Jacobian getJacobian(const ParameterVector &parameters);
		
		void calculateSolution();
	public:
		ParameterVector getSolution();
		Eigen::VectorXd getSolutionPredictions();
		Eigen::VectorXd getSolutionResiduals();
		double getDeviance();
		double getAIC();
		double getR2();
};

#endif
