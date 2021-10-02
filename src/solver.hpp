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
		double getGrowthRate(ParameterVector parameters, size_t growthGroup) const;
		Eigen::VectorXd getGrowthRates(ParameterVector parameters) const;
		double getCompetitionCoefficient(ParameterVector parameters, size_t rowGroup, size_t colGroup) const;
		Eigen::VectorXd getCompetitionCoefficientsRow(ParameterVector parameters, size_t rowGroup) const;
		
		Eigen::MatrixXd getColGroupedDesign() const;
		
		ParameterVector getInitialParameterValues() const;
		
		Eigen::VectorXd getPredictions(ParameterVector parameters) const;
		Eigen::VectorXd getResiduals(ParameterVector parameters) const;
};

#endif
