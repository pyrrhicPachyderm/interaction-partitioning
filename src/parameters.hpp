#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include "grouping.hpp"
#include "data.hpp"

enum GroupingType {GROWTH, ROW, COL, NUM_GROUPING_TYPES};
typedef std::array<Grouping, NUM_GROUPING_TYPES> GroupingSet;

class Parameters {
	protected:
		//The grouping sizes will be implicitly stored in the sizes of the growth rates vector and competition coefficients matrix.
		//The matrix will be stored in row-major order.
		Eigen::VectorXd growthRates;
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> competitionCoefficients;
	public:
		double getGrowthRate(size_t index) const;
		const Eigen::VectorXd &getGrowthRates() const;
		double getCompetitionCoefficient(size_t rowIndex, size_t colIndex) const;
		Eigen::VectorXd getCompetitionCoefficientsRow(size_t rowIndex) const;
		size_t getNumParameters() const;
	public:
		Parameters() = default;
		Parameters(Data data, GroupingSet groupings);
		
		//Some functions to allow converting back and forth as a pure vector, for systems (such as the NLS module) that need it that way.
		Parameters(Eigen::VectorXd parameters, GroupingSet groupings);
		Eigen::VectorXd getAsVector() const;
		size_t getAsVectorGrowthRateIndex(size_t index) const;
		size_t getAsVectorCompetitionCoefficientIndex(size_t rowIndex, size_t colIndex) const;
		
		//Also for the NLS module and the like, that want tolerances for each value in the pure vector.
		static Eigen::VectorXd getTolerances(Data data, GroupingSet groupings);
};

#endif
