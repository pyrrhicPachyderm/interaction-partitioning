#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include "grouping.hpp"
#include "groupingmove.hpp"
#include "data.hpp"
#include "distribution.hpp"

class Parameters {
	protected:
		size_t numRowSpecies = 0;
		size_t numColSpecies = 0;
		//The grouping sizes will be implicitly stored in the sizes of the growth rates vector and competition coefficients matrix.
		//The matrix will be stored in row-major order.
		Eigen::VectorXd growthRates;
		Eigen::MatrixXdRowMajor competitionCoefficients;
	public:
		double getNumRowSpecies() const;
		double getNumColSpecies() const;
		double getGrowthRate(size_t index) const;
		const Eigen::VectorXd &getGrowthRates() const;
		double getCompetitionCoefficient(size_t rowIndex, size_t colIndex) const;
		Eigen::VectorXd getCompetitionCoefficientsRow(size_t rowIndex) const;
		const Eigen::MatrixXdRowMajor &getCompetitionCoefficients() const;
		size_t getNumParameters() const;
	public:
		Parameters() = default;
		Parameters(Data data, GroupingSet groupings);
		Parameters(Parameters p, GroupingSet groupings); //Undoes the grouping, giving a full set of parameters for every species.
		
		//Some functions to allow converting back and forth as a pure vector, for systems (such as the NLS module) that need it that way.
		Parameters(Eigen::VectorXd parameters, GroupingSet groupings);
		Eigen::VectorXd getAsVector() const;
		size_t getAsVectorGrowthRateIndex(size_t index) const;
		size_t getAsVectorCompetitionCoefficientIndex(size_t rowIndex, size_t colIndex) const;
		
		//Also for the NLS module and the like, that want tolerances for each value in the pure vector.
		static Eigen::VectorXd getTolerances(Data data, GroupingSet groupings);
	public:
		//Splits or merges parameters as appropriate, returning the relevant component of the acceptance ratio.
		//That is, the Jacobian determinant, times the jumping density ratio of the additional random variables u^(k) and u^(k').
		double moveModel(GroupingType groupingType, MoveType moveType, const GroupingMove &groupingMove, Distribution<double> randomVariableDistribution, RandomGenerator &randomGenerator);
		
		//TODO: If we were to use a non-symmetric jumping density, this would need to return a component of the acceptance ratio.
		//But it would also need to work a bit differently in other ways.
		void moveParameters(Distribution<double> growthRateJump, Distribution<double> competitionCoefficientJump, RandomGenerator &randomGenerator);
};

//Augmented parameters, e.g. parameters plus a normal distribution variance parameter.
//Note that this does not overwrite functions from Parameters, so getNumParameters(), getAsVector(), and the like all ignore the additional parameters.
template<size_t nAug> class AugmentedParameters : public Parameters {
	public:
		typedef std::array<double, nAug> AdditionalParametersVector;
	protected:
		AdditionalParametersVector additionalParameters;
	public:
		using Parameters::Parameters;
		AugmentedParameters(Data data, GroupingSet groupings, AdditionalParametersVector additionalParameters):
			Parameters(data, groupings), additionalParameters(additionalParameters) {};
		AugmentedParameters(AugmentedParameters p, GroupingSet groupings): //Undoes the grouping, giving a full set of parameters for every species.
			Parameters((Parameters)p, groupings), additionalParameters(p.additionalParameters) {};
		const double &getAdditionalParameter(size_t index) const {
			return additionalParameters[index];
		}
		
		const AdditionalParametersVector &getAdditionalParameters() const {
			return additionalParameters;
		}
		
		void moveParameters(Distribution<double> growthRateJump, Distribution<double> competitionCoefficientJump, std::array<Distribution<double>, nAug> additionalParameterJumps, RandomGenerator &randomGenerator);
};

#endif
