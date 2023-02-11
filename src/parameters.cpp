#include "parameters.hpp"

#define RELATIVE_TOLERANCE 1e-6

double Parameters::getNumRowSpecies() const {
	return numRowSpecies;
}

double Parameters::getNumColSpecies() const {
	return numColSpecies;
}

double Parameters::getGrowthRate(size_t index) const {
	return growthRates[index];
}

const Eigen::VectorXd &Parameters::getGrowthRates() const {
	return growthRates;
}

double Parameters::getCompetitionCoefficient(size_t rowIndex, size_t colIndex) const {
	return competitionCoefficients(rowIndex, colIndex);
}

//Can't return a reference here, as MatrixXd::row() makes a temporary construct.
//TODO: Figure out how to do this.
Eigen::VectorXd Parameters::getCompetitionCoefficientsRow(size_t rowIndex) const {
	return competitionCoefficients.row(rowIndex);
}

const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &Parameters::getCompetitionCoefficients() const {
	return competitionCoefficients;
}

size_t Parameters::getNumParameters() const {
	return growthRates.size() + competitionCoefficients.size();
}

Parameters::Parameters(Data data, GroupingSet groupings):
	numRowSpecies(data.getNumRowSpecies()), numColSpecies(data.getNumColSpecies())
{
	//This is the standard constructor.
	//It will try to make somewhat reasonable initial guesses for all the parameter values.
	double growthRate = data.guessGrowthRate();
	double competitionCoefficient = data.guessCompetitionCoefficient();
	
	growthRates = Eigen::VectorXd::Constant(groupings[GROWTH].getNumGroups(), growthRate);
	competitionCoefficients = Eigen::MatrixXd::Constant(groupings[ROW].getNumGroups(), groupings[COL].getNumGroups(), competitionCoefficient);
}

Eigen::VectorXd Parameters::getTolerances(Data data, GroupingSet groupings) {
	//We need somewhat reasonable guesses for the magnitudes of the growth rates and the competition coefficients.
	//We will then multiply these by the RELATIVE_TOLERANCE.
	double growthRateTolerance = data.guessGrowthRate() * RELATIVE_TOLERANCE;
	double competitionCoefficientTolerance = data.guessCompetitionCoefficientMagnitude() * RELATIVE_TOLERANCE;
	
	size_t numGrowthRates = groupings[GROWTH].getNumGroups();
	size_t numCompetitionCoefficients = groupings[ROW].getNumGroups() * groupings[COL].getNumGroups();
	
	return (Eigen::VectorXd(numGrowthRates + numCompetitionCoefficients) <<
		Eigen::VectorXd::Constant(numGrowthRates, growthRateTolerance),
		Eigen::VectorXd::Constant(numCompetitionCoefficients, competitionCoefficientTolerance)
	).finished();
}

Parameters::Parameters(Eigen::VectorXd parameters, GroupingSet groupings):
	numRowSpecies(groupings[ROW].numSpecies), numColSpecies(groupings[COL].numSpecies)
{
	size_t numGrowthRates = groupings[GROWTH].getNumGroups();
	size_t numRowGroups = groupings[ROW].getNumGroups();
	size_t numColGroups = groupings[COL].getNumGroups();
	
	growthRates = parameters.segment(0, numGrowthRates);
	competitionCoefficients = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(&parameters[numGrowthRates], numRowGroups, numColGroups);
}

Parameters::Parameters(Parameters p, GroupingSet groupings):
	numRowSpecies(groupings[ROW].numSpecies), numColSpecies(groupings[COL].numSpecies)
{
	//This gives ungrouped parameters (one parameter per species/pair of species), based on grouped parameters and their groups.
	growthRates = Eigen::VectorXd(numRowSpecies);
	competitionCoefficients = Eigen::MatrixXd(numRowSpecies, numColSpecies);
	for(size_t i = 0; i < numRowSpecies; i++) {
		growthRates[i] = p.growthRates[groupings[GROWTH].getGroup(i)];
	}
	for(size_t i = 0; i < numRowSpecies; i++) {
		for(size_t j = 0; j < numColSpecies; j++) {
			competitionCoefficients(i,j) = p.competitionCoefficients(groupings[ROW].getGroup(i), groupings[COL].getGroup(j));
		}
	}
}

Eigen::VectorXd Parameters::getAsVector() const {
	return (Eigen::VectorXd(growthRates.size() + competitionCoefficients.size()) <<
		growthRates,
		Eigen::Map<const Eigen::VectorXd>(competitionCoefficients.data(), competitionCoefficients.size())
	).finished();
}

size_t Parameters::getAsVectorGrowthRateIndex(size_t index) const {
	return index;
}

size_t Parameters::getAsVectorCompetitionCoefficientIndex(size_t rowIndex, size_t colIndex) const {
	return growthRates.size() + competitionCoefficients.cols() * rowIndex + colIndex;
}

void Parameters::moveParameters(Distribution<double> growthRateJump, Distribution<double> competitionCoefficientJump) {
	for(size_t i = 0; i < (size_t)growthRates.size(); i++) {
		growthRates[i] += growthRateJump.getRandom();
	}
	for(size_t i = 0; i < (size_t)competitionCoefficients.rows(); i++) {
		for(size_t j = 0; j < (size_t)competitionCoefficients.cols(); j++) {
			competitionCoefficients(i, j) += competitionCoefficientJump.getRandom();
		}
	}
}

template<size_t nAug> void AugmentedParameters<nAug>::moveParameters(Distribution<double> growthRateJump, Distribution<double> competitionCoefficientJump, std::array<Distribution<double>, nAug> additionalParameterJumps) {
	Parameters::moveParameters(growthRateJump, competitionCoefficientJump); //Call the base class function.
	for(size_t i = 0; i < nAug; i++) {
		additionalParameters[i] += additionalParameterJumps[i].getRandom();
	}
}

//Explicitly instantiate.
template class AugmentedParameters<1>;
