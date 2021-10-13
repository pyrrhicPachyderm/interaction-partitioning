#include "parameters.hpp"

#define RELATIVE_TOLERANCE 1e-6

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

size_t Parameters::getNumParameters() const {
	return growthRates.size() + competitionCoefficients.size();
}

Parameters::Parameters(Data data, Grouping growthGrouping, Grouping rowGrouping, Grouping colGrouping) {
	//This is the standard constructor.
	//It will try to make somewhat reasonable initial guesses for all the parameter values.
	
	//It is reasonable to guess that all the competition coefficients are zero.
	double competitionCoefficient = 0.0;
	//As for the growth rates, we might assume that all the species are in one group, and that all competition coefficients are zero.
	//This gives us the average observed response.
	//If this is total, rather than per capita, we must divide by the average species density in the design.
	double growthRate = data.getResponse().mean();
	if(!data.isPerCapita) growthRate /= data.getDesign().mean();
	
	growthRates = Eigen::VectorXd::Constant(growthGrouping.getNumGroups(), growthRate);
	competitionCoefficients = Eigen::MatrixXd::Constant(rowGrouping.getNumGroups(), colGrouping.getNumGroups(), competitionCoefficient);
}

Eigen::VectorXd Parameters::getTolerances(Data data, Grouping growthGrouping, Grouping rowGrouping, Grouping colGrouping) {
	//We need somewhat reasonable guesses for the magnitudes of the growth rates and the competition coefficients.
	//We will then multiply these by the RELATIVE_TOLERANCE.
	//For the growth rates, we will use the same guess as for the initial values.
	double growthRateTolerance = data.getResponse().mean() * RELATIVE_TOLERANCE;
	if(!data.isPerCapita) growthRateTolerance /= data.getDesign().mean();
	//For the competition coefficients, we will assume that with all species present at average density, growth halts.
	//This gives us 1, divided by the square of average density, divided by the number of species.
	double competitionCoefficientTolerance = 1.0 / pow(data.getDesign().mean(), 2.0) / data.numSpecies * RELATIVE_TOLERANCE;
	
	size_t numGrowthRates = growthGrouping.getNumGroups();
	size_t numCompetitionCoefficients = rowGrouping.getNumGroups() * colGrouping.getNumGroups();
	
	return (Eigen::VectorXd(numGrowthRates + numCompetitionCoefficients) <<
		Eigen::VectorXd::Constant(numGrowthRates, growthRateTolerance),
		Eigen::VectorXd::Constant(numCompetitionCoefficients, competitionCoefficientTolerance)
	).finished();
}

Parameters::Parameters(Eigen::VectorXd parameters, Grouping growthGrouping, Grouping rowGrouping, Grouping colGrouping) {
	size_t numGrowthRates = growthGrouping.getNumGroups();
	size_t numRowGroups = rowGrouping.getNumGroups();
	size_t numColGroups = colGrouping.getNumGroups();
	
	growthRates = parameters.segment(0, numGrowthRates);
	competitionCoefficients = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(&parameters[numGrowthRates], numRowGroups, numColGroups);
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
