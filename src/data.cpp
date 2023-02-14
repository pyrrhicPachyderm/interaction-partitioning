#include <math.h>
#include <numeric>
#include <functional>
#include "parameters.hpp"
#include "ivp.hpp"
#include "data.hpp"

size_t Datasets::IndividualResponse::findNumFocals(std::vector<size_t> focals) {
	//We want the number of unique focals.
	//We sort and remove duplicates, then take the length.
	std::sort(focals.begin(), focals.end());
	auto last = std::unique(focals.begin(), focals.end());
	return last - focals.begin();
}

bool Datasets::IndividualResponse::areFocalsFirst() const {
	//Check whether the maximum numbered focal is lower than the number of focals.
	return *std::max_element(focal.begin(), focal.end()) < numRowSpecies;
}

Eigen::MatrixXdRowMajor Datasets::IndividualResponse::getColGroupedDesign(const Grouping &grouping) const {
	Eigen::MatrixXdRowMajor colGroupedDesign = Eigen::MatrixXdRowMajor::Zero(design.rows(), grouping.getNumGroups());
	
	for(size_t obs = 0; obs < (size_t)design.rows(); obs++) {
		for(size_t sp = 0; sp < (size_t)design.cols(); sp++) {
			colGroupedDesign(obs, grouping.getGroup(sp)) += design(obs, sp);
		}
	}
	
	return colGroupedDesign;
}

Eigen::VectorXd Datasets::IndividualResponse::getPredictions(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const {
	Eigen::VectorXd predictions = Eigen::VectorXd::Zero(response.size());
	
	Eigen::MatrixXdRowMajor colGroupedDesign = getColGroupedDesign(groupings[COL]);
	
	for(size_t obs = 0; obs < numObservations; obs++) {
		size_t focalGrowthGroup = groupings[GROWTH].getGroup(focal[obs]);
		size_t focalRowGroup = groupings[ROW].getGroup(focal[obs]);
		
		double focalDensity = 1.0;
		double focalGrowthRate = parameters.getGrowthRate(focalGrowthGroup);
		Eigen::VectorXd densities = colGroupedDesign.row(obs);
		Eigen::VectorXd competitionCoefficients = parameters.getCompetitionCoefficientsRow(focalRowGroup);
		
		predictions[obs] = model.getDerivative(focalDensity, focalGrowthRate, densities, competitionCoefficients);
	}
	
	return predictions;
}

Jacobian Datasets::IndividualResponse::getPredictionsJacobian(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const {
	Jacobian jacobian = Jacobian::Zero(numObservations, parameters.getNumParameters());
	
	Eigen::MatrixXdRowMajor colGroupedDesign = getColGroupedDesign(groupings[COL]);
	
	//Note that this is the Jacobian of the residuals, not of the predicted values.
	//As such, it is negated, compared to the predicted values.
	for(size_t obs = 0; obs < numObservations; obs++) {
		size_t focalGrowthGroup = groupings[GROWTH].getGroup(focal[obs]);
		size_t focalRowGroup = groupings[ROW].getGroup(focal[obs]);
		
		double focalDensity = 1.0;
		double focalGrowthRate = parameters.getGrowthRate(focalGrowthGroup);
		Eigen::VectorXd densities = colGroupedDesign.row(obs);
		Eigen::VectorXd competitionCoefficients = parameters.getCompetitionCoefficientsRow(focalRowGroup);
		
		//First, the derivatives with respect to the growth rates.
		//If it's not the growth rate of the observation's focal species, this is zero.
		//So there will be only one per observation.
		jacobian(obs, parameters.getAsVectorGrowthRateIndex(focalGrowthGroup)) = model.getGrowthRateJacobian(focalDensity, focalGrowthRate, densities, competitionCoefficients);
		
		//Second, the derivatives with respect to the competition coefficients.
		//If it's not a competition coefficient *on* the focal species, this is zero.
		//So there will be a number per row equal to the number of column groups.
		for(size_t colGroup = 0; colGroup < groupings[COL].getNumGroups(); colGroup++) {
			jacobian(obs, parameters.getAsVectorCompetitionCoefficientIndex(focalRowGroup, colGroup)) = model.getCompetitionCoefficientJacobian(focalDensity, focalGrowthRate, densities, competitionCoefficients, colGroup);
		}
	}
	
	return jacobian;
}

double Datasets::IndividualResponse::guessGrowthRate() const {
	//We might assume that all the species are in one group, and that all competition coefficients are zero.
	//This gives us the average observed response.
	return response.mean();
}

double Datasets::IndividualResponse::guessGrowthRateMagnitude() const {
	return guessGrowthRate();
}

double Datasets::IndividualResponse::guessCompetitionCoefficientMagnitude() const {
	//We will assume that with all species present at average density, growth halts.
	//This gives us 1, divided by the average density, divided by the number of species.
	return 1.0 / design.mean() / numColSpecies;
}

double Datasets::IndividualResponse::guessErrorVariance() const {
	//Simply return the variance of the response variable.
	Eigen::VectorXd residuals = response - Eigen::VectorXd::Constant(response.size(), response.mean());
	return residuals.dot(residuals) / (residuals.size() - 1);
}

size_t Datasets::TimeSeries::findNumExperiments(std::vector<size_t> id) {
	//We have one experiment per pair of consecutive time points with the same id.
	//That is, one experiment per row of the input, except for the first row with each id.
	//That is, one per row of input, minus one per unique id.
	size_t nRows = id.size();
	std::sort(id.begin(), id.end());
	size_t nIds = std::unique(id.begin(), id.end()) - id.begin();
	return nRows - nIds;
}

Datasets::TimeSeries::TimeSeries(const std::vector<size_t> &id, const Eigen::VectorXd &time, const Eigen::MatrixXdRowMajor &density):
	Datasets::Base(density.cols(), density.cols(), findNumExperiments(id) * density.cols()),
	initialDensity(numExperiments, numColSpecies),
	finalDensity(numExperiments, numColSpecies)
{
	assert(id.size() == (size_t)time.size());
	assert(id.size() == (size_t)density.rows());
	
	//First, we want to sort the rows, such that all entries with a given id are consecutive, and are sorted by time within ids.
	//Sorting all three (id, time, density) in place is terribly hard, so we'll just create and index list and sort that.
	std::vector<size_t> indexVec(id.size());
	std::iota(indexVec.begin(), indexVec.end(), 0); //Fill indexVec with 0, 1, 2, ...
	std::sort(indexVec.begin(), indexVec.end(), [id, time](size_t i1, size_t i2) -> bool {
		//Sort indices by id, or time on ties.
		if(id[i1] == id[i2]) return time[i1] < time[i2];
		else return id[i1] < id[i2];
	});
	
	//Populate the timeSpan, initialDensity and finalDensity vectors.
	//We're looking for consecutive pairs, so we'll inspect i and i-1, and start with i=1.
	size_t resultI = 0;
	for(size_t i = 1; i < indexVec.size(); i++) {
		size_t prevI = indexVec[i-1];
		size_t currI = indexVec[i];
		if(id[prevI] != id[currI]) continue; //We're looking for pairs with the same id.
		
		timeSpan.push_back(time[currI] - time[prevI]);
		initialDensity.row(resultI) = density.row(prevI);
		finalDensity.row(resultI) = density.row(currI);
		resultI++;
	}
}

size_t Datasets::TimeSeries::getNumSteps(double timeSpan) const {
	return (size_t)ceil(timeSpan / maxStepSize);
}

Eigen::VectorXd Datasets::TimeSeries::getPredictions(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const {
	Eigen::MatrixXdRowMajor predictions(initialDensity.rows(), initialDensity.cols());
	
	Parameters ungroupedParameters = Parameters(parameters, groupings);
	IVPDerivativeFunc derivativeFunc = std::bind(&Model::getDerivatives, model, std::placeholders::_2, ungroupedParameters.getGrowthRates(), ungroupedParameters.getCompetitionCoefficients());
	
	for(size_t i = 0; i < numExperiments; i++) {
		predictions.row(i) = solveIVP(derivativeFunc, initialDensity.row(i), 0.0, timeSpan[i], getNumSteps(timeSpan[i]), IVPStepFuncs::forwardEuler);
	}
	
	return predictions.reshaped<Eigen::RowMajor>();
}

double Datasets::TimeSeries::guessGrowthRate() const {
	//A positive growth rate (with initial gueses of 0 competition) leads to runaway exponential growth.
	//This gives enormous initial residuals, and hence a likelihood that rounds to 0.
	//This means likelihood ratios are always NaN.
	//To curb this, choose 0 as the initial growth rate.
	return 0.0;
}

double Datasets::TimeSeries::guessGrowthRateMagnitude() const {
	//Being time series data, the response is actual population growth, so 1 ought to be a decent guess for the growth rate magnitude.
	return 1.0;
}

double Datasets::TimeSeries::guessCompetitionCoefficientMagnitude() const {
	//We will assume that with all species present at average final density, growth halts.
	//This gives us 1, divided by the average density, divided by the number of species.
	return 1.0 / finalDensity.mean() / numColSpecies;
}

double Datasets::TimeSeries::guessErrorVariance() const {
	//Simply return the variance of the final density.
	Eigen::VectorXd residuals = finalDensity.reshaped<Eigen::RowMajor>() - Eigen::VectorXd::Constant(finalDensity.size(), finalDensity.mean());
	return residuals.dot(residuals) / (residuals.size() - 1);
}
