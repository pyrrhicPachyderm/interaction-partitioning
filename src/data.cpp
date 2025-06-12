#include <math.h>
#include <numeric>
#include <functional>
#include "parameters.hpp"
#include "ivp.hpp"
#include "data.hpp"

size_t Datasets::FocalResponse::findNumFocals(std::vector<size_t> focals) {
	//We want the number of unique focals.
	//We sort and remove duplicates, then take the length.
	std::sort(focals.begin(), focals.end());
	auto last = std::unique(focals.begin(), focals.end());
	return last - focals.begin();
}

bool Datasets::FocalResponse::areFocalsFirst() const {
	//Check whether the maximum numbered focal is lower than the number of focals.
	return *std::max_element(focal.begin(), focal.end()) < numRowSpecies;
}

Eigen::MatrixXdRowMajor Datasets::FocalResponse::getColGroupedDesign(const Grouping &grouping) const {
	Eigen::MatrixXdRowMajor colGroupedDesign = Eigen::MatrixXdRowMajor::Zero(design.rows(), grouping.getNumGroups());
	
	for(size_t obs = 0; obs < (size_t)design.rows(); obs++) {
		for(size_t sp = 0; sp < (size_t)design.cols(); sp++) {
			colGroupedDesign(obs, grouping.getGroup(sp)) += design(obs, sp);
		}
	}
	
	return colGroupedDesign;
}

Eigen::VectorXd Datasets::FocalResponse::getPredictions(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const {
	Eigen::VectorXd predictions = Eigen::VectorXd::Zero(response.size());
	
	Eigen::MatrixXdRowMajor colGroupedDesign = getColGroupedDesign(groupings[COL]);
	
	for(size_t obs = 0; obs < numObservations; obs++) {
		size_t focalGrowthGroup = groupings[GROWTH].getGroup(focal[obs]);
		size_t focalRowGroup = groupings[ROW].getGroup(focal[obs]);
		
		double focalDensity = isPerCapita ? 1.0 : design(obs, focal[obs]);
		double focalGrowthRate = parameters.getGrowthRate(focalGrowthGroup);
		Eigen::VectorXd densities = colGroupedDesign.row(obs);
		Eigen::VectorXd competitionCoefficients = parameters.getCompetitionCoefficientsRow(focalRowGroup);
		
		predictions[obs] = model.getResponse(focalDensity, focalGrowthRate, densities, competitionCoefficients);
	}
	
	return predictions;
}

Jacobian Datasets::FocalResponse::getPredictionsJacobian(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const {
	Jacobian jacobian = Jacobian::Zero(numObservations, parameters.getNumParameters());
	
	Eigen::MatrixXdRowMajor colGroupedDesign = getColGroupedDesign(groupings[COL]);
	
	//Note that this is the Jacobian of the residuals, not of the predicted values.
	//As such, it is negated, compared to the predicted values.
	for(size_t obs = 0; obs < numObservations; obs++) {
		size_t focalGrowthGroup = groupings[GROWTH].getGroup(focal[obs]);
		size_t focalRowGroup = groupings[ROW].getGroup(focal[obs]);
		
		double focalDensity = isPerCapita ? 1.0 : design(obs, focal[obs]);
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

double Datasets::FocalResponse::guessGrowthRate() const {
	//We might assume that all the species are in one group, and that all competition coefficients are zero.
	//This gives us the average observed response.
	return log(response.mean());
}

double Datasets::FocalResponse::guessGrowthRateMagnitude() const {
	return guessGrowthRate();
}

double Datasets::FocalResponse::guessCompetitionCoefficientMagnitude() const {
	//We will assume that with all species present at average density, growth halts.
	//This gives us 1, divided by the average density, divided by the number of species.
	return 1.0 / design.mean() / numColSpecies;
}

double Datasets::FocalResponse::guessErrorVariance() const {
	//Simply return the variance of the response variable.
	Eigen::VectorXd residuals = response - Eigen::VectorXd::Constant(response.size(), response.mean());
	return residuals.dot(residuals) / (residuals.size() - 1);
}

Datasets::TimeSeries::TimeSeries(const std::vector<size_t> &id, const Eigen::VectorXd &time, const Eigen::MatrixXdRowMajor &density) {
	assert(id.size() == (size_t)time.size());
	assert(id.size() == (size_t)density.rows());
	
	//Can't have non-focal species with time series data, so numRowSpecies = numColSpecies.
	numRowSpecies = density.cols();
	numColSpecies = density.cols();
	
	//First, we want to sort the rows, such that all entries with a given id are consecutive, and are sorted by time within ids.
	//Sorting all three (id, time, density) in place is terribly hard, so we'll just create and index list and sort that.
	std::vector<size_t> indexVec(id.size());
	std::iota(indexVec.begin(), indexVec.end(), 0); //Fill indexVec with 0, 1, 2, ...
	std::sort(indexVec.begin(), indexVec.end(), [id, time](size_t i1, size_t i2) -> bool {
		//Sort indices by id, or time on ties.
		if(id[i1] == id[i2]) return time[i1] < time[i2];
		else return id[i1] < id[i2];
	});
	
	//Populate timeSpan, includedSpecies, initialDensity and finalDensity.
	//Calculate numObservations while we're at it.
	//We're looking for consecutive pairs, so we'll inspect i and i-1, and start with i=1.
	numObservations = 0;
	for(size_t index = 1; index < indexVec.size(); index++) {
		size_t prevIndex = indexVec[index-1];
		size_t currIndex = indexVec[index];
		if(id[prevIndex] != id[currIndex]) continue; //We're looking for pairs with the same id.
		
		timeSpan.push_back(time[currIndex] - time[prevIndex]);
		includedSpecies.push_back(std::vector<size_t>());
		for(size_t i = 0; i < (size_t)density.cols(); i++) {
			if(density(prevIndex, i) > 0.0) {
				includedSpecies.back().push_back(i);
				numObservations++;
			}
		}
		
		initialDensity.push_back(Eigen::VectorXd(includedSpecies.back().size()));
		finalDensity.push_back(Eigen::VectorXd(includedSpecies.back().size()));
		for(size_t i = 0; i < includedSpecies.back().size(); i++) {
			size_t sp = includedSpecies.back()[i];
			initialDensity.back()[i] = density(prevIndex, sp);
			finalDensity.back()[i] = density(currIndex, sp);
		}
	}
	
	//Initialise the observations vector.
	observations = Eigen::VectorXd(numObservations);
	size_t obsIndex = 0;
	for(size_t i = 0; i < finalDensity.size(); i++) {
		for(size_t j = 0; j < (size_t)finalDensity[i].size(); j++) {
			observations[obsIndex++] = finalDensity[i][j];
		}
	}
}

size_t Datasets::TimeSeries::getNumSteps(double timeSpan) const {
	return (size_t)ceil(timeSpan / maxStepSize);
}

Eigen::VectorXd Datasets::TimeSeries::getGrowthRates(const Parameters &parameters, const GroupingSet &groupings, size_t experiment) const {
	const std::vector<size_t> &species = includedSpecies[experiment];
	size_t numSpecies = species.size();
	Eigen::VectorXd growthRates(numSpecies);
	
	for(size_t i = 0; i < numSpecies; i++) {
		growthRates[i] = parameters.getGrowthRate(groupings[GROWTH].getGroup(species[i]));
	}
	
	return growthRates;
}

Eigen::MatrixXdRowMajor Datasets::TimeSeries::getCompetitionCoefficients(const Parameters &parameters, const GroupingSet &groupings, size_t experiment) const {
	const std::vector<size_t> &species = includedSpecies[experiment];
	size_t numSpecies = species.size();
	Eigen::MatrixXdRowMajor competitionCoefficients(numSpecies, numSpecies);
	
	for(size_t i = 0; i < numSpecies; i++) {
		for(size_t j = 0; j < numSpecies; j++) {
			competitionCoefficients(i, j) = parameters.getCompetitionCoefficient(groupings[ROW].getGroup(species[i]), groupings[COL].getGroup(species[j]));
		}
	}
	
	return competitionCoefficients;
}

Eigen::VectorXd Datasets::TimeSeries::getPredictions(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const {
	Eigen::VectorXd predictions(observations.size());
	
	size_t resultIndex = 0;
	for(size_t i = 0; i < includedSpecies.size(); i++) {
		size_t numSpecies = includedSpecies[i].size();
		Eigen::VectorXd growthRates = getGrowthRates(parameters, groupings, i);
		Eigen::MatrixXdRowMajor competitionCoefficients = getCompetitionCoefficients(parameters, groupings, i);
		IVPDerivativeFunc derivativeFunc = std::bind(&Model::getResponses, model, std::placeholders::_2, growthRates, competitionCoefficients);
		
		predictions.segment(resultIndex, numSpecies) = solveIVP(derivativeFunc, initialDensity[i], 0.0, timeSpan[i], getNumSteps(timeSpan[i]), IVPStepFuncs::forwardEuler);
		
		resultIndex += numSpecies;
	}
	
	return predictions;
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
	//This gives us 1, divided by the average final density, divided by the average number of species per experiment.
	size_t averageNumSpecies = observations.size() / timeSpan.size();
	return 1.0 / observations.mean() / averageNumSpecies;
}

double Datasets::TimeSeries::guessErrorVariance() const {
	//Simply return the variance of the observations.
	Eigen::VectorXd residuals = observations - Eigen::VectorXd::Constant(observations.size(), observations.mean());
	return residuals.dot(residuals) / (residuals.size() - 1);
}
