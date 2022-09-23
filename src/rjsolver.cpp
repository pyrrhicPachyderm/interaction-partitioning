#include <numeric> //Gives std::accumulate.
#include <random>
#include "rjsolver.hpp"

#define RANDOM_SEED 42
#define MAX_TRANS_MODEL_JUMP_PROBABILITY 0.9
#define DESIRED_ACCEPTANCE_RATE 0.23 //For within-model jumps.
#define MAX_JUMP_VARIANCE_MULTIPLIER_CHANGE 2.0 //Controls how quickly dialIn changes the jump variance multiplier.

static std::default_random_engine randomNumberGenerator(RANDOM_SEED);

double ReversibleJumpSolver::getTransModelJumpProbability(GroupingIndexSet sourceGroupingIndices, GroupingIndexSet destGroupingIndices) const {
	return getTransModelJumpProbability(sourceGroupingIndices, destGroupingIndices, transModelJumpProbabilityMultiplier);
}

double ReversibleJumpSolver::getTransModelJumpProbability(GroupingIndexSet sourceGroupingIndices, GroupingIndexSet destGroupingIndices, double multiplier) const {
	double sourceHyperprior = hyperprior.getDensity(getGroupings(sourceGroupingIndices));
	double destHyperprior = hyperprior.getDensity(getGroupings(destGroupingIndices));
	return multiplier * std::min(1.0, destHyperprior/sourceHyperprior);
}

std::vector<double> ReversibleJumpSolver::getTransModelJumpProbabilities(GroupingType groupingType, MoveType moveType) const {
	return getTransModelJumpProbabilities(groupingType, moveType, currentGroupings, transModelJumpProbabilityMultiplier);
}

std::vector<double> ReversibleJumpSolver::getTransModelJumpProbabilities(GroupingType groupingType, MoveType moveType, GroupingIndexSet groupingIndices, double multiplier) const {
	const std::vector<size_t> &destIndices = groupingLattice.getMoveDests(moveType, groupingIndices[groupingType]);
	std::vector<double> probabilities;
	
	GroupingIndexSet newGroupingIndices = groupingIndices;
	for(size_t i = 0; i < destIndices.size(); i++) {
		newGroupingIndices[groupingType] = destIndices[i];
		probabilities.push_back(getTransModelJumpProbability(groupingIndices, newGroupingIndices, multiplier));
	}
	
	return probabilities;
}

double ReversibleJumpSolver::getUnscaledTotalTransModelJumpProbability(GroupingIndexSet groupingIndices) const {
	double result = 0.0;
	for(size_t groupingType = 0; groupingType < NUM_GROUPING_TYPES; groupingType++) {
		if(!isChangingGroupings[groupingType]) continue;
		for(size_t moveType = 0; moveType < NUM_MOVE_TYPES; moveType++) {
			std::vector<double> probabilities = getTransModelJumpProbabilities((GroupingType)groupingType, (MoveType)moveType, groupingIndices, 1.0);
			result += std::accumulate(probabilities.begin(), probabilities.end(), 0.0);
		}
	}
	return result;
}

double ReversibleJumpSolver::getUnscaledMaxTransModelJumpProbability(size_t recursionLevel, GroupingIndexSet groupingIndices) const {
	//We need to loop over each type of model that's changing; this is a NUM_GROUPING_TYPES-times nested loop.
	//However, whether each loop exists is conditional on isChangingGroupings.
	//So the neatest way to do this is recursion.
	//recursionLevel is an index into a GroupingSet; the index to loop over in this recursion.
	//groupingIndices is only determined for indices below recursionLevel.
	
	if(recursionLevel == NUM_GROUPING_TYPES) {
		return getUnscaledTotalTransModelJumpProbability(groupingIndices);
	}
	
	if(!isChangingGroupings[recursionLevel]) {
		//We're not looping at this level; just set the actual grouping we're using at this level.
		groupingIndices[recursionLevel] = currentGroupings[recursionLevel];
		return getUnscaledMaxTransModelJumpProbability(recursionLevel+1, groupingIndices);
	}
	
	double result = 0.0;
	for(size_t i = 0; i < groupingLattice.getNumGroupings(); i++) {
		groupingIndices[recursionLevel] = i;
		result = std::max(result, getUnscaledMaxTransModelJumpProbability(recursionLevel+1, groupingIndices));
	}
	return result;
}

double ReversibleJumpSolver::getTransModelJumpProbabilityMultiplier() const {
	//This is the value of c used in trans-model jump probabilities, as given defined in Green 1995 and the report.
	
	//The GroupingIndexSet to pass the recursive function doesn't matter.
	//It's only determined for indices below recursionLevel, which is zero.
	double unscaledMaxTransModelJumpProbability = getUnscaledMaxTransModelJumpProbability(0, GroupingIndexSet());
	
	return MAX_TRANS_MODEL_JUMP_PROBABILITY / unscaledMaxTransModelJumpProbability;
}

GroupingSet ReversibleJumpSolver::getGroupings(GroupingIndexSet groupingIndices) const {
	//TODO: This should be done based on NUM_GROUPING_TYPES, but constructing a std::array is hard.
	//Perhaps look at https://stackoverflow.com/a/32175958
	return GroupingSet({groupingLattice.getGrouping(groupingIndices[GROWTH]), groupingLattice.getGrouping(groupingIndices[ROW]), groupingLattice.getGrouping(groupingIndices[COL])});
}

ReversibleJumpSolver::GroupingIndexSet ReversibleJumpSolver::getGroupingIndices(GroupingSet groupings) const {
	//TODO: Should also be done based on NUM_GROUPING_TYPES.
	return GroupingIndexSet({groupingLattice.getIndex(groupings[GROWTH]), groupingLattice.getIndex(groupings[ROW]), groupingLattice.getIndex(groupings[COL])});
}

double ReversibleJumpSolver::getGrowthRateJumpVariance() const {
	return growthRateApproximatePosteriorVariance * jumpVarianceMultiplier;
}

double ReversibleJumpSolver::getCompetitionCoefficientJumpVariance() const {
	return competitionCoefficientApproximatePosteriorVariance * jumpVarianceMultiplier;
}

ReversibleJumpSolver::AdditionalParametersVector ReversibleJumpSolver::getAdditionalParametersJumpVariance() const {
	return additionalParametersApproximatePosteriorVariance * jumpVarianceMultiplier;
}

double ReversibleJumpSolver::getTransModelJumpVariance(GroupingType groupingType) const {
	if(groupingType == GROWTH) {
		return getGrowthRateJumpVariance();
	} else if(groupingType == ROW || groupingType == COL) {
		return getCompetitionCoefficientJumpVariance();
	} else __builtin_unreachable();
}

static double getVariance(std::vector<double> nums) {
	Eigen::Map<Eigen::VectorXd> vec(&nums[0], nums.size());
	Eigen::VectorXd residuals = vec - Eigen::VectorXd::Constant(vec.size(), vec.mean());
	return residuals.dot(residuals) / residuals.size();
}

static double getRandomProbability() {
	//A random double in [0,1).
	return std::uniform_real_distribution(0.0, 1.0)(randomNumberGenerator);
}

static double getRandomNormal(double variance) {
	return std::normal_distribution(0.0, sqrt(variance))(randomNumberGenerator);
}

static double getNormalDensity(double variance, double residual) {
	return 1.0 / sqrt(2 * M_PI * variance) * exp(-0.5 * residual*residual / variance);
}

//We need to use this, rather than direct setting, to ensure it also dirties the groupings.
//Further, this requires that this function is called after any changes to current or proposed groupings are performed.
void ReversibleJumpSolver::setIsProposing(bool b) {
	if(isProposing == b) return;
	isProposing = b;
	for(size_t i = 0; i < NUM_GROUPING_TYPES; i++) {
		if(currentGroupings[i] != proposedGroupings[i]) dirtyGrouping((GroupingType)i);
	}
}

double ReversibleJumpSolver::proposeTransModelJump(GroupingType groupingType, MoveType moveType, size_t adjIndex) {
	//adjIndex is the index into the adjacency list of the relevant part of groupingLattice.
	size_t newGroupingIndex = groupingLattice.getMoveDest(moveType, currentGroupings[groupingType], adjIndex);
	GroupingMove groupingMove = groupingLattice.getMove(moveType, currentGroupings[groupingType], adjIndex);
	
	proposedGroupings = currentGroupings;
	proposedGroupings[groupingType] = newGroupingIndex;
	
	double jumpVariance = getTransModelJumpVariance(groupingType);
	RandomVariableFunc getRandomVariable = std::bind(getRandomNormal, jumpVariance);
	RandomVariableDensityFunc getRandomVariableDensity = std::bind(getNormalDensity, jumpVariance, std::placeholders::_1);
	
	proposedParameters = currentParameters;
	double acceptanceRatio = proposedParameters.moveModel(groupingType, moveType, groupingMove, getRandomVariable, getRandomVariableDensity);
	
	acceptanceRatio *= getTransModelJumpProbability(currentGroupings, proposedGroupings) / getTransModelJumpProbability(proposedGroupings, currentGroupings);
	
	if(moveType == MERGE) proposedJumpType = MERGE_JUMP;
	else if(moveType == SPLIT) proposedJumpType = SPLIT_JUMP;
	else __builtin_unreachable();
	setIsProposing(true);
	
	return acceptanceRatio;
}

double ReversibleJumpSolver::proposeWithinModelJump() {
	proposedGroupings = currentGroupings;
	proposedParameters = currentParameters;
	
	RandomVariableFunc getGrowthRateJump = std::bind(getRandomNormal, getGrowthRateJumpVariance());
	RandomVariableFunc getCompetitionCoefficientJump = std::bind(getRandomNormal, getCompetitionCoefficientJumpVariance());
	std::array<RandomVariableFunc, NUM_ADDITIONAL_PARAMETERS> getAdditionalParameterJumps;
	for(size_t i = 0; i < NUM_ADDITIONAL_PARAMETERS; i++) {
		getAdditionalParameterJumps[i] = std::bind(getRandomNormal, getAdditionalParametersJumpVariance()[i]);
	}
	
	proposedParameters.moveParameters(getGrowthRateJump, getCompetitionCoefficientJump, getAdditionalParameterJumps);
	
	proposedJumpType = WITHIN_JUMP;
	setIsProposing(true);
	
	//TODO: If using a non-symmetric jumping density, the jumping density component of the acceptance ratio may not be 1.
	return 1.0;
}

double ReversibleJumpSolver::proposeJump() {
	double selector = getRandomProbability();
	
	//Iterate over each type of trans-model jump, and see if we do it.
	//Instead of accumulating some sum of previous probabilities, we'll just decrement selector.
	for(size_t groupingType = 0; groupingType < NUM_GROUPING_TYPES; groupingType++) {
		if(!isChangingGroupings[groupingType]) continue;
		for(size_t moveType = 0; moveType < NUM_MOVE_TYPES; moveType++) {
			std::vector<double> probabilities = getTransModelJumpProbabilities((GroupingType)groupingType, (MoveType)moveType);
			for(size_t adjIndex = 0; adjIndex < probabilities.size(); adjIndex++) {
				selector -= probabilities[adjIndex];
				if(selector < 0) {
					return proposeTransModelJump((GroupingType)groupingType, (MoveType)moveType, adjIndex);
				}
			}
		}
	}
	
	return proposeWithinModelJump();
}

void ReversibleJumpSolver::acceptJump() {
	currentGroupings = proposedGroupings;
	currentParameters = proposedParameters;
	setIsProposing(false);
}

void ReversibleJumpSolver::rejectJump() {
	setIsProposing(false);
}

double ReversibleJumpSolver::getErrorVariance() const {
	//This is the additional parameter, and is used in calculating the likelihood.
	return getParameters().getAdditionalParameter(0);
}

Eigen::VectorXd ReversibleJumpSolver::getResiduals() {
	return Solver::getResiduals(isProposing ? (Parameters)proposedParameters : (Parameters)currentParameters);
}

double ReversibleJumpSolver::getLikelihoodRatio(Eigen::VectorXd sourceResiduals, Eigen::VectorXd destResiduals, double sourceErrorVariance, double destErrorVariance) {
	//Having a getLikelihood() function and calling it for source and destination
	//just results in it returning 0 twice, as it multiplies several hundred small numbers together, and gets too small.
	//So we must calculate the likelihood ratio, getting the ratio as we go.
	double likelihoodRatio = 1.0;
	for(size_t i = 0; i < (size_t)sourceResiduals.size(); i++) {
		double sourceLikelihood = getNormalDensity(sourceErrorVariance, sourceResiduals[i]);
		double destLikelihood = getNormalDensity(destErrorVariance, destResiduals[i]);
		likelihoodRatio *= destLikelihood / sourceLikelihood;
	}
	return likelihoodRatio;
}

double ReversibleJumpSolver::getPriorDensity() const {
	double hyperpriorDensity = hyperprior.getDensity(getGroupings());
	
	//TODO: There should also be the prior of all the parameters here.
	//Currently, the competition coefficients have a prior.
	//But I'm trying some sort of improper flat prior for the growth rates; this is only really valid as long as their are no trans-model jumps changing the number of growth rates.
	//Lastly, some prior on the variance might be suitable.
	double parameterPriorDensity = 1.0;
	
	Eigen::MatrixXd competitionCoefficients = getParameters().getCompetitionCoefficients();
	for(size_t i = 0; i < (size_t)competitionCoefficients.rows(); i++) {
		for(size_t j = 0; j < (size_t)competitionCoefficients.cols(); j++) {
			parameterPriorDensity *= getNormalDensity(competitionCoefficientPriorVariance, competitionCoefficients(i,j));
		}
	}
	
	return hyperpriorDensity * parameterPriorDensity;
}

bool ReversibleJumpSolver::makeJump(bool canTransModelJump) {
	double sourcePrior = getPriorDensity();
	Eigen::VectorXd sourceResiduals = getResiduals();
	double sourceErrorVariance = getErrorVariance();
	
	double jumpingDensityRatio = canTransModelJump ? proposeJump() : proposeWithinModelJump();
	
	double destPrior = getPriorDensity();
	Eigen::VectorXd destResiduals = getResiduals();
	double destErrorVariance = getErrorVariance();
	
	double priorRatio = destPrior / sourcePrior;
	double likelihoodRatio = getLikelihoodRatio(sourceResiduals, destResiduals, sourceErrorVariance, destErrorVariance);
	double acceptanceRatio = priorRatio * likelihoodRatio * jumpingDensityRatio;
	
	
	double selector = getRandomProbability();
	if(selector < acceptanceRatio) {
		acceptJump();
		return true;
	} else {
		rejectJump();
		return false;
	}
}

void ReversibleJumpSolver::burnIn(size_t numJumps, bool canTransModelJump) {
	for(size_t i = 0; i < numJumps; i++) {
		makeJump(canTransModelJump);
	}
}

//A pair of functions to raise or lower a jump variance multiplier.
//They enact a maximum proportional change of MAX_JUMP_VARIANCE_MULTIPLIER_CHANGE.
//They enact no change with discrepantProportion = 0, up to full change with discrepantProportion = 1.
//These edit the jump variance in-place.
static void raiseJumpVarianceMultiplier(double &multiplier, double discrepantProportion) {
	multiplier *= 1.0 + discrepantProportion * (MAX_JUMP_VARIANCE_MULTIPLIER_CHANGE - 1.0);
}
static void lowerJumpVarianceMultiplier(double &multiplier, double discrepantProportion) {
	multiplier /= 1.0 + discrepantProportion * (MAX_JUMP_VARIANCE_MULTIPLIER_CHANGE - 1.0);
}


void ReversibleJumpSolver::dialIn(size_t jumpsPerDial, size_t numDials) {
	//Dials in the jumping variances.
	//Combining two different dialling in methods here.
	//The first to balance the sizes of jumps in different variables against one another.
	//The second to balance overall sizes of jumps.
	//The first tries to estimate the variance of each posterior.
	//The second aims for a given within-model acceptance ratio, and tries to equalise merge and split acceptance ratios.
	
	for(size_t i = 0; i < numDials; i++) {
		//For the first method.
		std::vector<double> growthRates;
		std::vector<double> competitionCoefficients;
		std::array<std::vector<double>, NUM_ADDITIONAL_PARAMETERS> additionalParameters;
		
		//For the second method.
		Eigen::Array<double, NUM_JUMP_TYPES, 1> numProposals = Eigen::Array<double, NUM_JUMP_TYPES, 1>::Zero();
		Eigen::Array<double, NUM_JUMP_TYPES, 1> numAccepts = Eigen::Array<double, NUM_JUMP_TYPES, 1>::Zero();
		
		for(size_t j = 0; j < jumpsPerDial; j++) {
			bool accepted = makeJump();
			numProposals[proposedJumpType] += 1.0;
			if(accepted) numAccepts[proposedJumpType] += 1.0;
			
			//We can't rely on there being a certain number of growth rates or competition coefficients.
			//And it's altogether too much work to tally multiple sets of growth rates or competition coefficients, and take the variance of each set.
			//So we just take the first of each.
			growthRates.push_back(getParameters().getGrowthRate(0));
			competitionCoefficients.push_back(getParameters().getCompetitionCoefficient(0,0));
			for(size_t i = 0; i < NUM_ADDITIONAL_PARAMETERS; i++) {
				additionalParameters[i].push_back(getParameters().getAdditionalParameter(i));
			}
		}
		growthRateApproximatePosteriorVariance = getVariance(growthRates);
		competitionCoefficientApproximatePosteriorVariance = getVariance(competitionCoefficients);
		for(size_t i = 0; i < NUM_ADDITIONAL_PARAMETERS; i++) {
			additionalParametersApproximatePosteriorVariance[i] = getVariance(additionalParameters[i]);
		}
		
		Eigen::Array<double, NUM_JUMP_TYPES, 1> acceptanceRates = numAccepts / numProposals;
		
		//If the acceptance rate is too high, we want to raise the withinModelJumpVarianceMultiplier, and vice versa.
		if(acceptanceRates[WITHIN_JUMP] > DESIRED_ACCEPTANCE_RATE) {
			double discrepantProportion = (acceptanceRates[WITHIN_JUMP] - DESIRED_ACCEPTANCE_RATE) / (1.0 - DESIRED_ACCEPTANCE_RATE);
			raiseJumpVarianceMultiplier(jumpVarianceMultiplier, discrepantProportion);
		} else {
			double discrepantProportion = (DESIRED_ACCEPTANCE_RATE - acceptanceRates[WITHIN_JUMP]) / DESIRED_ACCEPTANCE_RATE;
			lowerJumpVarianceMultiplier(jumpVarianceMultiplier, discrepantProportion);
		}
	}
}
