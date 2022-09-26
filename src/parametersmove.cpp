//This is a part of the parameters module, but is specifically the part that deals with a Reversible Jump MCMC trans-model move.
#include "parameters.hpp"

static Eigen::VectorXd mergeParameterVectors(const std::pair<Eigen::VectorXd, Eigen::VectorXd> &splitParameters, const GroupingMove &groupingMove) {
	size_t size1 = groupingMove.getSplitGroupSizes().first;
	size_t size2 = groupingMove.getSplitGroupSizes().second;
	return (splitParameters.first * size1 + splitParameters.second * size2) / (size1 + size2);
}

static std::pair<Eigen::VectorXd, Eigen::VectorXd> splitParameterVectors(const Eigen::VectorXd &mergedParameters, const GroupingMove &groupingMove, Distribution<double> randomVariableDistribution) {
	size_t size1 = groupingMove.getSplitGroupSizes().first;
	size_t size2 = groupingMove.getSplitGroupSizes().second;
	double scaling1 = (double)size2 / (size1 + size2);
	double scaling2 = (double)size1 / (size1 + size2);
	//Note that scaling1 uses size2.
	//This is such that parameter 1 moves more when group 1 is smaller (and group 2 is conversely larger), to preserve the weighted mean.
	
	Eigen::VectorXd randomVariable(mergedParameters.size());
	for(size_t i = 0; i < (size_t)randomVariable.size(); i++) {
		randomVariable[i] = randomVariableDistribution.getRandom();
	}
	
	Eigen::VectorXd splitParameters1 = mergedParameters + randomVariable * scaling1; //Raise the first parameter.
	Eigen::VectorXd splitParameters2 = mergedParameters - randomVariable * scaling2; //Lower the second.
	
	return std::make_pair(splitParameters1, splitParameters2);
}

static Eigen::VectorXd reverseEngineerRandomVariable(const Eigen::VectorXd &mergedParameters, const std::pair<Eigen::VectorXd, Eigen::VectorXd> &splitParameters) {
	return splitParameters.second - splitParameters.first;
}

static double getRandomVariableJumpingDensity(MoveType moveType, const Eigen::VectorXd &mergedParameters, const std::pair<Eigen::VectorXd, Eigen::VectorXd> &splitParameters, Distribution<double> randomVariableDistribution) {
	Eigen::VectorXd randomVariable = reverseEngineerRandomVariable(mergedParameters, splitParameters);
	
	double jumpingDensity = 1.0;
	for(size_t i = 0; i < (size_t)randomVariable.size(); i++) {
		jumpingDensity *= randomVariableDistribution.getDensity(randomVariable[i]);
	}
	
	//If we are performing a split, this is on the denominator of the ratio.
	if(moveType == SPLIT) jumpingDensity = 1 / jumpingDensity;
	
	return jumpingDensity;
}

static double getJacobianBlockDeterminant(MoveType moveType, double mergedParameter, double randomVariable, const std::pair<double, double> &splitParameters, const GroupingMove &groupingMove) {
	//See report for derivation.
	return 1.0;
}

static double getJacobianDeterminant(MoveType moveType, const Eigen::VectorXd &mergedParameters, const std::pair<Eigen::VectorXd, Eigen::VectorXd> &splitParameters, const GroupingMove &groupingMove) {
	Eigen::VectorXd randomVariable = reverseEngineerRandomVariable(mergedParameters, splitParameters);
	
	double determinant = 1.0;
	for(size_t i = 0; i < (size_t)randomVariable.size(); i++) {
		determinant *= getJacobianBlockDeterminant(moveType, mergedParameters[i], randomVariable[i], std::make_pair(splitParameters.first[i], splitParameters.second[i]), groupingMove);
	}
	return determinant;
}

double Parameters::moveModel(GroupingType groupingType, MoveType moveType, const GroupingMove &groupingMove, Distribution<double> randomVariableDistribution) {
	size_t mergedGroup = groupingMove.getMergedGroup();
	std::pair<size_t,size_t> splitGroups = groupingMove.getSplitGroups();
	
	//One of these will be taken from the input, and converted to the other, so both shall end up defined.
	Eigen::VectorXd mergedParameters;
	std::pair<Eigen::VectorXd, Eigen::VectorXd> splitParameters;
	
	//Only one of these will ever be intialised, and then written back.
	Eigen::VectorXd newGrowthRates;
	Eigen::MatrixXd newCompetitionCoefficients;
	
	//A number to use as the upper bounds of iteration later: the number of groups in the altered grouping pre-move.
	size_t oldGroupingSize;
	
	//Initialise newGrowthRates or newCompetitionCoefficients, as appropriate, to the correct size.
	//Also extract mergedParameters or splitParameters, as appropriate.
	//Also define oldGroupingSize.
	switch(groupingType) {
		case GROWTH:
			oldGroupingSize = growthRates.size();
			if(moveType == MERGE) {
				newGrowthRates = Eigen::VectorXd(growthRates.size() - 1);
				splitParameters = std::make_pair(growthRates.row(splitGroups.first), growthRates.row(splitGroups.second));
			} else if(moveType == SPLIT) {
				newGrowthRates = Eigen::VectorXd(growthRates.size() + 1);
				mergedParameters = growthRates.row(mergedGroup);
			} else __builtin_unreachable();
			break;
		case ROW:
			oldGroupingSize = competitionCoefficients.rows();
			if(moveType == MERGE) {
				newCompetitionCoefficients = Eigen::MatrixXd(competitionCoefficients.rows() - 1, competitionCoefficients.cols());
				splitParameters = std::make_pair(competitionCoefficients.row(splitGroups.first), competitionCoefficients.row(splitGroups.second));
			} else if(moveType == SPLIT) {
				newCompetitionCoefficients = Eigen::MatrixXd(competitionCoefficients.rows() + 1, competitionCoefficients.cols());
				mergedParameters = competitionCoefficients.row(mergedGroup);
			} else __builtin_unreachable();
			break;
		case COL:
			oldGroupingSize = competitionCoefficients.cols();
			if(moveType == MERGE) {
				newCompetitionCoefficients = Eigen::MatrixXd(competitionCoefficients.rows(), competitionCoefficients.cols() - 1);
				splitParameters = std::make_pair(competitionCoefficients.col(splitGroups.first), competitionCoefficients.col(splitGroups.second));
			} else if(moveType == SPLIT) {
				newCompetitionCoefficients = Eigen::MatrixXd(competitionCoefficients.rows(), competitionCoefficients.cols() + 1);
				mergedParameters = competitionCoefficients.col(mergedGroup);
			} else __builtin_unreachable();
			break;
		default:
			__builtin_unreachable();
	}
	
	//Convert mergedParameters to splitParameters, or vice versa, as appropriate.
	if(moveType == MERGE) {
		mergedParameters = mergeParameterVectors(splitParameters, groupingMove);
	} else if(moveType == SPLIT) {
		splitParameters = splitParameterVectors(mergedParameters, groupingMove, randomVariableDistribution);
	} else __builtin_unreachable();
	
	//Populate newGrowthRates or newCompetitionCoefficients, as appropriate.
	for(size_t i = 0; i < oldGroupingSize; i++) {
		switch(groupingType) {
			case GROWTH:
				if(moveType == SPLIT && i == mergedGroup) {
					newGrowthRates.row(splitGroups.first) = splitParameters.first;
					newGrowthRates.row(splitGroups.second) = splitParameters.second;
					break;
				}
				if(moveType == MERGE && i == splitGroups.first) {
					newGrowthRates.row(mergedGroup) = mergedParameters;
					break;
				}
				if(moveType == MERGE && i == splitGroups.second) break;
				newGrowthRates.row(groupingMove.getMapping(moveType, i)) = growthRates.row(i);
				break;
			case ROW:
				if(moveType == SPLIT && i == mergedGroup) {
					newCompetitionCoefficients.row(splitGroups.first) = splitParameters.first;
					newCompetitionCoefficients.row(splitGroups.second) = splitParameters.second;
					break;
				}
				if(moveType == MERGE && i == splitGroups.first) {
					newCompetitionCoefficients.row(mergedGroup) = mergedParameters;
					break;
				}
				if(moveType == MERGE && i == splitGroups.second) break;
				newCompetitionCoefficients.row(groupingMove.getMapping(moveType, i)) = competitionCoefficients.row(i);
				break;
			case COL:
				if(moveType == SPLIT && i == mergedGroup) {
					newCompetitionCoefficients.col(splitGroups.first) = splitParameters.first;
					newCompetitionCoefficients.col(splitGroups.second) = splitParameters.second;
					break;
				}
				if(moveType == MERGE && i == splitGroups.first) {
					newCompetitionCoefficients.col(mergedGroup) = mergedParameters;
					break;
				}
				if(moveType == MERGE && i == splitGroups.second) break;
				newCompetitionCoefficients.col(groupingMove.getMapping(moveType, i)) = competitionCoefficients.col(i);
				break;
			default:
				__builtin_unreachable();
		}
	}
	
	if(groupingType == GROWTH) growthRates = newGrowthRates;
	else if(groupingType == ROW || groupingType == COL) competitionCoefficients = newCompetitionCoefficients;
	else __builtin_unreachable();
	
	//Calculate and return the acceptance ratio.
	double acceptanceRatio =
		getRandomVariableJumpingDensity(moveType, mergedParameters, splitParameters, randomVariableDistribution) *
		getJacobianDeterminant(moveType, mergedParameters, splitParameters, groupingMove);
	return acceptanceRatio;
}
