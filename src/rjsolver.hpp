#ifndef REVERSIBLE_JUMP_SOLVER_HPP
#define REVERSIBLE_JUMP_SOLVER_HPP

#include "groupingmove.hpp"
#include "solver.hpp"
#include "priors.hpp"

#define DEFAULT_RANDOM_SEED 42

class ReversibleJumpSolver : public Solver {
	public:
		const static size_t NUM_ADDITIONAL_PARAMETERS = 1;
		typedef AugmentedParameters<NUM_ADDITIONAL_PARAMETERS>::AdditionalParametersVector AdditionalParametersVector;
		enum JumpType {MERGE_JUMP, SPLIT_JUMP, WITHIN_JUMP, NUM_JUMP_TYPES};
	protected:
		typedef std::array<bool, NUM_GROUPING_TYPES> GroupingBooleanSet;
		
		RandomGenerator randomGenerator = RandomGenerator(DEFAULT_RANDOM_SEED);
		
		Hyperprior hyperprior;
		AugmentedParametersPrior<NUM_ADDITIONAL_PARAMETERS> parametersPrior;
		
		GroupingSet initialGroupings; //For resetting to when starting a new chain.
		AugmentedParameters<NUM_ADDITIONAL_PARAMETERS> initialParameters = AugmentedParameters<NUM_ADDITIONAL_PARAMETERS>(data, initialGroupings, {{data.guessErrorVariance()}});
		
		GroupingSet currentGroupings = initialGroupings;
		AugmentedParameters<NUM_ADDITIONAL_PARAMETERS> currentParameters = initialParameters;
		
		GroupingBooleanSet isChangingGroupings;
		
		//These do not *need* default values, as a jump will be proposed before they are used, but need some value as Grouping does not have a default constructor.
		GroupingSet proposedGroupings = initialGroupings;
		AugmentedParameters<NUM_ADDITIONAL_PARAMETERS> proposedParameters = initialParameters;
		
		JumpType proposedJumpType;
		
		//Additional useful numbers.
		double transModelJumpProbabilityMultiplier = findTransModelJumpProbabilityMultiplier();
		double growthRateApproximatePosteriorVariance = guessInitialGrowthRateApproximatePosteriorVariance();
		double competitionCoefficientApproximatePosteriorVariance = guessInitialCompetitionCoefficientApproximatePosteriorVariance();
		AdditionalParametersVector additionalParametersApproximatePosteriorVariance = guessInitialAdditionalParametersApproximatePosteriorVariance();
		double withinModelJumpVarianceMultiplier = 1.0;
		double transModelJumpVarianceMultiplier = 1.0;
	public:
		ReversibleJumpSolver(Model model, Data data, Hyperprior hyperprior, AugmentedParametersPrior<NUM_ADDITIONAL_PARAMETERS> parametersPrior, GroupingSet groupings, GroupingBooleanSet isChangingGroupings):
			Solver(model, data),
			hyperprior(hyperprior),
			parametersPrior(parametersPrior),
			initialGroupings(groupings),
			isChangingGroupings(isChangingGroupings)
			{};
		
		void setSeed(size_t seed) {randomGenerator.seed(seed);};
		
		const GroupingSet &getGroupings() const {return currentGroupings;}
		const Grouping &getGrouping(GroupingType groupingType) const {return currentGroupings[groupingType];}
		const AugmentedParameters<NUM_ADDITIONAL_PARAMETERS> &getParameters() const {return currentParameters;}
	protected:
		double getTransModelJumpProbability(GroupingType groupingType, MoveType moveType) const;
		double getTransModelJumpProbability(GroupingType groupingType, MoveType moveType, bool reverse) const;
		size_t getNumTransModelJumps(GroupingType groupingType, MoveType moveType) const;
		
		//The function to get transModelJumpProbabilityMultiplier during initialisation, and some helper functions.
		double findTransModelJumpProbabilityMultiplier() const;
		double findUnscaledMaxTransModelJumpProbability(GroupingSizeSet groupingSizes, size_t recursionLevel) const;
		double findUnscaledMaxTransModelJumpProbability(GroupingSizeSet groupingSizes) const;
		double getUnscaledTransModelJumpProbability(GroupingSizeSet sourceGroupingSizes, GroupingType groupingType, MoveType moveType) const;
		double getUnscaledTransModelJumpProbability(GroupingSizeSet sourceGroupingSizes, GroupingType groupingType, MoveType moveType, bool reverse) const;
		double getUnscaledTransModelJumpProbability(GroupingSizeSet sourceGroupingSizes, GroupingSizeSet destGroupingSizes) const;
		
		double guessInitialGrowthRateApproximatePosteriorVariance() const;
		double guessInitialCompetitionCoefficientApproximatePosteriorVariance() const;
		AdditionalParametersVector guessInitialAdditionalParametersApproximatePosteriorVariance() const;
		
		double getJumpVarianceMultiplier(JumpType jumpType) const;
		Distribution<double> getGrowthRateJumpDistribution(JumpType jumpType) const;
		Distribution<double> getCompetitionCoefficientJumpDistribution(JumpType jumpType) const;
		std::array<Distribution<double>, NUM_ADDITIONAL_PARAMETERS> getAdditionalParametersJumpDistribution(JumpType jumpType) const;
		Distribution<double> getTransModelJumpDistribution(GroupingType groupingType, JumpType jumpType) const;
		
		double getRandomProbability();
		
		//The "propose" functions return the jumping density component of the acceptance ratio (including the Jacobian determinant).
		double proposeTransModelJump(GroupingType groupingType, MoveType moveType, size_t index);
		double proposeWithinModelJump();
		double proposeJump();
		
		void acceptJump();
		void rejectJump();
		
		double getLogLikelihood(double observation, double prediction, const AugmentedParameters<NUM_ADDITIONAL_PARAMETERS> &parameters) const;
		double getLogLikelihoodRatio();
		double getLogPriorRatio() const;
		
		bool makeJump(bool canTransModelJump);
		void burnIn(size_t numJumps, bool canTransModelJump);
	public:
		bool makeJump() {
			return makeJump(true);
		};
		
		void burnIn(size_t numJumps) {
			burnIn(numJumps, true);
		};
		
		void dialIn(size_t jumpsPerDial, size_t numDials);
		
		void resetChain();
};

#endif
