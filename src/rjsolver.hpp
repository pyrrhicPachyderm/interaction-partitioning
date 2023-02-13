#ifndef REVERSIBLE_JUMP_SOLVER_HPP
#define REVERSIBLE_JUMP_SOLVER_HPP

#include "groupingmove.hpp"
#include "solver.hpp"
#include "priors.hpp"

#define INITIAL_APPROXIMATE_POSTERIOR_VARIANCE_MULTIPLIER 0.01

class ReversibleJumpSolver : public Solver {
	public:
		const static size_t NUM_ADDITIONAL_PARAMETERS = 1;
		typedef AugmentedParameters<NUM_ADDITIONAL_PARAMETERS>::AdditionalParametersVector AdditionalParametersVector;
		enum JumpType {MERGE_JUMP, SPLIT_JUMP, WITHIN_JUMP, NUM_JUMP_TYPES};
	protected:
		typedef std::array<bool, NUM_GROUPING_TYPES> GroupingBooleanSet;
		
		Hyperprior hyperprior;
		AugmentedParametersPrior<NUM_ADDITIONAL_PARAMETERS> parametersPrior;
		
		GroupingSet initialGroupings; //For resetting to when starting a new chain.
		AugmentedParameters<NUM_ADDITIONAL_PARAMETERS> initialParameters;
		
		GroupingSet currentGroupings = initialGroupings;
		AugmentedParameters<NUM_ADDITIONAL_PARAMETERS> currentParameters = initialParameters;
		
		GroupingBooleanSet isChangingGroupings;
		
		//These do not *need* default values, as a jump will be proposed before they are used, but need some value as Grouping does not have a default constructor.
		GroupingSet proposedGroupings = initialGroupings;
		AugmentedParameters<NUM_ADDITIONAL_PARAMETERS> proposedParameters = initialParameters;
		
		JumpType proposedJumpType;
		
		//Additional useful numbers.
		double transModelJumpProbabilityMultiplier;
		double growthRateApproximatePosteriorVariance = data.guessGrowthRate() * INITIAL_APPROXIMATE_POSTERIOR_VARIANCE_MULTIPLIER;
		double competitionCoefficientApproximatePosteriorVariance = data.guessCompetitionCoefficientMagnitude() * INITIAL_APPROXIMATE_POSTERIOR_VARIANCE_MULTIPLIER;
		AdditionalParametersVector additionalParametersApproximatePosteriorVariance = {{data.guessErrorVariance() * INITIAL_APPROXIMATE_POSTERIOR_VARIANCE_MULTIPLIER}};
		double jumpVarianceMultiplier = 1.0;
	public:
		ReversibleJumpSolver(Model model, Data data, Hyperprior hyperprior, AugmentedParametersPrior<NUM_ADDITIONAL_PARAMETERS> parametersPrior, GroupingSet groupings, GroupingBooleanSet isChangingGroupings):
			Solver(model, data),
			hyperprior(hyperprior),
			parametersPrior(parametersPrior),
			initialGroupings(groupings),
			initialParameters(AugmentedParameters<NUM_ADDITIONAL_PARAMETERS>(data, groupings, {{data.guessErrorVariance()}})),
			isChangingGroupings(isChangingGroupings)
			{
				transModelJumpProbabilityMultiplier = getTransModelJumpProbabilityMultiplier();
			};
		
		const GroupingSet &getGroupings() const {return currentGroupings;}
		const Grouping &getGrouping(GroupingType groupingType) const {return currentGroupings[groupingType];}
		const AugmentedParameters<NUM_ADDITIONAL_PARAMETERS> &getParameters() const {return currentParameters;}
	protected:
		double getTransModelJumpProbability(GroupingType groupingType, MoveType moveType) const;
		double getTransModelJumpProbability(GroupingType groupingType, MoveType moveType, bool reverse) const;
		size_t getNumTransModelJumps(GroupingType groupingType, MoveType moveType) const;
		
		//The function to get transModelJumpProbabilityMultiplier during initialisation, and some helper functions.
		double getTransModelJumpProbabilityMultiplier() const;
		double getUnscaledMaxTransModelJumpProbability(GroupingSizeSet groupingSizes, size_t recursionLevel) const;
		double getUnscaledMaxTransModelJumpProbability(GroupingSizeSet groupingSizes) const;
		double getUnscaledTransModelJumpProbability(GroupingSizeSet sourceGroupingSizes, GroupingType groupingType, MoveType moveType) const;
		double getUnscaledTransModelJumpProbability(GroupingSizeSet sourceGroupingSizes, GroupingType groupingType, MoveType moveType, bool reverse) const;
		double getUnscaledTransModelJumpProbability(GroupingSizeSet sourceGroupingSizes, GroupingSizeSet destGroupingSizes) const;
		
		Distribution<double> getGrowthRateJumpDistribution() const;
		Distribution<double> getCompetitionCoefficientJumpDistribution() const;
		std::array<Distribution<double>, NUM_ADDITIONAL_PARAMETERS> getAdditionalParametersJumpDistribution() const;
		Distribution<double> getTransModelJumpDistribution(GroupingType groupingType) const;
		
		//The "propose" functions return the jumping density component of the acceptance ratio (including the Jacobian determinant).
		double proposeTransModelJump(GroupingType groupingType, MoveType moveType, size_t index);
		double proposeWithinModelJump();
		double proposeJump();
		
		void acceptJump();
		void rejectJump();
		
		Distribution<double> getErrorDistribution(const AugmentedParameters<NUM_ADDITIONAL_PARAMETERS> &parameters) const;
		
		double getLikelihoodRatio();
		double getPriorDensity(const AugmentedParameters<NUM_ADDITIONAL_PARAMETERS> &parameters, const GroupingSet &groupings) const;
		
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
