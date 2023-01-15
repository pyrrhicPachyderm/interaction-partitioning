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
		
		GroupingSet initialGroupings; //For resetting to when starting a new chain.
		AugmentedParameters<NUM_ADDITIONAL_PARAMETERS> initialParameters;
		
		GroupingSet currentGroupings = initialGroupings;
		AugmentedParameters<NUM_ADDITIONAL_PARAMETERS> currentParameters = initialParameters;
		
		GroupingBooleanSet isChangingGroupings;
		
		//If isProposing, then base class functions such as getResiduals() will use the propsed groupings.
		bool isProposing = false;
		
		//These do not *need* default values, as isProposing is initially false, but need some value as Grouping does not have a default constructor.
		GroupingSet proposedGroupings = initialGroupings;
		AugmentedParameters<NUM_ADDITIONAL_PARAMETERS> proposedParameters = initialParameters;
		
		JumpType proposedJumpType;
		
		//Additional useful numbers.
		double transModelJumpProbabilityMultiplier;
		double growthRateApproximatePosteriorVariance = data.guessGrowthRate() * INITIAL_APPROXIMATE_POSTERIOR_VARIANCE_MULTIPLIER;
		double competitionCoefficientApproximatePosteriorVariance = data.guessCompetitionCoefficientMagnitude() * INITIAL_APPROXIMATE_POSTERIOR_VARIANCE_MULTIPLIER;
		AdditionalParametersVector additionalParametersApproximatePosteriorVariance = {{data.getResponseVariance() * INITIAL_APPROXIMATE_POSTERIOR_VARIANCE_MULTIPLIER}};
		double jumpVarianceMultiplier = 1.0;
		double competitionCoefficientPriorVariance = pow(data.guessCompetitionCoefficientMagnitude(), 2);
	public:
		ReversibleJumpSolver(Data data, Hyperprior hyperprior, GroupingSet groupings, GroupingBooleanSet isChangingGroupings):
			Solver(data),
			hyperprior(hyperprior),
			initialGroupings(groupings),
			initialParameters(AugmentedParameters<NUM_ADDITIONAL_PARAMETERS>(data, groupings, {{data.getResponseVariance()}})),
			isChangingGroupings(isChangingGroupings)
			{
				transModelJumpProbabilityMultiplier = getTransModelJumpProbabilityMultiplier();
			};
	protected:
		//Functions to say that particular elements have changed and mark appropriate things as dirty.
		void dirtyDataSubclass() override {}
		void dirtyGroupingSubclass(GroupingType groupingType) override {}
		
		void setIsProposing(bool b);
	public:
		//Functions to retrieve groupings.
		const Grouping &getGrouping(GroupingType groupingType) const override {
			return isProposing ? proposedGroupings[groupingType] : currentGroupings[groupingType];
		}
		
		GroupingSet getGroupings() const {
			return isProposing ? proposedGroupings : currentGroupings;
		}
		
		const AugmentedParameters<NUM_ADDITIONAL_PARAMETERS> &getParameters() const {
			return isProposing ? proposedParameters : currentParameters;
		}
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
		
		Distribution<double> getErrorDistribution() const;
		Eigen::VectorXd getResiduals();
		
		double getLikelihoodRatio(Eigen::VectorXd sourceResiduals, Eigen::VectorXd destResiduals, Distribution<double> sourceErrorDistribution,  Distribution<double> destErrorDistribution);
		double getPriorDensity() const;
		
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
