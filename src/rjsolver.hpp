#ifndef REVERSIBLE_JUMP_SOLVER_HPP
#define REVERSIBLE_JUMP_SOLVER_HPP

#include "lattice.hpp"
#include "solver.hpp"

class ReversibleJumpSolver : public Solver {
	public:
		typedef std::function<double(GroupingSet groupings)> HyperpriorFunc;
		
		static double aicHyperprior(GroupingSet groupings);
	protected:
		typedef std::array<size_t, NUM_GROUPING_TYPES> GroupingIndexSet;
		typedef std::array<bool, NUM_GROUPING_TYPES> GroupingBooleanSet;
		
		GroupingLattice groupingLattice;
		HyperpriorFunc hyperpriorFunc;
		
		GroupingIndexSet currentGroupings;
		AugmentedParameters<1> currentParameters;
		
		GroupingBooleanSet isChangingGroupings;
		
		//If isProposing, then base class functions such as getResiduals() will use the propsed groupings.
		bool isProposing = false;
		
		GroupingIndexSet proposedGroupings;
		AugmentedParameters<1> proposedParameters;
		
		//Additional useful numbers.
		double transModelJumpProbabilityMultiplier;
		//TODO: Replace with approximatePosteriorVariance for each of the three.
		//TODO: Function to get variance, based on approximatePosteriorVariance and numDimensions.
		//TODO: Functions to get each of the three variances, plus trans-dimensional jump variance, the last given a GroupingType.
		//TODO: A dialIn() function, to set the approximatePosteriorVariances.
		double growthRateJumpVariance = 1e-3;
		double competitionCoefficientJumpVariance = 1e-5;
		double varianceJumpVariance = 1e-1;
	public:
		ReversibleJumpSolver(Data data, HyperpriorFunc hyperpriorFunc, GroupingSet groupings, GroupingBooleanSet isChangingGroupings):
			Solver(data), groupingLattice(GroupingLattice(data.numSpecies)), hyperpriorFunc(hyperpriorFunc), isChangingGroupings(isChangingGroupings) {
				currentGroupings = getGroupingIndices(groupings);
				currentParameters = AugmentedParameters<1>(data, groupings, {data.getResponseVariance()});
				//TODO: Set jump variances.
				transModelJumpProbabilityMultiplier = getTransModelJumpProbabilityMultiplier();
			};
		ReversibleJumpSolver(Data data, GroupingSet groupings, GroupingBooleanSet isChangingGroupings):
			ReversibleJumpSolver(data, aicHyperprior, groupings, isChangingGroupings) {};
	protected:
		//Functions to say that particular elements have changed and mark appropriate things as dirty.
		void dirtyDataSubclass() override {}
		void dirtyGroupingSubclass(GroupingType groupingType) override {}
		
		GroupingSet getGroupings(GroupingIndexSet groupingIndices) const;
		GroupingIndexSet getGroupingIndices(GroupingSet groupings) const;
		
		void setIsProposing(bool b);
	public:
		//Functions to retrieve groupings.
		const Grouping &getGrouping(GroupingType groupingType) const override {
			return groupingLattice.getGrouping(isProposing ? proposedGroupings[groupingType] : currentGroupings[groupingType]);
		}
		
		GroupingSet getGroupings() const {
			return getGroupings(isProposing ? proposedGroupings : currentGroupings);
		}
		
		const AugmentedParameters<1> &getParameters() const {
			return isProposing ? proposedParameters : currentParameters;
		}
	protected:
		//The following functions have an arbitrary version, which is used by getTransModelJumpProbabilityMultiplier,
		//and a canonical version, which uses the current parameters.
		//The type parameters decide which grouping the jump changes, and how it changes it.
		double getTransModelJumpProbability(GroupingIndexSet sourceGroupingIndices, GroupingIndexSet destGroupingIndices, double multiplier) const;
		double getTransModelJumpProbability(GroupingIndexSet sourceGroupingIndices, GroupingIndexSet destGroupingIndices) const;
		std::vector<double> getTransModelJumpProbabilities(GroupingType groupingType, MoveType moveType, GroupingIndexSet groupingIndices, double multiplier) const;
		std::vector<double> getTransModelJumpProbabilities(GroupingType groupingType, MoveType moveType) const;
		
		double getTransModelJumpProbabilityMultiplier() const;
		double getUnscaledMaxTransModelJumpProbability(size_t recursionLevel, GroupingIndexSet groupingIndices) const; //A helper function for the above.
		double getUnscaledTotalTransModelJumpProbability(GroupingIndexSet groupingIndices) const; //Another helper function.
		
		//The "propose" functions return the jumping density component of the acceptance ratio (including the Jacobian determinant).
		double proposeTransModelJump(GroupingType groupingType, MoveType moveType, size_t newGroupingIndex);
		double proposeWithinModelJump();
		double proposeJump();
		
		void acceptJump();
		void rejectJump();
		
		double getErrorVariance() const;
		Eigen::VectorXd getResiduals();
		
		double getLikelihoodRatio(Eigen::VectorXd sourceResiduals, Eigen::VectorXd destResiduals, double sourceErrorVariance, double destErrorVariance);
		double getPrior() const;
		
		void makeJump(bool canTransModelJump);
	public:
		void makeJump() {
			makeJump(true);
		};
		
		void burnIn(size_t numJumps);
};

#endif
