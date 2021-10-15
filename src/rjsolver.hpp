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
		Parameters currentParameters;
		
		GroupingBooleanSet isChangingGroupings;
		
		//If isProposing, then base class functions such as getResiduals() will use the propsed groupings.
		bool isProposing = false;
		
		GroupingIndexSet proposedGroupings;
		Parameters proposedParameters;
		
		//Additional useful numbers.
		double transModelJumpProbabilityMultiplier;
	public:
		ReversibleJumpSolver(Data data, HyperpriorFunc hyperpriorFunc, GroupingBooleanSet isChangingGroupings):
			Solver(data), groupingLattice(GroupingLattice(data.numSpecies)), hyperpriorFunc(hyperpriorFunc), isChangingGroupings(isChangingGroupings) {
				transModelJumpProbabilityMultiplier = getTransModelJumpProbabilityMultiplier();
			};
		ReversibleJumpSolver(Data data, GroupingBooleanSet isChangingGroupings):
			ReversibleJumpSolver(data, aicHyperprior, isChangingGroupings) {};
		
		void setGroupings(GroupingSet groupings) {
			currentGroupings = getGroupingIndices(groupings);
		};
	protected:
		//Functions to say that particular elements have changed and mark appropriate things as dirty.
		void dirtyDataSubclass() override {}
		void dirtyGroupingSubclass(GroupingType groupingType) override {}
		
		GroupingSet getGroupings(GroupingIndexSet groupingIndices) const;
		GroupingIndexSet getGroupingIndices(GroupingSet groupings) const;
	public:
		//Functions to retrieve groupings.
		const Grouping &getGroupingSubclass(GroupingType groupingType) const override {
			return groupingLattice.getGrouping(isProposing ? proposedGroupings[groupingType] : currentGroupings[groupingType]);
		}
		
		const Parameters &getParameters() const {
			return isProposing ? proposedParameters : currentParameters;
		}
	protected:
		//The following has an arbitrary version, which is used by getTransModelJumpProbabilityMultiplier,
		//and a canonical version, which uses the current parameters.
		//The template parameters decide which grouping the jump changes, and how it changes it.
		std::vector<double> getTransModelJumpProbabilities(GroupingType groupingType, MoveType moveType, GroupingIndexSet groupingIndices, double multiplier) const;
		std::vector<double> getTransModelJumpProbabilities(GroupingType groupingType, MoveType moveType) const;
		
		double getTransModelJumpProbabilityMultiplier() const;
		double getUnscaledMaxTransModelJumpProbability(size_t recursionLevel, GroupingIndexSet groupingIndices) const; //A helper function for the above.
		double getUnscaledTotalTransModelJumpProbability(GroupingIndexSet groupingIndices) const; //Another helper function.
};

#endif
