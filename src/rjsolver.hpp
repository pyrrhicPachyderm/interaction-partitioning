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
	public:
		ReversibleJumpSolver(Data data, HyperpriorFunc hyperpriorFunc, GroupingBooleanSet isChangingGroupings):
			Solver(data), groupingLattice(GroupingLattice(data.numSpecies)), hyperpriorFunc(hyperpriorFunc), isChangingGroupings(isChangingGroupings) {};
		ReversibleJumpSolver(Data data, GroupingBooleanSet isChangingGroupings):
			ReversibleJumpSolver(data, aicHyperprior, isChangingGroupings) {};
	protected:
		//Functions to say that particular elements have changed and mark appropriate things as dirty.
		void dirtyDataSubclass() override {}
		void dirtyGroupingSubclass(GroupingType groupingType) override {}
		
		GroupingSet getGroupings(GroupingIndexSet groupingIndices) const;
	public:
		//Functions to retrieve groupings.
		const Grouping &getGroupingSubclass(GroupingType groupingType) const override {
			return groupingLattice.getGrouping(isProposing ? proposedGroupings[groupingType] : currentGroupings[groupingType]);
		}
		
		const Parameters &getParameters() const {
			return isProposing ? proposedParameters : currentParameters;
		}
};

#endif
