#ifndef REVERSIBLE_JUMP_SOLVER_HPP
#define REVERSIBLE_JUMP_SOLVER_HPP

#include "lattice.hpp"
#include "solver.hpp"

class ReversibleJumpSolver : public Solver {
	public:
		typedef std::function<double(Grouping growthGrouping, Grouping rowGrouping, Grouping colGrouping)> HyperpriorFunc;
		
		static double aicHyperprior(Grouping growthGrouping, Grouping rowGrouping, Grouping colGrouping);
	protected:
		GroupingLattice groupingLattice;
		HyperpriorFunc hyperpriorFunc;
		
		size_t growthGrouping;
		size_t rowGrouping;
		size_t colGrouping;
		Parameters parameters;
		
		//If isProposing, then base class functions such as getResiduals() will use the propsed groupings.
		bool isProposing = false;
		
		size_t proposedGrowthGrouping;
		size_t proposedRowGrouping;
		size_t proposedColGrouping;
		Parameters proposedParameters;
	public:
		ReversibleJumpSolver(Data data, HyperpriorFunc hyperpriorFunc):
			Solver(data), groupingLattice(GroupingLattice(data.numSpecies)), hyperpriorFunc(hyperpriorFunc), growthGrouping(0), rowGrouping(0), colGrouping(0) {};
		ReversibleJumpSolver(Data data):
			ReversibleJumpSolver(data, aicHyperprior) {};
	protected:
		//Functions to say that particular elements have changed and mark appropriate things as dirty.
		void dirtyDataSubclass() {}
		void dirtyGrowthGroupingSubclass() {}
		void dirtyRowGroupingSubclass() {}
		void dirtyColGroupingSubclass() {}
	public:
		//Functions to retrieve groupings.
		const Grouping &getGrowthGrouping() const {
			return groupingLattice.getGrouping(isProposing ? proposedGrowthGrouping : growthGrouping);
		}
		const Grouping &getRowGrouping() const {
			return groupingLattice.getGrouping(isProposing ? proposedRowGrouping : rowGrouping);
		}
		const Grouping &getColGrouping() const {
			return groupingLattice.getGrouping(isProposing ? proposedColGrouping : colGrouping);
		}
		
		const Parameters &getParameters() const {
			return isProposing ? proposedParameters : parameters;
		}
};

#endif
