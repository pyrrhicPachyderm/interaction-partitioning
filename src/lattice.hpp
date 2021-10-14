#ifndef LATTICE_HPP
#define LATTICE_HPP

#include "grouping.hpp"

//A GroupingLattice is the set of all partitions of a given number of species, partially ordered by refinement.
//It forms a lattice, in the mathematical sense that there is a unique supremum and a unique infimum.
class GroupingLattice {
	public:
		const size_t numSpecies;
	protected:
		std::vector<Grouping> groupings;
		
		void constructAllGroupings();
	public:
		GroupingLattice(size_t numSpecies): numSpecies(numSpecies) {
			constructAllGroupings();
		};
};

#endif
