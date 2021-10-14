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
		
		//Adjacency lists for upward and downward moves.
		std::vector<std::vector<size_t>> merges;
		std::vector<std::vector<size_t>> splits;
		
		//Takes a grouping, and finds its index in groupings.
		size_t findMatch(std::vector<size_t> grouping);
		
		void constructAllGroupings();
		void constructAdjacencyLists();
	public:
		GroupingLattice(size_t numSpecies): numSpecies(numSpecies) {
			constructAllGroupings();
			constructAdjacencyLists();
		};
};

#endif
