#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <assert.h>
#include "grouping.hpp"

class GroupingMove;

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
		std::vector<std::vector<GroupingMove>> mergeMoves;
		std::vector<std::vector<GroupingMove>> splitMoves;
		
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

//A move between two partitions, by merging two groups or splitting one.
//The same move represents both directions.
class GroupingMove {
	protected:
		size_t mergedGroup;
		std::pair<size_t,size_t> splitGroups;
		
		//The maps are group-to-group mappings for the groups that stay the same.
		//Index into them using groups from the source grouping, and they return groups from the destination grouping.
		std::vector<size_t> mergeMap;
		std::vector<size_t> splitMap;
	public:
		//g1 and g2 are the two group numbers that are merged.
		GroupingMove(std::vector<size_t> splitGrouping, std::vector<size_t> mergedGrouping, size_t g1, size_t g2);
		
		size_t getMergeMapping(size_t group) {
			assert(group != splitGroups.first && group != splitGroups.second);
			return mergeMap[group];
		}
		
		size_t getSplitMapping(size_t group) {
			assert(group != mergedGroup);
			return splitMap[group];
		}
};

#endif
