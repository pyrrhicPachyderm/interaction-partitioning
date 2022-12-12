#ifndef GROUPING_MOVE_HPP
#define GROUPING_MOVE_HPP

#include <assert.h>
#include <array>
#include "grouping.hpp"

enum MoveType {MERGE, SPLIT, NUM_MOVE_TYPES};

class GroupingMove;

//A GroupingLattice is the set of all partitions of a given number of species, partially ordered by refinement.
//It forms a lattice, in the mathematical sense that there is a unique supremum and a unique infimum.
class GroupingLattice {
	public:
		const size_t numSpecies;
	protected:
		std::vector<Grouping> groupings;
		
		//Adjacency lists for upward and downward moves.
		std::array<std::vector<std::vector<size_t>>, NUM_MOVE_TYPES> moveDests;
		std::array<std::vector<std::vector<GroupingMove>>, NUM_MOVE_TYPES> moves;
		
		//Takes a grouping, and finds its index in groupings.
		size_t findMatch(const std::vector<size_t> &grouping) const;
		
		void constructAllGroupings();
		void constructAdjacencyLists();
	public:
		GroupingLattice(size_t numSpecies): numSpecies(numSpecies) {
			constructAllGroupings();
			constructAdjacencyLists();
		};
		
		size_t getNumGroupings() const {
			return groupings.size();
		}
		
		const Grouping &getGrouping(size_t index) const {
			return groupings[index];
		}
		
		size_t getIndex(const Grouping &grouping) const;
		
		size_t getMoveDest(MoveType moveType, size_t sourceIndex, size_t adjIndex) const {
			return moveDests[moveType][sourceIndex][adjIndex];
		}
		
		const std::vector<size_t> &getMoveDests(MoveType moveType, size_t sourceIndex) const {
			return moveDests[moveType][sourceIndex];
		}
		
		const GroupingMove &getMove(MoveType moveType, size_t sourceIndex, size_t adjIndex) const {
			return moves[moveType][sourceIndex][adjIndex];
		}
};

//A move between two partitions, by merging two groups or splitting one.
//The same move represents both directions.
class GroupingMove {
	protected:
		size_t mergedGroup;
		std::pair<size_t,size_t> splitGroups;
		
		size_t mergedGroupSize;
		std::pair<size_t,size_t> splitGroupSizes;
		
		//The maps are group-to-group mappings for the groups that stay the same.
		//Index into them using groups from the source grouping, and they return groups from the destination grouping.
		std::vector<size_t> mergeMap;
		std::vector<size_t> splitMap;
	public:
		//g1 and g2 are the two group numbers that are merged.
		GroupingMove(std::vector<size_t> splitGrouping, std::vector<size_t> mergedGrouping, size_t g1, size_t g2);
		
		//A move represents both directions, so grouping1 and grouping2 are interchangeable.
		GroupingMove(const Grouping &grouping1, const Grouping &grouping2);
		
		size_t getMergedGroup() const {
			return mergedGroup;
		}
		
		const std::pair<size_t,size_t> &getSplitGroups() const {
			return splitGroups;
		}
		
		size_t getMergedGroupSize() const {
			return mergedGroupSize;
		}
		
		const std::pair<size_t,size_t> &getSplitGroupSizes() const {
			return splitGroupSizes;
		}
		
		//Three functions to map old indices to new indices in a move.
		size_t getMergeMapping(size_t group) const {
			assert(group != splitGroups.first && group != splitGroups.second);
			return mergeMap[group];
		}
		
		size_t getSplitMapping(size_t group) const {
			assert(group != mergedGroup);
			return splitMap[group];
		}
		
		size_t getMapping(MoveType moveType, size_t group) const {
			if(moveType == MERGE) return getMergeMapping(group);
			if(moveType == SPLIT) return getSplitMapping(group);
			__builtin_unreachable();
		}
};

#endif
