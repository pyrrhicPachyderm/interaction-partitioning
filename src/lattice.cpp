#include <algorithm>
#include "lattice.hpp"

void GroupingLattice::constructAllGroupings() {
	Grouping grouping(numSpecies);
	do {
		groupings.push_back(grouping);
	} while(grouping.advance());
}

size_t GroupingLattice::findMatch(const std::vector<size_t> &grouping) const {
	//TODO:
	//An exhaustive search is horribly inefficient: O(n * Bell_n)
	//This can reduced to just O(n) with an appropriate tree structure,
	//with one level per species and the nodes of the next level indexed by the next group.
	//Implement this.
	//For now, the exhaustive search.
	
	for(size_t i = 0; i < groupings.size(); i++) {
		if(groupings[i].isMatch(grouping)) return i;
	}
	
	__builtin_unreachable();
}

size_t GroupingLattice::getIndex(const Grouping &grouping) const {
	return findMatch(grouping.getGroups());
}

//Takes a grouping and the numbers of two groups to merge, and returns the grouping with them merged.
static std::vector<size_t> mergeGroups(std::vector<size_t> grouping, size_t g1, size_t g2) {
	for(size_t i = 0; i < grouping.size(); i++) {
		if(grouping[i] == g1) grouping[i] = g2;
	}
	return Grouping::fixGrouping(grouping);
}

void GroupingLattice::constructAdjacencyLists() {
	//Set the adjacency lists to the right length.
	moveDests.fill(std::vector<std::vector<size_t>>(groupings.size()));
	moves.fill(std::vector<std::vector<GroupingMove>>(groupings.size()));
	
	for(size_t sourceIndex = 0; sourceIndex < groupings.size(); sourceIndex++) {
		const Grouping &source = groupings[sourceIndex];
		size_t numGroups = source.getNumGroups();
		
		//Iterate over each pair of groups, trying to merge them.
		for(size_t g1 = 0; g1 < numGroups; g1++) {
			for(size_t g2 = g1 + 1; g2 < numGroups; g2++) {
				std::vector<size_t> mergedGrouping = mergeGroups(source.getGroups(), g1, g2);
				size_t destIndex = findMatch(mergedGrouping);
				
				GroupingMove move(source.getGroups(), mergedGrouping, g1, g2);
				
				moveDests[MERGE][sourceIndex].push_back(destIndex);
				moveDests[SPLIT][destIndex].push_back(sourceIndex);
				moves[MERGE][sourceIndex].push_back(move);
				moves[SPLIT][destIndex].push_back(move);
			}
		}
	}
}

GroupingMove::GroupingMove(std::vector<size_t> splitGrouping, std::vector<size_t> mergedGrouping, size_t g1, size_t g2) {
	assert(splitGrouping.size() == mergedGrouping.size());
	
	splitGroups = std::make_pair(g1, g2);
	
	size_t numSplitGroups = *std::max_element(splitGrouping.begin(), splitGrouping.end()) + 1; //+1 for zero indexing.
	size_t numMergedGroups = *std::max_element(mergedGrouping.begin(), mergedGrouping.end()) + 1; //+1 for zero indexing.
	
	assert(numSplitGroups == numMergedGroups+1);
	
	mergeMap = std::vector<size_t>(numSplitGroups);
	splitMap = std::vector<size_t>(numMergedGroups);
	
	for(size_t i = 0; i < splitGrouping.size(); i++) {
		if(splitGrouping[i] == splitGroups.first || splitGrouping[i] == splitGroups.second) {
			mergedGroup = mergedGrouping[i];
		} else {
			mergeMap[splitGrouping[i]] = mergedGrouping[i];
			splitMap[mergedGrouping[i]] = splitGrouping[i];
		}
	}
}
