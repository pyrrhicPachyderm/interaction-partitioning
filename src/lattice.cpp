#include "lattice.hpp"

void GroupingLattice::constructAllGroupings() {
	Grouping grouping(numSpecies);
	do {
		groupings.push_back(grouping);
	} while(grouping.advance());
}

size_t GroupingLattice::findMatch(std::vector<size_t> grouping) {
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

//Takes a grouping and the numbers of two groups to merge, and returns the grouping with them merged.
static std::vector<size_t> mergeGroups(std::vector<size_t> grouping, size_t g1, size_t g2) {
	for(size_t i = 0; i < grouping.size(); i++) {
		if(grouping[i] == g1) grouping[i] = g2;
	}
	return Grouping::fixGrouping(grouping);
}

void GroupingLattice::constructAdjacencyLists() {
	//Set the adjacency lists to the right length.
	merges = std::vector<std::vector<size_t>>(groupings.size());
	splits = std::vector<std::vector<size_t>>(groupings.size());
	
	for(size_t sourceIndex = 0; sourceIndex < groupings.size(); sourceIndex++) {
		const Grouping &source = groupings[sourceIndex];
		size_t numGroups = source.getNumGroups();
		
		//Iterate over each pair of groups, trying to merge them.
		for(size_t g1 = 0; g1 < numGroups; g1++) {
			for(size_t g2 = g1 + 1; g2 < numGroups; g2++) {
				std::vector<size_t> mergedGrouping = mergeGroups(source.getGroups(), g1, g2);
				size_t destIndex = findMatch(mergedGrouping);
				
				merges[sourceIndex].push_back(destIndex);
				splits[destIndex].push_back(sourceIndex);
			}
		}
	}
}
