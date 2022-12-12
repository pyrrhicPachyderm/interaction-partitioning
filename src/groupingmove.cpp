#include <algorithm>
#include "groupingmove.hpp"

GroupingMove::GroupingMove(std::vector<size_t> splitGrouping, std::vector<size_t> mergedGrouping, size_t g1, size_t g2) {
	assert(splitGrouping.size() == mergedGrouping.size());
	
	splitGroups = std::make_pair(g1, g2);
	
	size_t numSplitGroups = *std::max_element(splitGrouping.begin(), splitGrouping.end()) + 1; //+1 for zero indexing.
	size_t numMergedGroups = *std::max_element(mergedGrouping.begin(), mergedGrouping.end()) + 1; //+1 for zero indexing.
	
	assert(numSplitGroups == numMergedGroups+1);
	
	mergeMap = std::vector<size_t>(numSplitGroups);
	splitMap = std::vector<size_t>(numMergedGroups);
	
	mergedGroupSize = 0; //A counter we will add to in the loop.
	splitGroupSizes = std::make_pair(0, 0); //More counters.
	
	for(size_t i = 0; i < splitGrouping.size(); i++) {
		//Manage the mapping.
		if(splitGrouping[i] == splitGroups.first || splitGrouping[i] == splitGroups.second) {
			mergedGroup = mergedGrouping[i];
		} else {
			mergeMap[splitGrouping[i]] = mergedGrouping[i];
			splitMap[mergedGrouping[i]] = splitGrouping[i];
		}
		
		//Manage the group sizes.
		if(splitGrouping[i] == splitGroups.first) {
			mergedGroupSize++;
			splitGroupSizes.first++;
		}
		if(splitGrouping[i] == splitGroups.second) {
			mergedGroupSize++;
			splitGroupSizes.second++;
		}
	}
}

GroupingMove::GroupingMove(const Grouping &grouping1, const Grouping &grouping2) {
	assert(grouping1.numSpecies == grouping2.numSpecies);
	
	size_t numGroups1 = grouping1.getNumGroups();
	size_t numGroups2 = grouping2.getNumGroups();
	
	assert(abs((long long int)numGroups1 - (long long int)numGroups2) == 1);
	
	const Grouping &mergedGrouping = numGroups1 < numGroups2 ? grouping1 : grouping2;
	const Grouping &splitGrouping = numGroups1 > numGroups2 ? grouping1 : grouping2;
	size_t numMergedGroups = std::min(numGroups1, numGroups2);
	size_t numSplitGroups = std::max(numGroups1, numGroups2);
	
	mergeMap = std::vector<size_t>(numSplitGroups);
	splitMap = std::vector<size_t>(numMergedGroups);
	bool foundSplit = false;
	std::vector<bool> splitMapped(numMergedGroups, false);
	
	for(size_t i = 0; i < splitGrouping.numSpecies; i++) {
		size_t merged = mergedGrouping.getGroup(i);
		size_t split = splitGrouping.getGroup(i);
		if(!foundSplit && splitMapped[merged] && splitMap[merged] != split) {
			foundSplit = true;
			mergedGroup = merged;
			splitGroups = std::make_pair(splitMap[merged], split);
		}
		
		//These will erroneously set even for the split and merged groupings,
		//but asserts on the getter functions ensure those vaues will never be used.
		mergeMap[split] = merged;
		splitMap[merged] = split;
		splitMapped[merged] = true;
	}
	
	//Set mergedGroupSize and splitGroupSizes.
	std::vector<size_t> mergedGroupSizeVector = mergedGrouping.getGroupSizes();
	std::vector<size_t> splitGroupSizeVector = splitGrouping.getGroupSizes();
	mergedGroupSize = mergedGroupSizeVector[mergedGroup];
	splitGroupSizes = std::make_pair(splitGroupSizeVector[splitGroups.first], splitGroupSizeVector[splitGroups.second]);
}
