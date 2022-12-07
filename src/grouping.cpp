#include <algorithm>
#include "grouping.hpp"

void Grouping::reset() {
	groups = std::vector<size_t>(numSpecies, 0);
}

void Grouping::separate() {
	for(size_t i = 0; i < numSpecies; i++) {
		groups[i] = i;
	}
}

bool Grouping::advance() {
	std::vector<size_t> maxGroups = getMaxGroups();
	return advanceIndex(maxGroups, numSpecies - 1);
}

std::vector<size_t> Grouping::getMaxGroups() const {
	std::vector<size_t> maxGroups(numSpecies, 0);
	//Start the loop from 1, as it must always be the case that groups[0] = 0,
	//hence also maxGroups[0] = 0.
	for(size_t i = 1; i < numSpecies; i++) {
		maxGroups[i] = std::max(groups[i], maxGroups[i-1]);
	}
	return maxGroups;
}

bool Grouping::advanceIndex(const std::vector<size_t> &maxGroups, size_t index) {
	//It is never valid to increment the first element, so reset.
	if(index == 0) {
		reset();
		return false;
	}
	
	//If we have seen an equal or higher group number to the left, this can safely be incremented.
	//We must reset all group numbers to the right of it to 0
	if(groups[index] <= maxGroups[index-1]) {
		groups[index]++;
		for(size_t i = index+1; i < numSpecies; i++) {
			groups[i] = 0;
		}
		return true;
	}
	//Otherwise, it is not valid to increment it, so we recurse.
	return advanceIndex(maxGroups, index-1);
}

size_t Grouping::getNumGroups() const {
	//The number of groups is simply the highest group number present.
	//Plus one, as the groups are zero indexed.
	return *std::max_element(groups.begin(), groups.end()) + 1;
}

std::vector<size_t> Grouping::getGroupSizes() const {
	std::vector<size_t> groupSizes;
	
	//Reserve at least enough capacity, so that no reallocs will be necessary.
	groupSizes.reserve(numSpecies);
	
	for(size_t i = 0; i < numSpecies; i++) {
		if(groups[i] >= groupSizes.size()) {
			//A new group must be added. As a grouping is rhyming scheme, we can assume this is the next group numerically.
			//So we just push a 1 to the end of the vector.
			groupSizes.push_back(1);
		} else {
			groupSizes[groups[i]] += 1;
		}
	}
	
	return groupSizes;
}

size_t Grouping::getNumMerges(size_t numGroups) {
	//Each pair of groups can be merged.
	//So this is (number of groups) choose 2.
	return numGroups * (numGroups - 1) / 2;
}

static inline size_t getNumSplitsSingleGroup(size_t groupSize) {
	//The number of ways of splitting a single group with k elements is:
	//2^(k-1) - 1 if k > 1
	//0 if k = 1
	//The latter is obvious: a group cannot be split in two if it has only one element.
	//For the former:
	//Assume, without loss of generality, that the first element of the group is assigned to the first new group.
	//The other k-1 elements of the group must be assigned to either the first new group or the second; there are 2^(k-1) ways of doing so.
	//However, one of these assigns every element to the first new group and none to the second; this is invalid, so we subtract that one.
	//Lastly, note that if k = 1, 2^(k-1) - 1 = 0, so the conditional is unnecessary.
	return (1 << (groupSize - 1)) - 1;
}

size_t Grouping::getMaxNumSplits(size_t numSpecies, size_t numGroups) {
	//The number of splits for a given grouping is documented in getNumSplitsSingleGroup().
	//It is clear that the greatest number of splits for a given numSpecies and numGroups
	//is achieved by making all but one group of size 1, and the last group as large as possible.
	//Proof:
	//Suppose two different groups both have more than one element.
	//By transferring an element from the smaller of these groups to the larger (or even in the case of equality)
	//the larger group now provides more than twice as many splits.
	//This must more than offset any splits lost by shrinking the smaller group.
	size_t bigGroupSize = numSpecies - numGroups + 1;
	//Only the big group can be split, not any of the singletons, so no sum is necessary here.
	return getNumSplitsSingleGroup(bigGroupSize);
}

size_t Grouping::getNumMerges() const {
	return getNumMerges(getNumGroups());
}

size_t Grouping::getNumSplits() const {
	//The number of possible split moves is the sum of the number of ways of splitting each group.
	std::vector<size_t> groupSizes = getGroupSizes();
	size_t numSplits = 0;
	for(size_t i = 0; i < groupSizes.size(); i++) {
		numSplits += getNumSplitsSingleGroup(groupSizes[i]);
	}
	return numSplits;
}

std::vector<size_t> Grouping::fixGrouping(std::vector<size_t> improperGrouping) {
	//We must make the improperGrouping a proper grouping.
	//Species are in the same group if and only if they have the same number in improperGrouping.
	//However, improperGrouping is not required to be a valid rhyming scheme.
	
	//We will use something akin to countingSort.
	//That is, we will assume the largest number in improperGrouping is not too large.
	//We will suffer issues with space complexity if that is not true.
	
	size_t maxImproperGroup = *std::max_element(improperGrouping.begin(), improperGrouping.end());
	
	//We will have a vector of mappings.
	//Index into this using the number of improperGrouping, and it will give the proper group number.
	//Must be +1 size due to zero indexing.
	std::vector<size_t> mappedGroup(maxImproperGroup+1);
	std::vector<size_t> isMapped(maxImproperGroup+1, false);
	
	std::vector<size_t> properGrouping(improperGrouping.size());
	
	size_t nextGroup = 0;
	for(size_t i = 0; i < improperGrouping.size(); i++) {
		if(!isMapped[improperGrouping[i]]) {
			mappedGroup[improperGrouping[i]] = nextGroup++;
			isMapped[improperGrouping[i]] = true;
		}
		properGrouping[i] = mappedGroup[improperGrouping[i]];
	}
	
	return properGrouping;
}

bool Grouping::isMatch(std::vector<size_t> grouping) const {
	for(size_t i = 0; i < numSpecies; i++) {
		if(grouping[i] != groups[i]) return false;
	}
	return true;
}

bool Grouping::operator==(const Grouping& g) const {
	if(numSpecies != g.numSpecies) return false;
	for(size_t i = 0; i < numSpecies; i++) {
		if(groups[i] != g.groups[i]) return false;
	}
	return true;
}

bool Grouping::operator!=(const Grouping& g) const {
	return !(*this == g);
}
