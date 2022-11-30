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
			goupSizes.push_back(1);
		} else {
			groupSizes[groups[i]] += 1;
		}
	}
	
	return groupSizes;
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
