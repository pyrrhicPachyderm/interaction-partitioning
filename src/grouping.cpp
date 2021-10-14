#include <algorithm>
#include "grouping.hpp"

void Grouping::reset() {
	groups = std::vector<size_t>(numSpecies, 0);
	maxGroups = std::vector<size_t>(numSpecies, 0);
}

void Grouping::separate() {
	for(size_t i = 0; i < numSpecies; i++) {
		groups[i] = i;
		maxGroups[i] = i;
	}
}

bool Grouping::advance() {
	return advanceIndex(numSpecies - 1);
}

bool Grouping::advanceIndex(size_t index) {
	//It is never valid to increment the first element, so reset.
	if(index == 0) {
		reset();
		return false;
	}
	
	//If we have seen an equal or higher group number to the left, this can safely be incremented.
	//We must reset all group numbers to the right of it to 0
	if(groups[index] <= maxGroups[index-1]) {
		groups[index]++;
		maxGroups[index] = std::max(maxGroups[index], groups[index]);
		for(size_t i = index+1; i < numSpecies; i++) {
			groups[i] = 0;
			maxGroups[i] = maxGroups[index];
		}
		return true;
	}
	//Otherwise, it is not valid to increment it, so we recurse.
	return advanceIndex(index-1);
}

size_t Grouping::getNumGroups() const {
	//The number of groups is simply the highest group number seen before or including the final element.
	//Plus one, as the groups are zero indexed.
	return maxGroups[numSpecies - 1] + 1;
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

bool Grouping::isMatch(std::vector<size_t> grouping) {
	for(size_t i = 0; i < numSpecies; i++) {
		if(grouping[i] != groups[i]) return false;
	}
	return true;
}
