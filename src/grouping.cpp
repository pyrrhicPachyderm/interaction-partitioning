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
