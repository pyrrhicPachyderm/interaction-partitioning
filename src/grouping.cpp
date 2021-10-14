#include "grouping.hpp"

void Grouping::reset() {
	group = std::vector<size_t>(numSpecies, 0);
	maxGroup = std::vector<size_t>(numSpecies, 0);
}

void Grouping::separate() {
	for(size_t i = 0; i < numSpecies; i++) {
		group[i] = i;
		maxGroup[i] = i;
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
	if(group[index] <= maxGroup[index-1]) {
		group[index]++;
		maxGroup[index] = std::max(maxGroup[index], group[index]);
		for(size_t i = index+1; i < numSpecies; i++) {
			group[i] = 0;
			maxGroup[i] = maxGroup[index];
		}
		return true;
	}
	//Otherwise, it is not valid to increment it, so we recurse.
	return advanceIndex(index-1);
}

size_t Grouping::getNumGroups() const {
	//The number of groups is simply the highest group number seen before or including the final element.
	//Plus one, as the groups are zero indexed.
	return maxGroup[numSpecies - 1] + 1;
}
