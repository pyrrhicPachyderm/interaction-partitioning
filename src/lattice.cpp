#include "lattice.hpp"

std::vector<Grouping> GroupingLattice::constructAllGroupings(size_t numSpecies) {
	std::vector<Grouping> groupings;
	Grouping grouping(numSpecies);
	
	do {
		groupings.push_back(grouping);
	} while(grouping.advance());
	
	return groupings;
}
