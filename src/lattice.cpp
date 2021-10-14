#include "lattice.hpp"

void GroupingLattice::constructAllGroupings() {
	Grouping grouping(numSpecies);
	do {
		groupings.push_back(grouping);
	} while(grouping.advance());
}
