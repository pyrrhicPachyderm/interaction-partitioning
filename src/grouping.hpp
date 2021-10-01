#ifndef GROUPING_HPP
#define GROUPING_HPP

#include <vector>

class Grouping {
	public:
		const size_t numSpecies;
	protected:
		//The vector group consists of one integer per species, indicating which group it belongs to.
		//It must form a valid rhyming scheme, 0-indexed.
		//That is, the first integer must be 0, and subsequent integers must be no more than 1 higher than the highest integer seen yet.
		std::vector<size_t> group;
	public:
		//If called without a provided group, put every species in the same group.
		Grouping(size_t numSpecies): numSpecies(numSpecies), group(std::vector<size_t>(numSpecies, 0)) {};
};

#endif
