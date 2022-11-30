#ifndef GROUPING_HPP
#define GROUPING_HPP

#include <stddef.h>
#include <vector>

class Grouping {
	public:
		const size_t numSpecies;
	protected:
		//The vector group consists of one integer per species, indicating which group it belongs to.
		//It must form a valid rhyming scheme, 0-indexed.
		//That is, the first integer must be 0, and subsequent integers must be no more than 1 higher than the highest integer seen yet.
		std::vector<size_t> groups;
	public:
		//Resets to the lexigraphically first grouping.
		void reset();
		
		//Resets to the lexigraphically last grouping.
		void separate();
		
		//Advances to the next grouping, per lexigraphic order.
		//Wraps back to the lexigraphically first grouping once it has been through them all.
		//Returns false if and only if it wrapped, facilitating a do while loop.
		bool advance();
	protected:
		//A helper function used by advance().
		//Returns the cumulative maximum of the group vector.
		std::vector<size_t> getMaxGroups() const;
		
		//A recursive function to be used in advancing to the lexigraphically next grouping.
		//Increments the specified index by one, if it is valid to do, and resets all species after it to group 0.
		//Recurses on the index to the left if this one cannot be validly incremented.
		//Calls reset() if this is the lexigraphically final grouping.
		//Returns false if and only if it reset, facilitating a do while loop.
		bool advanceIndex(const std::vector<size_t> &maxGroups, size_t index);
	public:
		//If called without a provided group, put every species in the same group.
		Grouping(size_t numSpecies): numSpecies(numSpecies) {
			reset();
		};
		
		const std::vector<size_t> &getGroups() const {return groups;}
		size_t getGroup(size_t species) const {return groups[species];}
		
		size_t getNumGroups() const;
		
		//Returns the number of elements in each group.
		std::vector<size_t> getGroupSizes() const;
		
		static std::vector<size_t> fixGrouping(std::vector<size_t> improperGrouping);
		
		//Checks whether the provided grouping matches this one.
		bool isMatch(std::vector<size_t> grouping) const;
};

#endif
