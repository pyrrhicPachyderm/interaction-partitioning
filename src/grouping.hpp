#ifndef GROUPING_HPP
#define GROUPING_HPP

#include <stddef.h>
#include <vector>
#include <array>

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
		//Always returns true.
		bool reset();
		
		//Resets to the lexigraphically last grouping.
		//Always returns true.
		bool separate();
		
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
		
		//A function to fix a grouping that is not a proper rhyming scheme.
		void fix();
	public:
		//If called without a provided group, put every species in the same group.
		Grouping(size_t numSpecies): numSpecies(numSpecies) {
			reset();
		};
		Grouping(const Grouping &g) = default;
		
		const std::vector<size_t> &getGroups() const {return groups;}
		size_t getGroup(size_t species) const {return groups[species];}
		
		size_t getNumGroups() const;
		
		//Returns the number of elements in each group.
		std::vector<size_t> getGroupSizes() const;
		
		//Returns the maximum number of merges/splits that come from a grouping with the given parameters.
		//Used for finding a constant when setting up an RJSolver.
		//For getNumMerges, this isn't just the max; it's always the same number.
		static size_t getNumMerges(size_t numGroups);
		static size_t getMaxNumSplits(size_t numSpecies, size_t numGroups);
		
		//Returns the actual number of merges/splits from a given grouping.
		size_t getNumMerges() const;
		size_t getNumSplits() const;
		
		//Gets the indexth merge/split.
		//The order is arbitrary, but consistent.
		//index must be at less than getNumMerges()/getNumSplits().
		Grouping getMerge(size_t index) const;
		Grouping getSplit(size_t index) const;
		
		//Operator overloading.
		void operator=(const Grouping& g);
		bool operator==(const Grouping& g) const;
		bool operator!=(const Grouping& g) const;
};

enum GroupingType {GROWTH, ROW, COL, NUM_GROUPING_TYPES};

typedef std::array<Grouping, NUM_GROUPING_TYPES> GroupingSet;
typedef std::array<size_t, NUM_GROUPING_TYPES> GroupingSizeSet;
//A GroupingSizeSet is the number of groups in each grouping of a GroupingSet.

extern GroupingSizeSet getGroupingSizeSet(const GroupingSet &groupingSet);

#endif
