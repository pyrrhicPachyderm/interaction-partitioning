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
		
		//The vector maxGroup consists of one integer per species.
		//It consists of a cumulative maximum of the group vector.
		//It is used in finding the next group.
		std::vector<size_t> maxGroup;
		
		//Resets to the lexigraphically first grouping.
		void reset() {
			group = std::vector<size_t>(numSpecies, 0);
			maxGroup = std::vector<size_t>(numSpecies, 0);
		}
		
		//Resets to the lexigraphically last grouping.
		void separate() {
			for(size_t i = 0; i < numSpecies; i++) {
				group[i] = i;
				maxGroup[i] = i;
			}
		}
		
		//A recursive function to be used in advancing to the lexigraphically nex grouping.
		//Increments the specified index by one, if it is valid to do, and resets all species after it to group 0.
		//Recurses on the index to the left if this one cannot be validly incremented.
		//Calls reset() if this is the lexigraphically final grouping.
		//Returns false if and only if it reset, facilitating a do while loop.
		bool advanceIndex(size_t index) {
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
	public:
		//If called without a provided group, put every species in the same group.
		Grouping(size_t numSpecies): numSpecies(numSpecies) {
			reset();
		};
		
		size_t getGroup(size_t species) const {
			return group[species];
		}
		
		size_t getNumGroups() const {
			//The number of groups is simply the highest group number seen before or including the final element.
			//Plus one, as the groups are zero indexed.
			return maxGroup[numSpecies - 1] + 1;
		}
		
		//Advances to the next grouping, per lexigraphic order.
		//Wraps back to the lexigraphically first grouping once it has been through them all.
		//Returns false if and only if it wrapped, facilitating a do while loop.
		bool advanceGrouping() {
			return advanceIndex(numSpecies - 1);
		}
};

#endif
