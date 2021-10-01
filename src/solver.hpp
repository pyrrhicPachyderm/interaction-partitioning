#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "grouping.hpp"
#include "data.hpp"

class Solver {
	protected:
		Data data;
		Grouping growthGrouping;
		Grouping rowGrouping;
		Grouping colGrouping;
	public:
		Solver(Data data, Grouping growthGrouping, Grouping rowGrouping, Grouping colGrouping):
			data(data), growthGrouping(growthGrouping), rowGrouping(rowGrouping), colGrouping(colGrouping)
		{
			assert(growthGrouping.numSpecies == data.numSpecies);
			assert(rowGrouping.numSpecies == data.numSpecies);
			assert(colGrouping.numSpecies == data.numSpecies);
		}
};

#endif
