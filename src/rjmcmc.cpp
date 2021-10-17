#include <stdio.h>
#include <unistd.h> //Gives getopt
#include "input.hpp"
#include "rjsolver.hpp"

#define BURN_IN 1e5
#define NUM_STEPS 1e5

int main(int argc, char **argv) {
	Input input = readInput(argc, argv);
	
	Grouping grouping(input.getData().numSpecies);
	grouping.separate();
	
	ReversibleJumpSolver solver(input.getData(), {grouping, grouping, grouping}, {false, true, true});
	
	solver.burnIn(BURN_IN);
	
	OutputColumn<Grouping> outputRowGroupings("row_group");
	OutputColumn<Grouping> outputColGroupings("col_group");
	for(size_t i = 0; i < NUM_STEPS; i++) {
		solver.makeJump();
		outputRowGroupings.insert(solver.getGrouping(ROW));
		outputColGroupings.insert(solver.getGrouping(COL));
	}
	
	outputTable(input.getOutputFile(), outputRowGroupings, outputColGroupings);
	
	return 0;
}
