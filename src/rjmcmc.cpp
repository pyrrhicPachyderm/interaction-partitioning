#include <stdio.h>
#include <unistd.h> //Gives getopt
#include "input.hpp"
#include "rjsolver.hpp"

#define JUMPS_PER_DIAL 1e3
#define NUM_DIALS 1e2
#define BURN_IN 1e5
#define NUM_STEPS 1e5

int main(int argc, char **argv) {
	Input input(argc, argv, {'a'});
	
	Grouping rowGrouping(input.getData().getNumRowSpecies());
	Grouping colGrouping(input.getData().getNumColSpecies());
	rowGrouping.separate();
	colGrouping.separate();
	
	Hyperprior hyperprior = input.getBoolOptResult('a') ?
		Hyperprior::aic() :
		Hyperprior::flat();
	
	ReversibleJumpSolver solver(input.getData(), hyperprior, {rowGrouping, rowGrouping, colGrouping}, {false, true, true});
	
	solver.dialIn(JUMPS_PER_DIAL, NUM_DIALS);
	solver.burnIn(BURN_IN);
	
	OutputColumn<Grouping> outputRowGroupings("row_group");
	OutputColumn<Grouping> outputColGroupings("col_group");
	OutputColumn<Parameters> outputParameters("parameters");
	OutputColumn<double> outputErrorVariance("parameters_error_variance");
	for(size_t i = 0; i < NUM_STEPS; i++) {
		solver.makeJump();
		outputRowGroupings.insert(solver.getGrouping(ROW));
		outputColGroupings.insert(solver.getGrouping(COL));
		outputParameters.insert(solver.getParameters());
		outputErrorVariance.insert(solver.getParameters().getAdditionalParameter(0));
	}
	
	outputTable(input.getOutputFile(), outputRowGroupings, outputColGroupings, outputParameters, outputErrorVariance);
	
	return 0;
}
