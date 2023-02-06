#include <stdio.h>
#include <unistd.h> //Gives getopt
#include "input.hpp"
#include "rjsolver.hpp"

#define DEFAULT_JUMPS_PER_DIAL 1000
#define DEFAULT_NUM_DIALS 100
#define DEFAULT_BURN_IN 100000
#define DEFAULT_NUM_STEPS 100000
#define DEFAULT_NUM_CHAINS 10

int main(int argc, char **argv) {
	Input input(argc, argv,
		{'a', 'g'},
		{'j', 'd', 'b', 's', 'c'},
		{DEFAULT_JUMPS_PER_DIAL, DEFAULT_NUM_DIALS, DEFAULT_BURN_IN, DEFAULT_NUM_STEPS, DEFAULT_NUM_CHAINS}
	);
	
	Grouping rowGrouping(input.getData().getNumRowSpecies());
	Grouping colGrouping(input.getData().getNumColSpecies());
	rowGrouping.separate();
	colGrouping.separate();
	
	Hyperprior hyperprior = input.getBoolOptResult('a') ?
		Hyperprior::aic() :
		Hyperprior::flat();
	
	bool isGrowthVarying = input.getBoolOptResult('g');
	
	size_t jumpsPerDial = input.getIntOptResult('j');
	size_t numDials = input.getIntOptResult('d');
	size_t burnIn = input.getIntOptResult('b');
	size_t numSteps = input.getIntOptResult('s');
	size_t numChains = input.getIntOptResult('c');
	
	ReversibleJumpSolver solver(input.getData(), hyperprior, {rowGrouping, rowGrouping, colGrouping}, {isGrowthVarying, true, true});
	
	OutputColumn<Grouping> outputGrowthGroupings("growth_group");
	OutputColumn<Grouping> outputRowGroupings("row_group");
	OutputColumn<Grouping> outputColGroupings("col_group");
	OutputColumn<Parameters> outputParameters("parameters");
	OutputColumn<double> outputErrorVariance("parameters_error_variance");
	OutputColumn<size_t> outputChainID("chain_id");
	
	solver.dialIn(jumpsPerDial, numDials);
	
	for(size_t chain = 0; chain < numChains; chain++) {
		solver.burnIn(burnIn);
		for(size_t i = 0; i < numSteps; i++) {
			solver.makeJump();
			outputGrowthGroupings.insert(solver.getGrouping(GROWTH));
			outputRowGroupings.insert(solver.getGrouping(ROW));
			outputColGroupings.insert(solver.getGrouping(COL));
			outputParameters.insert(Parameters(solver.getParameters(), solver.getGroupings()));
			outputErrorVariance.insert(solver.getParameters().getAdditionalParameter(0));
			outputChainID.insert(chain);
		}
	}
	
	outputTable(input.getOutputFile(), outputGrowthGroupings, outputRowGroupings, outputColGroupings, outputParameters, outputErrorVariance, outputChainID);
	
	return 0;
}
