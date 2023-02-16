#include <stdio.h>
#include <unistd.h> //Gives getopt
#include "input.hpp"
#include "rjsolver.hpp"

#define DEFAULT_JUMPS_PER_DIAL 50000
#define DEFAULT_NUM_DIALS 10
#define DEFAULT_BURN_IN 100000
#define DEFAULT_NUM_STEPS 100000
#define DEFAULT_NUM_CHAINS 10
#define DEFAULT_THINNING_FACTOR 1

#define RANDOM_SEED 42

int main(int argc, char **argv) {
	Input input(argc, argv, true,
		{'a', 'g'},
		{'j', 'd', 'b', 's', 'c', 't'},
		{DEFAULT_JUMPS_PER_DIAL, DEFAULT_NUM_DIALS, DEFAULT_BURN_IN, DEFAULT_NUM_STEPS, DEFAULT_NUM_CHAINS, DEFAULT_THINNING_FACTOR}
	);
	
	Grouping rowGrouping(input.getData().getNumRowSpecies());
	Grouping colGrouping(input.getData().getNumColSpecies());
	rowGrouping.separate();
	colGrouping.separate();
	
	Hyperprior hyperprior = input.getBoolOptResult('a') ?
		Hyperprior(new Hyperpriors::AIC()) :
		Hyperprior(new Hyperpriors::Flat());
	auto parametersPrior = input.getPriors<ReversibleJumpSolver::NUM_ADDITIONAL_PARAMETERS>();
	
	bool isGrowthVarying = input.getBoolOptResult('g');
	
	size_t jumpsPerDial = input.getIntOptResult('j');
	size_t numDials = input.getIntOptResult('d');
	size_t burnIn = input.getIntOptResult('b');
	size_t numSteps = input.getIntOptResult('s');
	size_t numChains = input.getIntOptResult('c');
	size_t thinningFactor = input.getIntOptResult('t');
	
	size_t numOutputsPerChain = numSteps / thinningFactor;
	size_t numOutputs = numOutputsPerChain * numChains;
	
	Model model = Model(new Models::LotkaVolterra());
	ReversibleJumpSolver masterSolver(model, input.getData(), hyperprior, parametersPrior, {rowGrouping, rowGrouping, colGrouping}, {isGrowthVarying, true, true});
	
	//Only groupings need actual default values; the rest have default constructors.
	OutputColumn<Grouping> outputGrowthGroupings("growth_group", numOutputs, masterSolver.getGrouping(GROWTH));
	OutputColumn<Grouping> outputRowGroupings("row_group", numOutputs, masterSolver.getGrouping(ROW));
	OutputColumn<Grouping> outputColGroupings("col_group", numOutputs, masterSolver.getGrouping(COL));
	OutputColumn<Parameters> outputParameters("parameters", numOutputs, Parameters());
	OutputColumn<double> outputErrorVariance("parameters_error_variance", numOutputs, 0.0);
	OutputColumn<size_t> outputChainID("chain_id", numOutputs, 0);
	
	masterSolver.setSeed(RANDOM_SEED);
	masterSolver.dialIn(jumpsPerDial, numDials);
	masterSolver.resetChain();
	
	for(size_t chain = 0; chain < numChains; chain++) {
		ReversibleJumpSolver solver = masterSolver;
		solver.setSeed(RANDOM_SEED + chain + 1); //+1 so that chain 0 doesn't use the same seed as dialing in.
		solver.burnIn(burnIn);
		
		size_t outputIndex = numOutputsPerChain * chain;
		for(size_t i = 0; i < numSteps; i++) {
			solver.makeJump();
			if(i%thinningFactor == 0) {
				outputGrowthGroupings.set(outputIndex, solver.getGrouping(GROWTH));
				outputRowGroupings.set(outputIndex, solver.getGrouping(ROW));
				outputColGroupings.set(outputIndex, solver.getGrouping(COL));
				outputParameters.set(outputIndex, Parameters(solver.getParameters(), solver.getGroupings()));
				outputErrorVariance.set(outputIndex, solver.getParameters().getAdditionalParameter(0));
				outputChainID.set(outputIndex, chain);
				
				outputIndex++;
			}
		}
	}
	
	outputTable(input.getOutputFile(), outputGrowthGroupings, outputRowGroupings, outputColGroupings, outputParameters, outputErrorVariance, outputChainID);
	
	return 0;
}
