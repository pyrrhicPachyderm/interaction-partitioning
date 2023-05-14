#include <stdio.h>
#include <unistd.h> //Gives getopt
#include <omp.h>
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
	
	GroupingSet groupings = {rowGrouping, rowGrouping, colGrouping};
	
	bool isGroupingGrowth = input.getBoolOptResult('g');
	std::array<bool, NUM_GROUPING_TYPES> isGrouping = {isGroupingGrowth, true, true};
	
	Hyperprior hyperprior = input.getBoolOptResult('a') ?
		Hyperprior(new Hyperpriors::AIC()) :
		Hyperprior(new Hyperpriors::Flat());
	
	std::string errorDistribution = input.getErrorDistribution();
	std::string additionalParameterName;
	
	ReversibleJumpSolverInterface *masterSolver = NULL;
	if(errorDistribution == "normal") {
		auto parametersPrior = input.getPriors<ReversibleJumpSolver<Distributions::Normal>::NUM_ADDITIONAL_PARAMETERS>();
		masterSolver = new ReversibleJumpSolver<Distributions::Normal>(input.getModel(), input.getData(), hyperprior, parametersPrior, groupings, isGrouping);
		additionalParameterName = "error_variance";
	} else if(errorDistribution == "gamma") {
		auto parametersPrior = input.getPriors<ReversibleJumpSolver<Distributions::Gamma2>::NUM_ADDITIONAL_PARAMETERS>();
		masterSolver = new ReversibleJumpSolver<Distributions::Gamma2>(input.getModel(), input.getData(), hyperprior, parametersPrior, groupings, isGrouping);
		additionalParameterName = "dispersion";
	} else if(errorDistribution == "negativebinomial") {
		auto parametersPrior = input.getPriors<ReversibleJumpSolver<Distributions::DiscreteWrapper<Distributions::NegativeBinomial2>>::NUM_ADDITIONAL_PARAMETERS>();
		masterSolver = new ReversibleJumpSolver<Distributions::DiscreteWrapper<Distributions::NegativeBinomial2>>(input.getModel(), input.getData(), hyperprior, parametersPrior, groupings, isGrouping);
		additionalParameterName = "dispersion";
	} else {
		fprintf(stderr, "Unrecognised error distribution.\n");
		exit(1);
	}
	
	size_t jumpsPerDial = input.getIntOptResult('j');
	size_t numDials = input.getIntOptResult('d');
	size_t burnIn = input.getIntOptResult('b');
	size_t numSteps = input.getIntOptResult('s');
	size_t numChains = input.getIntOptResult('c');
	size_t thinningFactor = input.getIntOptResult('t');
	
	size_t numOutputsPerChain = numSteps / thinningFactor;
	size_t numOutputs = numOutputsPerChain * numChains;
	
	//Only groupings need actual default values; the rest have default constructors.
	OutputColumn<Grouping> outputGrowthGroupings("growth_group", numOutputs, masterSolver->getGrouping(GROWTH));
	OutputColumn<Grouping> outputRowGroupings("row_group", numOutputs, masterSolver->getGrouping(ROW));
	OutputColumn<Grouping> outputColGroupings("col_group", numOutputs, masterSolver->getGrouping(COL));
	OutputColumn<Parameters> outputParameters("parameters", numOutputs, Parameters());
	OutputColumn<double> outputAdditionalParameter("parameters_" + additionalParameterName, numOutputs, 0.0);
	OutputColumn<size_t> outputChainID("chain_id", numOutputs, 0);
	
	masterSolver->setSeed(RANDOM_SEED);
	masterSolver->dialIn(jumpsPerDial, numDials);
	masterSolver->resetChain();
	
	#pragma omp parallel for
	for(size_t chain = 0; chain < numChains; chain++) {
		ReversibleJumpSolverInterface *solver = masterSolver->getCopy();
		solver->setSeed(RANDOM_SEED + chain + 1); //+1 so that chain 0 doesn't use the same seed as dialing in.
		solver->burnIn(burnIn);
		
		size_t outputIndex = numOutputsPerChain * chain;
		for(size_t i = 0; i < numSteps; i++) {
			solver->makeJump();
			if(i%thinningFactor == 0) {
				outputGrowthGroupings.set(outputIndex, solver->getGrouping(GROWTH));
				outputRowGroupings.set(outputIndex, solver->getGrouping(ROW));
				outputColGroupings.set(outputIndex, solver->getGrouping(COL));
				outputParameters.set(outputIndex, Parameters(solver->getParameters(), solver->getGroupings()));
				outputAdditionalParameter.set(outputIndex, solver->getAdditionalParameter(0));
				outputChainID.set(outputIndex, chain);
				
				outputIndex++;
			}
		}
		delete solver;
	}
	
	delete masterSolver;
	
	outputTable(input.getOutputFile(), outputGrowthGroupings, outputRowGroupings, outputColGroupings, outputParameters, outputAdditionalParameter, outputChainID);
	
	return 0;
}
