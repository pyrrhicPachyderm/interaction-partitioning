#include "input.hpp"
#include "mlsolver.hpp"

int main(int argc, char **argv) {
	Input input(argc, argv, false);
	
	std::string errorDistribution = input.getErrorDistribution();
	std::string additionalParameterName;
	
	MaximumLikelihoodSolverInterface *solver = NULL;
	if(errorDistribution == "normal") {
		solver = new GaussNewtonSolver(input.getModel(), input.getData());
	} else if(errorDistribution == "gamma") {
		solver = new NLoptSolver<Distributions::Gamma2>(input.getModel(), input.getData());
	} else if(errorDistribution == "negativebinomial") {
		solver = new NLoptSolver<Distributions::DiscreteWrapper<Distributions::NegativeBinomial2>>(input.getModel(), input.getData());
	} else {
		fprintf(stderr, "Unrecognised error distribution.\n");
		exit(1);
	}
	
	solver->updateGrouping(GROWTH, &Grouping::separate);
	
	OutputColumn<Grouping> outputRowGroupings("row_group");
	OutputColumn<Grouping> outputColGroupings("col_group");
	OutputColumn<double> outputAICs("aic");
	OutputColumn<double> outputAICcs("aicc");
	OutputColumn<double> outputR2s("R2");
	OutputColumn<Parameters> outputParameters("parameters");
	//TODO: Output additional parameters for generalised models.
	
	do {
		do {
			outputRowGroupings.insert(solver->getGrouping(ROW));
			outputColGroupings.insert(solver->getGrouping(COL));
			outputAICs.insert(solver->getAIC());
			outputAICcs.insert(solver->getAICc());
			outputR2s.insert(solver->getR2());
			outputParameters.insert(Parameters(solver->getSolutionParameters(), solver->getGroupings()));
		} while(solver->updateGrouping(ROW, &Grouping::advance));
	} while(solver->updateGrouping(COL, &Grouping::advance));
	
	delete solver;
	
	outputTable(input.getOutputFile(), outputRowGroupings, outputColGroupings, outputParameters, outputAICs, outputAICcs, outputR2s);
	
	return 0;
}
