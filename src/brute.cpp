#include "input.hpp"
#include "mlsolver.hpp"

int main(int argc, char **argv) {
	Input input(argc, argv, false, {'n'}, {}, {}, {'f', 'x'}, {0, 1e-3});
	
	bool isNull = input.getBoolOptResult('n');
	double ftol = input.getDoubleOptResult('f');
	double xtol = input.getDoubleOptResult('x');
	std::string errorDistribution = input.getErrorDistribution();
	std::string additionalParameterName;
	
	MaximumLikelihoodSolverInterface *masterSolver = NULL;
	if(errorDistribution == "normal") {
		masterSolver = new GaussNewtonSolver(input.getModel(), input.getData());
	} else if(errorDistribution == "gamma") {
		masterSolver = new NLoptSolver<Distributions::Gamma2>(input.getModel(), input.getData(), ftol, xtol);
	} else if(errorDistribution == "negativebinomial") {
		masterSolver = new NLoptSolver<Distributions::DiscreteWrapper<Distributions::NegativeBinomial2>>(input.getModel(), input.getData(), ftol, xtol);
	} else {
		fprintf(stderr, "Unrecognised error distribution.\n");
		exit(1);
	}
	
	//Build a vector of all the GroupingSets to evaluate.
	GroupingSet groupingSet = masterSolver->getGroupings();
	groupingSet[GROWTH].separate();
	std::vector<GroupingSet> groupingSets;
	if(isNull) {
		groupingSets.push_back(groupingSet);
	} else {
		do {
			do {
				groupingSets.push_back(groupingSet);
			} while(groupingSet[ROW].advance());
		} while(groupingSet[COL].advance());
	}
	
	//Only groupings need actual default values; the rest have default constructors.
	OutputColumn<Grouping> outputRowGroupings("row_group", groupingSets.size(), masterSolver->getGrouping(ROW));
	OutputColumn<Grouping> outputColGroupings("col_group", groupingSets.size(), masterSolver->getGrouping(COL));
	OutputColumn<double> outputAICs("aic", groupingSets.size(), 0.0);
	OutputColumn<double> outputAICcs("aicc", groupingSets.size(), 0.0);
	OutputColumn<double> outputR2s("R2", groupingSets.size(), 0.0);
	OutputColumn<Parameters> outputParameters("parameters", groupingSets.size(), Parameters());
	//TODO: Output additional parameters for generalised models.
	
	#pragma omp parallel for schedule(dynamic)
	for(size_t i = 0; i < groupingSets.size(); i++) {
		MaximumLikelihoodSolverInterface *solver = masterSolver->getCopy();
		solver->setGroupings(groupingSets[i]);
		
		outputRowGroupings.set(i, solver->getGrouping(ROW));
		outputColGroupings.set(i, solver->getGrouping(COL));
		outputAICs.set(i, solver->getAIC(isNull));
		outputAICcs.set(i, solver->getAICc(isNull));
		outputR2s.set(i, solver->getR2(isNull));
		outputParameters.set(i, Parameters(solver->getSolutionParameters(isNull), solver->getGroupings()));
		
		delete solver;
	}
	
	delete masterSolver;
	
	outputTable(input.getOutputFile(), outputRowGroupings, outputColGroupings, outputParameters, outputAICs, outputAICcs, outputR2s);
	
	return 0;
}
