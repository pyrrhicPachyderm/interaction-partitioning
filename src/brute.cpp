#include "input.hpp"
#include "mlsolver.hpp"

int main(int argc, char **argv) {
	Input input = readInput(argc, argv);
	
	MaximumLikelihoodSolver solver = MaximumLikelihoodSolver(input.getData());
	solver.updateGrowthGrouping(&Grouping::separate);
	
	OutputColumn<Grouping> outputRowGroupings("row_group");
	OutputColumn<Grouping> outputColGroupings("col_group");
	OutputColumn<double> outputAICs("aic");
	OutputColumn<double> outputR2s("R2");
	
	do {
		do {
			outputRowGroupings.insert(solver.getRowGrouping());
			outputColGroupings.insert(solver.getColGrouping());
			outputAICs.insert(solver.getAIC());
			outputR2s.insert(solver.getR2());
		} while(solver.updateRowGrouping(&Grouping::advance));
	} while(solver.updateColGrouping(&Grouping::advance));
	
	outputTable(input.getOutputFile(), outputRowGroupings, outputColGroupings, outputAICs, outputR2s);
	
	return 0;
}
