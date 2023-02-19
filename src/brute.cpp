#include "input.hpp"
#include "mlsolver.hpp"

int main(int argc, char **argv) {
	Input input(argc, argv, false);
	
	if(input.getErrorDistribution() != "normal") {
		fprintf(stderr, "Maximum likelihood solver does not support non-normal error distributions.\n");
		exit(1);
	}
	
	Model model = Model(new Models::LotkaVolterra());
	MaximumLikelihoodSolver solver = MaximumLikelihoodSolver(model, input.getData());
	solver.updateGrouping(GROWTH, &Grouping::separate);
	
	OutputColumn<Grouping> outputRowGroupings("row_group");
	OutputColumn<Grouping> outputColGroupings("col_group");
	OutputColumn<double> outputAICs("aic");
	OutputColumn<double> outputAICcs("aicc");
	OutputColumn<double> outputR2s("R2");
	OutputColumn<Parameters> outputParameters("parameters");
	
	do {
		do {
			outputRowGroupings.insert(solver.getGrouping(ROW));
			outputColGroupings.insert(solver.getGrouping(COL));
			outputAICs.insert(solver.getAIC());
			outputAICcs.insert(solver.getAICc());
			outputR2s.insert(solver.getR2());
			outputParameters.insert(Parameters(solver.getSolution(), solver.getGroupings()));
		} while(solver.updateGrouping(ROW, &Grouping::advance));
	} while(solver.updateGrouping(COL, &Grouping::advance));
	
	outputTable(input.getOutputFile(), outputRowGroupings, outputColGroupings, outputParameters, outputAICs, outputAICcs, outputR2s);
	
	return 0;
}
