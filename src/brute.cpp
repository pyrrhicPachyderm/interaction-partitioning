#include <stdio.h>
#include <unistd.h> //Gives getopt
#include "io.hpp"
#include "solver.hpp"

#define NUM_MANDATORY_ARGS 4

//Prints the usage spiel to stderr.
static void printUsage(int argc, const char **argv) {
	//Making of use of concatenated string literals in the following:
	fprintf(stderr, "Usage: %s FOCAL_VECTOR_FILE RESPONSE_VECTOR_FILE DESIGN_MATRIX_FILE OUTPUT_FILE [options]\n\n"
		"\tThe focal vector is an integer vector giving the index of the focal species for each observation, 0-indexed.\n"
		"\tThe response vector is a numeric vector giving the growth, fecundity, or other response variable for each observation.\n"
		"\tThe design matrix is a numeric matrix with one row per observation and one column per species, giving the design densities.\n"
		"\tVectors should be written as column vectors.\n"
		"\tA file may be given as -, to use stdin/stdout.\n"
		"\n"
		"Options:\n"
		
		"\t-?\n"
		"\t\tShow this help message and exit.\n"
		
		"\t-p\n"
		"\t\tResponse variables are per capita growth/fecundity.\n"
		"\t\tOtherwise, the response is assumed to be a measure of the total size of the next generation.\n"
		,
		argv[0]
	);
}

int main(int argc, char **argv) {
	bool isPerCapita = false;
	
	//Parse the options. getopt should shuffle all the mandatory arguments to the end by itself.
	int opt;
	while((opt = getopt(argc, argv, "p")) != -1) { //No initial colon means leaving on automatic error reporting.
		switch(opt) {
			case 'p':
				isPerCapita = true;
				break;
			case '?': //This covers an explicit ? as well as other errors.
				//TODO: Some additional error parsing. Perhaps see https://stackoverflow.com/a/44371579
				printUsage(argc, (const char**)argv);
				return 1;
			default: //This should never happen.
				return 1;
		}
	}
	
	//Get the mandatory arguments.
	int remainingArgs = argc - optind;
	if(remainingArgs != NUM_MANDATORY_ARGS) {
		fprintf(stderr, "Wrong number of mandatory arguments\n");
		printUsage(argc, (const char**)argv);
		return 1;
	}
	std::vector<size_t> focal = readIndexVector(argv[optind]);
	Eigen::VectorXd response = readDoubleVector(argv[optind+1]);
	Eigen::MatrixXd design = readDoubleMatrix(argv[optind+2]);
	const char *outputFile = strdupa(argv[optind+3]);
	
	Data data = Data(
		focal,
		response,
		design,
		isPerCapita
	);
	
	Solver solver = Solver(data);
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
	
	outputTable(outputFile, outputRowGroupings, outputColGroupings, outputAICs, outputR2s);
	
	return 0;
}
