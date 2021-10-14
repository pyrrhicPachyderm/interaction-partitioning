#include <stdio.h>
#include <unistd.h> //Gives getopt
#include "input.hpp"

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

Input readInput(int argc, char** argv) {
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
				exit(1);
			default: //This should never happen.
				exit(1);
		}
	}
	
	//Get the mandatory arguments.
	int remainingArgs = argc - optind;
	if(remainingArgs != NUM_MANDATORY_ARGS) {
		fprintf(stderr, "Wrong number of mandatory arguments\n");
		printUsage(argc, (const char**)argv);
		exit(1);
	}
	std::vector<size_t> focal = readIndexVector(argv[optind]);
	Eigen::VectorXd response = readDoubleVector(argv[optind+1]);
	Eigen::MatrixXd design = readDoubleMatrix(argv[optind+2]);
	const char *outputFilePtr = argv[optind+3];
	
	Data data(
		focal,
		response,
		design,
		isPerCapita
	);
	
	return Input(data, outputFilePtr);
}
