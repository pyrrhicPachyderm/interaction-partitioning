#include <stdio.h>
#include <unistd.h> //Gives getopt
#include "input.hpp"

#define NUM_MANDATORY_ARGS 4

const static std::string DEFAULT_OPTS_STRING = "p";

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
		"\t\tOtherwise, the response is assumed to be a measure of the total change in population size from one generation to the next.\n"
		,
		argv[0]
	);
	//TODO: Take boolOpts and some sort of boolOptsHelpStrings as an argument, and include those as well.
}

Input::Input(int argc, char** argv, std::vector<char> boolOpts):
	boolOpts(boolOpts)
{
	setOptsString();
	
	bool isPerCapita = false;
	
	//Parse the options. getopt should shuffle all the mandatory arguments to the end by itself.
	int opt;
	while((opt = getopt(argc, argv, optsString.c_str())) != -1) { //No initial colon means leaving on automatic error reporting.
		switch(opt) {
			case 'p':
				isPerCapita = true;
				break;
			case '?': //This covers an explicit ? as well as other errors.
				//TODO: Some additional error parsing. Perhaps see https://stackoverflow.com/a/44371579
				printUsage(argc, (const char**)argv);
				exit(1);
			default:
				parseOptInput(opt, optarg);
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
	outputFile = argv[optind+3];
	
	data = Data(
		focal,
		response,
		design,
		isPerCapita
	);
}

void Input::setOptsString() {
	optsString = DEFAULT_OPTS_STRING;
	
	for(char opt : boolOpts) {
		optsString.push_back(opt);
	}
}

//Returns the index of the first occurrence of opt in opts.
//Horrible algorithmically, but we don't expect more than a handful of opts and a handful of calls.
//Returns -1 on an error.
static int getOptIndex(std::vector<char> opts, char opt) {
	for(size_t i = 0; i < opts.size(); i++) {
		if(opts[i] == opt) return i;
	}
	return -1;
}

bool Input::getBoolOptResult(char opt) const {
	int boolIndex = getOptIndex(boolOpts, opt);
	assert(boolIndex >= 0);
	return boolOptResults[boolIndex];
}

void Input::parseOptInput(char opt, const char *optarg) {
	int boolIndex = getOptIndex(boolOpts, opt);
	if(boolIndex >= 0) boolOptResults[boolIndex] = true;
}
