#include <stdio.h>
#include <unistd.h> //Gives getopt
#include "input.hpp"

#define NUM_MANDATORY_ARGS 4

const static std::string DEFAULT_OPTS_STRING = "p";

//Prints the usage spiel to stderr.
static void printUsage(int argc, const char **argv) {
	//Making of use of concatenated string literals in the following:
	fprintf(stderr, "Usage: %s OUTPUT_FILE FOCAL_VECTOR_FILE RESPONSE_VECTOR_FILE DESIGN_MATRIX_FILE [options]\n\n"
		"\tThe focal vector is an integer vector giving the index of the focal species for each observation, 0-indexed.\n"
		"\tThe response vector is a numeric vector giving the growth, fecundity, or other response variable for each observation.\n"
		"\tThe design matrix is a numeric matrix with one row per observation and one column per species, giving the design densities.\n"
		"\tVectors should be written as column vectors.\n"
		"\tA file may be given as -, to use stdin/stdout.\n"
		"\n"
		"Options:\n"
		
		"\t-?\n"
		"\t\tShow this help message and exit.\n"
		,
		argv[0]
	);
	//TODO: Take boolOpts and some sort of boolOptsHelpStrings as an argument, and include those as well.
}

Input::Input(int argc, char** argv, std::vector<char> boolOpts, std::vector<char> intOpts, std::vector<size_t> intOptDefaults):
	boolOpts(boolOpts),
	intOpts(intOpts),
	intOptResults(intOptDefaults)
{
	assert(intOpts.size() == intOptResults.size());
	
	std::string optsString = getOptsString();
	
	//Parse the options. getopt should shuffle all the mandatory arguments to the end by itself.
	int opt;
	while((opt = getopt(argc, argv, optsString.c_str())) != -1) { //No initial colon means leaving on automatic error reporting.
		switch(opt) {
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
	outputFile = argv[optind];
	std::vector<size_t> focal = readIndexVector(argv[optind+1]);
	Eigen::VectorXd response = readDoubleVector(argv[optind+2]);
	Eigen::MatrixXd design = readDoubleMatrix(argv[optind+3]);
	
	data = Data(
		focal,
		response,
		design
	);
}

std::string Input::getOptsString() {
	std::string optsString = DEFAULT_OPTS_STRING;
	
	for(char opt : boolOpts) {
		optsString.push_back(opt);
	}
	for(char opt : intOpts) {
		optsString.push_back(opt);
		optsString.push_back(':');
	}
	
	return optsString;
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
	int index = getOptIndex(boolOpts, opt);
	assert(index >= 0);
	return boolOptResults[index];
}

size_t Input::getIntOptResult(char opt) const {
	int index = getOptIndex(intOpts, opt);
	assert(index >= 0);
	return intOptResults[index];
}

void Input::parseOptInput(char opt, const char *optarg) {
	int index = -1;
	index = getOptIndex(boolOpts, opt);
	if(index >= 0) boolOptResults[index] = true;
	index = getOptIndex(intOpts, opt);
	if(index >= 0) intOptResults[index] = strtoull(optarg, NULL, 0); //TODO: Error handling for strtoull.
}
