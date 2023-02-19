#include <stdio.h>
#include <unistd.h> //Gives getopt
#include "utils/string.hpp"
#include "utils/array.hpp"
#include "input.hpp"

#define NUM_MANDATORY_ARGS 7
#define NUM_UNAUGMENTED_PRIORS 2

const static std::string DEFAULT_OPTS_STRING = "p";

//Prints the usage spiel to stderr.
static void printUsage(int argc, const char **argv, bool needsPriors) {
	//Making of use of concatenated string literals in the following:
	fprintf(stderr, "Usage:\n"
		"\t%s [options] OUTPUT_FILE MODEL ERROR_DISTRIBUTION indv|pop FOCAL_VECTOR_FILE RESPONSE_VECTOR_FILE DESIGN_MATRIX_FILE%s\n"
		"\t%s [options] OUTPUT_FILE MODEL ERROR_DISTRIBUTION time ID_VECTOR_FILE TIME_VECTOR_FILE DENSITY_MATRIX_FILE%s\n"
		"\n"
		"\tThe model should be, e.g. LotkaVolterra, BevertonHolt. Case insensitive.\n"
		"\tThe error distribution should be, e.g. Normal, Gamma. Case insensitive. Maximum Likelihood methods only support Normal.\n"
		"\tThe argument 'indv', 'pop' or 'time' distinguishes between individual response data, pop size response data, and time series data.\n"
		"\tAll vector and matrix files should be whitespace-separated tables.\n"
		"\tVectors should be written as column vectors.\n"
		"%s"
		"\tA file may be given as -, to use stdin/stdout.\n"
		"\tIndividual response data:\n"
		"\t\tThe focal vector is an integer vector giving the index of the focal species for each observation, 0-indexed.\n"
		"\t\tThe response vector is a numeric vector giving the growth, fecundity, or other response variable for each observation.\n"
		"\t\tThe design matrix is a numeric matrix with one row per observation and one column per species, giving the design densities.\n"
		"\tPopulation size response data:\n"
		"\t\tAs individual response data, except the response vector gives the total final population size of the focal species.\n"
		"\tTime series data:\n"
		"\t\tThe id vector is an integer vector, identifying which experiment each observation was a part of.\n"
		"\t\tThe time vector is a numeric vector giving the time stamp at which densities were measured.\n"
		"\t\tThe density matrix is a numeric matrix with one row per observation and one column per species, giving the measured densities.\n"
		"\n"
		"Options:\n"
		"\t-?\n"
		"\t\tShow this help message and exit.\n"
		,
		argv[0],
		!needsPriors ? "" : " PRIORS_FILE",
		argv[0],
		!needsPriors ? "" : " PRIORS_FILE",
		!needsPriors ? "" :
			"\tA priors file has three lines, giving the priors for the growth rates, competition coefficients, and error variance respectively.\n"
			"\tEach line begins with the name of a distribution (case-insensitive, without spaces), then has a whitespace separated list of parameters.\n"
	);
	//TODO: Take boolOpts and some sort of boolOptsHelpStrings as an argument, and include those as well.
}

static Model matchModelString(std::string modelString) {
	if(modelString == "lotkavolterra") {
		return Model(new Models::LotkaVolterra());
	} else if (modelString == "bevertonholt") {
		return Model(new Models::BevertonHolt());
	} else if (modelString == "ricker") {
		return Model(new Models::Ricker());
	} else {
		fprintf(stderr, "Unrecognised model.\n");
		exit(1);
	}
}

Input::Input(int argc, char** argv, bool needsPriors, std::vector<char> boolOpts, std::vector<char> intOpts, std::vector<size_t> intOptDefaults):
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
				printUsage(argc, (const char**)argv, needsPriors);
				exit(1);
			default:
				parseOptInput(opt, optarg);
		}
	}
	
	//Get the mandatory arguments.
	int numRemainingArgs = argc - optind;
	int numRequiredArgs = NUM_MANDATORY_ARGS + (int)needsPriors; //The priors file is an additional argument.
	if(numRemainingArgs != numRequiredArgs) {
		fprintf(stderr, "Wrong number of mandatory arguments\n");
		printUsage(argc, (const char**)argv, needsPriors);
		exit(1);
	}
	
	outputFile = argv[optind++];
	model = matchModelString(argv[optind++]);
	errorDistribution = argv[optind++];
	strToLowerCase(errorDistribution);
	
	std::string dataType(argv[optind++]);
	if(dataType == "indv") {
		std::vector<size_t> focal = readIndexVector(argv[optind++]);
		Eigen::VectorXd response = readDoubleVector(argv[optind++]);
		Eigen::MatrixXdRowMajor design = readDoubleMatrix(argv[optind++]);
		data = Data(new Datasets::IndividualResponse(
			focal,
			response,
			design
		));
	} else if(dataType == "pop") {
		std::vector<size_t> focal = readIndexVector(argv[optind++]);
		Eigen::VectorXd response = readDoubleVector(argv[optind++]);
		Eigen::MatrixXdRowMajor design = readDoubleMatrix(argv[optind++]);
		data = Data(new Datasets::PopSizeResponse(
			focal,
			response,
			design
		));
	} else if(dataType == "time") {
		std::vector<size_t> id = readIndexVector(argv[optind++]);
		Eigen::VectorXd time = readDoubleVector(argv[optind++]);
		Eigen::MatrixXdRowMajor density = readDoubleMatrix(argv[optind++]);
		data = Data(new Datasets::TimeSeries(
			id,
			time,
			density
		));
	} else {
		fprintf(stderr, "Invalid data type\n");
		printUsage(argc, (const char**)argv, needsPriors);
		exit(1);
	}
	
	if(needsPriors) {
		priors = readDistributionList(argv[optind++]);
	}
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

template<size_t nAug> AugmentedParametersPrior<nAug> Input::getPriors() const {
	if(priors.size() != NUM_UNAUGMENTED_PRIORS + nAug) {
		fprintf(stderr, "Wrong number of priors\n");
		//TODO: Integrate this into the constructor so that it *can* print usage.
		//printUsage(argc, (const char**)argv, needsPriors);
		exit(1);
	}
	
	return AugmentedParametersPrior<nAug>(priors[0], priors[1], array_map(
		[this] (size_t index) -> decltype(priors)::value_type {
			return priors[NUM_UNAUGMENTED_PRIORS + index];
		},
		make_index_array<nAug>()
	));
}

//Explicitly instantiate.
template AugmentedParametersPrior<1> Input::getPriors() const;
