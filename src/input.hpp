#ifndef INPUT_HPP
#define INPUT_HPP

#include "io.hpp"
#include "data.hpp"

class Input {
	protected:
		std::vector<char> boolOpts;
		std::vector<bool> boolOptResults = std::vector<bool>(boolOpts.size(), false);
		std::vector<char> intOpts;
		std::vector<size_t> intOptResults;
		
		Data data;
		std::string outputFile;
		std::vector<Distribution<double>> priors;
		
		std::string getOptsString();
		void parseOptInput(char opt, const char *optarg);
	public:
		Input(int argc, char** argv, bool needsPriors, std::vector<char> boolOpts, std::vector<char> intOpts, std::vector<size_t> intOptDefaults);
		Input(int argc, char** argv, bool needsPriors, std::vector<char> boolOpts): Input(argc, argv, needsPriors, boolOpts, {}, {}) {};
		Input(int argc, char** argv, bool needsPriors): Input(argc, argv, needsPriors, {}) {};
		
		const Data &getData() const {return data;};
		const std::string &getOutputFile() const {return outputFile;};
		
		bool getBoolOptResult(char opt) const;
		size_t getIntOptResult(char opt) const;
};

#endif
