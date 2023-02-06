#ifndef INPUT_HPP
#define INPUT_HPP

#include "io.hpp"
#include "data.hpp"

class Input {
	protected:
		std::vector<char> boolOpts;
		std::vector<bool> boolOptResults = std::vector<bool>(boolOpts.size(), false);
		
		std::string optsString;
		
		Data data;
		std::string outputFile;
		
		void setOptsString();
		void parseOptInput(char opt, const char *optarg);
	public:
		Input(int argc, char** argv, std::vector<char> boolOpts);
		Input(int argc, char** argv): Input(argc, argv, {}) {};
		
		const Data &getData() {return data;};
		const std::string &getOutputFile() {return outputFile;};
	public:
		bool getBoolOptResult(char opt) const;
};

#endif
