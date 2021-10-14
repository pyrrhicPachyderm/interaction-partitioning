#ifndef INPUT_HPP
#define INPUT_HPP

#include "io.hpp"
#include "data.hpp"

class Input {
	protected:
		Data data;
		char *outputFile;
	public:
		Input(Data data, const char *outputFilePtr): data(data) {
			outputFile = strdup(outputFilePtr);
		};
		
		~Input() {
			free(outputFile);
		}
		
		const Data &getData() {return data;};
		const char *getOutputFile() {return outputFile;};
};

extern Input readInput(int argc, char** argv);

#endif
