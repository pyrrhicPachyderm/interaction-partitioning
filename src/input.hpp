#ifndef INPUT_HPP
#define INPUT_HPP

#include "io.hpp"
#include "data.hpp"

class Input {
	protected:
		Data data;
		std::string outputFile;
	public:
		Input(Data data, std::string outputFile): data(data), outputFile(outputFile) {};
		
		const Data &getData() {return data;};
		const std::string &getOutputFile() {return outputFile;};
};

extern Input readInput(int argc, char** argv);

#endif
