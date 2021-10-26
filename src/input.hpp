#ifndef INPUT_HPP
#define INPUT_HPP

#include "io.hpp"
#include "data.hpp"

class Input {
	protected:
		Data data;
		std::string outputFile;
	public:
		Input(int argc, char** argv);
		
		const Data &getData() {return data;};
		const std::string &getOutputFile() {return outputFile;};
};

#endif
