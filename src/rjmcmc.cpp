#include <stdio.h>
#include <unistd.h> //Gives getopt
#include "input.hpp"
#include "rjsolver.hpp"

int main(int argc, char **argv) {
	Input input = readInput(argc, argv);
	
	Grouping grouping(input.getData().numSpecies);
	grouping.separate();
	
	ReversibleJumpSolver solver(input.getData(), {grouping, grouping, grouping}, {false, true, true});
	
	return 0;
}
