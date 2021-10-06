#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iterator>
#include "grouping.hpp"
#include "io.hpp"

const char *OUTPUT_TABLE_SEPARATOR = "\t";

template<typename StreamT> static StreamT &openFile(const char *filename) {
	StreamT &file = *(new StreamT);
	file.open(filename);
	if(!file) {
		fprintf(stderr, "Could not open file %s\n", filename);
		exit(1);
	}
	return file;
}

std::istream &openInput(const char *filename) {
	if(strcmp(filename, "-") == 0) return std::cin;
	return openFile<std::ifstream>(filename);
}

std::ostream &openOutput(const char *filename) {
	if(strcmp(filename, "-") == 0) return std::cout;
	return openFile<std::ofstream>(filename);
}

//////////////////
//Input functions.
//////////////////

//Set numCols to 0 for a dynamically sized matrix.
//In this case, the number of columns is returned through numCols.
template<typename T> static std::vector<T> readMatrix(const char *filename, size_t *numCols) {
	std::istream &file = openInput(filename);
	
	std::vector<T> result;
	
	std::string lineString;
	while(std::getline(file, lineString)) {
		std::istringstream buffer = std::istringstream(lineString);
		std::vector<T> line = std::vector<T>{std::istream_iterator<T>(buffer), std::istream_iterator<T>()};
		
		if(*numCols != 0 && *numCols != line.size()) {
			fprintf(stderr, "%s has the wrong number of columns, or an inconsistent number of columns\n", filename);
			exit(1);
		} else {
			*numCols = line.size();
			result.insert(result.end(), line.begin(), line.end());
		}
	}
	
	return result;
}

std::vector<size_t> readIndexVector(const char *filename) {
	size_t numCols = 1;
	return readMatrix<size_t>(filename, &numCols);
}

Eigen::VectorXd readDoubleVector(const char *filename) {
	size_t numCols = 1;
	std::vector<double> raw = readMatrix<double>(filename, &numCols);
	
	Eigen::VectorXd result = Eigen::VectorXd::Map(&raw[0], raw.size());
	return result;
}

Eigen::MatrixXd readDoubleMatrix(const char *filename) {
	size_t numCols = 0;
	std::vector<double> raw = readMatrix<double>(filename, &numCols);
	
	//Ensure the mapping is done row-major.
	Eigen::MatrixXd result = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Map(&raw[0], raw.size() / numCols, numCols);
	return result;
}

///////////////////
//Output functions.
///////////////////

template<typename T> void OutputColumn<T>::printHeader(std::ostream &stream) const {
	stream << name;
}

static void printMultiHeader(std::ostream &stream, const std::string name, size_t numCols) {
	for(size_t i = 0; i < numCols; i++) {
		stream << name << "_" << i;
		if(i != numCols-1) stream << OUTPUT_TABLE_SEPARATOR;
	}
}

//When using printMultiHeader, we will need to assume every entry in the column has the same length.
//The calling function checks that there is at least one entry in the column, sow e can just grab the first.

template<> void OutputColumn<Grouping>::printHeader(std::ostream &stream) const {
	size_t numSpecies = column[0].numSpecies;
	printMultiHeader(stream, name, numSpecies);
}

template<typename T> void OutputColumn<T>::printElement(std::ostream &stream, size_t index) const {
	stream << column[index];
}

template<> void OutputColumn<Grouping>::printElement(std::ostream &stream, size_t index) const {
	size_t numSpecies = column[index].numSpecies;
	for(size_t i = 0; i < numSpecies; i++) {
		stream << column[index].getGroup(i);
		if(i != numSpecies-1) stream << OUTPUT_TABLE_SEPARATOR;
	}
}

//Explicitly instantiate the templates.
//Templated functions existing only in the .cpp file require explicit instantiation.
template class OutputColumn<double>;
template class OutputColumn<Grouping>;
