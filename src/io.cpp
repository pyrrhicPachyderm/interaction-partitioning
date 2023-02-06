#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iterator>
#include "grouping.hpp"
#include "parameters.hpp"
#include "io.hpp"

const char *OUTPUT_TABLE_SEPARATOR = "\t";
const char *OUTPUT_GROWTH_RATE_STRING = "r";
const char *OUTPUT_COMPETITION_COEFFICIENT_STRING = "alpha";
const char *OUTPUT_NULL_VALUE_STRING = "NA";

template<typename StreamT> static StreamT &openFile(std::string filename) {
	StreamT &file = *(new StreamT);
	file.open(filename);
	if(!file) {
		fprintf(stderr, "Could not open file %s\n", filename.c_str());
		exit(1);
	}
	return file;
}

std::istream &openInput(std::string filename) {
	if(filename == "-") return std::cin;
	return openFile<std::ifstream>(filename);
}

std::ostream &openOutput(std::string filename) {
	if(filename == "-") return std::cout;
	return openFile<std::ofstream>(filename);
}

//////////////////
//Input functions.
//////////////////

//Set numCols to 0 for a dynamically sized matrix.
//In this case, the number of columns is returned through numCols.
template<typename T> static std::vector<T> readMatrix(std::string filename, size_t *numCols) {
	std::istream &file = openInput(filename);
	
	std::vector<T> result;
	
	std::string lineString;
	while(std::getline(file, lineString)) {
		std::istringstream buffer = std::istringstream(lineString);
		std::vector<T> line = std::vector<T>{std::istream_iterator<T>(buffer), std::istream_iterator<T>()};
		
		if(*numCols != 0 && *numCols != line.size()) {
			fprintf(stderr, "%s has the wrong number of columns, or an inconsistent number of columns\n", filename.c_str());
			exit(1);
		} else {
			*numCols = line.size();
			result.insert(result.end(), line.begin(), line.end());
		}
	}
	
	return result;
}

std::vector<size_t> readIndexVector(std::string filename) {
	size_t numCols = 1;
	return readMatrix<size_t>(filename, &numCols);
}

Eigen::VectorXd readDoubleVector(std::string filename) {
	size_t numCols = 1;
	std::vector<double> raw = readMatrix<double>(filename, &numCols);
	
	Eigen::VectorXd result = Eigen::VectorXd::Map(&raw[0], raw.size());
	return result;
}

Eigen::MatrixXd readDoubleMatrix(std::string filename) {
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
//The calling function checks that there is at least one entry in the column, so we can just grab the first.

template<> void OutputColumn<Grouping>::printHeader(std::ostream &stream) const {
	size_t numSpecies = column[0].numSpecies;
	printMultiHeader(stream, name, numSpecies);
}

template<> void OutputColumn<Eigen::VectorXd>::printHeader(std::ostream &stream) const {
	size_t numCols = column[0].size();
	printMultiHeader(stream, name, numCols);
}

template<> void OutputColumn<Parameters>::printHeader(std::ostream &stream) const {
	size_t numRowSpecies = column[0].getNumRowSpecies();
	size_t numColSpecies = column[0].getNumColSpecies();
	printMultiHeader(stream, name + "_" + OUTPUT_GROWTH_RATE_STRING, numRowSpecies);
	stream << OUTPUT_TABLE_SEPARATOR;
	
	for(size_t i = 0; i < numRowSpecies; i++) {
		printMultiHeader(stream, name + "_" + OUTPUT_COMPETITION_COEFFICIENT_STRING + "_" + std::to_string(i), numColSpecies);
		if(i != numRowSpecies-1) stream << OUTPUT_TABLE_SEPARATOR;
	}
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

template<> void OutputColumn<Eigen::VectorXd>::printElement(std::ostream &stream, size_t index) const {
	size_t numCols = column[index].size();
	for(size_t i = 0; i < numCols; i++) {
		stream << column[index][i];
		if(i != numCols-1) stream << OUTPUT_TABLE_SEPARATOR;
	}
}

template<> void OutputColumn<Parameters>::printElement(std::ostream &stream, size_t index) const {
	size_t numRowSpecies = column[index].getNumRowSpecies();
	size_t numColSpecies = column[index].getNumColSpecies();
	
	Eigen::VectorXd growthRates = column[index].getGrowthRates();
	for(size_t i = 0; i < numRowSpecies; i++) {
		if(i < (size_t)growthRates.size()) {
			stream << growthRates[i];
		} else {
			stream << OUTPUT_NULL_VALUE_STRING;
		}
		stream << OUTPUT_TABLE_SEPARATOR;
	}
	
	Eigen::MatrixXd competitionCoefficients = column[index].getCompetitionCoefficients();
	for(size_t i = 0; i < numRowSpecies; i++) {
		for(size_t j = 0; j < numColSpecies; j++) {
			if(i < (size_t)competitionCoefficients.rows() && j < (size_t)competitionCoefficients.cols()) {
				stream << competitionCoefficients(i, j);
			} else {
				stream << OUTPUT_NULL_VALUE_STRING;
			}
			if(i != numRowSpecies-1 || j != numColSpecies-1) stream << OUTPUT_TABLE_SEPARATOR;
		}
	}
}

//Explicitly instantiate the templates.
//Templated functions existing only in the .cpp file require explicit instantiation.
template class OutputColumn<double>;
template class OutputColumn<size_t>;
template class OutputColumn<Grouping>;
template class OutputColumn<Eigen::VectorXd>;
template class OutputColumn<Parameters>;
