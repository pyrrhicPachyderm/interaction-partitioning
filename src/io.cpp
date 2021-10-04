#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iterator>
#include "grouping.hpp"
#include "io.hpp"

const char *OUTPUT_TABLE_SEPARATOR = "\t";

template<typename StreamT> StreamT openFile(const char *filename) {
	StreamT file;
	file.open(filename);
	if(!file) {
		fprintf(stderr, "Could not open file %s\n", filename);
		exit(1);
	}
	return file;
}

//////////////////
//Input functions.
//////////////////

//Set numCols to 0 for a dynamically sized matrix.
//In this case, the number of columns is returned through numCols.
template<typename T> static std::vector<T> readMatrix(const char *filename, size_t *numCols) {
	std::ifstream file = openFile<std::ifstream>(filename);
	
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
	
	file.close();
	
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

template<typename T> void OutputColumn<T>::printHeader(std::ostream stream) const {
	stream << name;
}

template<> void OutputColumn<Grouping>::printHeader(std::ostream stream) const {
	//We need to know the number of species in a grouping to print the appropriate number of headers.
	//We will assume all groupings in the vector have the same number of species.
	//So we will just check the first.
	//This requires us to be sure that there is a first.
	if(column.empty()) {
		stream << name;
		return;
	}
	
	size_t numSpecies = column[0].numSpecies;
	for(size_t i = 0; i < numSpecies; i++) {
		stream << name << "_" << i;
		if(i != numSpecies-1) stream << OUTPUT_TABLE_SEPARATOR;
	}
}

template<typename T> void OutputColumn<T>::printElement(std::ostream stream, size_t index) const {
	stream << column[index];
}

template<> void OutputColumn<Grouping>::printElement(std::ostream stream, size_t index) const {
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
