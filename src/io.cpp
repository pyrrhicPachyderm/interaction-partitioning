#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iterator>
#include "io.hpp"

template<typename StreamT> static StreamT openFile(const char *filename) {
	StreamT file;
	file.open(filename);
	if(!file) {
		fprintf(stderr, "Could not open file %s\n", filename);
		exit(1);
	}
	return file;
}

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
