#ifndef IO_HPP
#define IO_HPP

#include <stddef.h>
#include <vector>
#include <iostream>
#include <Eigen/Core>
#include "distribution.hpp"

/////////////////////////////////
//"Internal" variables/functions.
/////////////////////////////////

//These should only be needed by the io module.
//However, they are used by outputTable, which needs to be defined in the header file.
//It is far too complicated a template to explicitly instantiate.

extern const char *OUTPUT_TABLE_SEPARATOR;

std::istream &openInput(std::string filename);
std::ostream &openOutput(std::string filename);

//////////////////
//Input functions.
//////////////////

extern std::vector<size_t> readIndexVector(std::string filename);
extern Eigen::VectorXd readDoubleVector(std::string filename);
extern Eigen::MatrixXd readDoubleMatrix(std::string filename);

///////////////////
//Output functions.
///////////////////

template<typename T> class OutputColumn {
	protected:
		const std::string name;
		std::vector<T> column;
	public:
		OutputColumn(std::string name): name(name) {};
		
		void insert(T element) {
			column.push_back(element);
		}
		
		size_t getLength() const {
			return column.size();
		}
		
		void printHeader(std::ostream &stream) const;
		
		void printElement(std::ostream &stream, size_t index) const;
};

//The following functions use a template parameter pack.
//It's impractical to explicitly instantiate this, so alas it must live in the header file.

template<class T, class... OutputColumnTs> void outputTableHeader(std::ostream &stream, const OutputColumn<T>& column, const OutputColumnTs&... columns) {
	if(column.getLength() > 0) {
		column.printHeader(stream);
		if constexpr(sizeof...(columns) > 0) {
			stream << OUTPUT_TABLE_SEPARATOR;
			outputTableHeader(stream, columns...);
		}
	} else if constexpr(sizeof...(columns) > 0) {
		outputTableHeader(stream, columns...);
	}
}

template<class T, class... OutputColumnTs> void outputTableRow(std::ostream &stream, size_t row, const OutputColumn<T>& column, const OutputColumnTs&... columns) {
	if(column.getLength() > 0) {
		column.printElement(stream, row);
		if constexpr(sizeof...(columns) > 0) {
			stream << OUTPUT_TABLE_SEPARATOR;
			outputTableRow(stream, row, columns...);
		}
	} else if constexpr(sizeof...(columns) > 0) {
		outputTableRow(stream, row, columns...);
	}
}

template<class T, class... OutputColumnTs> void outputTable(std::string filename, const OutputColumn<T>& column, const OutputColumnTs&... columns) {
	//We can't really iterate over the columns to check they're all the same length, so we're just going to assume.
	//TODO: Improve this.
	size_t numRows = column.getLength();
	
	std::ostream &stream = openOutput(filename);
	
	outputTableHeader(stream, column, columns...);
	stream << "\n";
	
	for(size_t row = 0; row < numRows; row++) {
		outputTableRow(stream, row, column, columns...);
		stream << "\n";
	}
	
	stream.flush();
}

#endif
