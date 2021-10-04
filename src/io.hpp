#ifndef IO_HPP
#define IO_HPP

#include <vector>
#include <fstream>
#include <Eigen/Core>

extern const char *OUTPUT_TABLE_SEPARATOR;

template<typename StreamT> StreamT openFile(const char *filename);

extern std::vector<size_t> readIndexVector(const char *filename);
extern Eigen::VectorXd readDoubleVector(const char *filename);
extern Eigen::MatrixXd readDoubleMatrix(const char *filename);

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
		
		void printHeader(std::ostream stream) const;
		
		void printElement(std::ostream stream, size_t index) const;
};

//The following functions use a template parameter pack.
//It's impractical to explicitly instantiate this, so alas it must live in the header file.

template<class T, class... OutputColumnTs> void outputTableHeader(std::ostream stream, OutputColumn<T> column, OutputColumnTs... columns) {
	column.printHeader(stream);
	if(sizeof...(columns) > 0) {
		stream << OUTPUT_TABLE_SEPARATOR;
		outputTableHeader(stream, columns...);
	}
}

template<class T, class... OutputColumnTs> void outputTableRow(std::ostream stream, size_t row, OutputColumn<T> column, OutputColumnTs... columns) {
	column.printElement(stream, row);
	if(sizeof...(columns) > 0) {
		stream << OUTPUT_TABLE_SEPARATOR;
		outputTableRow(stream, row, columns...);
	}
}

template<class T, class... OutputColumnTs> void outputTable(const char *filename, OutputColumn<T> column, OutputColumnTs... columns) {
	//We can't really iterate over the columns to check they're all the same length, so we're just going to assume.
	//TODO: Improve this.
	size_t numRows = column.getLength();
	
	std::ofstream stream = openFile<std::ofstream>(filename);
	
	outputTableHeader(stream, column, columns...);
	stream << std::endl;
	
	for(size_t row = 0; row < numRows; row++) {
		outputTableRow(stream, row, column, columns...);
		stream << std::endl;
	}
	
	stream.close();
}

#endif