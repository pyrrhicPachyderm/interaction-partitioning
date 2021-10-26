#ifndef DATA_HPP
#define DATA_HPP

#include <assert.h>
#include <vector>
#include <Eigen/Core>

class Data {
	protected:
		size_t numSpecies;
		size_t numObservations;
		bool isPerCapita; //Is the response variable total, or per capita?
		
		//The index of the species used as the focal in each observation.
		std::vector<size_t> focal;
		
		//The vector of response variables, one for each observation.
		Eigen::VectorXd response;
		
		//The design density matrix, with numSpecies columns and numObservations rows.
		Eigen::MatrixXd design;
	public:
		Data() = default;
		Data(const std::vector<size_t> &focal, const Eigen::VectorXd &response, const Eigen::MatrixXd &design, bool isPerCapita):
			numSpecies(design.cols()), numObservations(design.rows()), isPerCapita(isPerCapita), focal(focal), response(response), design(design)
		{
			assert(focal.size() == numObservations);
			assert((size_t)response.size() == numObservations);
		};
		
		size_t getNumSpecies() const {return numSpecies;};
		size_t getNumObservations() const {return numObservations;};
		bool getIsPerCapita() const {return isPerCapita;};
		
		const std::vector<size_t> &getFocal() const {return focal;};
		const Eigen::VectorXd &getResponse() const {return response;};
		const Eigen::MatrixXd &getDesign() const {return design;};
		
		double getResponseMean() const;
		double getResponseVariance() const;
		
		//Functions to get rough guesses of parameter values from the data.
		//Useful for initial values of iterative processes, or tolerances, but not much else.
		//0 is a good guess for the values of competition coefficients,
		//so there is separate a function for guessing the approximate maximum magnitude of them.
		double guessGrowthRate() const;
		double guessCompetitionCoefficient() const {return 0.0;}
		double guessCompetitionCoefficientMagnitude() const;
};

#endif
