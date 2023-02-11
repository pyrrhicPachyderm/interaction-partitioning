#ifndef DATA_HPP
#define DATA_HPP

#include <stddef.h>
#include <assert.h>
#include <vector>
#include <Eigen/Core>
#include "grouping.hpp"

enum GroupingType {GROWTH, ROW, COL, NUM_GROUPING_TYPES};

class Data {
	protected:
		size_t numRowSpecies;
		size_t numColSpecies;
		size_t numObservations;
		
		//The index of the species used as the focal in each observation.
		std::vector<size_t> focal;
		
		//The vector of response variables, one for each observation.
		Eigen::VectorXd response;
		
		//The design density matrix, with numColSpecies columns and numObservations rows.
		Eigen::MatrixXd design;
	private:
		static size_t findNumFocals(std::vector<size_t> focals);
		
		//There may be fewer focals (rowSpecies) than total species (columnSpecies).
		//To keep things simple, we will require that all focal species precede all non-focal species.
		//This is checked by areFocalsFirst().
		bool areFocalsFirst() const;
	public:
		Data() = default;
		Data(const std::vector<size_t> &focal, const Eigen::VectorXd &response, const Eigen::MatrixXd &design):
			numRowSpecies(findNumFocals(focal)),
			numColSpecies(design.cols()),
			numObservations(design.rows()),
			focal(focal),
			response(response),
			design(design)
		{
			assert(focal.size() == numObservations);
			assert((size_t)response.size() == numObservations);
			assert(areFocalsFirst());
		};
		
		size_t getNumSpecies(GroupingType type) const {
			switch(type) {
				case GROWTH:
					return numRowSpecies;
				case ROW:
					return numRowSpecies;
				case COL:
					return numColSpecies;
				default:
					__builtin_unreachable();
			}
		}
		size_t getNumRowSpecies() const {return numRowSpecies;};
		size_t getNumColSpecies() const {return numColSpecies;};
		size_t getNumObservations() const {return numObservations;};
		
		const std::vector<size_t> &getFocal() const {return focal;};
		const Eigen::VectorXd &getResponse() const {return response;};
		const Eigen::MatrixXd &getDesign() const {return design;};
		Eigen::MatrixXd getColGroupedDesign(const Grouping &grouping) const;
		
		//Functions to get rough guesses of parameter values from the data.
		//Useful for initial values of iterative processes, or tolerances, but not much else.
		//0 is a good guess for the values of competition coefficients,
		//so there is a separate function for guessing the approximate maximum magnitude of them.
		double guessGrowthRate() const;
		double guessCompetitionCoefficient() const {return 0.0;}
		double guessCompetitionCoefficientMagnitude() const;
		double guessErrorVariance() const;
};

#endif
