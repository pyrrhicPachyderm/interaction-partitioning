#ifndef DATA_HPP
#define DATA_HPP

#include <stddef.h>
#include <assert.h>
#include <vector>
#include <Eigen/Core>
#include "model.hpp"
#include "grouping.hpp"

class Parameters;

namespace Datasets {
	//A Jacobian will be represented simply as a matrix.
	//It has a number of rows equal to the number of observations, and
	//a number of columns equal to the number of parameters.
	//It treats parameters in the order given by Parameters::getAsVector().
	//TODO: Do this more nicely.
	typedef Eigen::MatrixXd Jacobian;
	
	class Base {
		protected:
			size_t numRowSpecies;
			size_t numColSpecies;
			size_t numObservations;
		public:
			Base() = default;
			Base(size_t numRowSpecies, size_t numColSpecies, size_t numObservations):
				numRowSpecies(numRowSpecies), numColSpecies(numColSpecies), numObservations(numObservations) {};
			
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
			
			virtual Eigen::VectorXd getObservations() const = 0;
			virtual Eigen::VectorXd getPredictions(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const = 0;
			Eigen::VectorXd getResiduals(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const {
				return getObservations() - getPredictions(model, parameters, groupings);
			}
			
			//NLS requires the Jacobian of the residuals, which is the negative of the Jacobian of the predicted values.
			virtual Jacobian getPredictionsJacobian(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const = 0;
			Jacobian getResidualsJacobian(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const {
				return -getPredictionsJacobian(model, parameters, groupings);
			}
			
			//Functions to get rough guesses of parameter values from the data.
			//Useful for initial values of iterative processes, or tolerances, but not much else.
			//0 is a good guess for the values of competition coefficients,
			//so there is a separate function for guessing the approximate maximum magnitude of them.
			virtual double guessGrowthRate() const = 0;
			double guessCompetitionCoefficient() const {return 0.0;}
			virtual double guessCompetitionCoefficientMagnitude() const = 0;
			virtual double guessErrorVariance() const = 0;
	};
	
	class IndividualResponse : public Base {
		protected:
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
			IndividualResponse() = default;
			IndividualResponse(const std::vector<size_t> &focal, const Eigen::VectorXd &response, const Eigen::MatrixXd &design):
				Base(findNumFocals(focal), design.cols(), design.rows()),
				focal(focal),
				response(response),
				design(design)
			{
				assert(focal.size() == numObservations);
				assert((size_t)response.size() == numObservations);
				assert(areFocalsFirst());
			};
		protected:
			Eigen::MatrixXd getColGroupedDesign(const Grouping &grouping) const;
		public:
			Eigen::VectorXd getObservations() const override {return response;};
			Eigen::VectorXd getPredictions(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const override;
			
			Jacobian getPredictionsJacobian(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const override;
			
			double guessGrowthRate() const override;
			double guessCompetitionCoefficientMagnitude() const override;
			double guessErrorVariance() const override;
	};
}

class Data {
	protected:
		std::shared_ptr<const Datasets::Base> d;
	public:
		Data() = default;
		Data(const Datasets::Base *d):
			d(d) {}
		
		size_t getNumSpecies(GroupingType type) const {return d->getNumSpecies(type);}
		size_t getNumRowSpecies() const {return d->getNumRowSpecies();};
		size_t getNumColSpecies() const {return d->getNumColSpecies();};
		size_t getNumObservations() const {return d->getNumObservations();};
		
		Eigen::VectorXd getObservations() const {return d->getObservations();}
		Eigen::VectorXd getPredictions(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const {return d->getPredictions(model, parameters, groupings);}
		Eigen::VectorXd getResiduals(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const {return d->getResiduals(model, parameters, groupings);}
		
		Datasets::Jacobian getPredictionsJacobian(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const {return d->getPredictionsJacobian(model, parameters, groupings);}
		Datasets::Jacobian getResidualsJacobian(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const {return d->getResidualsJacobian(model, parameters, groupings);}
		
		double guessGrowthRate() const {return d->guessGrowthRate();}
		double guessCompetitionCoefficient() const {return d->guessCompetitionCoefficient();}
		double guessCompetitionCoefficientMagnitude() {return d->guessCompetitionCoefficientMagnitude();}
		double guessErrorVariance() const {return d->guessErrorVariance();}
};

#endif
