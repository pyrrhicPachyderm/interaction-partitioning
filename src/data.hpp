#ifndef DATA_HPP
#define DATA_HPP

#include <stddef.h>
#include <assert.h>
#include <vector>
#include "types.hpp"
#include "model.hpp"
#include "grouping.hpp"

class Parameters;

namespace Datasets {
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
			
			virtual const Eigen::VectorXd &getObservations() const = 0;
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
			//The same *may* be true for the growth rates, depending on the data type (see time series),
			//So that is also split.
			virtual double guessGrowthRate() const = 0;
			virtual double guessGrowthRateMagnitude() const = 0;
			double guessCompetitionCoefficient() const {return 0.0;}
			virtual double guessCompetitionCoefficientMagnitude() const = 0;
			virtual double guessErrorVariance() const = 0;
			
			//Virtual destructor, as this is an abstract class.
			virtual ~Base() {};
	};
	
	class FocalResponse : public Base {
		protected:
			const bool isPerCapita;
			
			//The index of the species used as the focal in each observation.
			std::vector<size_t> focal;
			
			//The vector of response variables, one for each observation.
			Eigen::VectorXd response;
			
			//The design density matrix, with numColSpecies columns and numObservations rows.
			Eigen::MatrixXdRowMajor design;
		private:
			static size_t findNumFocals(std::vector<size_t> focals);
			
			//There may be fewer focals (rowSpecies) than total species (columnSpecies).
			//To keep things simple, we will require that all focal species precede all non-focal species.
			//This is checked by areFocalsFirst().
			bool areFocalsFirst() const;
		public:
			FocalResponse(bool isPerCapita): isPerCapita(isPerCapita) {};
			FocalResponse(bool isPerCapita, const std::vector<size_t> &focal, const Eigen::VectorXd &response, const Eigen::MatrixXdRowMajor &design):
				Base(findNumFocals(focal), design.cols(), design.rows()),
				isPerCapita(isPerCapita),
				focal(focal),
				response(response),
				design(design)
			{
				assert(focal.size() == numObservations);
				assert((size_t)response.size() == numObservations);
				assert(areFocalsFirst());
			};
		protected:
			Eigen::MatrixXdRowMajor getColGroupedDesign(const Grouping &grouping) const;
		public:
			const Eigen::VectorXd &getObservations() const override {return response;};
			Eigen::VectorXd getPredictions(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const override;
			
			Jacobian getPredictionsJacobian(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const override;
			
			double guessGrowthRate() const override;
			double guessGrowthRateMagnitude() const override;
			double guessCompetitionCoefficientMagnitude() const override;
			double guessErrorVariance() const override;
	};
	
	class IndividualResponse : public FocalResponse {
		public:
			IndividualResponse(): FocalResponse(true) {};
			IndividualResponse(const std::vector<size_t> &focal, const Eigen::VectorXd &response, const Eigen::MatrixXdRowMajor &design):
				FocalResponse(true, focal, response, design) {};
	};
	
	class PopSizeResponse : public FocalResponse {
		public:
			PopSizeResponse(): FocalResponse(false) {};
			PopSizeResponse(const std::vector<size_t> &focal, const Eigen::VectorXd &response, const Eigen::MatrixXdRowMajor &design):
				FocalResponse(false, focal, response, design) {};
	};
	
	class TimeSeries : public Base {
		protected:
			//Time series data is loaded in with an arbitrary number of time points per experiment.
			//But we will analyse it as one experiment per pair of consecutive time points.
			
			//The time span of each experiment.
			std::vector<double> timeSpan;
			
			//The density of each species at the start and end of each experiment.
			//Species with zero density at the start of an experiment are not included.
			//To account for this, the indices of all the species included in each experiment are also stored.
			std::vector<std::vector<size_t>> includedSpecies;
			std::vector<Eigen::VectorXd> initialDensity;
			std::vector<Eigen::VectorXd> finalDensity;
			
			//Simply finalDensity, but stored as a sinlge linear vector for convenient return.
			Eigen::VectorXd observations;
			
			double maxStepSize = 1; //TODO: make this a command line option.
			size_t getNumSteps(double timeSpan) const;
		public:
			TimeSeries() = default;
			TimeSeries(const std::vector<size_t> &id, const Eigen::VectorXd &time, const Eigen::MatrixXdRowMajor &density);
		protected:
			//A pair of helper functions to get growth rate vectors and alpha matrices,
			//ungrouped and subsetted for the included species of a particular experiment.
			Eigen::VectorXd getGrowthRates(const Parameters &parameters, const GroupingSet &groupings, size_t experiment) const;
			Eigen::MatrixXdRowMajor getCompetitionCoefficients(const Parameters &parameters, const GroupingSet &groupings, size_t experiment) const;
		public:
			const Eigen::VectorXd &getObservations() const override {return observations;};
			Eigen::VectorXd getPredictions(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const override;
			
			Jacobian getPredictionsJacobian(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const override {
				//This Jacobian would have to be found numerically, not analytically.
				//For now, that's too hard, so we won't be able to use a maximum likelihood solver with time series data.
				assert(false);
			}
			
			double guessGrowthRate() const override;
			double guessGrowthRateMagnitude() const override;
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
		
		const Eigen::VectorXd &getObservations() const {return d->getObservations();}
		Eigen::VectorXd getPredictions(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const {return d->getPredictions(model, parameters, groupings);}
		Eigen::VectorXd getResiduals(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const {return d->getResiduals(model, parameters, groupings);}
		
		Jacobian getPredictionsJacobian(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const {return d->getPredictionsJacobian(model, parameters, groupings);}
		Jacobian getResidualsJacobian(const Model &model, const Parameters &parameters, const GroupingSet &groupings) const {return d->getResidualsJacobian(model, parameters, groupings);}
		
		double guessGrowthRate() const {return d->guessGrowthRate();}
		double guessGrowthRateMagnitude() const {return d->guessGrowthRateMagnitude();}
		double guessCompetitionCoefficient() const {return d->guessCompetitionCoefficient();}
		double guessCompetitionCoefficientMagnitude() const {return d->guessCompetitionCoefficientMagnitude();}
		double guessErrorVariance() const {return d->guessErrorVariance();}
};

#endif
