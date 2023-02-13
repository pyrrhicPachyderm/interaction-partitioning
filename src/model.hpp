#ifndef MODEL_HPP
#define MODEL_HPP

#include <memory>
#include "types.hpp"

namespace Models {
	//An abstract class, for other models to inherit from.
	class Base {
		public:
			//The core function of a model: dN/dt or (N_{t+1} - N_t).
			virtual double getDerivative(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients) const = 0;
			
			//A wrapper around getDerivatives that handles a population.
			//It takes one set of densities, treating each as focal (in turn) and as competitors.
			Eigen::VectorXd getDerivatives(const Eigen::VectorXd &densities, const Eigen::VectorXd &growthRates, const Eigen::MatrixXdRowMajor &competitionCoefficients) const;
			
			//These get single elements of the Jacobian matrix (an element for a single parameter and observation).
			//They are the rate of change of dN/dt with respect to the parameter of interest.
			//Note that the competition coefficient one requires one extra parameter: the index of which competition coefficient we're interested in.
			//Also note that these functions cannot be called for every element of the Jacobian matrix: they cannot be called for parameters which do not affect the focal species of the observation.
			//All elements for which they cannot be called are zero.
			//This makes them useless for time series data.
			virtual double getGrowthRateJacobian(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients) const = 0;
			virtual double getCompetitionCoefficientJacobian(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients, size_t index) const = 0;
			
			//Virtual destructor, as this is an abstract class.
			virtual ~Base() {};
	};
	
	class LotkaVolterra : public Base {
		public:
			double getDerivative(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients) const override;
			double getGrowthRateJacobian(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients) const override;
			double getCompetitionCoefficientJacobian(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients, size_t index) const override;
	};
}

//A wrapper class to deal with maintaining the pointer required for runtime polymorphism.
class Model {
	protected:
		std::shared_ptr<const Models::Base> m;
	public:
		Model(const Models::Base *m):
			m(m) {}
		
		double getDerivative(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients) const {
			return m->getDerivative(focalDensity, focalGrowthRate, densities, competitionCoefficients);
		};
		Eigen::VectorXd getDerivatives(const Eigen::VectorXd &densities, const Eigen::VectorXd &growthRates, const Eigen::MatrixXdRowMajor &competitionCoefficients) const {
			return m->getDerivatives(densities, growthRates, competitionCoefficients);
		};
		double getGrowthRateJacobian(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients) const {
			return m->getGrowthRateJacobian(focalDensity, focalGrowthRate, densities, competitionCoefficients);
		};
		double getCompetitionCoefficientJacobian(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients, size_t index) const {
			return m->getCompetitionCoefficientJacobian(focalDensity, focalGrowthRate, densities, competitionCoefficients, index);
		};
};

#endif
