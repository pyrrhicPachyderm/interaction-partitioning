#ifndef MODEL_HPP
#define MODEL_HPP

#include <memory>
#include <Eigen/Core>

namespace Models {
	//An abstract class, for other models to inherit from.
	class Base {
		public:
			//The core function of a model: dN/dt or (N_{t+1} - N_t).
			virtual double getDerivative(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients) const = 0;
			
			//TODO: Eigen::VectorXd getDerivatives(const Eigen::VectorXd &focalDensities, const Eigen::VectorXd &focalGrowthRates, const Eigen::VectorXd &densities, const Eigen::MatrixXd &competitionCoefficients) const;
			
			//TODO: getGrowthRateJacobian()
			//TODO: getCompetitionCoefficientJacobian()
			
			//Virtual destructor, as this is an abstract class.
			virtual ~Base() {};
	};
	
	class LotkaVolterra : public Base {
		public:
			double getDerivative(double focalDensity, double focalGrowthRate, const Eigen::VectorXd &densities, const Eigen::VectorXd &competitionCoefficients) const override;
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
};

#endif
