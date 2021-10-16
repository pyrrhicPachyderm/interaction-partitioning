#ifndef DATA_HPP
#define DATA_HPP

#include <assert.h>
#include <vector>
#include <Eigen/Core>

class Data {
	public:
		const size_t numSpecies;
		const size_t numObservations;
		const bool isPerCapita; //Is the response variable total, or per capita?
	protected:
		//The index of the species used as the focal in each observation.
		std::vector<size_t> focal;
		
		//The vector of response variables, one for each observation.
		Eigen::VectorXd response;
		
		//The design density matrix, with numSpecies columns and numObservations rows.
		Eigen::MatrixXd design;
	public:
		Data(const std::vector<size_t> &focal, const Eigen::VectorXd &response, const Eigen::MatrixXd &design, bool isPerCapita):
			numSpecies(design.cols()), numObservations(design.rows()), isPerCapita(isPerCapita), focal(focal), response(response), design(design)
		{
			assert(focal.size() == numObservations);
			assert((size_t)response.size() == numObservations);
		};
		
		const std::vector<size_t> &getFocal() const {
			return focal;
		};
		
		const Eigen::VectorXd &getResponse() const {
			return response;
		};
		
		const Eigen::MatrixXd &getDesign() const {
			return design;
		};
		
		double getResponseMean() const {
			return response.mean();
		};
		
		double getResponseVariance() const {
			Eigen::VectorXd residuals = response - Eigen::VectorXd::Constant(response.size(), response.mean());
			return residuals.dot(residuals) / residuals.size();
		};
};

#endif
