#ifndef DATA_HPP
#define DATA_HPP

#include <assert.h>
#include <vector>
#include <Eigen/Core>

class Data {
	public:
		const size_t numSpecies;
		const size_t numObservations;
	protected:
		//The index of the species used as the focal in each observation.
		std::vector<size_t> focal;
		
		//The vector of response variables, one for each observation.
		Eigen::VectorXd response;
		
		//The design density matrix, with numSpecies columns and numObservations rows.
		Eigen::MatrixXd design;
	public:
		Data(std::vector<size_t> focal, Eigen::VectorXd response, Eigen::MatrixXd design):
			numSpecies(design.cols()), numObservations(design.rows()), focal(focal), response(response), design(design)
		{
			assert(focal.size() == numObservations);
			assert((size_t)response.rows() == numObservations);
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
};

#endif
