#pragma once
#include <vector>
#include <Eigen/Dense>
#include "hpolyhedron.hpp"


namespace cis2m {
	class BrunovskyForm {
		public:
			/// Constructor
			BrunovskyForm(const Eigen::MatrixXd& Ad, const Eigen::MatrixXd& Bd);

			/// Destructor
			~BrunovskyForm();

			/**
			 * Get the Brunovsky Form A and B matrices
			 */
			std::pair<Eigen::MatrixXd, Eigen::MatrixXd>  GetDynSystem();

			/** 
			 * Transform the dynamics constraints expressed with 
			 * the Controllability Form basis 
			 */
			HPolyhedron GetDynConstraints(const HPolyhedron& dyncnstr);

			/**
			 * Transform the input constraints considering the transformation
			 * used to get the Brunovksy form.
			 */
			HPolyhedron GetInputConstraints(const HPolyhedron& InputCnstr);

			/**
			 * Transform the disturbance matrix with the Controllability Form
			 * basis
			 */
			Eigen::MatrixXd GetDisturbanceMatrix(const Eigen::MatrixXd& DisturbanceMatrix);

			/**
			 *  Get the Transformation Matrix (Change of basis to get the 
			 *  Controllability Form)
			 */
			Eigen::MatrixXd GetTransformationMatrix();


			/**
			 * Return the controllability indexes
			 */
			std::vector<int> GetControllabilityIndexes();


			/**
			 * Return the max controllability index
			 */
			int GetMaxControllabilityIndex();


		private:
			/// Dynamics A matrix in controllability form
			Eigen::MatrixXd A_cntrl_form_;

			/// Dynamics B matrix in controllability form
			Eigen::MatrixXd B_cntrl_form_;

			/// Dynamics A matrix in Brunovsky form
			Eigen::MatrixXd A_brunovsky_form_;
			Eigen::MatrixXd B_brunovsky_form_;

			/// Controllability Indexes
			std::vector<int> controllability_indexes_;

			/// Max Controllability Index
			int max_controllability_index_;

			/// Transformation Matrix
			Eigen::MatrixXd TransformationMatrix_;

			/// Brunovsky Form Bm Matrix
			Eigen::MatrixXd Bm_;
	};
}
