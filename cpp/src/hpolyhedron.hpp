#pragma once

#include <Eigen/Dense>

namespace cis2m {
	class HPolyhedron {
		public:
			/**
			 * Constructors
			 */
			HPolyhedron();
			HPolyhedron(const Eigen::MatrixXd& A, const Eigen::MatrixXd& b);
			HPolyhedron(
					const Eigen::MatrixXd& Ai,
					const Eigen::MatrixXd& bi,
					const Eigen::MatrixXd& Ae,
					const Eigen::MatrixXd& be);
			~HPolyhedron();

			Eigen::VectorXd ComputeSupport(const Eigen::MatrixXd& A_other) const;

			HPolyhedron Projection(int el_index, double tol=1e-6);

			HPolyhedron affineT(const Eigen::MatrixXd& T);

			// Operators
			HPolyhedron& operator=(const HPolyhedron& other);

			HPolyhedron& operator+=(const HPolyhedron& P);
			HPolyhedron& operator-=(const HPolyhedron& P);


			// Getters
			Eigen::MatrixXd Ai() const;
			Eigen::VectorXd bi() const; 
			Eigen::MatrixXd Ae() const;
			Eigen::VectorXd be() const; 

			int GetNumInequalities();
			int GetSpaceDim();

		private:
			Eigen::MatrixXd Ai_;
			Eigen::VectorXd bi_;

			Eigen::MatrixXd Ae_;
			Eigen::VectorXd be_;

			int SpaceDim_;
			int NumIneqs_;
	};

	// Operators  outside che class
	HPolyhedron operator+(HPolyhedron Pl, HPolyhedron Pr);
	HPolyhedron operator-(HPolyhedron Pl, HPolyhedron Pr);
}
