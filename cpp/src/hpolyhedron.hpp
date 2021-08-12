#pragma once

#include <Eigen/Dense>

namespace cis2m {
	class HPolyhedron {
		public:
			HPolyhedron();
			HPolyhedron(const Eigen::MatrixXd& A, const Eigen::MatrixXd& b);
			~HPolyhedron();


			Eigen::VectorXd ComputeSupport(const Eigen::MatrixXd& A_other) const;

			// Operators
			HPolyhedron& operator=(const HPolyhedron& other);

			HPolyhedron& operator+=(const HPolyhedron& P);
			HPolyhedron& operator-=(const HPolyhedron& P);

			HPolyhedron& operator+(const HPolyhedron& P);
			HPolyhedron& operator-(const HPolyhedron& P);

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
}
