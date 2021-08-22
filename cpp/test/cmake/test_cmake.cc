#include <iostream>
#include <cis2m/hpolyhedron.hpp>
#include <Eigen/Dense>

int main() {
	Eigen::MatrixXd A(Eigen::MatrixXd(Eigen::MatrixXd::Identity(3,3)));
	Eigen::VectorXd B(Eigen::VectorXd(Eigen::VectorXd::Zero(3,1)));
	cis2m::HPolyhedron poly(A, B);
	return 0;
}
