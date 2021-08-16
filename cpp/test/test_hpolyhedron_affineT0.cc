#include "hpolyhedron.hpp"
#include <Eigen/Dense>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

int main() {
	MatrixXd A(4, 2);
	A << 1, 0, 1, 3, -1, 0, -1, -1;
	VectorXd b(4);
	b << 4, 4, 4, 4;
	cis2m::HPolyhedron p1(A, b);
	std::cout << "P1: " << std::endl;
	std::cout << "A: " << std::endl << p1.Ai() << std::endl;
	std::cout << "b: " << std::endl << p1.bi().transpose() << std::endl << std::endl;	

	// Affine Transformation (Full Rank)
	MatrixXd R(2,2);
	R << 1, 2, 1, 0;

	cis2m::HPolyhedron p2 = p1.affineT(R);
	std::cout << " ==== Affine Transformation ==== " << std::endl;
	std::cout << "P4 = R * P1: " << std::endl;
	std::cout << "R: " << std::endl << R << std::endl;
	std::cout << "A: " << std::endl << p2.Ai() << std::endl;
	std::cout << "b: " << std::endl << p2.bi().transpose() << std::endl << std::endl;
	std::cout << "Ae: " << std::endl << p2.Ae() << std::endl;
	std::cout << "be: " << std::endl << p2.be().transpose() << std::endl;

	double expectedA[][2] {0, 1, 1.5, -0.5, 0, -1, -0.5, -0.5};
	double expectedb[] {4, 4, 4, 4};

	bool badflag = 0;
	MatrixXd A_tr = p2.Ai();
	VectorXd b_tr = p2.bi();
	for (int i = 0; i < A_tr.rows(); i++) {
		for (int j = 0; j < A_tr.cols(); j++) {
			if (A_tr(i, j) != expectedA[i][j])  {
				badflag = 1;
				std::cout << "Something wrong here!" << std::endl;
				return badflag;
			}
		}
		if (b_tr(i) != expectedb[i]) {
			badflag = 1;
			std::cout << "Something wrong here!" << std::endl;
			return badflag;
		}
	}

	if (!badflag) {
		std::cout << "Looks good!" << std::endl;
	} 

	return badflag;
}
