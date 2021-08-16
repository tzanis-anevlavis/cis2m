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


	// Affine Transformation
	MatrixXd R(2,2);
	R << 1, 2, 1, 0;

	cis2m::HPolyhedron p4 = p1.affineT(R);
	std::cout << " ==== Affine Transformation ==== " << std::endl;
	std::cout << "P4 = R * P1: " << std::endl;
	std::cout << "R: " << std::endl << R << std::endl;
	std::cout << "A: " << std::endl << p4.Ai() << std::endl;
	std::cout << "b: " << std::endl << p4.bi().transpose() << std::endl << std::endl;
	std::cout << "Ae: " << std::endl << p4.Ae() << std::endl;
	std::cout << "be: " << std::endl << p4.be().transpose() << std::endl;

	bool badflag = 0;
	double expectedA[][2] {0, 1, 1.5, -0.5, 0, -1, -0.5, -0.5};
	double expectedb[] {4, 4, 4, 4};

	for (int i = 0; i < 4; i++) {
		MatrixXd A = p4.Ai();
		VectorXd b = p4.bi();
		for (int j = 0; j < 2; j++) {
			if (A(i, j) != expectedA[i][j])  {
				badflag = 1;
				std::cout << "Something wrong here!" << std::endl;
				return badflag;
			}
		}
		if (b(i) != expectedb[i]) {
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
