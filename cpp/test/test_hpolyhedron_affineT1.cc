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


	// Affine Transformation (Not full rank)
	MatrixXd R2(2,2);
	R2 << 1, 2, 0, 0;

	cis2m::HPolyhedron p5 = p1.affineT(R2);
	std::cout << " ==== Affine Transformation ==== " << std::endl;
	std::cout << "P4 = R2 * P1: " << std::endl;
	std::cout << "R2: " << std::endl << R2 << std::endl;
	std::cout << "A: " << std::endl << p5.Ai() << std::endl;
	std::cout << "b: " << std::endl << p5.bi().transpose() << std::endl << std::endl;
	std::cout << "Ae: " << std::endl << p5.Ae() << std::endl;
	std::cout << "be: " << std::endl << p5.be().transpose() << std::endl << std::endl;

	double expectedA2[][2] {3.0, 0, 0, 0, -1.0, 0.0};
	double expectedb2[] {12, 8, 12};

	bool badflag = 0;
	MatrixXd A_tr = p5.Ai();
	VectorXd b_tr = p5.bi();
	for (int i = 0; i < A_tr.rows(); i++) {
		for (int j = 0; j < A_tr.cols(); j++) {
			if (A_tr(i, j) != expectedA2[i][j])  {
				badflag = 1;
				std::cout << "Something wrong here!" << std::endl;
				std::cout << std::endl << A_tr << std::endl;
				std::cout << std::endl << b_tr << std::endl;
				return badflag;
			}
		}
		if (b_tr(i) != expectedb2[i]) {
			badflag = 1;
			std::cout << "Something wrong here!" << std::endl;
			std::cout << std::endl << b_tr << std::endl;
			return badflag;
		}
	}


	if (!badflag) {
		std::cout << "Looks good!" << std::endl;
	} 

	return badflag;
}
