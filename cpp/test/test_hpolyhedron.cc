#include "hpolyhedron.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <vector>

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


	MatrixXd A2(4, 2);
	A2 << 2, 0, 0, 2, -1, 0, 0, -1;
	VectorXd b2(4);
	b2 << 1, 1, 1, 1;
	cis2m::HPolyhedron p2(A2, b2);
	std::cout << "P2: " << std::endl;
	std::cout << "A: " << std::endl << p2.Ai() << std::endl;
	std::cout << "b: " << std::endl << p2.bi().transpose() << std::endl << std::endl;	

	// Minkowski 1
	cis2m::HPolyhedron p3 = (p1 - p2);
	std::cout << " ===== Minkowski Difference ===== " << std::endl;
	std::cout << "P3 = P1 - P2: " << std::endl;
	std::cout << "A: " << std::endl << p3.Ai() << std::endl;
	std::cout << "b: " << std::endl << p3.bi().transpose() << std::endl;

	bool badflag = 0;

	VectorXd newbi = p3.bi();;
	std::vector<double> expected {3.5, 2.0, 3.0, 2.0};
	for (int i = 0; i < 4; i++) {
		if (newbi(i) != expected[i])  {
			badflag = 1;
			std::cout << "Something wrong here!" << std::endl;
			return badflag;
		}
	}
	std::cout << "Looks good!" << std::endl;
	std::cout << " =====  ===== " << std::endl;
	std::cout << std::endl;


	// Minkowski 2
	MatrixXd AA1(5, 2);
	AA1 << 1, 0, 0, 1, 0, -1.1, -1, 0, -1, -1;
	VectorXd bb1(5);
	bb1 << 6, 6, 0, 0, -4;
	cis2m::HPolyhedron pp1(AA1, bb1);
	std::cout << "P1: " << std::endl;
	std::cout << "A: " << std::endl << pp1.Ai() << std::endl;
	std::cout << "b: " << std::endl << pp1.bi().transpose() << std::endl << std::endl;	


	MatrixXd AA2(4, 2);
	AA2 << 1, 0, 0, 1, -1, 0, 0, -1;
	VectorXd bb2(4);
	bb2 << 1, 1, 0, 0;
	cis2m::HPolyhedron pp2(AA2, bb2);
	std::cout << "P2: " << std::endl;
	std::cout << "A: " << std::endl << pp2.Ai() << std::endl;
	std::cout << "b: " << std::endl << pp2.bi().transpose() << std::endl << std::endl;	

	cis2m::HPolyhedron pp3 = (pp1 - pp2);
	std::cout << " ===== Minkowski Difference ===== " << std::endl;
	std::cout << "P3 = P1 - P2: " << std::endl;
	std::cout << "A: " << std::endl << pp3.Ai() << std::endl;
	std::cout << "b: " << std::endl << pp3.bi().transpose() << std::endl;

	std::vector<double> expected2 {5.0, 5.0, 0, 0, -4.0};
	for (int i = 0; i < 5; i++) {
		if (pp3.bi()(i) != expected2[i])  {
			badflag = 1;
			std::cout << "Something wrong here!" << std::endl;
			return badflag;
		}
	}
	std::cout << "Looks good!" << std::endl;
	std::cout << " =====  ===== " << std::endl;
	std::cout << std::endl;




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
	std::cout << "Looks good!" << std::endl;
	std::cout << " =====  ===== " << std::endl;
	std::cout << std::endl;

	// Affine Transformation
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

	double expectedA2[][2] {-1.0, 0, 1.0, 0.0};
	double expectedb2[] {12, 4};

	double expectedAe2[][2] {0, -1.0};
	double expectedbe2[] {0};

	for (int i = 0; i < 2; i++) {
		MatrixXd A = p5.Ai();
		VectorXd b = p5.bi();
		for (int j = 0; j < 2; j++) {
			if (A(i, j) != expectedA2[i][j])  {
				badflag = 1;
				std::cout << "Something wrong here!" << std::endl;
				return badflag;
			}
		}
		if (b(i) != expectedb2[i]) {
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
