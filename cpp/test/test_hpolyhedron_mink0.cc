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

	if (!badflag) {
		std::cout << "Looks good!" << std::endl;
	} 

	return badflag;
}
