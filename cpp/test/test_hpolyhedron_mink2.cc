#include "hpolyhedron.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;

int main() {
	// Minkowski 2
	MatrixXd AA1(5, 2);
	AA1 << 1, 0, 0, 1, 0, -1, -1, 0, -1, -1;
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

	bool badflag = 0;
	std::vector<double> expected2 {5.0, 5.0, 0, 0, -4.0};
	for (int i = 0; i < 5; i++) {
		if (pp3.bi()(i) != expected2[i])  {
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
