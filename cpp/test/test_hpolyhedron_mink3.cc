#include "hpolyhedron.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;

int main() {
	MatrixXd A(18, 9);
	A << 
		-0.0072,   0,    0.0072,         0,         0,         0,         0,         0,         0,
		0.0072,    0,   -0.0072,         0,         0,         0,         0,         0,         0,
		0,         0,         0,   -0.0072,         0,    0.0072,         0,         0,         0,
		0,         0,         0,    0.0072,         0,   -0.0072,         0,         0,         0,
		0,         0,         0,         0,         0,         0,   -0.0072,         0,    0.0072,
		0,         0,         0,         0,         0,         0,    0.0072,         0,   -0.0072,
		0.1200,   -0.2400,    0.1200,    0,         0,         0,         0,         0,         0,
		-0.1200,    0.2400,   -0.1200,   0,         0,         0,         0,         0,         0,
		0,         0,         0,    0.1200,   -0.2400,    0.1200,         0,         0,         0, 
		0,         0,         0,   -0.1200,    0.2400,   -0.1200,         0,         0,         0,
		0,         0,         0,         0,         0,         0,    0.1200,   -0.2400,    0.1200,
		0,         0,         0,         0,         0,         0,   -0.1200,    0.2400,   -0.1200,
		0,         0,         0,         0,         0,         0,    0.0001,    0.0004,    0.0001,
		0,         0,         0,    0.0001,    0.0004,    0.0001,         0,         0,         0,
		-0.0001,   -0.0004,   -0.0001,   0,         0,         0,         0,         0,         0,
		0.0001,    0.0004,    0.0001,    0,         0,         0,         0,         0,         0,
		0,         0,         0,   -0.0002,   -0.0006,   -0.0002,         0,         0,         0,
		0,         0,         0,         0,         0,         0,   -0.0001,   -0.0004,   -0.0001;


	VectorXd b(18);
	b << 0.8660, 0.8660, 0.8660, 0.8660, 0.8660, 0.8660, 2.8319, 2.8319, 2.8319, 2.8319, 2.8319, 2.8319, 0.9487, 0.9487, 0.9487, 0.9487, -0.8321, 0.9487;

	cis2m::HPolyhedron p1(A, b);
	std::cout << "P1: " << std::endl;
	std::cout << "A: " << std::endl << p1.Ai() << std::endl;
	std::cout << "b: " << std::endl << p1.bi().transpose() << std::endl << std::endl;	


	MatrixXd A2(18, 9);
	A2 <<
		0.0003,    0.0012,    0.0003,         0,         0,         0,         0,         0,         0,
		0,         0,         0,    0.0003,    0.0012,    0.0003,         0,         0,         0,
		0,         0,         0,         0,         0,         0,    0.0003,    0.0012,    0.0003,
		-0.0072,         0,    0.0072,         0,         0,         0,         0,         0,         0,
		0,         0,         0,   -0.0072,         0,    0.0072,         0,         0,         0,
		0,         0,         0,         0,         0,         0,   -0.0072,         0,    0.0072,
		0.1200,   -0.2400,    0.1200,         0,         0,         0,         0,         0,         0,
		0,         0,         0,    0.1200,   -0.2400,    0.1200,         0,         0,         0,
		0,         0,         0,         0,         0,         0,    0.1200,   -0.2400,    0.1200,
		-0.0003,   -0.0012,   -0.0003,         0,         0,         0,         0,         0,         0,
		0,         0,         0,   -0.0003,   -0.0012,   -0.0003,         0,         0,         0,
		0,         0,         0,         0,         0,         0,   -0.0003,   -0.0012,   -0.0003,
		0.0072,         0,   -0.0072,         0,         0,         0,         0,         0,         0,
		0,         0,         0,    0.0072,         0,   -0.0072,         0,         0,         0,
		0,         0,         0,         0,         0,         0,    0.0072,         0,   -0.0072,
		-0.1200,    0.2400,   -0.1200,         0,         0,         0,         0,         0,         0,
		0,         0,         0,   -0.1200,    0.2400,   -0.1200,         0,         0,         0,
		0,         0,         0,         0,         0,         0,   -0.1200,    0.2400,   -0.1200;


	VectorXd b2(VectorXd::Identity(18,1));
	for (int i = 0; i < 18; i++) {
		b2(i) = 0.01;
	}

	cis2m::HPolyhedron p2(A2, b2);
	std::cout << "P2: " << std::endl;
	std::cout << "A: " << std::endl << p2.Ai() << std::endl;
	std::cout << "b: " << std::endl << p2.bi().transpose() << std::endl << std::endl;	

	std::cout << std::endl;
	std::cout << " ===== Minkowski Difference ===== " << std::endl;
	cis2m::HPolyhedron p3 = (p1 - p2);
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
