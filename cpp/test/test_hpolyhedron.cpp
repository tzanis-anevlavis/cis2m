#include "hpolyhedron.hpp"
#include <Eigen/Dense>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

int main() {
	MatrixXd A(4, 2);
	A << 1, 0, 0, 1, -1, 0, 0, -1;
	VectorXd b(4);
	b << 4, 4, 4, 4;
	cis2m::HPolyhedron p1(A, b);

	MatrixXd A2(4, 2);
	A2 << 1, 0, 0, 1, -1, 0, 0, -1;
	VectorXd b2(4);
	b2 << 1, 1, 1, 1;
	cis2m::HPolyhedron p2(A2, b2);

	p1 = p1 - p2;

	bool flag = 0;

	VectorXd newbi = p1.bi();;
	for (int i = 0; i < 4; i++) {
		if (newbi(i) != 3)  {
			flag = 1;
			break;
		}
	}

	std::cout << std::endl << p1.Ai() << std::endl;
	std::cout << std::endl << p1.bi() << std::endl;

	return flag;
}
