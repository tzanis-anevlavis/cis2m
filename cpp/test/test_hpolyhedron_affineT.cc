#include "hpolyhedron.hpp"
#include <Eigen/Dense>
#include <iostream>

#include <assert.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

int main() {
    double tol = 1e-13;
    MatrixXd A(4, 2);
    A << 1, 0, 1, 3, -1, 0, -1, -1;
    VectorXd b(4);
    b << 4, 4, 4, 4;
    cis2m::HPolyhedron P1(A, b);

    // Affine Transformation (Full Rank)
    MatrixXd R1(2, 2);
    R1 << 1, 2, 1, 0;
    cis2m::HPolyhedron P2 = P1.affineT(R1);
    // Expected result
    MatrixXd expectedA2(4, 2);
    expectedA2 << 0, 1, 1.5, -0.5, 0, -1, -0.5, -0.5;
    VectorXd expectedb2(4);
    expectedb2 << 4, 4, 4, 4;
    // Check
    for (int i = 0; i < P2.GetNumInequalities(); i++) {
        for (int j = 0; j < P2.GetSpaceDim(); j++) {
            assert(abs(P2.Ai()(i, j) - expectedA2(i, j)) < tol);
        }
        assert(abs(P2.bi()(i) - expectedb2(i)) < tol);
    }

    // Affine Transformation (Not full rank)
    MatrixXd R2(2, 2);
    R2 << 1, 2, 0, 0;
    cis2m::HPolyhedron P3 = P1.affineT(R1);
    // Expected result
    MatrixXd expectedA3(3, 2);
    expectedA3 << 3.0, 0, 0, 0, -1.0, 0.0;
    VectorXd expectedb3(3);
    expectedb3 << 12, 8, 12;
    // Check
    for (int i = 0; i < P3.GetNumInequalities(); i++) {
        for (int j = 0; j < P3.GetSpaceDim(); j++) {
            assert(abs(P3.Ai()(i, j) - expectedA3(i, j)) < tol);
        }
        assert(abs(P3.bi()(i) - expectedb3(i)) < tol);
    }

    return 0;
}
