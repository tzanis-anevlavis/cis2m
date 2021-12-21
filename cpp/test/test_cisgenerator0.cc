/*****************************************************************************
 * Copyright (c) 2021, University of California Los Angeles.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * * Neither the name of the copyright holder nor the names of its contributors
 *   may be used to endorse or promote products derived from this software
 *   without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *****************************************************************************/

#include "cis_generator.hpp"

#include <chrono>
#include <fstream>
#include <iostream>

// #include "gurobiRoutines.h"

using namespace std::chrono;

using Eigen::MatrixXd;

int main() {
    // System Dynamics
    MatrixXd At(9, 9);
    At << 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0,
        0, 0, 0, 0, 0, 0, 0.1200, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0.1200, 0, 0,
        1.0, 0, 0, 0, 0, 0, 0, 0.1200, 0, 0, 1.0, 0, 0, 0, 0.0072, 0, 0, 0.1200,
        0, 0, 1.0, 0, 0, 0, 0.0072, 0, 0, 0.1200, 0, 0, 1.0, 0, 0, 0, 0.0072, 0,
        0, 0.1200, 0, 0, 1.0;
    MatrixXd A = At.transpose();

    MatrixXd B(9, 3);
    B << 0.0003, 0, 0, 0, 0.0003, 0, 0, 0, 0.0003, 0.0072, 0, 0, 0, 0.0072, 0,
        0, 0, 0.0072, 0.1200, 0, 0, 0, 0.1200, 0, 0, 0, 0.1200;

    // Safe Set
    MatrixXd Gx(30, 9);
    Gx << 1.0, 0, 0, 0, 0, 0, 0, 0, 0, -1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0,
        0, 0, 0, 0, 0, 0, 0, -1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0,
        0, 0, 0, 0, -1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0,
        0, -1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, -1.0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, -1.0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, -1.0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, -1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1.0, 0, 0, 0, 0, 0, 0, 0, 0, -1.0, 0.0, 0.0, 0.3162, 0, 0, 0, 0, 0, 0,
        0.0, 0.0, 0.3162, 0, 0, 0, 0, 0, 0, 0, 0.3162, 0.0, 0, 0, 0, 0, 0, 0,
        0.0, 0.3162, 0.0, 0, 0, 0, 0, 0, 0, -0.3162, 0.0, 0.0, 0, 0, 0, 0, 0, 0,
        -0.3162, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0.3162, 0.0, 0.0, 0, 0, 0, 0, 0, 0,
        0.3162, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0.0, -0.5547, 0.0, 0, 0, 0, 0, 0, 0,
        0.0, -0.5547, 0.0, 0, 0, 0, 0, 0, 0, 0.0, 0.0, -0.3162, 0, 0, 0, 0, 0,
        0, 0.0, 0.0, -0.3162, 0, 0, 0, 0, 0, 0;
    MatrixXd Fx(30, 1);
    Fx << 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 0.866, 0.866, 0.866, 0.866, 0.866,
        0.866, 2.8319, 2.8319, 2.8319, 2.8319, 2.8319, 2.8319, 0.9487, 0.9487,
        0.9487, 0.9487, 0.9487, 0.9487, 0.9487, 0.9487, -0.8321, -0.8321,
        0.9487, 0.9487;
    cis2m::HPolyhedron SafeSet(Gx, Fx);

    // // Print out input data:
    // std::cout << "Original A: " << std::endl << A << std::endl;
    // std::cout << "Original B: " << std::endl << B << std::endl;
    // std::cout << "Original State Gx = " << std::endl;
    // std::cout << Gx << std::endl;
    // std::cout << "Original State Fx = " << std::endl;
    // std::cout << Fx.transpose() << std::endl;

    // ==============================================================================================
    // CIS (No disturbance)
    cis2m::CISGenerator cisg(2, 0, A, B);

    std::cout << "Computing the CIS..." << std::endl;
    // Computing the CIS given a Safe set
    cisg.computeImplicitCIS(SafeSet, 2, 0);
    cis2m::HPolyhedron CIS = cisg.Fetch_CIS();

    // std::cout << "CIS Size: " << CIS.Ai().rows() << " X " << CIS.Ai().cols()
    // << std::endl;
    //   std::cout << "CIS: " << std::endl << cisg.Fetch_A_State() << std::endl;
    //   std::cout << "CIS A:  " << std::endl << CIS.Ai() << std::endl;
    // std::cout << "CIS b:  " << std::endl << CIS.bi() << std::endl;
    std::cout << "CIS Size: " << CIS.Ai().rows() << " X " << CIS.Ai().cols()
              << std::endl;

    //   // Save in file:
    //   std::ofstream myfile;
    //   myfile.open("cis_out_nodist.csv");
    //   myfile << "Original A: " << std::endl << A << std::endl;
    //   myfile << "Original B: " << std::endl << B << std::endl;
    //   myfile << "Original State Gx: " << std::endl << Gx << std::endl;
    //   myfile << "Original State Fx: " << std::endl << Fx << std::endl;
    //   myfile << "CIS A:  " << std::endl << CIS.Ai() << std::endl;
    //   myfile << "CIS b:  " << std::endl << CIS.bi() << std::endl;
    //   myfile.close();

    // Check emptiness:
    bool isEmpty_Safe = SafeSet.isEmpty();
    std::cout << "SafeSet Not-Empty: " << std::endl
              << isEmpty_Safe << std::endl;

    //   bool isEmpty_CIS = CIS.isEmpty();
    //   std::cout << "CIS Empty: " << std::endl << isEmpty_CIS << std::endl;

    cis2m::HPolyhedron CISTEST(
        CIS.Ai(),
        CIS.bi()); //  - 20 * Eigen::VectorXd::Ones(270));
    bool isEmpty_CIST = CISTEST.isEmpty();
    std::cout << "CISTEST Empty: " << std::endl << isEmpty_CIST << std::endl;

    //   // Use Gurobi:
    //   // RCIS inequalitites:
    //   Eigen::MatrixXd Aineq = SafeSet.Ai();
    //   Eigen::VectorXd bineq = SafeSet.bi();
    //   // Original cost: || u - udes ||^2.
    //   MatrixXd H =
    //       Eigen::MatrixXd::Zero(SafeSet.GetSpaceDim(),
    //       SafeSet.GetSpaceDim());
    //   VectorXd c = Eigen::VectorXd::Zero(SafeSet.GetSpaceDim());
    //   // Call solver:
    //   opt_result res = solveGurobi(H, c, Aineq, bineq, 0);
    //   std::cout << "SafeSet Not-Empty: " << std::endl << res.solved <<
    //   std::endl;

    bool isInv_CIST = CISTEST.isPositivelyInvariant(cisg.Fetch_A_lifted());
    std::cout << "CISTEST Invariant: " << std::endl << isInv_CIST << std::endl;

    return 0;
}
