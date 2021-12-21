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

#include <iomanip> // std::setprecision
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

int main() {
    // Safe set:
    MatrixXd Gx(6, 3);
    Gx << Eigen::MatrixXd::Identity(3, 3), -Eigen::MatrixXd::Identity(3, 3);
    VectorXd Fx(6, 1);
    Fx << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;
    cis2m::HPolyhedron SafeSet(Gx, Fx);

    // Disturbance set:
    MatrixXd Gw(6, 3);
    Gw << MatrixXd::Identity(3, 3), -MatrixXd::Identity(3, 3);
    VectorXd Fw(6, 1);
    Fw << 0.01, 0.01, 0.01, 0.01, 0.01, 0.01;
    cis2m::HPolyhedron W(Gw, Fw);

    // System Dynamics
    MatrixXd A(3, 3);
    A << 0, 1.0, 0, 0, 0, 1.0, 0, 0, 0;
    MatrixXd B(3, 1);
    B << 0, 0, 1.0;
    MatrixXd E = MatrixXd::Identity(3, 3);

    // RCIS
    cis2m::CISGenerator cisg(2, 0, A, B, E);
    // Add disturbance information
    cisg.AddDisturbanceSet(W);

    // Computing the CIS given a Safe set
    cisg.computeImplicitCIS(SafeSet, 2, 0);
    cis2m::HPolyhedron CIS = cisg.Fetch_CIS();

    //   std::cout << "A: " << std::endl << CIS.Ai() << std::endl;
    //   std::cout << "b: " << std::endl << CIS.bi() << std::endl;

    //   std::cout << "A: " << std::endl << CIS.Ai() << std::endl;

    //   std::cout << "CIS Size: " << CIS.Ai().rows() << " X " <<
    //   CIS.Ai().cols()
    //             << std::endl;
    //   std::cout << "CIS: " << std::endl << cisg.Fetch_A_State() << std::endl;

    cis2m::HPolyhedron CISTT(CIS.Ai(), CIS.bi());
    std::cout << CISTT.isEmpty() << std::endl;

    //   MatrixXd Alifted = cisg.Fetch_A_lifted();
    //   std::cout << "Alifted: " << std::endl << Alifted << std::endl;

    //   std::cout << CISTT.isPositivelyInvariant(Alifted) << std::endl;
    //   std::cout << CISTT.isPositivelyInvariant(A, E, W) << std::endl;

    return 0;
}
