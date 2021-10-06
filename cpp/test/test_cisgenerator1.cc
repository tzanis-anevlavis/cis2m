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
  //   MatrixXd Gx(30, 9);
  //   Gx << 1.0, 0, 0, 0, 0, 0, 0, 0, 0, -1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0,
  //   0, 0,
  //       0, 0, 0, 0, 0, 0, -1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0,
  //       0, 0, 0, 0, -1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0,
  //       0, 0, -1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0,
  //       -1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, -1.0,
  //       0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, -1.0, 0, 0,
  //       0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, -1.0, 0, 0, 0, 0,
  //       0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, -1.0, 0.0, 0.0, 0.3162,
  //       0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.3162, 0, 0, 0, 0, 0, 0, 0, 0.3162, 0.0,
  //       0, 0, 0, 0, 0, 0, 0.0, 0.3162, 0.0, 0, 0, 0, 0, 0, 0, -0.3162, 0.0,
  //       0.0, 0, 0, 0, 0, 0, 0, -0.3162, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0.3162,
  //       0.0, 0.0, 0, 0, 0, 0, 0, 0, 0.3162, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0.0,
  //       -0.5547, 0.0, 0, 0, 0, 0, 0, 0, 0.0, -0.5547, 0.0, 0, 0, 0, 0, 0, 0,
  //       0.0, 0.0, -0.3162, 0, 0, 0, 0, 0, 0, 0.0, 0.0, -0.3162, 0, 0, 0, 0,
  //       0, 0;
  //   MatrixXd Fx(30, 1);
  //   Fx << 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 0.866, 0.866, 0.866, 0.866, 0.866,
  //   0.866,
  //       2.8319, 2.8319, 2.8319, 2.8319, 2.8319, 2.8319, 0.9487, 0.9487,
  //       0.9487, 0.9487, 0.9487, 0.9487, 0.9487, 0.9487, -0.8321, -0.8321,
  //       0.9487, 0.9487;

  // Safe set (Irredundant):
  MatrixXd Gx(18, 9);
  Gx << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0,
      -0.554700000000000, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
  VectorXd Fx(18, 1);
  Fx << 3.0000, 3.0000, 3.0000, 3.0000, 3.0000, 0.8660, 0.8660, 0.8660, 0.8660,
      0.8660, 0.8660, 2.8319, 2.8319, 2.8319, 2.8319, 2.8319, 2.8319, -0.8321;

  cis2m::HPolyhedron SafeSet(Gx, Fx);

  // Disturbance set:
  MatrixXd Gw(18, 9);
  Gw << MatrixXd::Identity(9, 9), -MatrixXd::Identity(9, 9);
  VectorXd Fw(18, 1);
  Fw << 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.05, 0.05, 0.05, 0.01, 0.01, 0.01,
      0.01, 0.01, 0.01, 0.05, 0.05, 0.05;
  cis2m::HPolyhedron W(Gw, Fw);

  // System Dynamics
  MatrixXd At(9, 9);
  At << 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0,
      0, 0, 0, 0, 0, 0.1200, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0.1200, 0, 0, 1.0, 0,
      0, 0, 0, 0, 0, 0.1200, 0, 0, 1.0, 0, 0, 0, 0.0072, 0, 0, 0.1200, 0, 0,
      1.0, 0, 0, 0, 0.0072, 0, 0, 0.1200, 0, 0, 1.0, 0, 0, 0, 0.0072, 0, 0,
      0.1200, 0, 0, 1.0;
  MatrixXd A = At.transpose();
  MatrixXd B(9, 3);
  B << 0.000288, 0, 0, 0, 0.000288, 0, 0, 0, 0.000288, 0.0072, 0, 0, 0, 0.0072,
      0, 0, 0, 0.0072, 0.1200, 0, 0, 0, 0.1200, 0, 0, 0, 0.1200;
  MatrixXd E = MatrixXd::Identity(9, 9);

  //   std::cout << "A: " << std::endl << A << std::endl;
  //   std::cout << "B: " << std::endl << B << std::endl;
  //   std::cout << "E: " << std::endl << E << std::endl;

  // RCIS
  cis2m::CISGenerator cisg(2, 0, A, B, E);
  // Add disturbance information
  cisg.AddDisturbanceSet(W);

  // CIS
  //   cis2m::CISGenerator cisg(A, B);

  // Computing the CIS given a Safe set
  cisg.computeCIS(SafeSet, 2, 0);
  cis2m::HPolyhedron CIS = cisg.Fetch_CIS();

  //   std::cout << "A: " << std::endl
  //             << std::scientific << std::setprecision(20) << CIS.Ai()
  //             << std::endl;
  //   std::cout << "b: " << std::endl
  //             << std::scientific << std::setprecision(20) << CIS.bi()
  //             << std::endl;

  //   std::cout << "A: " << std::endl << CIS.Ai() << std::endl;

  //   std::cout << "CIS Size: " << CIS.Ai().rows() << " X " << CIS.Ai().cols()
  //             << std::endl;
  //   std::cout << "CIS: " << std::endl << cisg.Fetch_A_State() << std::endl;

  //   cis2m::HPolyhedron CISTT(CIS.Ai(), CIS.bi());
  //   std::cout << CISTT.isEmpty() << std::endl;

  //   MatrixXd Alifted = cisg.Fetch_A_lifted();
  //   std::cout << "Alifted: " << std::endl << Alifted << std::endl;

  //   std::cout << CISTT.isPositivelyInvariant(Alifted) << std::endl;
  //   std::cout << CISTT.isPositivelyInvariant(A, E, W) << std::endl;

  return 0;
}
