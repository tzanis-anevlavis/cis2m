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

/* This test constructs an input system (A,B,E), state safe set (Gx,Fx), and
 * transforms it to Brunovsky form (Abru, Bbru, Ebru) and (Gxbru, Fx) */
#include "cis_generator.hpp"

#include <chrono>
#include <fstream>
#include <iomanip> // std::setprecision
#include <iostream>

using namespace std::chrono;

using Eigen::MatrixXd;
using Eigen::VectorXd;

int main() {
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

  // Irredundant safe set:
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
      0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, -0.55470,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
  VectorXd Fx(18, 1);
  Fx << 3.0000, 3.0000, 3.0000, 3.0000, 3.0000, 0.8660, 0.8660, 0.8660, 0.8660,
      0.8660, 0.8660, 2.8319, 2.8319, 2.8319, 2.8319, 2.8319, 2.8319, -0.8321;
  cis2m::HPolyhedron SafeSet(Gx, Fx);

  //   std::cout << "========= Original Space =========" << std::endl;
  //   std::cout << "A: " << std::endl
  //             << std::scientific << std::setprecision(20) << A << std::endl;
  //   std::cout << "B: " << std::endl
  //             << std::scientific << std::setprecision(20) << B << std::endl;
  //   std::cout << "E: " << std::endl
  //             << std::scientific << std::setprecision(20) << E << std::endl;
  //   std::cout << "Gx: " << std::endl
  //             << std::scientific << std::setprecision(20) << Gx << std::endl;

  // ==============================================================================================
  //  Transformation to Brunovsky space
  // ==============================================================================================
  cis2m::CISGenerator cisg(2, 0, A, B);
  cis2m::BrunovskyForm *bru_form = cisg.getBrunovskyForm();

  std::pair<MatrixXd, MatrixXd> sysBru = bru_form->GetSystem();
  MatrixXd T = bru_form->GetTransformationMatrix();

  std::cout << "========= Brunovsky Space =========" << std::endl;
  //   std::cout << "Abru: " << std::endl
  //             << std::scientific << std::setprecision(20) << sysBru.first
  //             << std::endl;
  //   std::cout << "Bbru: " << std::endl
  //             << std::scientific << std::setprecision(20) << sysBru.second
  //             << std::endl;
  //   std::cout << "Ebru: " << std::endl
  //             << std::scientific << std::setprecision(20) << T * E <<
  //             std::endl;
  //   std::cout << "T: " << std::endl
  //             << std::scientific << std::setprecision(20) << T << std::endl;
  //   std::cout << "Gxbru: " << std::endl
  //             << std::scientific << std::setprecision(20) << Gx * T.inverse()
  //             << std::endl;

  MatrixXd A_BF_extended(10, 3);
  A_BF_extended.block(0, 0, 3, 3) << Eigen::MatrixXd::Identity(3, 3);
  std::cout << "A_BF_extended: " << std::endl << A_BF_extended << std::endl;

  return 0;
}
