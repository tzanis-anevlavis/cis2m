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
  At << 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0,
      0, 0, 0, 0, 0, 0.1200, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0.1200, 0, 0, 1.0, 0,
      0, 0, 0, 0, 0, 0.1200, 0, 0, 1.0, 0, 0, 0, 0.0072, 0, 0, 0.1200, 0, 0,
      1.0, 0, 0, 0, 0.0072, 0, 0, 0.1200, 0, 0, 1.0, 0, 0, 0, 0.0072, 0, 0,
      0.1200, 0, 0, 1.0;
  MatrixXd A = At.transpose();

  MatrixXd B(9, 3);
  B << 0.0003, 0, 0, 0, 0.0003, 0, 0, 0, 0.0003, 0.0072, 0, 0, 0, 0.0072, 0, 0,
      0, 0.0072, 0.1200, 0, 0, 0, 0.1200, 0, 0, 0, 0.1200;

  // ==============================================================================================
  // CIS (No disturbance)
  cis2m::CISGenerator cisg(A, B);

  cis2m::BrunovskyForm *bru_form = cisg.getBrunovskyForm();

  std::pair<MatrixXd, MatrixXd> dynSys = bru_form->GetDynSystem();

  std::cout << "A: " << std::endl << dynSys.first << std::endl;
  std::cout << "B: " << std::endl << dynSys.second << std::endl;

  return 0;
}
