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
 **************************************************************2***************/

/* This test constructs an input system (A,B,E), state safe set (Gx,Fx), and
 * transforms it to Brunovsky form (Abru, Bbru, Ebru) and (Gxbru, Fx) */
#include "cis_generator.hpp"

#include <chrono>
#include <fstream>
#include <iomanip> // std::setprecision
#include <iostream>

#include <assert.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

int main() {
    // ==============================================================================================
    // System Dynamics
    // ==============================================================================================
    MatrixXd A(3, 3);
    A << 1.0, 0.12, 0.0072, 0, 1.0, 0.12, 0, 0, 1.0;
    MatrixXd B(3, 1);
    B << 0.000288, 0.0072, 0.12;
    MatrixXd E = MatrixXd::Identity(3, 3);

    // ==============================================================================================
    //  Transformation to Brunovsky space
    // ==============================================================================================
    cis2m::BrunovskyForm bru_form(A, B, E);
    std::pair<MatrixXd, MatrixXd> sysBru = bru_form.GetBrunovskySystem();

    // ==============================================================================================
    // Expected Brunovsky System
    // ==============================================================================================
    MatrixXd Abru(3, 3);
    Abru << 0, 1.0, 0, 0, 0, 1.0, 0, 0, 0;
    MatrixXd Bbru(3, 1);
    Bbru << 0, 0, 1;
    MatrixXd Ebru(3, 3);
    Ebru << 578.7037037037037037, -69.4444444444444444, 2.7777777777777778,
        578.7037037037037037, 0, -1.3888888888888888, 578.7037037037037037,
        69.4444444444444444, 2.7777777777777778;
    MatrixXd Tbru = bru_form.GetTransformationMatrix();

    // ==============================================================================================
    // Check equalities
    // ==============================================================================================
    double tol = 1e-13;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            assert(abs(Abru(i, j) - sysBru.first(i, j)) < tol);
            assert(abs((Tbru * E)(i, j) - Ebru(i, j)) < tol);
        }
        assert(abs(Bbru(i) - sysBru.second(i)) < tol);
    }
}
