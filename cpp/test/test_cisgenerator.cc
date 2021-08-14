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

/*
#include <drake/systems/framework/framework_common.h>
#include <drake/systems/framework/vector_base.h>

#include "drake/geometry/optimization/iris.h"

#include "drake/geometry/optimization/point.h"
#include "drake/geometry/optimization/vpolytope.h"
#include "drake/geometry/scene_graph.h"
#include "drake/math/rigid_transform.h"
*/
#include "cis_generator.hpp"

#include <chrono>
#include <iostream>

using namespace std::chrono;
/*
using namespace drake;
using namespace geometry;
using namespace optimization;
*/

using Eigen::MatrixXd;


int main() {
	auto start = high_resolution_clock::now();

	/*
	ConvexSets obstacles;
	const HPolyhedron domain = HPolyhedron::MakeBox(Eigen::Vector3d(-1.0, -1.0, 0), Eigen::Vector3d(1.0, 1.0, 1.0));

	// Obstacle in the center
	obstacles.emplace_back(VPolytope::MakeBox(Eigen::Vector3d(-0.2, -0.2, 0.0), Eigen::Vector3d(0.2, 0.2, 0.5)));
	obstacles.emplace_back(VPolytope::MakeBox(Eigen::Vector3d(-0.5, -0.5, 0.0), Eigen::Vector3d(-0.3, -0.3, 0.5)));


	std::vector<Eigen::Vector3d> points;
	points.push_back(Eigen::Vector3d(-0.7, 0.0, 0.5));
	points.push_back(Eigen::Vector3d(0.0, -0.7, 0.5));
	points.push_back(Eigen::Vector3d(0.7, 0.0, 0.5));
	points.push_back(Eigen::Vector3d(0.0, 0.7, 0.5));

	std::vector<HPolyhedron> volumes;

//	std::cout << "Computing Volumes..." << std::endl;
	for (auto& el : points) {
		HPolyhedron hp = Iris(obstacles, el, domain);
		volumes.push_back(hp);
	}
//	std::cout << "Done!" << std::endl;
//	*/




	// System Dynamics 
	MatrixXd At(9, 9);
	At << 1.0000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0000, 0, 0, 0, 0, 0, 0, 0.1200, 0, 0, 1.0000, 0, 0, 0, 0, 0, 0, 0.1200, 0, 0, 1.0000, 0, 0, 0, 0, 0, 0, 0.1200, 0, 0, 1.0000, 0, 0, 0, 0.0072, 0, 0, 0.1200, 0, 0, 1.0000, 0, 0, 0, 0.0072, 0, 0, 0.1200, 0, 0, 1.0000, 0, 0, 0, 0.0072, 0, 0, 0.1200, 0, 0, 1.0000;
	MatrixXd A = At.transpose();

	MatrixXd B(9,3);
	B << 0.0003, 0, 0, 0, 0.0003, 0, 0, 0, 0.0003, 0.0072, 0, 0, 0, 0.0072, 0, 0, 0, 0.0072, 0.1200, 0, 0, 0, 0.1200, 0, 0, 0, 0.1200;

	MatrixXd Gx(30, 9);
	Gx <<  1.0000, 0, 0, 0, 0, 0, 0, 0, 0, -1.0000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0000, 0, 0, 0, 0, 0, 0, 0, 0, -1.0000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0000, 0, 0, 0, 0, 0, 0, 0, 0, -1.0000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0000, 0, 0, 0, 0, 0, 0, 0, 0, -1.0000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0000, 0, 0, 0, 0, 0, 0, 0, 0, -1.0000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0000, 0, 0, 0, 0, 0, 0, 0, 0, -1.0000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0000, 0, 0, 0, 0, 0, 0, 0, 0, -1.0000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0000, 0, 0, 0, 0, 0, 0, 0, 0, -1.0000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0000, 0, 0, 0, 0, 0, 0, 0, 0, -1.0000, 0.0000, 0.0000, 0.3162, 0, 0, 0, 0, 0, 0, 0.0000, 0.0000, 0.3162, 0, 0, 0, 0, 0, 0, 0, 0.3162, 0.0000, 0, 0, 0, 0, 0, 0, 0.0, 0.3162, 0.0000, 0, 0, 0, 0, 0, 0, -0.3162, 0.0, 0.0, 0, 0, 0, 0, 0, 0, -0.3162, 0.0, 0.0000, 0, 0, 0, 0, 0, 0, 0.3162, 0.0, 0.0000, 0, 0, 0, 0, 0, 0, 0.3162, 0.0, 0.0000, 0, 0, 0, 0, 0, 0, 0.0, -0.5547, 0.0, 0, 0, 0, 0, 0, 0, 0.0000, -0.5547, 0.0000, 0, 0, 0, 0, 0, 0, 0.0, 0.0000, -0.3162, 0, 0, 0, 0, 0, 0, 0.0000, 0.0, -0.3162, 0, 0, 0, 0, 0, 0; 

	MatrixXd Fx(30, 1);
	Fx << 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 0.866, 0.866, 0.866, 0.866, 0.866, 0.866, 2.8319, 2.8319, 2.8319, 2.8319, 2.8319, 2.8319, 0.9487, 0.9487, 0.9487, 0.9487, 0.9487, 0.9487, 0.9487, 0.9487, -0.8321, -0.8321, 0.9487, 0.9487;


	std::cout << "Original A: " << std::endl << A << std::endl;
	std::cout << "Original B: " << std::endl << B << std::endl;

	std::cout << "Original State Gx = " << std::endl;
	std::cout << Gx << std::endl;
	std::cout << "Original State Fx = " << std::endl;
	std::cout << Fx << std::endl;


	// ================================
	// Disturbance 
	MatrixXd Wa(MatrixXd::Zero(2 * 9, 9));
	Wa.block(0, 0, 9, 9) = MatrixXd::Identity(9, 9);
	Wa.block(9, 0, 9, 9) = -MatrixXd::Identity(9, 9);
	MatrixXd Wb(MatrixXd::Zero(2 * 9, 1));
	for (int i = 0; i < 2 * 9; i++) {
		Wb(i) = 0.01;
	}
	cis2m::HPolyhedron PP(Wa, Wb);

	// ==============================================================================================
	// CIS
	std::cout << "TESTING CISGenerator..." << std::endl;
	std::cout << "Creating class..." << std::endl;
	//cis2m::CISGenerator cisg(A, B, MatrixXd::Identity(9, 9));
	cis2m::CISGenerator cisg(A, B);


	cisg.AddDisturbanceSet(PP);
	cis2m::HPolyhedron SafeSet(Gx, Fx);

	std::cout << "Computing CIS..." << std::endl;
	cis2m::HPolyhedron CIS = cisg.computeCIS(SafeSet, 6, 0);

	std::cout << "CIS Size: " << CIS.Ai().rows() << " X " << CIS.Ai().cols() << std::endl;

	/*
	for (int i = 0; i < CIS.Ai().cols() - 9; i++) {
			int deleted_col = CIS.Ai().cols() - 1;
			CIS = CIS.Projection(deleted_col);
			std::cout << "CIS Size: " << CIS.Ai().rows() << " X " << CIS.Ai().cols() << std::endl;
	}
	*/

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	std::cout << "Time: " << duration.count() << std::endl;


	//bool isContained = region.PointInSet(Eigen::Vector2d(0.3, 0));
	// Setting the option keeps the sample in the set (but misses above the box).
	/*
	IrisOptions options;
	options.require_sample_point_is_contained = true;
	region = Iris(obstacles, sample, domain, options);
	*/

	return 0;
}
