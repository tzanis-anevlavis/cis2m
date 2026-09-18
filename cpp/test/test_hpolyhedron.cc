/*****************************************************************************
 % Copyright (C) 2026, T.Anevlavis.
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

#include "hpolyhedron.hpp"

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <utility>

namespace {

using cis2m::HPolyhedron;
using Eigen::Index;
constexpr double kTol = 1e-7;
const double kInf = std::numeric_limits<double>::infinity();
const double kNaN = std::numeric_limits<double>::quiet_NaN();

VectorXd Vec(std::initializer_list<double> values) {
    VectorXd result(values.size());
    Index i = 0;
    for (double value : values) {
        result(i++) = value;
    }
    return result;
}

HPolyhedron Box(const VectorXd& lower, const VectorXd& upper) {
    const Index n = lower.size();
    MatrixXd A(2 * n, n);
    A << MatrixXd::Identity(n, n), -MatrixXd::Identity(n, n);
    VectorXd b(2 * n);
    b << upper, -lower;
    return HPolyhedron(A, b);
}

HPolyhedron UnitBox(Index n) {
    return Box(-VectorXd::Ones(n), VectorXd::Ones(n));
}

HPolyhedron Point(const VectorXd& x) {
    return HPolyhedron(MatrixXd(0, x.size()), VectorXd(0), MatrixXd::Identity(x.size(), x.size()), x);
}

HPolyhedron Simplex(Index n) {
    MatrixXd A(n + 1, n);
    A << -MatrixXd::Identity(n, n), MatrixXd::Ones(1, n);
    VectorXd b = VectorXd::Zero(n + 1);
    b(n) = 1;
    return HPolyhedron(A, b);
}

HPolyhedron LegacyPolygon() {
    MatrixXd A(4, 2);
    A << 1, 0, 1, 3, -1, 0, -1, -1;
    return HPolyhedron(A, VectorXd::Constant(4, 4));
}

MatrixXd RandomMatrix(std::mt19937& rng, Index rows, Index cols) {
    std::uniform_real_distribution<double> dist(-1, 1);
    MatrixXd A(rows, cols);
    for (Index i = 0; i < rows; ++i) {
        for (Index j = 0; j < cols; ++j) {
            A(i, j) = dist(rng);
        }
    }
    return A;
}

MatrixXd RandomInvertibleMap(std::mt19937& rng, Index n) {
    // Strict diagonal dominance keeps the map invertible and reasonably
    // conditioned while retaining dense, coupled coordinates.
    MatrixXd T = 0.05 * RandomMatrix(rng, n, n);
    T.diagonal().array() += 1.0;
    return T;
}

HPolyhedron TransformedUnitBox(const VectorXd& center, const MatrixXd& T) {
    // x = center + T*z, |z| <= 1, represented directly in H-form.
    const MatrixXd inverse = T.inverse();
    MatrixXd A(2 * center.size(), center.size());
    A << inverse, -inverse;
    const VectorXd offset = inverse * center;
    VectorXd b(2 * center.size());
    b << VectorXd::Ones(center.size()) + offset,
         VectorXd::Ones(center.size()) - offset;
    return HPolyhedron(A, b);
}

HPolyhedron TransformedSimplex(const VectorXd& center, const MatrixXd& T) {
    // x = center + T*z, where z >= 0 and sum(z) <= 1.
    const Index n = center.size();
    const MatrixXd inverse = T.inverse();
    MatrixXd A(n + 1, n);
    A.topRows(n) = -inverse;
    A.bottomRows(1) = MatrixXd::Ones(1, n) * inverse;
    VectorXd b(n + 1);
    b.head(n) = -inverse * center;
    b(n) = 1.0 + (inverse * center).sum();
    return HPolyhedron(A, b);
}

HPolyhedron AffineCube(
    const VectorXd& center,
    const MatrixXd& coordinate_system,
    Index intrinsic_dimension) {
    // The first intrinsic_dimension transformed coordinates range over a unit
    // box. The remaining coordinates are fixed, giving exact equalities.
    const Index n = center.size();
    const Index k = intrinsic_dimension;
    const MatrixXd inverse = coordinate_system.inverse();
    MatrixXd Ai(2 * k, n);
    Ai << inverse.topRows(k), -inverse.topRows(k);
    const VectorXd offset = inverse.topRows(k) * center;
    VectorXd bi(2 * k);
    bi << VectorXd::Ones(k) + offset, VectorXd::Ones(k) - offset;
    const MatrixXd Ae = inverse.bottomRows(n - k);
    const VectorXd be = Ae * center;
    return HPolyhedron(Ai, bi, Ae, be);
}

void ExpectVectorNear(const VectorXd& actual, const VectorXd& expected, double tol = kTol) {
    ASSERT_EQ(actual.size(), expected.size());
    ASSERT_TRUE(actual.allFinite());
    for (Index i = 0; i < actual.size(); ++i) {
        EXPECT_NEAR(actual(i), expected(i), tol) << "row " << i;
    }
}

// Equality of sets does not require identical normal scaling, ordering, or
// redundant rows. Test every defining halfspace in both directions using LPs.
void ExpectSameSet(const HPolyhedron& actual, const HPolyhedron& expected, double tol = kTol) {
    ASSERT_TRUE(actual.IsValid());
    ASSERT_TRUE(expected.IsValid());
    ASSERT_EQ(actual.Dimension(), expected.Dimension());
    ASSERT_EQ(actual.IsFeasible(), expected.IsFeasible());
    if (!actual.IsFeasible()) {
        return;
    }
    for (int side = 0; side < 2; ++side) {
        const auto& P = side == 0 ? actual : expected;
        const auto& Q = side == 0 ? expected : actual;
        MatrixXd normals(Q.NumInequalities() + 2 * Q.NumEqualities(), Q.Dimension());
        normals << Q.Ai(), Q.Ae(), -Q.Ae();
        VectorXd bounds(normals.rows());
        bounds << Q.bi(), Q.be(), -Q.be();
        const VectorXd support = P.ComputeSupport(normals);
        ASSERT_TRUE(support.allFinite()) << "side " << side;
        for (Index i = 0; i < bounds.size(); ++i) {
            EXPECT_LE(support(i), bounds(i) + tol) << "side " << side << ", row " << i;
        }
    }
}

HPolyhedron Image(const HPolyhedron& P, const MatrixXd& T, int method) {
    if (method == 0) {
        return P.AffineTransform(T);
    }
    if (method == 1) {
        return P.AffineTransform_SVD(T);
    }
    return P.AffineTransform_QR(T);
}

TEST(HPolyhedronConstruction, InvalidEmptyAndFullSpaceAreDistinct) {
    HPolyhedron invalid;
    EXPECT_FALSE(invalid.IsValid());
    EXPECT_FALSE(invalid.IsFeasible());
    EXPECT_FALSE(invalid.IsFullSpace());
    EXPECT_EQ(invalid.Dimension(), 0u);
    EXPECT_FALSE(HPolyhedron::FullSpace(0).IsValid());
    EXPECT_FALSE(HPolyhedron::EmptySet(0).IsValid());
    EXPECT_FALSE(HPolyhedron::FullSpace(std::numeric_limits<size_t>::max()).IsValid());
    EXPECT_FALSE(HPolyhedron::EmptySet(std::numeric_limits<size_t>::max()).IsValid());

    const auto full = HPolyhedron::FullSpace(3);
    ASSERT_TRUE(full.IsValid());
    EXPECT_TRUE(full.IsFullSpace());
    EXPECT_TRUE(full.IsFeasible());
    EXPECT_TRUE(full.Contains(Vec({-1e6, 0, 1e6})));
    EXPECT_EQ(full.Ai().cols(), 3);
    EXPECT_EQ(full.Ae().cols(), 3);
    EXPECT_EQ(full.NumInequalities(), 0u);
    EXPECT_EQ(full.NumEqualities(), 0u);

    const auto empty = HPolyhedron::EmptySet(3);
    ASSERT_TRUE(empty.IsValid());
    EXPECT_FALSE(empty.IsFeasible());
    EXPECT_FALSE(empty.IsFullSpace());
    EXPECT_FALSE(empty.Contains(VectorXd::Zero(3)));
}

TEST(HPolyhedronConstruction, EmptyBlocksPreserveAmbientDimension) {
    const auto explicit_full = HPolyhedron(MatrixXd(0, 4), VectorXd(0));
    EXPECT_TRUE(explicit_full.IsFullSpace());
    EXPECT_EQ(explicit_full.Ae().cols(), 4);

    // An omitted block can infer its width from the other block.
    const HPolyhedron equality_only(MatrixXd(), VectorXd(),
                                    MatrixXd::Identity(2, 2), Vec({2, -3}));
    ASSERT_TRUE(equality_only.IsValid());
    EXPECT_EQ(equality_only.Ai().cols(), 2);
    EXPECT_TRUE(equality_only.Contains(Vec({2, -3})));
    const HPolyhedron inequality_only(MatrixXd::Identity(2, 2), Vec({1, 2}),
                                      MatrixXd(), VectorXd());
    EXPECT_TRUE(inequality_only.IsValid());
    EXPECT_EQ(inequality_only.Ae().cols(), 2);
    EXPECT_FALSE(HPolyhedron(MatrixXd(), VectorXd()).IsValid());
}

TEST(HPolyhedronConstruction, RejectsMismatchedRowsAndColumns) {
    const MatrixXd A = MatrixXd::Identity(2, 2);
    const VectorXd b = VectorXd::Ones(2);
    EXPECT_FALSE(HPolyhedron(A, VectorXd::Zero(1)).IsValid());
    EXPECT_FALSE(HPolyhedron(A, b, A, VectorXd::Zero(1)).IsValid());
    EXPECT_FALSE(HPolyhedron(A, b, MatrixXd::Ones(1, 3), Vec({0})).IsValid());
    EXPECT_FALSE(HPolyhedron(MatrixXd(0, 3), VectorXd(0), A, b).IsValid());
    EXPECT_FALSE(HPolyhedron(A, b, MatrixXd(0, 3), VectorXd(0)).IsValid());
    EXPECT_FALSE(HPolyhedron(MatrixXd::Zero(1, 0), Vec({0})).IsValid());
}

TEST(HPolyhedronConstruction, RejectsNonfiniteDataInEveryBlock) {
    for (double invalid : {kNaN, kInf, -kInf}) {
        MatrixXd A = MatrixXd::Identity(2, 2);
        VectorXd b = VectorXd::Ones(2);
        A(0, 0) = invalid;
        EXPECT_FALSE(HPolyhedron(A, b).IsValid());
        EXPECT_FALSE(HPolyhedron(MatrixXd::Identity(2, 2), b, A, b).IsValid());
        A.setIdentity();
        b(0) = invalid;
        EXPECT_FALSE(HPolyhedron(A, b).IsValid());
        EXPECT_FALSE(HPolyhedron(A, VectorXd::Ones(2), A, b).IsValid());
    }
}

TEST(HPolyhedronConstruction, CopyAssignmentAndSelfAssignmentPreserveAllConstraints) {
    const HPolyhedron original = Point(Vec({2, -1, 3}));
    HPolyhedron copy(original);
    HPolyhedron assigned;
    assigned = copy;
    assigned = assigned;
    ExpectSameSet(assigned, original);
    copy = HPolyhedron::FullSpace(3);
    EXPECT_FALSE(assigned.IsFullSpace());
    assigned = HPolyhedron();
    EXPECT_FALSE(assigned.IsValid());
}

TEST(HPolyhedronConstruction, MovingPreservesDestinationAndInvalidatesSource) {
    for (const auto& original : {UnitBox(2), Point(Vec({1, -2})),
                                 HPolyhedron::FullSpace(2), HPolyhedron::EmptySet(2)}) {
        HPolyhedron source = original;
        HPolyhedron destination(std::move(source));
        ExpectSameSet(destination, original);
        EXPECT_FALSE(source.IsValid());
        EXPECT_FALSE(source.IsFullSpace());
        EXPECT_FALSE(source.IsFeasible());
        source = UnitBox(3);
        source = std::move(destination);
        ExpectSameSet(source, original);
        EXPECT_FALSE(destination.IsValid());
        EXPECT_FALSE(destination.IsFullSpace());
        HPolyhedron& alias = source;
        source = std::move(alias);
        ExpectSameSet(source, original);
    }
}

TEST(HPolyhedronConstruction, ShowReportsTheRepresentationWithoutChangingIt) {
    std::ostringstream output;
    auto* previous = std::cout.rdbuf(output.rdbuf());
    const auto P = Point(Vec({1, -2}));
    P.Show();
    std::cout.rdbuf(previous);
    EXPECT_NE(output.str().find("dimension=2"), std::string::npos);
    EXPECT_NE(output.str().find("equalities=2"), std::string::npos);
    EXPECT_TRUE(P.Contains(Vec({1, -2})));
}

TEST(HPolyhedronSupport, ScalingFailuresDoNotBecomeFeasibilityOrSupportCertificates) {
    // Finite input can still exceed double precision during LP equilibration.
    // Failure is distinct from a proven empty set or a finite support value.
    for (const auto& P : {
            HPolyhedron(MatrixXd::Constant(1, 1, 1e300), Vec({1e-300})),
            HPolyhedron(MatrixXd::Constant(1, 1, 1e-300), Vec({1e300})),
            HPolyhedron(MatrixXd::Ones(1, 1), Vec({1e-320})),
            HPolyhedron(MatrixXd::Constant(1, 1, 1e-320), Vec({1}))}) {
        ASSERT_TRUE(P.IsValid());
        EXPECT_TRUE(std::isnan(P.ComputeSupport(MatrixXd::Ones(1, 1))(0)));
        EXPECT_THROW(P.IsFeasible(), std::runtime_error);
        auto unchanged = P;
        unchanged.RemoveRedundantConstraints();
        EXPECT_TRUE(unchanged.Ai().isApprox(P.Ai()));
        EXPECT_TRUE(unchanged.bi().isApprox(P.bi()));
    }
    // Scaling an objective can overflow or underflow even with a usable model.
    EXPECT_TRUE(std::isnan(Box(Vec({-1e12}), Vec({1e12}))
        .ComputeSupport(MatrixXd::Constant(1, 1, 1e300))(0)));
    EXPECT_TRUE(std::isnan(Box(Vec({-1e-12}), Vec({1e-12}))
        .ComputeSupport(MatrixXd::Constant(1, 1, 1e-320))(0)));
}

TEST(HPolyhedronMembership, IncludesBoundaryAndHonorsEqualityTolerance) {
    const HPolyhedron P = UnitBox(2).Intersection(
        HPolyhedron(MatrixXd(0, 2), VectorXd(0), Vec({1, 1}).transpose(), Vec({0.5})));
    EXPECT_TRUE(P.Contains(Vec({1, -0.5})));
    EXPECT_FALSE(P.Contains(Vec({1.01, -0.51})));
    EXPECT_FALSE(P.Contains(Vec({0, 0})));
    EXPECT_TRUE(P.Contains(Vec({0.25, 0.2500001}), 1e-6));
    EXPECT_FALSE(P.Contains(Vec({0.25, 0.2500001}), 1e-9));
    EXPECT_TRUE(Point(Vec({1})).Contains(Vec({1}), 0));
}

TEST(HPolyhedronMembership, RejectsInvalidPointsAndTolerances) {
    for (const auto& P : {HPolyhedron(), UnitBox(2), HPolyhedron::FullSpace(2)}) {
        EXPECT_FALSE(P.Contains(Vec({0})));
        for (double invalid : {kNaN, kInf, -kInf}) {
            EXPECT_FALSE(P.Contains(Vec({invalid, 0})));
        }
        for (double invalid : {-1.0, kNaN, kInf}) {
            EXPECT_FALSE(P.Contains(Vec({0, 0}), invalid));
        }
    }
}

TEST(HPolyhedronMembership, RejectsFiniteInputsWhoseResidualsOverflow) {
    const VectorXd point = Vec({1e308, 0});
    const MatrixXd normal = (MatrixXd(1, 2) << 1e308, 0).finished();
    EXPECT_FALSE(HPolyhedron(normal, Vec({0})).Contains(point));
    EXPECT_FALSE(HPolyhedron(
        MatrixXd(0, 2), VectorXd(0), normal, Vec({0})).Contains(point));
}

TEST(HPolyhedronFeasibility, DetectsInequalityAndEqualityContradictions) {
    EXPECT_TRUE(Simplex(10).IsFeasible());
    EXPECT_TRUE(HPolyhedron(MatrixXd::Identity(3, 3), Vec({-1, 2, 3})).IsFeasible());
    EXPECT_FALSE(Box(Vec({1}), Vec({0})).IsFeasible());
    EXPECT_FALSE(HPolyhedron(MatrixXd(0, 1), VectorXd(0),
        MatrixXd::Ones(2, 1), Vec({1, 2})).IsFeasible());
    // A constant equality must distinguish 0=0 from 0=nonzero.
    EXPECT_TRUE(HPolyhedron(MatrixXd(0, 2), VectorXd(0),
        MatrixXd::Zero(1, 2), Vec({0})).IsFeasible());
    EXPECT_FALSE(HPolyhedron(MatrixXd(0, 2), VectorXd(0),
        MatrixXd::Zero(1, 2), Vec({1e-12})).IsFeasible());
    EXPECT_FALSE(HPolyhedron(MatrixXd::Zero(1, 2), Vec({-1e-12})).IsFeasible());
}

TEST(HPolyhedronSupport, BoxesAndSimplicesHaveAnalyticSupportsUpToTenDimensions) {
    std::mt19937 rng(421);
    for (Index n = 1; n <= 10; ++n) {
        SCOPED_TRACE(n);
        const MatrixXd directions = RandomMatrix(rng, 40, n);
        const VectorXd center = RandomMatrix(rng, n, 1);
        const VectorXd radii = VectorXd::LinSpaced(n, 0.5, 2.0);
        const auto box = Box(center - radii, center + radii);
        ExpectVectorNear(box.ComputeSupport(directions),
                         directions * center + directions.cwiseAbs() * radii);
        VectorXd expected(directions.rows());
        for (Index i = 0; i < directions.rows(); ++i) {
            expected(i) = std::max(0.0, directions.row(i).maxCoeff());
        }
        ExpectVectorNear(Simplex(n).ComputeSupport(directions), expected);
    }
}

TEST(HPolyhedronSupport, OneHundredRandomParallelotopesMatchClosedFormSupport) {
    std::mt19937 rng(98341);
    for (int trial = 0; trial < 100; ++trial) {
        SCOPED_TRACE(trial);
        const Index n = 1 + trial % 10;
        const VectorXd center = RandomMatrix(rng, n, 1);
        const MatrixXd T = RandomInvertibleMap(rng, n);
        const MatrixXd directions = RandomMatrix(rng, 20, n);
        const HPolyhedron P = TransformedUnitBox(center, T);

        // For P = center + T[-1,1]^n, h_P(d) has this closed form.
        const VectorXd expected = directions * center + (directions * T).cwiseAbs().rowwise().sum();
        ExpectVectorNear(P.ComputeSupport(directions), expected, 1e-6);
    }
}

TEST(HPolyhedronLinearProgram, DenseRowScaledAffineCubesMatchClosedFormSupport) {
    std::mt19937 rng(73019);
    std::uniform_real_distribution<double> exponent_dist(-10.0, 10.0);
    for (int trial = 0; trial < 60; trial++) {
        SCOPED_TRACE(trial);
        const Index n = 2 + trial % 9;
        const Index k = 1 + (3 * trial) % (n - 1);
        const VectorXd center = RandomMatrix(rng, n, 1);
        const MatrixXd coordinates = RandomInvertibleMap(rng, n);
        const MatrixXd basis = coordinates.leftCols(k);
        const HPolyhedron affine_cube = AffineCube(center, coordinates, k);

        // Independently scale every row over twenty orders of magnitude. The
        // represented set and its support function must remain unchanged.
        MatrixXd Ai = affine_cube.Ai();
        VectorXd bi = affine_cube.bi();
        MatrixXd Ae = affine_cube.Ae();
        VectorXd be = affine_cube.be();
        for (Index i = 0; i < Ai.rows(); i++) {
            const double scale = std::pow(10.0, exponent_dist(rng));
            Ai.row(i) *= scale;
            bi(i) *= scale;
        }
        for (Index i = 0; i < Ae.rows(); i++) {
            const double scale = std::pow(10.0, exponent_dist(rng));
            Ae.row(i) *= scale;
            be(i) *= scale;
        }
        const HPolyhedron scaled(Ai, bi, Ae, be);
        ASSERT_TRUE(scaled.IsFeasible());
        EXPECT_TRUE(scaled.IsBounded());

        const MatrixXd directions = RandomMatrix(rng, 30, n);
        const VectorXd expected = directions * center
            + (directions * basis).cwiseAbs().rowwise().sum();
        ExpectVectorNear(scaled.ComputeSupport(directions), expected, 2e-6);
    }
}

TEST(HPolyhedronSupport, PointsAffineSlicesAndNegativeSupport) {
    const MatrixXd directions = (MatrixXd(4, 2) << 1, 0, 0, 1, 1, -2, -1, 0).finished();
    ExpectVectorNear(Point(Vec({2, -3})).ComputeSupport(directions),
                     directions * Vec({2, -3}));
    const HPolyhedron slice((MatrixXd(2, 2) << 0, 1, 0, -1).finished(), Vec({3, 1}),
                             Vec({1, 0}).transpose(), Vec({2}));
    ExpectVectorNear(slice.ComputeSupport(directions), Vec({2, 3, 4, -2}));
    // Represent the same affine slice using paired inequalities.
    MatrixXd A(4, 2);
    A << 0, 1, 0, -1, 1, 0, -1, 0;
    ExpectSameSet(slice, HPolyhedron(A, Vec({3, 1, 2, -2})));
}

TEST(HPolyhedronSupport, DistinguishesUnboundedFiniteAndZeroDirections) {
    const HPolyhedron P(Vec({1, 0}).transpose(), Vec({2}));
    MatrixXd directions(5, 2);
    directions << 1, 0, -1, 0, 0, 1, 0, 0, 2, 0;
    const VectorXd support = P.ComputeSupport(directions);
    EXPECT_DOUBLE_EQ(support(0), 2);
    EXPECT_EQ(support(1), kInf);
    EXPECT_EQ(support(2), kInf);
    EXPECT_DOUBLE_EQ(support(3), 0);
    EXPECT_DOUBLE_EQ(support(4), 4); // A bounded objective after an unbounded one.
    const VectorXd full = HPolyhedron::FullSpace(2).ComputeSupport(directions);
    EXPECT_EQ(full(0), kInf);
    EXPECT_EQ(full(3), 0);
}

TEST(HPolyhedronSupport, ScalesDirectionsAndConstraintsWithoutLosingSmallValues) {
    const MatrixXd A = (MatrixXd(4, 2) << 1e-12, 0, -1e12, 0, 0, 1e-8, 0, -1e8).finished();
    const auto P = HPolyhedron(A, Vec({2e-12, -1e12, 3e-8, 1e8}));
    MatrixXd directions(4, 2);
    directions << 1e-12, -2e-12, -1e-10, 0, 1e6, 1, 0, 0;
    const VectorXd actual = P.ComputeSupport(directions);
    EXPECT_NEAR(actual(0), 4e-12, 1e-20);
    EXPECT_NEAR(actual(1), -1e-10, 1e-18);
    EXPECT_NEAR(actual(2), 2000003, kTol);
    EXPECT_EQ(actual(3), 0);
}

TEST(HPolyhedronSupport, InfeasibleAndInvalidResultsAreNaN) {
    const MatrixXd directions = (MatrixXd(3, 1) << 1, 0, -1).finished();
    for (const auto& P : {HPolyhedron(), HPolyhedron::EmptySet(1), Box(Vec({1}), Vec({0}))}) {
        const VectorXd support = P.ComputeSupport(directions);
        ASSERT_EQ(support.size(), 3);
        EXPECT_TRUE(support.array().isNaN().all());
    }
    EXPECT_TRUE(UnitBox(2).ComputeSupport(directions).array().isNaN().all());
    for (double bad : {kNaN, kInf, -kInf}) {
        EXPECT_TRUE(UnitBox(1).ComputeSupport(MatrixXd::Constant(1, 1, bad)).array().isNaN().all());
    }
    EXPECT_EQ(UnitBox(2).ComputeSupport(MatrixXd(0, 2)).size(), 0);
}

TEST(HPolyhedronBoundedness, BoundedUnboundedEmptyAndLowerDimensionalSets) {
    EXPECT_TRUE(UnitBox(10).IsBounded());
    EXPECT_TRUE(Point(Vec({1, 2, 3})).IsBounded());
    EXPECT_TRUE(HPolyhedron::EmptySet(2).IsBounded());
    EXPECT_FALSE(HPolyhedron::FullSpace(2).IsBounded());
    EXPECT_FALSE(HPolyhedron().IsBounded());
    EXPECT_FALSE(HPolyhedron(Vec({1, 0}).transpose(), Vec({1})).IsBounded());
    EXPECT_FALSE(HPolyhedron(MatrixXd(0, 2), VectorXd(0),
        Vec({1, 0}).transpose(), Vec({1})).IsBounded());
}

TEST(HPolyhedronSubset, HandlesEqualityEmptyFullAndUnboundedSets) {
    EXPECT_TRUE(Simplex(2).IsSubsetOf(UnitBox(2)));
    EXPECT_FALSE(UnitBox(2).IsSubsetOf(Simplex(2)));
    EXPECT_TRUE(Point(Vec({0.5, 0.5})).IsSubsetOf(Simplex(2)));
    EXPECT_FALSE(UnitBox(2).IsSubsetOf(Point(Vec({0, 0}))));
    EXPECT_TRUE(HPolyhedron::EmptySet(2).IsSubsetOf(UnitBox(2)));
    EXPECT_FALSE(UnitBox(2).IsSubsetOf(HPolyhedron::EmptySet(2)));
    EXPECT_TRUE(UnitBox(2).IsSubsetOf(HPolyhedron::FullSpace(2)));
    EXPECT_FALSE(HPolyhedron::FullSpace(2).IsSubsetOf(UnitBox(2)));
    const HPolyhedron halfspace(Vec({1, 0}).transpose(), Vec({2}));
    EXPECT_TRUE(halfspace.IsSubsetOf(halfspace));
    EXPECT_FALSE(UnitBox(1).IsSubsetOf(UnitBox(2)));
    EXPECT_FALSE(HPolyhedron().IsSubsetOf(UnitBox(2)));
    EXPECT_FALSE(UnitBox(2).IsSubsetOf(HPolyhedron()));
    for (double bad : {-1.0, kInf, kNaN}) {
        EXPECT_FALSE(halfspace.IsSubsetOf(halfspace, bad));
    }
    EXPECT_TRUE(Point(Vec({1 + 1e-8})).IsSubsetOf(UnitBox(1), 1e-6));
    EXPECT_FALSE(Point(Vec({1 + 1e-8})).IsSubsetOf(UnitBox(1), 1e-10));
}

TEST(HPolyhedronSubset, DenseNestedParallelotopesHaveTheExpectedOrdering) {
    std::mt19937 rng(29611);
    for (int trial = 0; trial < 50; trial++) {
        SCOPED_TRACE(trial);
        const Index n = 1 + trial % 10;
        const VectorXd center = RandomMatrix(rng, n, 1);
        const MatrixXd T = RandomInvertibleMap(rng, n);
        const HPolyhedron inner = TransformedUnitBox(center, T);
        const HPolyhedron outer = TransformedUnitBox(center, 1.5 * T);
        EXPECT_TRUE(inner.IsSubsetOf(outer, 1e-9));
        EXPECT_FALSE(outer.IsSubsetOf(inner, 1e-9));

        // Translation by a sufficiently large generator moves the entire set
        // beyond one of the original parallelotope's facets.
        const HPolyhedron translated = TransformedUnitBox(center + 3.0 * T.col(0), T);
        EXPECT_FALSE(translated.IsSubsetOf(outer, 1e-9));
    }
}

TEST(HPolyhedronNormalization, PreservesSignsZeroBoundsAndAffineEqualities) {
    MatrixXd A(4, 2);
    A << -2, 0, 2, 0, 0, -1, 0, 1;
    VectorXd b = Vec({-2, 6, 0, 4});
    const HPolyhedron original(A, b, Vec({2, 2}).transpose(), Vec({4}));
    HPolyhedron normalized = original;
    normalized.NormalizeConstraints();
    ExpectVectorNear(normalized.bi(), Vec({-1, 1, 0, 1}));
    ExpectVectorNear(normalized.be(), Vec({1}));
    ExpectSameSet(normalized, original);
    auto negative = Point(Vec({-2}));
    negative.NormalizeConstraints(0);
    ExpectVectorNear(negative.be(), Vec({-1}));
    EXPECT_TRUE(negative.Contains(Vec({-2})));
    HPolyhedron zero_rhs = Point(Vec({0}));
    zero_rhs.NormalizeConstraints(0);
    EXPECT_TRUE(zero_rhs.Ae().allFinite());
    EXPECT_TRUE(zero_rhs.Contains(Vec({0})));

    for (double bad : {-1.0, kNaN, kInf}) {
        auto unchanged = original;
        unchanged.NormalizeConstraints(bad);
        EXPECT_TRUE(unchanged.Ai().isApprox(original.Ai(), 0));
    }
    HPolyhedron invalid;
    invalid.NormalizeConstraints();
    EXPECT_FALSE(invalid.IsValid());
    auto full = HPolyhedron::FullSpace(2);
    full.NormalizeConstraints();
    EXPECT_TRUE(full.IsFullSpace());
}

TEST(HPolyhedronRedundancy, RemovesDuplicatesSequentiallyWithoutLosingBounds) {
    // Two copies of each essential bound must never delete one another.
    MatrixXd A(6, 1);
    A << 1, 2, -1, -3, 1, 0;
    const VectorXd b = Vec({1, 2, 1, 3, 2, 0});
    HPolyhedron P(A, b);
    P.RemoveRedundantConstraints(1e-9);
    EXPECT_EQ(P.NumInequalities(), 2u);
    ExpectSameSet(P, UnitBox(1));
    P.RemoveRedundantConstraints(1e-9);
    EXPECT_EQ(P.NumInequalities(), 2u);
}

TEST(HPolyhedronRedundancy, EqualityCanMakeEveryInequalityRedundant) {
    HPolyhedron P(UnitBox(2).Ai(), UnitBox(2).bi(),
                  MatrixXd::Identity(2, 2), Vec({0, 0}));
    P.RemoveRedundantConstraints(1e-9);
    EXPECT_EQ(P.NumInequalities(), 0u);
    EXPECT_EQ(P.NumEqualities(), 2u);
    EXPECT_FALSE(P.IsFullSpace());
    ExpectSameSet(P, Point(Vec({0, 0})));
}

TEST(HPolyhedronRedundancy, ConstantRowsInfeasibilityAndInvalidParameters) {
    HPolyhedron tautology(MatrixXd::Zero(1, 2), Vec({1}));
    tautology.RemoveRedundantConstraints(kTol);
    EXPECT_TRUE(tautology.IsFullSpace());
    HPolyhedron contradiction(MatrixXd::Zero(1, 2), Vec({-1}));
    contradiction.RemoveRedundantConstraints(kTol);
    EXPECT_FALSE(contradiction.IsFeasible());
    auto inconsistent = Box(Vec({2}), Vec({1}));
    inconsistent.RemoveRedundantConstraints(kTol);
    EXPECT_FALSE(inconsistent.IsFeasible());
    auto essential = UnitBox(2);
    essential.RemoveRedundantConstraints(kTol);
    EXPECT_EQ(essential.NumInequalities(), 4u);
    const auto original = essential;
    for (double bad : {-1.0, kNaN, kInf}) {
        essential.RemoveRedundantConstraints(bad);
        EXPECT_TRUE(essential.Ai().isApprox(original.Ai(), 0));
    }
    HPolyhedron invalid;
    invalid.RemoveRedundantConstraints(kTol);
    EXPECT_FALSE(invalid.IsValid());
}

TEST(HPolyhedronRedundancy, OneHundredRandomRedundantHalfspacesReduceToTheKnownBox) {
    std::mt19937 rng(19037);
    std::uniform_real_distribution<double> exponent_dist(-8.0, 8.0);
    for (int trial = 0; trial < 100; ++trial) {
        SCOPED_TRACE(trial);
        const Index n = 1 + trial % 10;
        const HPolyhedron box = UnitBox(n);
        constexpr Index kAdditionalRows = 8;
        MatrixXd A(2 * n + kAdditionalRows, n);
        VectorXd b(A.rows());
        A.topRows(2 * n) = box.Ai();
        b.head(2 * n) = box.bi();

        // Every added row strictly contains the unit box because its bound is
        // larger than max_{|x|<=1} g*x = ||g||_1. Random row scales exercise
        // normalization without changing that mathematical certificate.
        for (Index i = 0; i < kAdditionalRows; ++i) {
            const Eigen::RowVectorXd normal = RandomMatrix(rng, 1, n);
            const double row_scale = std::pow(10.0, exponent_dist(rng));
            A.row(2 * n + i) = row_scale * normal;
            b(2 * n + i) = row_scale * (normal.cwiseAbs().sum() + 0.25);
        }

        HPolyhedron P(A, b);
        P.RemoveRedundantConstraints(1e-9);
        ASSERT_EQ(P.NumInequalities(), static_cast<size_t>(2 * n));
        EXPECT_TRUE(P.Ai().isApprox(box.Ai(), 0.0));
        EXPECT_TRUE(P.bi().isApprox(box.bi(), 0.0));
    }
}

TEST(HPolyhedronRedundancy, DenseShuffledRepresentationsRetainEveryEssentialFacet) {
    std::mt19937 rng(87103);
    std::uniform_real_distribution<double> exponent_dist(-8.0, 8.0);
    for (int trial = 0; trial < 40; trial++) {
        SCOPED_TRACE(trial);
        const Index n = 2 + trial % 6;
        const VectorXd center = RandomMatrix(rng, n, 1);
        const MatrixXd T = RandomInvertibleMap(rng, n);
        const HPolyhedron expected = TransformedUnitBox(center, T);
        const Index duplicate_count = n;
        constexpr Index kAdditionalRows = 6;
        const Index row_count = expected.Ai().rows() + duplicate_count + kAdditionalRows;
        MatrixXd unshuffled_A(row_count, n);
        VectorXd unshuffled_b(row_count);

        // Preserve one independently scaled copy of every essential facet.
        for (Index i = 0; i < expected.Ai().rows(); i++) {
            const double scale = std::pow(10.0, exponent_dist(rng));
            unshuffled_A.row(i) = scale * expected.Ai().row(i);
            unshuffled_b(i) = scale * expected.bi()(i);
        }
        // Add differently scaled duplicates. Sequential removal must retain one
        // representative of each facet instead of deleting both copies.
        for (Index i = 0; i < duplicate_count; i++) {
            const double scale = std::pow(10.0, exponent_dist(rng));
            unshuffled_A.row(expected.Ai().rows() + i) = scale * expected.Ai().row(i);
            unshuffled_b(expected.Ai().rows() + i) = scale * expected.bi()(i);
        }
        // Each remaining row lies strictly outside the parallelotope because
        // h_P(g) = g*c + ||g*T||_1 is available in closed form.
        for (Index i = 0; i < kAdditionalRows; i++) {
            const Eigen::RowVectorXd normal = RandomMatrix(rng, 1, n);
            const double support = normal.dot(center)
                + (normal * T).cwiseAbs().sum();
            const double scale = std::pow(10.0, exponent_dist(rng));
            const Index row = expected.Ai().rows() + duplicate_count + i;
            unshuffled_A.row(row) = scale * normal;
            unshuffled_b(row) = scale * (support + 0.5);
        }

        std::vector<Index> order(static_cast<size_t>(row_count));
        std::iota(order.begin(), order.end(), 0);
        std::shuffle(order.begin(), order.end(), rng);
        MatrixXd A(row_count, n);
        VectorXd b(row_count);
        for (Index i = 0; i < row_count; i++) {
            A.row(i) = unshuffled_A.row(order[static_cast<size_t>(i)]);
            b(i) = unshuffled_b(order[static_cast<size_t>(i)]);
        }

        HPolyhedron actual(A, b);
        actual.RemoveRedundantConstraints(1e-9);
        ASSERT_EQ(actual.NumInequalities(), static_cast<size_t>(2 * n));
        ExpectSameSet(actual, expected, 2e-6);
    }
}

TEST(HPolyhedronProjection, OneSidedFamiliesDisappearAndIndependentRowsRemain) {
    for (double sign : {1.0, -1.0}) {
        HPolyhedron P((MatrixXd(2, 2) << sign, -1, 0, 1).finished(), Vec({0, 2}));
        ExpectSameSet(P.Projection(0), HPolyhedron(MatrixXd::Ones(1, 1), Vec({2})));
        const HPolyhedron one_row((MatrixXd(1, 2) << sign, -1).finished(), Vec({0}));
        EXPECT_TRUE(one_row.Projection(0).IsFullSpace());
    }
}

TEST(HPolyhedronProjection, EqualitySubstitutionRetainsInducedBounds) {
    // x+y=1, x<=2, y<=2 projects to -1<=y<=2. Leaving x's coefficient
    // active after substitution used to discard the lower bound on y.
    const HPolyhedron P(MatrixXd::Identity(2, 2), Vec({2, 2}),
                         Vec({1, 1}).transpose(), Vec({1}));
    ExpectSameSet(P.Projection(0), Box(Vec({-1}), Vec({2})));
    EXPECT_FALSE(P.Projection(0).Contains(Vec({-1.1})));

    // The pivot is not the first equation; substitution must update every
    // other equality and preserve nonzero right-hand sides.
    const HPolyhedron Q(MatrixXd(0, 3), VectorXd(0),
        (MatrixXd(2, 3) << 0, 1, 1, 2, 1, 0).finished(), Vec({3, 2}));
    const auto projected = Q.Projection(0);
    EXPECT_TRUE(projected.Contains(Vec({1, 2})));
    EXPECT_FALSE(projected.Contains(Vec({1, 3})));
    ExpectSameSet(projected, HPolyhedron(MatrixXd(0, 2), VectorXd(0),
        Vec({1, 1}).transpose(), Vec({3})));
}

TEST(HPolyhedronProjection, ContradictoryEqualitiesAndInequalitiesStayEmpty) {
    const HPolyhedron equality_empty(MatrixXd(0, 2), VectorXd(0),
        (MatrixXd(2, 2) << 1, 0, 2, 0).finished(), Vec({1, 3}));
    EXPECT_FALSE(equality_empty.Projection(0).IsFeasible());
    const HPolyhedron inequality_empty((MatrixXd(2, 2) << 1, 0, -1, 0).finished(), Vec({0, -1}));
    EXPECT_FALSE(inequality_empty.Projection(0).IsFeasible());
    ExpectSameSet(Point(Vec({2, -3})).Projection(0), Point(Vec({-3})));
    EXPECT_TRUE(HPolyhedron::FullSpace(4).Projection(2).IsFullSpace());
}

TEST(HPolyhedronProjection, NonzeroSmallCoefficientsAreNotDropped) {
    // Eliminate x from x<=0 and -epsilon*x+y<=0. For any positive epsilon,
    // y<=0 is necessary. Treating a small coefficient as zero or comparing
    // against the norm of another row changes the projection.
    for (double scale : {1e-12, 1.0, 1e12}) {
        MatrixXd A(2, 2);
        A << scale, 0, -1e-12, 1;
        ExpectSameSet(HPolyhedron(A, Vec({0, 0})).Projection(0),
                      HPolyhedron(MatrixXd::Ones(1, 1), Vec({0})));
    }
    const HPolyhedron P(Vec({1, 0}).transpose(), Vec({0}),
                         Vec({1e-12, 1}).transpose(), Vec({0}));
    ExpectSameSet(P.Projection(0), HPolyhedron(MatrixXd::Constant(1, 1, -1), Vec({0})));
}

TEST(HPolyhedronProjection, CoordinateSelectionOrderingAndLargeDimensions) {
    for (Index n : {2, 5, 10}) {
        for (size_t removed = 0; removed < static_cast<size_t>(n); ++removed) {
            SCOPED_TRACE(n);
            SCOPED_TRACE(removed);
            ExpectSameSet(UnitBox(n).Projection(removed), UnitBox(n - 1));
            ExpectSameSet(Simplex(n).Projection(removed), Simplex(n - 1));
        }
    }
    const auto P = Box(Vec({1, 2, 3, 4}), Vec({2, 4, 6, 8}));
    ExpectSameSet(P.ProjectOnto({3, 1}), Box(Vec({4, 2}), Vec({8, 4})));
    ExpectSameSet(P.ProjectOnto({2}), Box(Vec({3}), Vec({6})));
    ExpectSameSet(P.ProjectOnto({0, 1, 2, 3}), P);
}

TEST(HPolyhedronProjection, RejectsInvalidRequestsAndExcessiveGrowth) {
    EXPECT_FALSE(HPolyhedron().Projection(0).IsValid());
    EXPECT_FALSE(UnitBox(1).Projection(0).IsValid());
    EXPECT_FALSE(UnitBox(2).Projection(2).IsValid());
    EXPECT_FALSE(UnitBox(2).ProjectOnto({}).IsValid());
    EXPECT_FALSE(UnitBox(2).ProjectOnto({0, 0}).IsValid());
    EXPECT_FALSE(UnitBox(2).ProjectOnto({2}).IsValid());
    EXPECT_FALSE(HPolyhedron().ProjectOnto({0}).IsValid());
    for (double bad : {-1.0, kNaN, kInf}) {
        EXPECT_FALSE(UnitBox(2).Projection(0, bad).IsValid());
        EXPECT_FALSE(UnitBox(2).ProjectOnto({0}, bad).IsValid());
    }
    // 317 upper and 317 lower rows require >100000 FM combinations.
    MatrixXd A = MatrixXd::Ones(634, 2);
    A.bottomRows(317).col(0).setConstant(-1);
    const HPolyhedron excessive_growth(A, VectorXd::Ones(634));
    EXPECT_FALSE(excessive_growth.Projection(0).IsValid());
    EXPECT_FALSE(excessive_growth.ProjectOnto({1}).IsValid());

    // Independent rows alone are also subject to the projection row cap.
    MatrixXd independent_A = MatrixXd::Zero(100001, 2);
    independent_A.col(1).setOnes();
    EXPECT_FALSE(HPolyhedron(independent_A, VectorXd::Ones(100001))
        .Projection(0).IsValid());

    // Row normalization must reject overflow in either the inequality or
    // equality right-hand side instead of projecting a changed set.
    const HPolyhedron inequality_overflow(
        (MatrixXd(1, 2) << 1e-320, 0).finished(), Vec({1}));
    EXPECT_FALSE(inequality_overflow.Projection(0).IsValid());
    const HPolyhedron equality_overflow(
        MatrixXd(0, 2), VectorXd(0),
        (MatrixXd(1, 2) << 1e-320, 0).finished(), Vec({1}));
    EXPECT_FALSE(equality_overflow.Projection(0).IsValid());
}

TEST(HPolyhedronProjection, RandomBoundedSetsMatchLiftedDirectionSupports) {
    std::mt19937 rng(1067);
    for (int trial = 0; trial < 20; ++trial) {
        SCOPED_TRACE(trial);
        const Index n = 3 + trial % 3;
        MatrixXd A(2 * n + 4, n);
        A << UnitBox(n).Ai(), RandomMatrix(rng, 4, n);
        VectorXd b = VectorXd::Ones(A.rows());
        const HPolyhedron P(A, b);
        const Index removed = trial % n;
        const auto projected = P.Projection(removed, 1e-10);
        ASSERT_TRUE(projected.IsValid());
        const MatrixXd directions = RandomMatrix(rng, 12, n - 1);
        MatrixXd lifted = MatrixXd::Zero(directions.rows(), n);
        lifted.leftCols(removed) = directions.leftCols(removed);
        lifted.rightCols(n - removed - 1) = directions.rightCols(n - removed - 1);
        ExpectVectorNear(projected.ComputeSupport(directions), P.ComputeSupport(lifted));
    }
}

TEST(HPolyhedronProjection, RandomParallelotopeProjectionsMatchClosedFormSupport) {
    std::mt19937 rng(70163);
    for (int trial = 0; trial < 50; ++trial) {
        SCOPED_TRACE(trial);
        const Index n = 2 + trial % 5;
        const VectorXd center = RandomMatrix(rng, n, 1);
        const MatrixXd T = RandomInvertibleMap(rng, n);
        const HPolyhedron P = TransformedUnitBox(center, T);
        const Index removed = trial % n;
        const HPolyhedron projected = P.Projection(removed, 1e-9);
        ASSERT_TRUE(projected.IsValid());

        const MatrixXd directions = RandomMatrix(rng, 20, n - 1);
        MatrixXd lifted = MatrixXd::Zero(directions.rows(), n);
        lifted.leftCols(removed) = directions.leftCols(removed);
        lifted.rightCols(n - removed - 1) = directions.rightCols(n - removed - 1);
        const VectorXd expected = lifted * center
            + (lifted * T).cwiseAbs().rowwise().sum();
        ExpectVectorNear(projected.ComputeSupport(directions), expected, 1e-6);
    }
}

TEST(HPolyhedronProjection, RandomCoordinateSelectionsPreserveRequestedOrder) {
    std::mt19937 rng(46021);
    for (int trial = 0; trial < 20; ++trial) {
        SCOPED_TRACE(trial);
        const Index n = 3 + trial % 2;
        const VectorXd center = RandomMatrix(rng, n, 1);
        const MatrixXd T = RandomInvertibleMap(rng, n);
        const HPolyhedron P = TransformedUnitBox(center, T);
        const std::vector<size_t> coordinates = {
            static_cast<size_t>(n - 1), 0};
        const HPolyhedron projected = P.ProjectOnto(coordinates, 1e-9);
        ASSERT_TRUE(projected.IsValid());

        const MatrixXd directions = RandomMatrix(rng, 20, 2);
        MatrixXd lifted = MatrixXd::Zero(directions.rows(), n);
        lifted.col(n - 1) = directions.col(0);
        lifted.col(0) = directions.col(1);
        const VectorXd expected = lifted * center
            + (lifted * T).cwiseAbs().rowwise().sum();
        ExpectVectorNear(projected.ComputeSupport(directions), expected, 1e-6);
    }
}

TEST(HPolyhedronProjection, LowerDimensionalAffineCubesMatchClosedFormProjections) {
    std::mt19937 rng(62011);
    for (int trial = 0; trial < 40; trial++) {
        SCOPED_TRACE(trial);
        const Index n = 3 + trial % 5;
        const Index k = 1 + (2 * trial) % (n - 1);
        const Index projected_dimension = 1 + (3 * trial) % (n - 1);
        const VectorXd center = RandomMatrix(rng, n, 1);
        const MatrixXd coordinates = RandomInvertibleMap(rng, n);
        const MatrixXd basis = coordinates.leftCols(k);
        const HPolyhedron P = AffineCube(center, coordinates, k);

        std::vector<size_t> permutation(static_cast<size_t>(n));
        std::iota(permutation.begin(), permutation.end(), 0);
        std::shuffle(permutation.begin(), permutation.end(), rng);
        const std::vector<size_t> retained(
            permutation.begin(), permutation.begin() + projected_dimension);
        const HPolyhedron projected = P.ProjectOnto(retained, 1e-9);
        ASSERT_TRUE(projected.IsValid());
        ASSERT_TRUE(projected.IsFeasible());
        EXPECT_TRUE(projected.IsBounded());

        VectorXd projected_center(projected_dimension);
        MatrixXd projected_basis(projected_dimension, k);
        for (Index i = 0; i < projected_dimension; i++) {
            const Index coordinate = static_cast<Index>(retained[static_cast<size_t>(i)]);
            projected_center(i) = center(coordinate);
            projected_basis.row(i) = basis.row(coordinate);
        }
        const MatrixXd directions = RandomMatrix(rng, 25, projected_dimension);
        const VectorXd expected = directions * projected_center
            + (directions * projected_basis).cwiseAbs().rowwise().sum();
        MatrixXd lifted_directions = MatrixXd::Zero(directions.rows(), n);
        for (Index i = 0; i < projected_dimension; i++) {
            const Index coordinate = static_cast<Index>(retained[static_cast<size_t>(i)]);
            lifted_directions.col(coordinate) = directions.col(i);
        }
        {
            SCOPED_TRACE("source affine cube");
            ExpectVectorNear(P.ComputeSupport(lifted_directions), expected, 2e-6);
        }
        {
            SCOPED_TRACE("projected affine cube");
            ExpectVectorNear(projected.ComputeSupport(directions), expected, 2e-6);
        }
    }
}

TEST(HPolyhedronLinearImage, LegacyInvertiblePolygonExample) {
    MatrixXd T(2, 2);
    T << 1, 2, 1, 0;
    MatrixXd expected_A(4, 2);
    expected_A << 0, 1, 1.5, -0.5, 0, -1, -0.5, -0.5;
    const HPolyhedron expected(expected_A, VectorXd::Constant(4, 4));
    for (int method = 0; method < 3; ++method) {
        ExpectSameSet(Image(LegacyPolygon(), T, method), expected);
    }
}

TEST(HPolyhedronLinearImage, LegacySingularPolygonExampleIncludesImageEquality) {
    MatrixXd T(2, 2);
    T << 1, 2, 0, 0;
    // The original polygon has y1 range [-12,4]; y2 is identically zero.
    const HPolyhedron expected((MatrixXd(2, 2) << 1, 0, -1, 0).finished(),
                               Vec({4, 12}), Vec({0, 1}).transpose(), Vec({0}));
    for (int method = 0; method < 3; ++method) {
        const auto image = Image(LegacyPolygon(), T, method);
        ExpectSameSet(image, expected);
        EXPECT_FALSE(image.Contains(Vec({0, 1})));
    }
}

TEST(HPolyhedronLinearImage, AllRanksAndRectangularShapesRespectSupportAndRange) {
    std::vector<MatrixXd> maps;
    maps.push_back(MatrixXd::Identity(2, 2));
    maps.push_back((MatrixXd(2, 2) << 1, 0.5, 0, -1).finished());
    maps.push_back((MatrixXd(2, 2) << 0, -1, 1, 0).finished());
    maps.push_back((MatrixXd(1, 2) << 1, 2).finished());
    maps.push_back((MatrixXd(3, 2) << 1, 0, 0, 1, 1, 1).finished());
    maps.push_back((MatrixXd(3, 2) << 1, 2, 2, 4, 0, 0).finished());
    maps.push_back(MatrixXd::Zero(3, 2));
    std::mt19937 rng(748);
    for (size_t i = 0; i < maps.size(); ++i) {
        const auto& T = maps[i];
        const MatrixXd directions = RandomMatrix(rng, 20, T.rows());
        // Image of the unit box has h(d) = ||d*T||_1, independent of the
        // implementation of images, projection, or support of the source.
        const VectorXd expected = (directions * T).cwiseAbs().rowwise().sum();
        for (int method = 0; method < 3; ++method) {
            SCOPED_TRACE(i);
            SCOPED_TRACE(method);
            const auto image = Image(UnitBox(2), T, method);
            ASSERT_TRUE(image.IsValid());
            EXPECT_EQ(image.Dimension(), static_cast<size_t>(T.rows()));
            ExpectVectorNear(image.ComputeSupport(directions), expected);
            EXPECT_TRUE(image.Contains(T * Vec({0.25, -0.5})));
            if (i == 5) {
                EXPECT_FALSE(image.Contains(Vec({0, 0, 1})));
            }
            if (i == 6) {
                ExpectSameSet(image, Point(VectorXd::Zero(3)));
            }
        }
    }
}

TEST(HPolyhedronLinearImage, AffineHullsFullSpaceAndEmptyInputs) {
    MatrixXd T(3, 2);
    T << 1, 2, 2, 4, 0, 0;
    const auto image = HPolyhedron::FullSpace(2).AffineTransform(T);
    EXPECT_TRUE(image.Contains(Vec({3, 6, 0})));
    EXPECT_FALSE(image.Contains(Vec({3, 5, 0})));
    EXPECT_FALSE(image.Contains(Vec({3, 6, 1})));
    EXPECT_FALSE(image.IsBounded());
    ExpectSameSet(Point(Vec({1, 2})).AffineTransform(T), Point(Vec({5, 10, 0})));
    EXPECT_FALSE(HPolyhedron::EmptySet(2).AffineTransform(T).IsFeasible());
    EXPECT_FALSE(HPolyhedron::EmptySet(2).AffineTransform(MatrixXd::Zero(1, 2)).IsFeasible());
    ExpectSameSet(UnitBox(2).AffineTransform(MatrixXd::Zero(1, 2)), Point(Vec({0})));
    ExpectSameSet(UnitBox(2).AffineTransform(MatrixXd::Zero(2, 2)),
        Point(VectorXd::Zero(2)));
}

TEST(HPolyhedronLinearImage, PropagatesProjectionGrowthFailure) {
    // A rank-one map from R^3 introduces two kernel coordinates. The first
    // elimination exceeds the Fourier-Motzkin cap; the image must stay invalid
    // and the remaining elimination must not run on it.
    MatrixXd A = MatrixXd::Ones(634, 3);
    A.bottomRows(317) *= -1;
    const HPolyhedron P(A, VectorXd::Ones(634));
    const MatrixXd T = (MatrixXd(1, 3) << 1, 0, 0).finished();
    EXPECT_FALSE(P.AffineTransform_SVD(T).IsValid());
}

TEST(HPolyhedronLinearImage, ScalingIllConditioningAndTenDimensionalMaps) {
    std::mt19937 rng(975);
    for (double scale : {1e-8, 1.0, 1e8}) {
        const MatrixXd T = scale * MatrixXd::Identity(3, 3);
        for (int method = 0; method < 3; ++method) {
            const auto image = Image(UnitBox(3), T, method);
            const VectorXd support = image.ComputeSupport(MatrixXd::Identity(3, 3));
            ExpectVectorNear(support / scale, VectorXd::Ones(3));
        }
    }
    MatrixXd ill = MatrixXd::Identity(3, 3);
    ill(2, 2) = 1e-12;
    for (int method = 0; method < 3; ++method) {
        const auto image = Image(UnitBox(3), ill, method);
        const VectorXd support = image.ComputeSupport(MatrixXd::Identity(3, 3));
        EXPECT_NEAR(support(2), 1e-12, 1e-18);
    }
    // Invertible and sparse coordinate maps avoid exponential elimination
    // while exercising the dimensions required by the invariant-set workflow.
    MatrixXd T = MatrixXd::Identity(10, 10) + 0.05 * RandomMatrix(rng, 10, 10);
    const MatrixXd directions = RandomMatrix(rng, 20, 10);
    ExpectVectorNear(Simplex(10).AffineTransform(T).ComputeSupport(directions),
        (directions * T).rowwise().maxCoeff().cwiseMax(0.0));
    T = MatrixXd::Identity(10, 10).topRows(4);
    ExpectSameSet(UnitBox(10).AffineTransform(T), UnitBox(4));
}

TEST(HPolyhedronLinearImage, RejectsInvalidDataWithoutEigenAssertions) {
    for (int method = 0; method < 3; ++method) {
        EXPECT_FALSE(Image(HPolyhedron(), MatrixXd::Identity(2, 2), method).IsValid());
        EXPECT_FALSE(Image(UnitBox(2), MatrixXd::Identity(3, 3), method).IsValid());
        EXPECT_FALSE(Image(UnitBox(2), MatrixXd(0, 2), method).IsValid());
        for (double bad : {kNaN, kInf}) {
            EXPECT_FALSE(Image(UnitBox(2), MatrixXd::Constant(2, 2, bad), method).IsValid());
        }
    }
    for (double bad : {0.0, -1.0, kNaN, kInf}) {
        EXPECT_FALSE(UnitBox(2).AffineTransform(MatrixXd::Identity(2, 2), bad).IsValid());
    }
    // Explicitly exercise the condition-threshold SVD selection.
    ExpectSameSet(UnitBox(2).AffineTransform(
        (MatrixXd(2, 2) << 2, 0, 0, 1).finished(), 1),
        Box(Vec({-2, -1}), Vec({2, 1})));
}

TEST(HPolyhedronSetOperations, IntersectionRetainsBothConstraintTypes) {
    const auto P = UnitBox(2).Intersection(
        HPolyhedron(MatrixXd(0, 2), VectorXd(0), Vec({1, 1}).transpose(), Vec({0.5})));
    EXPECT_TRUE(P.Contains(Vec({0.25, 0.25})));
    EXPECT_FALSE(P.Contains(Vec({1, 1})));
    EXPECT_FALSE(P.Contains(Vec({-2, 2.5})));
    ExpectSameSet(UnitBox(2).Intersection(HPolyhedron::FullSpace(2)), UnitBox(2));
    EXPECT_FALSE(UnitBox(2).Intersection(HPolyhedron::EmptySet(2)).IsFeasible());
    EXPECT_FALSE(UnitBox(2).Intersection(UnitBox(3)).IsValid());
    EXPECT_FALSE(UnitBox(2).Intersection(HPolyhedron()).IsValid());
    EXPECT_FALSE(HPolyhedron().Intersection(UnitBox(2)).IsValid());
}

TEST(HPolyhedronSetOperations, PreimageSupportsSingularAndRectangularMaps) {
    const MatrixXd T = (MatrixXd(2, 3) << 1, 1, 0, 0, 0, 2).finished();
    const auto P = UnitBox(2).Preimage(T);
    EXPECT_TRUE(P.Contains(Vec({100, -100, 0.5})));
    EXPECT_FALSE(P.Contains(Vec({1, 1, 0})));
    EXPECT_FALSE(P.Contains(Vec({0, 0, 0.6})));
    const auto equality_preimage = Point(Vec({1, 2})).Preimage(T);
    EXPECT_TRUE(equality_preimage.Contains(Vec({0.5, 0.5, 1})));
    EXPECT_FALSE(equality_preimage.Contains(Vec({0.5, 0.5, 0.9})));
    EXPECT_TRUE(UnitBox(2).Preimage(MatrixXd::Zero(2, 3)).IsFeasible());
    EXPECT_FALSE(Point(Vec({1, 2})).Preimage(MatrixXd::Zero(2, 3)).IsFeasible());
    EXPECT_FALSE(UnitBox(2).Preimage(MatrixXd::Identity(3, 3)).IsValid());
    EXPECT_FALSE(UnitBox(2).Preimage(MatrixXd(2, 0)).IsValid());
    EXPECT_FALSE(UnitBox(2).Preimage(MatrixXd::Constant(2, 2, kNaN)).IsValid());
    EXPECT_FALSE(HPolyhedron().Preimage(T).IsValid());
}

TEST(HPolyhedronSetOperations, RandomRectangularPreimagesMatchPointwiseDefinition) {
    std::mt19937 rng(53381);
    for (int trial = 0; trial < 50; trial++) {
        SCOPED_TRACE(trial);
        const Index target_dimension = 2 + trial % 7;
        const Index source_dimension = target_dimension + 2;
        const Index intrinsic_dimension = (trial % 2 == 0)
            ? target_dimension
            : target_dimension - 1;
        const VectorXd center = RandomMatrix(rng, target_dimension, 1);
        const MatrixXd coordinates = RandomInvertibleMap(rng, target_dimension);
        const MatrixXd basis = coordinates.leftCols(intrinsic_dimension);
        const HPolyhedron target = AffineCube(center, coordinates, intrinsic_dimension);
        MatrixXd T(target_dimension, source_dimension);
        T << MatrixXd::Identity(target_dimension, target_dimension),
               0.2 * RandomMatrix(rng, target_dimension, 2);
        const HPolyhedron preimage = target.Preimage(T);
        ASSERT_TRUE(preimage.IsValid());

        // Construct an in-set source point exactly, including when the target
        // is lower-dimensional and represented by equalities.
        const VectorXd coefficients = RandomMatrix(rng, intrinsic_dimension, 1);
        const VectorXd image_point = center + basis * coefficients;
        VectorXd source_point = VectorXd::Zero(source_dimension);
        source_point.head(target_dimension) = image_point;
        EXPECT_TRUE(target.Contains(T * source_point));
        EXPECT_TRUE(preimage.Contains(source_point));

        // Move outside either an affine-hull equality or a box facet.
        VectorXd outside_image;
        if (intrinsic_dimension < target_dimension) {
            outside_image = image_point + 0.1 * coordinates.col(intrinsic_dimension);
        }
        else {
            VectorXd outside_coefficients = coefficients;
            outside_coefficients(0) = 1.5;
            outside_image = center + basis * outside_coefficients;
        }
        source_point.head(target_dimension) = outside_image;
        EXPECT_FALSE(target.Contains(T * source_point));
        EXPECT_FALSE(preimage.Contains(source_point));

        // The defining identity must hold pointwise for arbitrary source data.
        for (int sample = 0; sample < 20; sample++) {
            const VectorXd x = RandomMatrix(rng, source_dimension, 1);
            EXPECT_EQ(preimage.Contains(x), target.Contains(T * x));
        }
    }
}

TEST(HPolyhedronSum, ExactSumIntroducesMissingFacetDirections) {
    // Box + diamond is an octagon. Updating box bounds by support alone
    // would return [-2,2]^2, incorrectly including (2,2).
    const HPolyhedron diamond((MatrixXd(4, 2) << 1, 1, 1, -1, -1, 1, -1, -1).finished(),
                              VectorXd::Ones(4));
    MatrixXd A(8, 2);
    A << UnitBox(2).Ai(), diamond.Ai();
    const HPolyhedron expected(A, Vec({2, 2, 2, 2, 3, 3, 3, 3}));
    const auto result = UnitBox(2) + diamond;
    ExpectSameSet(result, expected);
    ExpectSameSet(diamond + UnitBox(2), expected);
    EXPECT_TRUE(result.Contains(Vec({2, 1})));
    EXPECT_FALSE(result.Contains(Vec({2, 2})));
}

TEST(HPolyhedronSum, RedundantNormalsDoNotIntroduceSpuriousPoints) {
    MatrixXd A(3, 1);
    A << 1, -1, 1;
    const HPolyhedron P(A, Vec({1, 1, 9}));
    ExpectSameSet(P + Box(Vec({2}), Vec({3})), Box(Vec({1}), Vec({4})));
}

TEST(HPolyhedronSum, AffineSlicesSingletonsSelfSumAndUnboundedSets) {
    const HPolyhedron P((MatrixXd(2, 2) << 0, 1, 0, -1).finished(), Vec({1, 1}),
                         Vec({1, 0}).transpose(), Vec({0}));
    const HPolyhedron Q((MatrixXd(2, 2) << 0, 1, 0, -1).finished(), Vec({2, 2}),
                         Vec({1, 0}).transpose(), Vec({1}));
    ExpectSameSet(P + Q, HPolyhedron(P.Ai(), Vec({3, 3}), P.Ae(), Vec({1})));
    ExpectSameSet(UnitBox(2) + Point(Vec({2, -3})), Box(Vec({1, -4}), Vec({3, -2})));
    auto self = UnitBox(3);
    self += self;
    ExpectSameSet(self, Box(VectorXd::Constant(3, -2), VectorXd::Constant(3, 2)));
    const HPolyhedron upper(MatrixXd::Ones(1, 1), Vec({1}));
    const HPolyhedron lower(MatrixXd::Constant(1, 1, -1), Vec({0}));
    EXPECT_TRUE((upper + lower).IsFullSpace());
    ExpectSameSet(upper + upper, HPolyhedron(MatrixXd::Ones(1, 1), Vec({2})));
    EXPECT_TRUE((UnitBox(2) + HPolyhedron::FullSpace(2)).IsFullSpace());
    EXPECT_TRUE((HPolyhedron::FullSpace(2) + UnitBox(2)).IsFullSpace());
}

TEST(HPolyhedronSum, RandomTrianglesSatisfySupportAddition) {
    std::mt19937 rng(245);
    for (int trial = 0; trial < 20; ++trial) {
        SCOPED_TRACE(trial);
        const MatrixXd T = MatrixXd::Identity(2, 2) + 0.2 * RandomMatrix(rng, 2, 2);
        const auto Q = Simplex(2).AffineTransform(T);
        const auto sum = Simplex(2) + Q;
        const MatrixXd d = RandomMatrix(rng, 25, 2);
        // Each simplex's vertices are 0 and the columns of I or T.
        const VectorXd expected = d.rowwise().maxCoeff().cwiseMax(0.0)
            + (d * T).rowwise().maxCoeff().cwiseMax(0.0);
        ExpectVectorNear(sum.ComputeSupport(d), expected);
    }
}

TEST(HPolyhedronSum, DenseTranslatedSimplicesSatisfySupportAdditionInHigherDimensions) {
    std::mt19937 rng(90127);
    for (int trial = 0; trial < 6; trial++) {
        SCOPED_TRACE(trial);
        const Index n = 2 + trial % 3;
        const VectorXd center_P = RandomMatrix(rng, n, 1);
        const VectorXd center_Q = RandomMatrix(rng, n, 1);
        const MatrixXd T_P = RandomInvertibleMap(rng, n);
        const MatrixXd T_Q = RandomInvertibleMap(rng, n);
        const HPolyhedron P = TransformedSimplex(center_P, T_P);
        const HPolyhedron Q = TransformedSimplex(center_Q, T_Q);
        const HPolyhedron sum = P + Q;
        ASSERT_TRUE(sum.IsValid());
        ASSERT_TRUE(sum.IsFeasible());
        EXPECT_TRUE(sum.IsBounded());

        const MatrixXd directions = RandomMatrix(rng, 30, n);
        const VectorXd support_P = directions * center_P
            + (directions * T_P).rowwise().maxCoeff().cwiseMax(0.0);
        const VectorXd support_Q = directions * center_Q
            + (directions * T_Q).rowwise().maxCoeff().cwiseMax(0.0);
        {
            SCOPED_TRACE("first simplex");
            ExpectVectorNear(P.ComputeSupport(directions), support_P, 2e-6);
        }
        {
            SCOPED_TRACE("second simplex");
            ExpectVectorNear(Q.ComputeSupport(directions), support_Q, 2e-6);
        }
        {
            SCOPED_TRACE("Minkowski sum");
            ExpectVectorNear(sum.ComputeSupport(directions), support_P + support_Q, 2e-6);
        }

        // Independently sample feasible points from both simplices. Their sums
        // must belong to the computed Minkowski sum.
        for (int sample = 0; sample < 10; sample++) {
            VectorXd z_P = RandomMatrix(rng, n, 1).cwiseAbs();
            VectorXd z_Q = RandomMatrix(rng, n, 1).cwiseAbs();
            z_P /= std::max(1.0, z_P.sum());
            z_Q /= std::max(1.0, z_Q.sum());
            const VectorXd point = center_P + T_P * z_P + center_Q + T_Q * z_Q;
            EXPECT_TRUE(sum.Contains(point, 2e-6));
        }
    }
}

TEST(HPolyhedronSum, PropagatesProjectionGrowthFailure) {
    // Exact Minkowski addition projects one auxiliary variable per dimension.
    // A failed first projection must invalidate the result and stop the loop.
    MatrixXd A = MatrixXd::Ones(634, 2);
    A.bottomRows(317) *= -1;
    const HPolyhedron P(A, VectorXd::Ones(634));
    EXPECT_FALSE((P + UnitBox(2)).IsValid());
}

TEST(HPolyhedronDifference, AsymmetricIntervalsAndLegacyPolygon) {
    ExpectSameSet(Box(Vec({0}), Vec({2})) - Box(Vec({0}), Vec({1})),
                  Box(Vec({0}), Vec({1})));
    ExpectSameSet(Box(Vec({0}), Vec({2})) - Box(Vec({2}), Vec({3})),
                  Box(Vec({-2}), Vec({-1})));
    const auto Q = Box(Vec({-1, -1}), Vec({0.5, 0.5}));
    const auto result = LegacyPolygon() - Q;
    ExpectVectorNear(result.bi(), Vec({3.5, 2, 3, 2}));
    EXPECT_TRUE(result.Ai().isApprox(LegacyPolygon().Ai(), 0));

    const HPolyhedron P((MatrixXd(5, 2) << 1, 0, 0, 1, 0, -1, -1, 0, -1, -1).finished(),
                         Vec({6, 6, 0, 0, -4}));
    ExpectVectorNear((P - Box(Vec({0, 0}), Vec({1, 1}))).bi(), Vec({5, 5, 0, 0, -4}));
}

TEST(HPolyhedronDifference, ShrinkingCanProduceASingletonOrAnEmptySet) {
    ExpectSameSet(UnitBox(2) - UnitBox(2), Point(Vec({0, 0})));
    EXPECT_FALSE((UnitBox(2) - Box(Vec({-2, -2}), Vec({2, 2}))).IsFeasible());
    auto self = Box(Vec({1, 2}), Vec({3, 4}));
    self -= self;
    ExpectSameSet(self, Point(Vec({0, 0})));
}

TEST(HPolyhedronDifference, AffineEqualitiesAndIncompatibleDisturbances) {
    const HPolyhedron P((MatrixXd(2, 2) << 0, 1, 0, -1).finished(), Vec({5, 5}),
                         Vec({1, 0}).transpose(), Vec({3}));
    const HPolyhedron Q(P.Ai(), Vec({2, 2}), P.Ae(), Vec({1}));
    ExpectSameSet(P - Q, HPolyhedron(P.Ai(), Vec({3, 3}), P.Ae(), Vec({2})));
    EXPECT_FALSE((P - UnitBox(2)).IsFeasible());
    // Normal scaling must not allow variation along an equality to disappear.
    const HPolyhedron tiny(MatrixXd(0, 2), VectorXd(0),
                            Vec({1e-12, 0}).transpose(), Vec({3e-12}));
    EXPECT_FALSE((tiny - UnitBox(2)).IsFeasible());
    ExpectSameSet(tiny - Point(Vec({1, 0})), HPolyhedron(MatrixXd(0, 2), VectorXd(0),
        Vec({1, 0}).transpose(), Vec({2})));
}

TEST(HPolyhedronDifference, UnboundedSupportsDoNotProduceInfiniteConstraintData) {
    const HPolyhedron halfspace(Vec({1, 0}).transpose(), Vec({2}));
    const HPolyhedron ray((MatrixXd(3, 2) << -1, 0, 0, 1, 0, -1).finished(), Vec({0, 0, 0}));
    const auto empty = halfspace - ray;
    ASSERT_TRUE(empty.IsValid());
    EXPECT_TRUE(empty.Ai().allFinite());
    EXPECT_TRUE(empty.bi().allFinite());
    EXPECT_FALSE(empty.IsFeasible());

    const HPolyhedron strip((MatrixXd(2, 2) << 1, 0, -1, 0).finished(), Vec({1, 1}));
    ExpectSameSet(strip - strip, HPolyhedron(MatrixXd(0, 2), VectorXd(0),
        Vec({1, 0}).transpose(), Vec({0})));
    EXPECT_TRUE((HPolyhedron::FullSpace(2) - strip).IsFullSpace());
    EXPECT_TRUE((HPolyhedron::FullSpace(2) - HPolyhedron::FullSpace(2)).IsFullSpace());
    EXPECT_FALSE((strip - HPolyhedron::FullSpace(2)).IsFeasible());
    // A noncanonical representation of full space must not shrink to empty.
    const HPolyhedron tautology(MatrixXd::Zero(1, 2), Vec({1}));
    EXPECT_TRUE((tautology - HPolyhedron::FullSpace(2)).IsFeasible());
}

TEST(HPolyhedronDifference, MappedDisturbanceShrinksOnlyTheStatePartOfJointConstraints) {
    // Joint variables [z1,z2,r]. The last two rows are pure-input bounds;
    // the other rows couple z and r. W is the translated rectangle below.
    MatrixXd G(6, 3);
    G << 1, 2, 1, -1, 1, -2, 0, -1, 3, 1, 0, 0, 0, 0, 1, 0, 0, -1;
    const VectorXd F = Vec({10, 10, 10, 10, 4, 4});
    MatrixXd T(3, 2);
    T << 1, 2, -1, 1, 0, 0;
    const auto W = Box(Vec({-1, 0}), Vec({2, 1}));
    const auto shrunk = HPolyhedron(G, F).PontryaginDifferenceOfLinearImage(W, T);
    // Enumerating rectangle vertices is an independent support oracle.
    MatrixXd vertices(2, 4);
    vertices << -1, -1, 2, 2, 0, 1, 0, 1;
    const VectorXd expected = F - (G * T * vertices).rowwise().maxCoeff();
    ExpectVectorNear(shrunk.bi(), expected);
    EXPECT_TRUE(shrunk.Ai().isApprox(G, 0));
    EXPECT_EQ(shrunk.bi()(4), F(4));
    EXPECT_EQ(shrunk.bi()(5), F(5));
}

TEST(HPolyhedronDifference, RandomMappedSupportsMatchVertexOraclesUpToTenDimensions) {
    std::mt19937 rng(81023);
    for (int trial = 0; trial < 100; ++trial) {
        SCOPED_TRACE(trial);
        const Index n = 1 + trial % 10;
        const Index d = 1 + (trial / 10) % 4;
        MatrixXd G(2 * n + 4, n);
        G << UnitBox(n).Ai(), RandomMatrix(rng, 4, n);
        const MatrixXd T = RandomMatrix(rng, n, d);
        const VectorXd lower = RandomMatrix(rng, d, 1);
        const VectorXd upper = lower + VectorXd::Constant(d, 0.2);
        MatrixXd vertices(d, 1 << d);
        for (Index v = 0; v < vertices.cols(); ++v) {
            for (Index j = 0; j < d; ++j) {
                vertices(j, v) = ((v & (1 << j)) != 0) ? upper(j) : lower(j);
            }
        }
        const VectorXd F = VectorXd::Constant(G.rows(), 20);
        const auto result = HPolyhedron(G, F).PontryaginDifferenceOfLinearImage(Box(lower, upper), T);
        ASSERT_TRUE(result.IsValid());
        ExpectVectorNear(result.bi(), F - (G * T * vertices).rowwise().maxCoeff());
    }
}

TEST(HPolyhedronDifference, RandomMappedNonBoxSimplicesMatchVertexOracles) {
    std::mt19937 rng(34061);
    for (int trial = 0; trial < 80; trial++) {
        SCOPED_TRACE(trial);
        const Index n = 1 + trial % 10;
        const Index disturbance_dimension = 2 + trial % 5;
        MatrixXd G(2 * n + 6, n);
        G << UnitBox(n).Ai(), RandomMatrix(rng, 6, n);
        const VectorXd F = VectorXd::Constant(G.rows(), 50.0);
        const MatrixXd T_disturbance = RandomMatrix(rng, n, disturbance_dimension);
        const VectorXd center = RandomMatrix(rng, disturbance_dimension, 1);
        const MatrixXd T_simplex = RandomInvertibleMap(rng, disturbance_dimension);
        const HPolyhedron W = TransformedSimplex(center, T_simplex);

        const HPolyhedron result = HPolyhedron(G, F).PontryaginDifferenceOfLinearImage(W, T_disturbance);
        ASSERT_TRUE(result.IsValid());
        ASSERT_TRUE(result.IsFeasible());
        const MatrixXd directions = G * T_disturbance;
        // The simplex vertices are center and center plus each column of T_simplex.
        const VectorXd support =
            directions * center +
            (directions * T_simplex).rowwise().maxCoeff().cwiseMax(0.0);
        ExpectVectorNear(result.bi(), F - support, 2e-6);
        EXPECT_TRUE(result.Ai().isApprox(G, 0.0));
    }
}

TEST(HPolyhedronDifference, EqualityRobustificationUsesConstantAffineHullSupport) {
    std::mt19937 rng(48109);
    for (int trial = 0; trial < 30; trial++) {
        SCOPED_TRACE(trial);
        const Index disturbance_dimension = 2 + trial % 5;
        const Index n = disturbance_dimension + 2;
        const VectorXd center = RandomMatrix(rng, disturbance_dimension, 1);
        const MatrixXd coordinates = RandomInvertibleMap(rng, disturbance_dimension);
        const MatrixXd basis = coordinates.leftCols(disturbance_dimension - 1);
        const HPolyhedron W = AffineCube(
            center, coordinates, disturbance_dimension - 1);
        MatrixXd T = MatrixXd::Zero(n, disturbance_dimension);
        T.topRows(disturbance_dimension).setIdentity();

        MatrixXd G(2 * n + 3, n);
        G << UnitBox(n).Ai(), RandomMatrix(rng, 3, n);
        const VectorXd F = VectorXd::Constant(G.rows(), 30.0);
        Eigen::RowVectorXd equality = Eigen::RowVectorXd::Zero(n);
        equality.head(disturbance_dimension) = coordinates.inverse().bottomRows(1);
        const double equality_bound = 2.5;
        const HPolyhedron P(G, F, equality, Vec({equality_bound}));

        const HPolyhedron result = P.PontryaginDifferenceOfLinearImage(W, T);
        ASSERT_TRUE(result.IsValid());
        ASSERT_TRUE(result.IsFeasible());
        const MatrixXd directions = G * T;
        const VectorXd support = directions * center
            + (directions * basis).cwiseAbs().rowwise().sum();
        ExpectVectorNear(result.bi(), F - support, 2e-6);
        ASSERT_EQ(result.NumEqualities(), 1u);
        EXPECT_TRUE(result.Ae().isApprox(equality, 0.0));
        EXPECT_NEAR(result.be()(0), equality_bound - equality.head(disturbance_dimension).dot(center), 2e-7);

        // The first transformed coordinate varies over W. Using that normal as
        // a robust equality must therefore make the difference empty.
        Eigen::RowVectorXd varying_equality = Eigen::RowVectorXd::Zero(n);
        varying_equality.head(disturbance_dimension) = coordinates.inverse().topRows(1);
        const HPolyhedron varying_P(G, F, varying_equality, Vec({equality_bound}));
        EXPECT_FALSE(varying_P.PontryaginDifferenceOfLinearImage(W, T).IsFeasible());
    }
}

TEST(HPolyhedronDifference, MappedNonBoxAndLowerDimensionalDisturbances) {
    // The translated diagonal segment has vertices (1,-1) and (2,1).
    const HPolyhedron W((MatrixXd(2, 2) << 1, 0, -1, 0).finished(), Vec({2, -1}),
                         Vec({2, -1}).transpose(), Vec({3}));
    const MatrixXd T = (MatrixXd(2, 2) << 1, 2, -1, 1).finished();
    const MatrixXd vertices = (MatrixXd(2, 2) << 1, 2, -1, 1).finished();
    ExpectVectorNear(UnitBox(2).PontryaginDifferenceOfLinearImage(W, T).bi(),
                     UnitBox(2).bi() - (UnitBox(2).Ai() * T * vertices).rowwise().maxCoeff());
    // A simplex is non-box; its vertices are zero and the canonical axes.
    const MatrixXd mixed = (MatrixXd(2, 3) << 1, -2, 3, -1, 4, 2).finished();
    ExpectVectorNear(UnitBox(2).PontryaginDifferenceOfLinearImage(Simplex(3), mixed).bi(),
        UnitBox(2).bi() - (UnitBox(2).Ai() * mixed).rowwise().maxCoeff().cwiseMax(0.0));
    ExpectSameSet(UnitBox(2).PontryaginDifferenceOfLinearImage(W, MatrixXd::Zero(2, 2)), UnitBox(2));
    ExpectSameSet(UnitBox(2).PontryaginDifferenceOfLinearImage(HPolyhedron::FullSpace(3), MatrixXd::Zero(2, 3)), UnitBox(2));
}

TEST(HPolyhedronDifference, NumericalSupportFailureReturnsInvalidSet) {
    // W is valid as an H-representation, but its coordinate equilibration
    // overflows. A failed support query must not be interpreted as an empty
    // set or as zero disturbance.
    const HPolyhedron W(MatrixXd::Ones(1, 1), Vec({1e-320}));
    const HPolyhedron result = UnitBox(1).PontryaginDifferenceOfLinearImage(
        W, MatrixXd::Identity(1, 1));
    EXPECT_FALSE(result.IsValid());
}

TEST(HPolyhedronSetOperations, InvalidOperandsInvalidateResultsWithoutMutatingRhs) {
    const auto P = UnitBox(2);
    for (const auto& bad : {HPolyhedron(), UnitBox(1)}) {
        EXPECT_FALSE((P + bad).IsValid());
        EXPECT_FALSE((bad + P).IsValid());
        EXPECT_FALSE((P - bad).IsValid());
        EXPECT_FALSE((bad - P).IsValid());
    }
    EXPECT_TRUE(P.Contains(Vec({1, 1})));
    EXPECT_FALSE(P.PontryaginDifferenceOfLinearImage(HPolyhedron(), MatrixXd::Identity(2, 2)).IsValid());
    EXPECT_FALSE(P.PontryaginDifferenceOfLinearImage(UnitBox(3), MatrixXd::Identity(2, 2)).IsValid());
    EXPECT_FALSE(P.PontryaginDifferenceOfLinearImage(UnitBox(2), MatrixXd::Identity(3, 2)).IsValid());
    EXPECT_FALSE(P.PontryaginDifferenceOfLinearImage(UnitBox(2), MatrixXd::Constant(2, 2, kNaN)).IsValid());
    EXPECT_FALSE(HPolyhedron().PontryaginDifferenceOfLinearImage(P, MatrixXd::Identity(2, 2)).IsValid());
}

// These tests mirror the polyhedral operations used by Eq. (15).
TEST(HPolyhedronAlgorithmOperations, NilpotentDisturbanceSequenceSaturatesAtNu) {
    MatrixXd A = MatrixXd::Zero(3, 3);
    A(0, 1) = A(1, 2) = 1;
    const VectorXd E = Vec({0, 0, 1});
    MatrixXd G(5, 4);
    G << 1, 1, 1, 2, -1, -2, -3, 1, 2, -1, 4, 0, 0, 0, 0, 1, 0, 0, 0, -1;
    const VectorXd F = Vec({20, 20, 20, 4, 4});
    HPolyhedron shrunk(G, F);
    const auto W = Box(Vec({-0.1}), Vec({0.2}));
    VectorXd propagated = E;
    VectorXd expected = F;
    for (int t = 0; t < 6; ++t) {
        MatrixXd T = MatrixXd::Zero(4, 1);
        T.topRows(3) = propagated;
        const VectorXd d = G * T;
        for (Index i = 0; i < G.rows(); ++i) {
            expected(i) -= d(i) >= 0 ? 0.2 * d(i) : -0.1 * d(i);
        }
        shrunk = shrunk.PontryaginDifferenceOfLinearImage(W, T);
        ExpectVectorNear(shrunk.bi(), expected);
        propagated = A * propagated;
        if (t >= 2) {
            ExpectVectorNear(shrunk.bi(), Vec({19.4, 19.4, 18.7, 4, 4}));
        }
    }
}

TEST(HPolyhedronAlgorithmOperations, LiftedPreimageIntersectionAndProjection) {
    // x+ = r, r = v0, v0+ = v1, v1+ = v0, with |x|<=2 and |r|<=1.
    MatrixXd lifted(3, 3);
    lifted << 0, 1, 0, 0, 0, 1, 0, 1, 0;
    MatrixXd selector(2, 3);
    selector << 1, 0, 0, 0, 1, 0;
    const auto safe = Box(Vec({-2, -1}), Vec({2, 1})).Preimage(selector);
    auto C = safe.Intersection(safe.Preimage(lifted));
    ExpectSameSet(C, Box(Vec({-2, -1, -1}), Vec({2, 1, 1})));
    EXPECT_TRUE(C.IsSubsetOf(C.Preimage(lifted)));
    ExpectSameSet(C.ProjectOnto({0}), Box(Vec({-2}), Vec({2})));
}

}
