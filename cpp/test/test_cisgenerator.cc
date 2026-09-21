#include "cis_generator.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

using Eigen::Index;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using cis2m::CISComponent;
using cis2m::ControlledInvariantSetGenerator;
using cis2m::CISOptions;
using cis2m::HPolyhedron;

namespace {

HPolyhedron Box(const VectorXd& lower, const VectorXd& upper) {
    const Index n = lower.size();
    MatrixXd G(2 * n, n);
    G.topRows(n).setIdentity();
    G.bottomRows(n) = -MatrixXd::Identity(n, n);
    VectorXd F(2 * n);
    F << upper, -lower;
    return HPolyhedron(G, F);
}

HPolyhedron SymmetricBox(Index dimension, double radius) {
    return Box(VectorXd::Constant(dimension, -radius),
               VectorXd::Constant(dimension, radius));
}

CISOptions Single(std::size_t tau, std::size_t lambda, bool implicit = true) {
    CISOptions options;
    options.tau = tau;
    options.lambda = lambda;
    options.is_implicit = implicit;
    return options;
}

struct TransformedProblem {
    MatrixXd A;
    MatrixXd B;
    MatrixXd E;
    HPolyhedron safe_set;
    HPolyhedron disturbance_set;
};

TransformedProblem DenseTransformedProblem(Index n, Index m) {
    MatrixXd Ac = MatrixXd::Zero(n, n);
    MatrixXd Bc = MatrixXd::Zero(n, m);
    Index first = 0;
    for (Index channel = 0; channel < m; channel++) {
        const Index length = n / m + (channel < n % m ? 1 : 0);
        for (Index state = first; state + 1 < first + length; state++) {
            Ac(state, state + 1) = 1.0;
        }
        Bc(first + length - 1, channel) = 1.0;
        first += length;
    }

    MatrixXd T = MatrixXd::Identity(n, n);
    for (Index row = 0; row < n; row++) {
        for (Index column = 0; column < n; column++) {
            if (row != column) {
                const int pattern = static_cast<int>((3 * row + 5 * column + 1) % 7) - 3;
                T(row, column) = 0.012 * pattern;
            }
        }
    }
    MatrixXd Am(m, n);
    for (Index row = 0; row < m; row++) {
        for (Index column = 0; column < n; column++) {
            const int pattern = static_cast<int>((3 * row + 2 * column + 1) % 5) - 2;
            Am(row, column) = 0.025 * pattern;
        }
    }
    MatrixXd Bm = MatrixXd::Identity(m, m);
    for (Index row = 0; row < m; row++) {
        for (Index column = 0; column < m; column++) {
            const int pattern = static_cast<int>((2 * row + 3 * column + 1) % 5) - 2;
            Bm(row, column) += 0.04 * pattern;
        }
    }
    const MatrixXd inverse_T = T.inverse();
    const MatrixXd A = inverse_T * (Ac + Bc * Am) * T;
    const MatrixXd B = inverse_T * Bc * Bm;

    MatrixXd Ec(n, m);
    for (Index row = 0; row < n; row++) {
        for (Index column = 0; column < m; column++) {
            const int pattern = static_cast<int>((3 * row + 2 * column + 2) % 5) - 2;
            Ec(row, column) = 0.012 * pattern;
        }
    }
    const MatrixXd E = inverse_T * Ec;

    const HPolyhedron canonical_box = SymmetricBox(n + m, 2.5);
    MatrixXd Gc(canonical_box.Ai().rows() + 3, n + m);
    Gc.topRows(canonical_box.Ai().rows()) = canonical_box.Ai();
    for (Index row = canonical_box.Ai().rows(); row < Gc.rows(); row++) {
        for (Index column = 0; column < Gc.cols(); column++) {
            const int pattern = static_cast<int>((4 * row + 3 * column + 2) % 9) - 4;
            Gc(row, column) = 0.04 * pattern;
        }
    }
    VectorXd Fc(canonical_box.bi().size() + 3);
    Fc.head(canonical_box.bi().size()) = canonical_box.bi();
    Fc.tail(3).setConstant(1.4);

    // Pull a dense set in (z,r) coordinates back to the physical (x,u) space.
    MatrixXd joint_map = MatrixXd::Zero(n + m, n + m);
    joint_map.topLeftCorner(n, n) = T;
    joint_map.bottomLeftCorner(m, n) = Am * T;
    joint_map.bottomRightCorner(m, m) = Bm;
    const HPolyhedron safe_set = HPolyhedron(Gc, Fc).Preimage(joint_map);
    const VectorXd lower = VectorXd::Constant(m, -0.05);
    const VectorXd upper = VectorXd::Constant(m, 0.08);
    return {A, B, E, safe_set, Box(lower, upper)};
}

Index CountSignificantEntries(const MatrixXd& matrix, double tolerance = 1e-12) {
    return (matrix.cwiseAbs().array() > tolerance).count();
}

void ExpectMatrixNear(const MatrixXd& actual, const MatrixXd& expected, double tol = 1e-9) {
    ASSERT_EQ(actual.rows(), expected.rows());
    ASSERT_EQ(actual.cols(), expected.cols());
    EXPECT_LE((actual - expected).norm(), tol * std::max(1.0, expected.norm()));
}

void ExpectSameRepresentation(const HPolyhedron& actual, const HPolyhedron& expected) {
    ASSERT_TRUE(actual.IsValid());
    ASSERT_TRUE(expected.IsValid());
    ExpectMatrixNear(actual.Ai(), expected.Ai());
    ExpectMatrixNear(actual.bi(), expected.bi());
    ExpectMatrixNear(actual.Ae(), expected.Ae());
    ExpectMatrixNear(actual.be(), expected.be());
}

// Check every constraint of A_lifted*C + E_lifted*W <= C by support LPs.
void VerifyLiftedInvariance(const CISComponent& component,
                            const HPolyhedron& W = HPolyhedron::EmptySet(1)) {
    ASSERT_TRUE(component.is_implicit);
    ASSERT_TRUE(component.set.IsValid());
    ASSERT_TRUE(component.set.IsFeasible());
    const HPolyhedron& C = component.set;
    const Index ni = C.Ai().rows();
    const Index ne = C.Ae().rows();
    MatrixXd normals(ni + 2 * ne, C.Ai().cols());
    normals << C.Ai(), C.Ae(), -C.Ae();
    VectorXd bounds(ni + 2 * ne);
    bounds << C.bi(), C.be(), -C.be();
    const VectorXd support = C.ComputeSupport(normals * component.lifted_dynamics);
    const VectorXd disturbance_support = !W.IsFeasible() ?
        VectorXd::Zero(normals.rows()).eval() :
        W.ComputeSupport(normals * component.lifted_disturbance);
    ASSERT_TRUE(support.allFinite());
    ASSERT_TRUE(disturbance_support.allFinite());
    for (Index row = 0; row < normals.rows(); row++) {
        EXPECT_LE(support(row) + disturbance_support(row), bounds(row) + 2e-6)
            << "Lifted invariance constraint " << row;
    }
}

// Compute each disturbance support directly for an axis-aligned W, without
// using the generator's shrinking sequence or Pontryagin-difference helper.
double BoxSupport(const Eigen::RowVectorXd& direction,
                  const VectorXd& center, const VectorXd& radius) {
    return direction.dot(center) + direction.cwiseAbs().dot(radius);
}

bool FiniteTrajectoryOracle(const CISComponent& component,
                            const HPolyhedron& safe_set,
                            const VectorXd& initial,
                            Index state_dimension,
                            Index horizon,
                            const VectorXd& center = VectorXd(),
                            const VectorXd& radius = VectorXd()) {
    const Index n = state_dimension;
    const Index m = component.input_from_state.rows();
    const Index lifted_dimension = initial.size();
    MatrixXd joint_map = MatrixXd::Zero(n + m, lifted_dimension);
    joint_map.topLeftCorner(n, n).setIdentity();
    joint_map.bottomLeftCorner(m, n) = component.input_from_state;
    joint_map.bottomRightCorner(m, lifted_dimension - n) = component.input_from_virtual;

    VectorXd nominal = initial;
    MatrixXd disturbance_power = MatrixXd::Identity(lifted_dimension, lifted_dimension);
    std::vector<MatrixXd> disturbance_maps;
    for (Index t = 0; t < horizon; t++) {
        const VectorXd joint = joint_map * nominal;
        for (Index row = 0; row < safe_set.Ai().rows(); row++) {
            const Eigen::RowVectorXd normal = safe_set.Ai().row(row);
            double worst = normal.dot(joint);
            for (const MatrixXd& map : disturbance_maps) {
                worst += BoxSupport(normal * joint_map * map, center, radius);
            }
            if (worst > safe_set.bi()(row) + 1e-7) {
                return false;
            }
        }
        nominal = component.lifted_dynamics * nominal;
        disturbance_maps.push_back(disturbance_power * component.lifted_disturbance);
        disturbance_power = disturbance_power * component.lifted_dynamics;
    }
    return true;
}

// Characterize admissible (x,u) pairs and project them to x.
void VerifyProjectedControlledInvariance(const HPolyhedron& projected,
                                         const HPolyhedron& safe_set,
                                         const MatrixXd& A, const MatrixXd& B,
                                         const MatrixXd& E = MatrixXd(),
                                         const HPolyhedron& W = HPolyhedron::EmptySet(1)) {
    ASSERT_TRUE(projected.IsValid());
    ASSERT_TRUE(projected.IsFeasible());
    HPolyhedron next = projected;
    if (W.IsFeasible()) {
        next = next.PontryaginDifferenceOfLinearImage(W, E);
    }
    ASSERT_TRUE(next.IsValid());
    MatrixXd dynamics(A.rows(), A.cols() + B.cols());
    dynamics << A, B;
    const HPolyhedron predecessor = safe_set.Intersection(next.Preimage(dynamics));
    ASSERT_TRUE(predecessor.IsValid());
    std::vector<std::size_t> state_coordinates(static_cast<std::size_t>(A.rows()));
    for (Index i = 0; i < A.rows(); i++) {
        state_coordinates[static_cast<std::size_t>(i)] = static_cast<std::size_t>(i);
    }
    const HPolyhedron admissible_states = predecessor.ProjectOnto(state_coordinates);
    ASSERT_TRUE(admissible_states.IsValid());
    EXPECT_TRUE(projected.IsSubsetOf(admissible_states, 1e-6));
}

TEST(ControlledInvariantSetGeneratorTest, RejectsInvalidSystemsAndDisturbanceMatrices) {
    MatrixXd A(2, 2);
    A << 0, 1, 0, 0;
    MatrixXd B(2, 1);
    B << 0, 1;
    EXPECT_THROW(ControlledInvariantSetGenerator(MatrixXd::Zero(2, 3), B), std::invalid_argument);
    EXPECT_THROW(ControlledInvariantSetGenerator(A, MatrixXd::Zero(3, 1)), std::invalid_argument);
    EXPECT_THROW(ControlledInvariantSetGenerator(A, MatrixXd::Zero(2, 1)), std::invalid_argument);
    EXPECT_THROW(ControlledInvariantSetGenerator(MatrixXd::Zero(2, 2), B), std::invalid_argument);
    MatrixXd nonfinite = A;
    nonfinite(0, 0) = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(ControlledInvariantSetGenerator(nonfinite, B), std::invalid_argument);
    EXPECT_THROW(ControlledInvariantSetGenerator(A, B, MatrixXd::Zero(1, 1)), std::invalid_argument);
    MatrixXd E = MatrixXd::Ones(2, 1);
    E(0, 0) = std::numeric_limits<double>::infinity();
    EXPECT_THROW(ControlledInvariantSetGenerator(A, B, E), std::invalid_argument);
}

TEST(ControlledInvariantSetGeneratorTest, RejectsInvalidSetsAndParameters) {
    MatrixXd A = MatrixXd::Zero(1, 1);
    MatrixXd B = MatrixXd::Identity(1, 1);
    MatrixXd E = MatrixXd::Identity(1, 1);
    const ControlledInvariantSetGenerator nominal(A, B);
    const ControlledInvariantSetGenerator disturbed(A, B, E);
    const HPolyhedron safe = SymmetricBox(2, 1.0);
    const HPolyhedron W = SymmetricBox(1, 0.1);
    EXPECT_THROW(nominal.Compute(safe, CISOptions()), std::invalid_argument);
    EXPECT_THROW(nominal.Compute(HPolyhedron(), Single(0, 1)), std::invalid_argument);
    EXPECT_THROW(nominal.Compute(SymmetricBox(3, 1.0), Single(0, 1)), std::invalid_argument);
    EXPECT_THROW(nominal.Compute(MatrixXd::Ones(2, 2), VectorXd::Ones(1), Single(0, 1)),
                 std::invalid_argument);
    EXPECT_THROW(disturbed.Compute(safe, HPolyhedron(), Single(0, 1)), std::invalid_argument);
    EXPECT_THROW(disturbed.Compute(safe, HPolyhedron::EmptySet(2), Single(0, 1)),
                 std::invalid_argument);
    EXPECT_THROW(disturbed.Compute(safe, HPolyhedron::FullSpace(1), Single(0, 1)),
                 std::invalid_argument);
    EXPECT_THROW(disturbed.Compute(safe, SymmetricBox(2, 1.0), Single(0, 1)),
                 std::invalid_argument);
    EXPECT_THROW(nominal.Compute(safe, Single(0, std::numeric_limits<std::size_t>::max())),
                 std::invalid_argument);
    CISOptions oversized_hierarchy;
    oversized_hierarchy.hierarchy_level = std::numeric_limits<std::size_t>::max();
    EXPECT_THROW(nominal.Compute(safe, oversized_hierarchy), std::invalid_argument);
}

TEST(ControlledInvariantSetGeneratorTest, OmittedOrEmptyDisturbanceWithEUsesNominalDynamics) {
    MatrixXd A(2, 2);
    A << 0.3, 1.0, 0.0, 0.2;
    MatrixXd B(2, 1);
    B << 0.0, 1.0;
    MatrixXd E(2, 2);
    E << 0.1, -0.3, 0.2, 0.4;
    const HPolyhedron safe = SymmetricBox(3, 2.0);
    const ControlledInvariantSetGenerator nominal(A, B);
    const ControlledInvariantSetGenerator with_E(A, B, E);

    for (bool implicit : {true, false}) {
        const CISOptions options = Single(1, 2, implicit);
        const auto expected = nominal.Compute(safe, options);
        const auto omitted = with_E.Compute(safe, options);
        const auto empty = with_E.Compute(safe, HPolyhedron::EmptySet(2), options);
        const auto matrix = with_E.Compute(safe.Ai(), safe.bi(), options);
        ASSERT_EQ(expected.size(), 1u);
        ASSERT_EQ(omitted.size(), 1u);
        ASSERT_EQ(empty.size(), 1u);
        ASSERT_EQ(matrix.size(), 1u);
        for (const auto* actual : {&omitted.front(), &empty.front(), &matrix.front()}) {
            ExpectSameRepresentation(actual->set, expected.front().set);
            ExpectMatrixNear(actual->lifted_dynamics, expected.front().lifted_dynamics);
            EXPECT_EQ(actual->lifted_disturbance.rows(), 5);
            EXPECT_EQ(actual->lifted_disturbance.cols(), 0);
            if (implicit) {
                VerifyLiftedInvariance(*actual);
            }
        }
    }
}

TEST(ControlledInvariantSetGeneratorTest, DisturbanceSetWithoutEIsIgnored) {
    MatrixXd A(2, 2);
    A << 0.0, 1.0, 0.0, 0.0;
    MatrixXd B(2, 1);
    B << 0.0, 1.0;
    const ControlledInvariantSetGenerator generator(A, B);
    const HPolyhedron safe = SymmetricBox(3, 1.0);
    const HPolyhedron W = SymmetricBox(2, 0.2);
    const CISOptions options = Single(0, 2);
    const auto expected = generator.Compute(safe, options);
    const auto ignored = generator.Compute(safe, W, options);
    const auto matrix = generator.Compute(safe.Ai(), safe.bi(), W.Ai(), W.bi(), options);
    ASSERT_EQ(expected.size(), 1u);
    ASSERT_EQ(ignored.size(), 1u);
    ASSERT_EQ(matrix.size(), 1u);
    ExpectSameRepresentation(ignored.front().set, expected.front().set);
    ExpectSameRepresentation(matrix.front().set, expected.front().set);
    EXPECT_EQ(ignored.front().lifted_disturbance.cols(), 0);
    VerifyLiftedInvariance(ignored.front());
}

TEST(ControlledInvariantSetGeneratorTest, ScalarNominalSetAndInputMapsAreExact) {
    MatrixXd A = MatrixXd::Constant(1, 1, 0.0);
    MatrixXd B = MatrixXd::Constant(1, 1, 1.0);
    const ControlledInvariantSetGenerator generator(A, B);
    const HPolyhedron safe = SymmetricBox(2, 1.0);
    const auto components = generator.Compute(safe, Single(0, 1));
    ASSERT_EQ(components.size(), 1u);
    const CISComponent& result = components.front();
    EXPECT_EQ(result.set.Dimension(), 2u);
    EXPECT_EQ(result.tau, 0u);
    EXPECT_EQ(result.lambda, 1u);
    EXPECT_TRUE(result.set.IsSubsetOf(safe));
    EXPECT_TRUE(safe.IsSubsetOf(result.set));
    ExpectMatrixNear(result.lifted_dynamics, (MatrixXd(2, 2) << 0, 1, 0, 1).finished());
    ExpectMatrixNear(result.input_from_state, MatrixXd::Zero(1, 1));
    ExpectMatrixNear(result.input_from_virtual, MatrixXd::Identity(1, 1));
    EXPECT_EQ(result.lifted_disturbance.cols(), 0);
    VerifyLiftedInvariance(result);

    const auto explicit_result = generator.Compute(safe, Single(0, 1, false));
    ASSERT_EQ(explicit_result.size(), 1u);
    EXPECT_EQ(explicit_result.front().set.Dimension(), 1u);
    EXPECT_TRUE(explicit_result.front().set.IsSubsetOf(SymmetricBox(1, 1.0)));
    EXPECT_TRUE(SymmetricBox(1, 1.0).IsSubsetOf(explicit_result.front().set));
    VerifyProjectedControlledInvariance(explicit_result.front().set, safe, A, B);
}

TEST(ControlledInvariantSetGeneratorTest, ScalarAsymmetricDisturbanceShrinksAtTimeOne) {
    MatrixXd A = MatrixXd::Zero(1, 1);
    MatrixXd B = MatrixXd::Identity(1, 1);
    MatrixXd E = MatrixXd::Identity(1, 1);
    const ControlledInvariantSetGenerator generator(A, B, E);
    const HPolyhedron safe = SymmetricBox(2, 1.0);
    VectorXd lower(1), upper(1);
    lower << -0.2;
    upper << 0.3;
    const HPolyhedron W = Box(lower, upper);
    const auto result = generator.Compute(safe, W, Single(0, 1));
    ASSERT_EQ(result.size(), 1u);
    VectorXd expected_lower(2), expected_upper(2);
    expected_lower << -1.0, -0.8;
    expected_upper << 1.0, 0.7;
    const HPolyhedron expected = Box(expected_lower, expected_upper);
    EXPECT_TRUE(result.front().set.IsSubsetOf(expected, 1e-7));
    EXPECT_TRUE(expected.IsSubsetOf(result.front().set, 1e-7));
    VerifyLiftedInvariance(result.front(), W);

    const auto projected = generator.Compute(safe.Ai(), safe.bi(), W.Ai(), W.bi(),
                                             Single(0, 1, false));
    ASSERT_EQ(projected.size(), 1u);
    EXPECT_TRUE(projected.front().set.IsSubsetOf(SymmetricBox(1, 1.0)));
    EXPECT_TRUE(SymmetricBox(1, 1.0).IsSubsetOf(projected.front().set));
    VerifyProjectedControlledInvariance(projected.front().set, safe, A, B, E, W);
}

TEST(ControlledInvariantSetGeneratorTest, TranslatedDisturbanceWithoutZeroIsHandled) {
    const MatrixXd A = MatrixXd::Zero(1, 1);
    const MatrixXd B = MatrixXd::Identity(1, 1);
    const ControlledInvariantSetGenerator generator(A, B, B);
    const HPolyhedron safe = SymmetricBox(2, 1.0);
    VectorXd lower(1), upper(1);
    lower << 0.2;
    upper << 0.3;
    const HPolyhedron W = Box(lower, upper);
    const auto results = generator.Compute(safe, W, Single(0, 1));
    ASSERT_EQ(results.size(), 1u);
    VectorXd expected_lower(2), expected_upper(2);
    expected_lower << -1.0, -1.0;
    expected_upper << 1.0, 0.7;
    const HPolyhedron expected = Box(expected_lower, expected_upper);
    EXPECT_TRUE(results.front().set.IsSubsetOf(expected));
    EXPECT_TRUE(expected.IsSubsetOf(results.front().set));
    VerifyLiftedInvariance(results.front(), W);
}

TEST(ControlledInvariantSetGeneratorTest, TwoStepDisturbanceShrinkReachesAndRepeatsItsLimit) {
    MatrixXd A(2, 2);
    A << 0, 1, 0, 0;
    MatrixXd B(2, 1);
    B << 0, 1;
    const MatrixXd E = B;
    const HPolyhedron box = SymmetricBox(3, 2.0);
    MatrixXd G(box.Ai().rows() + 1, 3);
    G.topRows(box.Ai().rows()) = box.Ai();
    G.bottomRows(1) << 1, 1, 0;
    VectorXd F(box.bi().size() + 1);
    F.head(box.bi().size()) = box.bi();
    F.tail(1) << 2;
    const HPolyhedron safe(G, F);
    VectorXd lower(1), upper(1);
    lower << -0.1;
    upper << 0.2;
    const HPolyhedron W = Box(lower, upper);
    const ControlledInvariantSetGenerator generator(A, B, E);
    const auto results = generator.Compute(safe, W, Single(0, 2));
    ASSERT_EQ(results.size(), 1u);
    const CISComponent& result = results.front();
    ASSERT_EQ(result.set.NumInequalities(), 28u);
    const Index mixed_row = G.rows() - 1;
    EXPECT_NEAR(result.set.bi()(mixed_row), 2.0, 1e-10); // S_0.
    EXPECT_NEAR(result.set.bi()(G.rows() + mixed_row), 1.8, 1e-10); // S_1.
    EXPECT_NEAR(result.set.bi()(2 * G.rows() + mixed_row), 1.6, 1e-10); // S_2 = S_infinity.
    EXPECT_NEAR(result.set.bi()(3 * G.rows() + mixed_row), 1.6, 1e-10); // Plateau.
    VerifyLiftedInvariance(result, W);
}

TEST(ControlledInvariantSetGeneratorTest, StateOnlyConstraintsHaveFreeInputAndMatchJointForm) {
    MatrixXd A(2, 2);
    A << 0, 1, 0, 0;
    MatrixXd B(2, 1);
    B << 0, 1;
    const ControlledInvariantSetGenerator generator(A, B);
    const HPolyhedron state_safe = SymmetricBox(2, 1.0);
    MatrixXd Gjoint = MatrixXd::Zero(state_safe.Ai().rows(), 3);
    Gjoint.leftCols(2) = state_safe.Ai();
    const HPolyhedron joint_safe(Gjoint, state_safe.bi());
    const CISOptions options = Single(1, 2);
    const auto state_result = generator.Compute(state_safe, options);
    const auto joint_result = generator.Compute(joint_safe, options);
    const auto matrix_result = generator.Compute(state_safe.Ai(), state_safe.bi(), options);
    ASSERT_EQ(state_result.size(), 1u);
    ASSERT_EQ(joint_result.size(), 1u);
    ASSERT_EQ(matrix_result.size(), 1u);
    EXPECT_EQ(state_result.front().set.Dimension(), 5u);
    ExpectSameRepresentation(state_result.front().set, joint_result.front().set);
    ExpectSameRepresentation(state_result.front().set, matrix_result.front().set);
    VerifyLiftedInvariance(state_result.front());
}

TEST(ControlledInvariantSetGeneratorTest, StateOnlyDisturbedSafeSetMatchesJointForm) {
    MatrixXd A(2, 2);
    A << 0, 1, 0, 0;
    MatrixXd B(2, 1);
    B << 0, 1;
    const ControlledInvariantSetGenerator generator(A, B, B);
    const HPolyhedron state_safe = SymmetricBox(2, 1.0);
    MatrixXd joint_G = MatrixXd::Zero(state_safe.Ai().rows(), 3);
    joint_G.leftCols(2) = state_safe.Ai();
    const HPolyhedron joint_safe(joint_G, state_safe.bi());
    const HPolyhedron W = SymmetricBox(1, 0.1);
    const auto state_result = generator.Compute(state_safe, W, Single(0, 2));
    const auto joint_result = generator.Compute(joint_safe, W, Single(0, 2));
    ASSERT_EQ(state_result.size(), 1u);
    ASSERT_EQ(joint_result.size(), 1u);
    ExpectSameRepresentation(state_result.front().set, joint_result.front().set);
    VerifyLiftedInvariance(state_result.front(), W);
}

TEST(ControlledInvariantSetGeneratorTest, TwoStateDisturbedProjectionIsControlledInvariant) {
    MatrixXd A(2, 2);
    A << 0, 1, 0, 0;
    MatrixXd B(2, 1);
    B << 0, 1;
    const MatrixXd E = B;
    const ControlledInvariantSetGenerator generator(A, B, E);
    const HPolyhedron state_safe = SymmetricBox(2, 1.0);
    MatrixXd joint_G = MatrixXd::Zero(state_safe.Ai().rows(), 3);
    joint_G.leftCols(2) = state_safe.Ai();
    const HPolyhedron joint_safe(joint_G, state_safe.bi());
    const HPolyhedron W = SymmetricBox(1, 0.1);
    const auto results = generator.Compute(state_safe, W, Single(0, 2, false));
    ASSERT_EQ(results.size(), 1u);
    EXPECT_EQ(results.front().set.Dimension(), 2u);
    EXPECT_TRUE(results.front().set.Contains(VectorXd::Zero(2)));
    VerifyProjectedControlledInvariance(results.front().set, joint_safe, A, B, E, W);
}

TEST(ControlledInvariantSetGeneratorTest, MixedTwoStateLiftedAndProjectedSetsAreInvariant) {
    MatrixXd A(2, 2);
    A << 0.4, 1.0, 0.2, -0.3;
    MatrixXd B(2, 1);
    B << 0.5, 1.0;
    MatrixXd E(2, 2);
    E << 0.1, 0.04,
        -0.15, 0.08;

    const HPolyhedron box = SymmetricBox(3, 1.5);
    MatrixXd G(box.Ai().rows() + 2, 3);
    G.topRows(box.Ai().rows()) = box.Ai();
    G.bottomRows(2) << 0.6, -0.4, 0.8,
                       -0.3, 0.7, -0.5;
    VectorXd F(box.bi().size() + 2);
    F.head(box.bi().size()) = box.bi();
    F.tail(2) << 1.0, 1.1;
    const HPolyhedron safe(G, F);
    const HPolyhedron disturbance_box = SymmetricBox(2, 0.1);
    MatrixXd Gw(disturbance_box.Ai().rows() + 1, 2);
    Gw.topRows(disturbance_box.Ai().rows()) = disturbance_box.Ai();
    Gw.bottomRows(1) << 1.0, 1.0;
    VectorXd Fw(disturbance_box.bi().size() + 1);
    Fw.head(disturbance_box.bi().size()) = disturbance_box.bi();
    Fw.tail(1) << 0.1;
    const HPolyhedron W(Gw, Fw);
    EXPECT_FALSE(W.Contains(VectorXd::Constant(2, 0.1)));
    EXPECT_TRUE(W.Contains(VectorXd::Zero(2)));

    const ControlledInvariantSetGenerator nominal(A, B);
    const ControlledInvariantSetGenerator disturbed(A, B, E);
    CISOptions implicit_options;
    implicit_options.hierarchy_level = 2;
    CISOptions explicit_options = implicit_options;
    explicit_options.is_implicit = false;
    const auto nominal_lifted = nominal.Compute(safe, implicit_options);
    const auto nominal_projected = nominal.Compute(safe, explicit_options);
    const auto robust_lifted = disturbed.Compute(safe, W, implicit_options);
    const auto robust_projected = disturbed.Compute(safe, W, explicit_options);
    ASSERT_EQ(nominal_lifted.size(), 2u);
    ASSERT_EQ(nominal_projected.size(), 2u);
    ASSERT_EQ(robust_lifted.size(), 2u);
    ASSERT_EQ(robust_projected.size(), 2u);

    const VectorXd interior = (VectorXd(2) << 0.1, -0.1).finished();
    for (std::size_t i = 0; i < 2; i++) {
        SCOPED_TRACE(i);
        for (const auto* component : {&nominal_projected[i], &robust_projected[i]}) {
            EXPECT_EQ(component->set.Dimension(), 2u);
            EXPECT_EQ(component->tau, 1u - i);
            EXPECT_EQ(component->lambda, i + 1u);
            EXPECT_TRUE(component->set.IsFeasible());
            EXPECT_TRUE(component->set.Contains(interior));
            EXPECT_FALSE(component->set.Contains(VectorXd::Constant(2, 2.0)));
        }
        EXPECT_EQ(nominal_lifted[i].set.Dimension(), 4u);
        EXPECT_EQ(robust_lifted[i].set.Dimension(), 4u);
        VerifyLiftedInvariance(nominal_lifted[i]);
        VerifyProjectedControlledInvariance(nominal_projected[i].set, safe, A, B);
        VerifyLiftedInvariance(robust_lifted[i], W);
        VerifyProjectedControlledInvariance(robust_projected[i].set, safe, A, B, E, W);
        EXPECT_TRUE(robust_lifted[i].set.IsSubsetOf(nominal_lifted[i].set, 1e-6));
        EXPECT_TRUE(robust_projected[i].set.IsSubsetOf(nominal_projected[i].set, 1e-6));
    }
}

TEST(ControlledInvariantSetGeneratorTest, FullAndEmptySafeSetsStayDistinct) {
    const ControlledInvariantSetGenerator generator(MatrixXd::Zero(1, 1), MatrixXd::Identity(1, 1));
    const auto full = generator.Compute(HPolyhedron::FullSpace(1), Single(0, 1));
    ASSERT_EQ(full.size(), 1u);
    EXPECT_TRUE(full.front().set.IsFullSpace());
    EXPECT_EQ(full.front().set.Dimension(), 2u);
    const auto empty = generator.Compute(HPolyhedron::EmptySet(2), Single(0, 1));
    ASSERT_EQ(empty.size(), 1u);
    EXPECT_TRUE(empty.front().set.IsValid());
    EXPECT_FALSE(empty.front().set.IsFeasible());
}

TEST(ControlledInvariantSetGeneratorTest, FullSpaceSafeSetReturnsFullSpaceAtEveryHierarchyLevel) {
    MatrixXd A(2, 2);
    A << 0.0, 1.0, 0.0, 0.0;
    MatrixXd B(2, 1);
    B << 0.0, 1.0;
    const ControlledInvariantSetGenerator generator(A, B, B);
    const HPolyhedron W = SymmetricBox(1, 0.1);
    CISOptions options;
    options.hierarchy_level = 3;
    for (bool implicit : {true, false}) {
        options.is_implicit = implicit;
        for (std::size_t safe_dimension : {2u, 3u}) {
            const auto components = generator.Compute(HPolyhedron::FullSpace(safe_dimension),
                                                      W, options);
            ASSERT_EQ(components.size(), 3u);
            for (std::size_t i = 0; i < components.size(); i++) {
                const CISComponent& component = components[i];
                EXPECT_TRUE(component.set.IsFullSpace());
                EXPECT_EQ(component.set.Dimension(), implicit ? 5u : 2u);
                EXPECT_EQ(component.tau, 2u - i);
                EXPECT_EQ(component.lambda, i + 1u);
                EXPECT_EQ(component.lifted_dynamics.rows(), 5);
                EXPECT_EQ(component.lifted_disturbance.rows(), 5);
                EXPECT_EQ(component.lifted_disturbance.cols(), 1);
                ExpectMatrixNear(component.lifted_disturbance.topRows(2), B);
            }
        }
    }
}

TEST(ControlledInvariantSetGeneratorTest, EqualityConstraintsCanMakeRobustSetEmpty) {
    const MatrixXd A = MatrixXd::Zero(1, 1);
    const MatrixXd B = MatrixXd::Identity(1, 1);
    MatrixXd Ae = MatrixXd::Identity(2, 2);
    const HPolyhedron singleton(MatrixXd(0, 2), VectorXd(0), Ae, VectorXd::Zero(2));
    const ControlledInvariantSetGenerator nominal(A, B);
    const auto nominal_result = nominal.Compute(singleton, Single(0, 1));
    ASSERT_EQ(nominal_result.size(), 1u);
    EXPECT_TRUE(nominal_result.front().set.IsFeasible());
    EXPECT_EQ(nominal_result.front().set.NumEqualities(), 4u);
    VerifyLiftedInvariance(nominal_result.front());

    const ControlledInvariantSetGenerator disturbed(A, B, MatrixXd::Identity(1, 1));
    const auto robust_result = disturbed.Compute(singleton, SymmetricBox(1, 0.1), Single(0, 1));
    ASSERT_EQ(robust_result.size(), 1u);
    EXPECT_TRUE(robust_result.front().set.IsValid());
    EXPECT_FALSE(robust_result.front().set.IsFeasible());
    const auto projected = disturbed.Compute(singleton, SymmetricBox(1, 0.1),
                                             Single(0, 1, false));
    ASSERT_EQ(projected.size(), 1u);
    EXPECT_EQ(projected.front().set.Dimension(), 1u);
    EXPECT_TRUE(projected.front().set.IsValid());
    EXPECT_FALSE(projected.front().set.IsFeasible());
}

TEST(ControlledInvariantSetGeneratorTest, SingletonDisturbancePreservesRobustEqualities) {
    const MatrixXd A = MatrixXd::Zero(1, 1);
    const MatrixXd B = MatrixXd::Identity(1, 1);
    MatrixXd Ae = MatrixXd::Identity(2, 2);
    VectorXd be(2);
    be << 0.0, -0.2;
    const HPolyhedron safe(MatrixXd(0, 2), VectorXd(0), Ae, be);
    MatrixXd disturbance_equality = MatrixXd::Identity(1, 1);
    VectorXd disturbance_value(1);
    disturbance_value << 0.2;
    const HPolyhedron W(MatrixXd(0, 1), VectorXd(0),
                        disturbance_equality, disturbance_value);
    const ControlledInvariantSetGenerator generator(A, B, B);
    const auto results = generator.Compute(safe, W, Single(0, 1));
    ASSERT_EQ(results.size(), 1u);
    const CISComponent& result = results.front();
    ASSERT_TRUE(result.set.IsFeasible());
    VectorXd expected(2);
    expected << 0.0, -0.2;
    EXPECT_TRUE(result.set.Contains(expected));
    EXPECT_FALSE(result.set.Contains(VectorXd::Zero(2)));
    VerifyLiftedInvariance(result, W);
}

TEST(ControlledInvariantSetGeneratorTest, ZeroDisturbanceSetMatchesNominalConstruction) {
    MatrixXd A(2, 2);
    A << 0.2, 1.0, 0.1, 0.3;
    MatrixXd B(2, 1);
    B << 0, 1;
    MatrixXd E(2, 1);
    E << 0.4, -0.2;
    const HPolyhedron safe = SymmetricBox(3, 1.5);
    const HPolyhedron zero_disturbance(
        MatrixXd(0, 1), VectorXd(0), MatrixXd::Identity(1, 1), VectorXd::Zero(1));
    const ControlledInvariantSetGenerator nominal(A, B);
    const ControlledInvariantSetGenerator disturbed(A, B, E);
    const CISOptions options = Single(1, 2);
    const auto nominal_result = nominal.Compute(safe, options);
    const auto disturbed_result = disturbed.Compute(safe, zero_disturbance, options);
    ASSERT_EQ(nominal_result.size(), 1u);
    ASSERT_EQ(disturbed_result.size(), 1u);
    ExpectSameRepresentation(disturbed_result.front().set, nominal_result.front().set);
    ExpectMatrixNear(disturbed_result.front().lifted_dynamics,
                     nominal_result.front().lifted_dynamics);
    VerifyLiftedInvariance(disturbed_result.front(), zero_disturbance);
}

TEST(ControlledInvariantSetGeneratorTest, HierarchyReturnsEachLassoComponentInOrder) {
    const ControlledInvariantSetGenerator generator(MatrixXd::Zero(1, 1), MatrixXd::Identity(1, 1));
    const HPolyhedron safe = SymmetricBox(2, 1.0);
    CISOptions options;
    options.hierarchy_level = 3;
    options.lambda = 0; // Ignored in hierarchy mode.
    options.tau = std::numeric_limits<std::size_t>::max();
    const auto hierarchy = generator.Compute(safe, options);
    ASSERT_EQ(hierarchy.size(), 3u);
    for (std::size_t i = 0; i < hierarchy.size(); i++) {
        EXPECT_EQ(hierarchy[i].lambda, i + 1);
        EXPECT_EQ(hierarchy[i].tau, 2 - i);
        EXPECT_EQ(hierarchy[i].set.Dimension(), 4u);
        const auto individual = generator.Compute(safe, Single(2 - i, i + 1));
        ASSERT_EQ(individual.size(), 1u);
        ExpectSameRepresentation(hierarchy[i].set, individual.front().set);
        ExpectMatrixNear(hierarchy[i].lifted_dynamics, individual.front().lifted_dynamics);
        VerifyLiftedInvariance(hierarchy[i]);
    }
}

TEST(ControlledInvariantSetGeneratorTest, ExplicitHierarchyMatchesIndividualProjections) {
    const ControlledInvariantSetGenerator generator(MatrixXd::Zero(1, 1), MatrixXd::Identity(1, 1));
    const HPolyhedron safe = SymmetricBox(2, 1.0);
    CISOptions options;
    options.hierarchy_level = 3;
    options.is_implicit = false;
    const auto hierarchy = generator.Compute(safe, options);
    ASSERT_EQ(hierarchy.size(), 3u);
    for (std::size_t i = 0; i < hierarchy.size(); i++) {
        ASSERT_FALSE(hierarchy[i].is_implicit);
        EXPECT_EQ(hierarchy[i].set.Dimension(), 1u);
        const auto individual = generator.Compute(safe, Single(2 - i, i + 1, false));
        ASSERT_EQ(individual.size(), 1u);
        EXPECT_TRUE(hierarchy[i].set.IsSubsetOf(individual.front().set));
        EXPECT_TRUE(individual.front().set.IsSubsetOf(hierarchy[i].set));
    }
}

TEST(ControlledInvariantSetGeneratorTest, ProjectedSetsAreControlledInvariantAndMonotone) {
    MatrixXd A = MatrixXd::Constant(1, 1, 0.7);
    MatrixXd B = MatrixXd::Constant(1, 1, 1.0);
    const ControlledInvariantSetGenerator generator(A, B);
    const HPolyhedron safe = SymmetricBox(2, 1.0);
    const auto short_lasso = generator.Compute(safe, Single(0, 1, false));
    const auto longer_transient = generator.Compute(safe, Single(1, 1, false));
    const auto longer_period = generator.Compute(safe, Single(0, 2, false));
    ASSERT_EQ(short_lasso.size(), 1u);
    ASSERT_EQ(longer_transient.size(), 1u);
    ASSERT_EQ(longer_period.size(), 1u);
    VerifyProjectedControlledInvariance(short_lasso.front().set, safe, A, B);
    VerifyProjectedControlledInvariance(longer_transient.front().set, safe, A, B);
    VerifyProjectedControlledInvariance(longer_period.front().set, safe, A, B);
    EXPECT_TRUE(short_lasso.front().set.IsSubsetOf(longer_transient.front().set, 1e-6));
    EXPECT_TRUE(short_lasso.front().set.IsSubsetOf(longer_period.front().set, 1e-6));
}

TEST(ControlledInvariantSetGeneratorTest, MixedConstraintsAndFeedbackMatchFiniteTrajectoryOracle) {
    MatrixXd A = MatrixXd::Constant(1, 1, 0.7);
    MatrixXd B = MatrixXd::Constant(1, 1, 1.5);
    MatrixXd E = MatrixXd::Constant(1, 1, 0.5);
    MatrixXd G(6, 2);
    G << 1, 0, -1, 0, 0, 1, 0, -1, 1, 0.6, -0.4, 1;
    VectorXd F(6);
    F << 1.2, 1.2, 1.0, 1.0, 1.3, 1.1;
    const HPolyhedron safe(G, F);
    VectorXd lower(1), upper(1), center(1), radius(1);
    lower << -0.1;
    upper << 0.2;
    center << 0.05;
    radius << 0.15;
    const HPolyhedron W = Box(lower, upper);
    const ControlledInvariantSetGenerator generator(A, B, E);
    const auto results = generator.Compute(safe, W, Single(1, 2));
    ASSERT_EQ(results.size(), 1u);
    const CISComponent& result = results.front();
    ASSERT_TRUE(result.set.IsFeasible());
    EXPECT_EQ(result.set.Dimension(), 4u);
    VerifyLiftedInvariance(result, W);
    ExpectMatrixNear(result.lifted_dynamics.topLeftCorner(1, 1) +
                         B * (-result.input_from_state), A);

    std::mt19937 rng(71);
    std::uniform_real_distribution<double> distribution(-1.5, 1.5);
    int accepted = 0;
    int rejected = 0;
    for (int trial = 0; trial < 120; trial++) {
        VectorXd initial(4);
        for (Index i = 0; i < initial.size(); i++) {
            initial(i) = distribution(rng);
        }
        const bool expected = FiniteTrajectoryOracle(result, safe, initial, 1, 4, center, radius);
        accepted += expected;
        rejected += !expected;
        EXPECT_EQ(result.set.Contains(initial, 2e-7), expected) << "trial=" << trial;
    }
    EXPECT_GT(accepted, 0);
    EXPECT_GT(rejected, 0);

    const auto projected = generator.Compute(safe, W, Single(1, 2, false));
    ASSERT_EQ(projected.size(), 1u);
    VerifyProjectedControlledInvariance(projected.front().set, safe, A, B, E, W);
}

TEST(ControlledInvariantSetGeneratorTest, NontrivialMultiInputTransformationAndRobustInvariance) {
    MatrixXd Ac = MatrixXd::Zero(3, 3);
    Ac(0, 1) = 1.0;
    MatrixXd Bc = MatrixXd::Zero(3, 2);
    Bc(1, 0) = 1.0;
    Bc(2, 1) = 1.0;
    MatrixXd T(3, 3);
    T << 1, 0.2, -0.1, 0.1, 1, 0.3, -0.2, 0.1, 1;
    MatrixXd Am(2, 3);
    Am << 0.3, -0.2, 0.1, -0.1, 0.2, 0.4;
    MatrixXd Bm(2, 2);
    Bm << 1.2, 0.2, -0.1, 0.9;
    const MatrixXd A = T.inverse() * (Ac + Bc * Am) * T;
    const MatrixXd B = T.inverse() * Bc * Bm;
    MatrixXd E(3, 2);
    E << 0.1, 0, 0, 0.2, -0.1, 0.1;

    const HPolyhedron base_safe = SymmetricBox(5, 2.5);
    MatrixXd G(base_safe.Ai().rows() + 2, 5);
    G.topRows(base_safe.Ai().rows()) = base_safe.Ai();
    G.bottomRows(2) << 0.5, -0.2, 0.3, 0.8, -0.7,
                       -0.4, 0.2, 0.5, -0.3, 0.6;
    VectorXd F(base_safe.bi().size() + 2);
    F.head(base_safe.bi().size()) = base_safe.bi();
    F.tail(2).setConstant(2.0);
    const HPolyhedron safe(G, F);
    const HPolyhedron W = SymmetricBox(2, 0.2);
    const ControlledInvariantSetGenerator generator(A, B, E);
    const auto results = generator.Compute(safe, W, Single(1, 2));
    ASSERT_EQ(results.size(), 1u);
    const CISComponent& result = results.front();
    ASSERT_TRUE(result.set.IsFeasible());
    EXPECT_EQ(generator.Transformation().MaxControllabilityIndex(), 2u);
    EXPECT_EQ(result.set.Dimension(), 9u);
    EXPECT_EQ(result.lifted_disturbance.rows(), 9);
    ExpectMatrixNear(result.lifted_dynamics.topLeftCorner(3, 3),
                     A + B * result.input_from_state);
    ExpectMatrixNear(result.lifted_dynamics.topRightCorner(3, 6),
                     B * result.input_from_virtual);
    MatrixXd expected_H = MatrixXd::Zero(2, 6);
    expected_H(0, 0) = 1.0;
    expected_H(1, 3) = 1.0;
    MatrixXd expected_P = MatrixXd::Zero(6, 6);
    for (Index channel = 0; channel < 2; channel++) {
        const Index first = 3 * channel;
        expected_P(first, first + 1) = 1.0;
        expected_P(first + 1, first + 2) = 1.0;
        expected_P(first + 2, first + 1) = 1.0;
    }
    ExpectMatrixNear(result.lifted_dynamics.bottomRightCorner(6, 6), expected_P);
    ExpectMatrixNear(generator.Transformation().InputTransformationMatrix() *
                         result.input_from_virtual, expected_H);
    VerifyLiftedInvariance(result, W);

    std::mt19937 rng(93);
    std::uniform_real_distribution<double> small_distribution(-0.3, 0.3);
    std::uniform_real_distribution<double> wide_distribution(-3.0, 3.0);
    const VectorXd center = VectorXd::Zero(2);
    const VectorXd radius = VectorXd::Constant(2, 0.2);
    int accepted = 0;
    int rejected = 0;
    for (int trial = 0; trial < 100; trial++) {
        VectorXd initial(9);
        for (Index i = 0; i < initial.size(); i++) {
            initial(i) = trial < 50 ? small_distribution(rng) : wide_distribution(rng);
        }
        const bool expected = FiniteTrajectoryOracle(result, safe, initial, 3, 5, center, radius);
        accepted += expected;
        rejected += !expected;
        EXPECT_EQ(result.set.Contains(initial, 2e-7), expected) << "trial=" << trial;
    }
    EXPECT_GT(accepted, 0);
    EXPECT_GT(rejected, 0);

    // Two physical input channels remain coupled through the mixed safe set
    // after projecting the one-sample lifted controller state.
    const auto projected = generator.Compute(safe, W, Single(0, 1, false));
    ASSERT_EQ(projected.size(), 1u);
    EXPECT_EQ(projected.front().set.Dimension(), 3u);
    EXPECT_TRUE(projected.front().set.Contains(VectorXd::Zero(3)));
    VerifyProjectedControlledInvariance(projected.front().set, safe, A, B, E, W);
}

TEST(ControlledInvariantSetGeneratorTest, DenseHighDimensionalSingleInputImplicitSetsAreInvariant) {
    const std::vector<Index> dimensions = {6, 8, 10, 15, 20};
    const std::vector<std::size_t> transients = {0, 1, 1, 0, 0};
    const std::vector<std::size_t> periods = {1, 2, 3, 1, 1};
    for (Index trial = 0; trial < static_cast<Index>(dimensions.size()); trial++) {
        const Index n = dimensions[static_cast<std::size_t>(trial)];
        const auto problem = DenseTransformedProblem(n, 1);
        const std::size_t tau = transients[static_cast<std::size_t>(trial)];
        const std::size_t lambda = periods[static_cast<std::size_t>(trial)];
        const std::size_t q = tau + lambda;

        ASSERT_GT(CountSignificantEntries(problem.A), n * n / 2) << "n=" << n;
        ASSERT_GT(CountSignificantEntries(problem.safe_set.Ai()),
                  problem.safe_set.Ai().size() / 2) << "n=" << n;
        const ControlledInvariantSetGenerator nominal(problem.A, problem.B);
        const ControlledInvariantSetGenerator disturbed(problem.A, problem.B, problem.E);
        ASSERT_EQ(nominal.Transformation().MaxControllabilityIndex(),
                  static_cast<std::size_t>(n)) << "n=" << n;
        ASSERT_EQ(disturbed.Transformation().MaxControllabilityIndex(),
                  static_cast<std::size_t>(n)) << "n=" << n;

        const auto nominal_results = nominal.Compute(
            problem.safe_set, Single(tau, lambda));
        const auto robust_results = disturbed.Compute(
            problem.safe_set, problem.disturbance_set, Single(tau, lambda));
        ASSERT_EQ(nominal_results.size(), 1u) << "n=" << n;
        ASSERT_EQ(robust_results.size(), 1u) << "n=" << n;
        const CISComponent& nominal_result = nominal_results.front();
        const CISComponent& robust_result = robust_results.front();
        EXPECT_TRUE(nominal_result.is_implicit) << "n=" << n;
        EXPECT_TRUE(robust_result.is_implicit) << "n=" << n;
        EXPECT_TRUE(nominal_result.set.IsValid()) << "n=" << n;
        EXPECT_TRUE(robust_result.set.IsValid()) << "n=" << n;
        ASSERT_TRUE(nominal_result.set.IsFeasible()) << "n=" << n;
        ASSERT_TRUE(robust_result.set.IsFeasible()) << "n=" << n;
        EXPECT_EQ(nominal_result.set.Dimension(), static_cast<std::size_t>(n) + q);
        EXPECT_EQ(robust_result.set.Dimension(), static_cast<std::size_t>(n) + q);
        VerifyLiftedInvariance(nominal_result);
        VerifyLiftedInvariance(robust_result, problem.disturbance_set);
    }
}

TEST(ControlledInvariantSetGeneratorTest, DenseHighDimensionalTwoInputImplicitSetsAreInvariant) {
    const std::vector<Index> dimensions = {15, 20};
    for (Index trial = 0; trial < static_cast<Index>(dimensions.size()); trial++) {
        const Index n = dimensions[static_cast<std::size_t>(trial)];
        constexpr Index m = 2;
        const auto problem = DenseTransformedProblem(n, m);
        const std::size_t tau = static_cast<std::size_t>(trial);
        constexpr std::size_t lambda = 1;
        const std::size_t q = tau + lambda;

        ASSERT_GT(CountSignificantEntries(problem.A), n * n / 2) << "n=" << n;
        ASSERT_GT(CountSignificantEntries(problem.safe_set.Ai()),
                  problem.safe_set.Ai().size() / 2) << "n=" << n;
        const ControlledInvariantSetGenerator nominal(problem.A, problem.B);
        const ControlledInvariantSetGenerator disturbed(problem.A, problem.B, problem.E);
        const std::size_t expected_max_index = static_cast<std::size_t>((n + m - 1) / m);
        ASSERT_EQ(nominal.Transformation().MaxControllabilityIndex(), expected_max_index)
            << "n=" << n;
        ASSERT_EQ(disturbed.Transformation().MaxControllabilityIndex(), expected_max_index)
            << "n=" << n;

        const auto nominal_results = nominal.Compute(
            problem.safe_set, Single(tau, lambda));
        const auto robust_results = disturbed.Compute(
            problem.safe_set, problem.disturbance_set, Single(tau, lambda));
        ASSERT_EQ(nominal_results.size(), 1u) << "n=" << n;
        ASSERT_EQ(robust_results.size(), 1u) << "n=" << n;
        const CISComponent& nominal_result = nominal_results.front();
        const CISComponent& robust_result = robust_results.front();
        EXPECT_TRUE(nominal_result.is_implicit) << "n=" << n;
        EXPECT_TRUE(robust_result.is_implicit) << "n=" << n;
        EXPECT_TRUE(nominal_result.set.IsValid()) << "n=" << n;
        EXPECT_TRUE(robust_result.set.IsValid()) << "n=" << n;
        ASSERT_TRUE(nominal_result.set.IsFeasible()) << "n=" << n;
        ASSERT_TRUE(robust_result.set.IsFeasible()) << "n=" << n;
        EXPECT_EQ(nominal_result.set.Dimension(),
                  static_cast<std::size_t>(n) + static_cast<std::size_t>(m) * q);
        EXPECT_EQ(robust_result.set.Dimension(),
                  static_cast<std::size_t>(n) + static_cast<std::size_t>(m) * q);
        VerifyLiftedInvariance(nominal_result);
        VerifyLiftedInvariance(robust_result, problem.disturbance_set);
    }
}

TEST(ControlledInvariantSetGeneratorTest, DISABLED_DenseTransformedExplicitProjectionNumericalRegression) {
    // Disabled until Projection classifies roundoff-scale coefficients using a
    // scale-aware zero test. Exact sign checks currently treat coefficients as
    // small as 1e-19 as genuine Fourier-Motzkin bounds, create thousands of
    // ill-conditioned row pairs, and can leave SCIP reporting unresolved
    // numerical trouble for several minutes.
    constexpr Index n = 6;
    const auto problem = DenseTransformedProblem(n, 1);
    const ControlledInvariantSetGenerator generator(problem.A, problem.B, problem.E);
    const auto projected = generator.Compute(
        problem.safe_set, problem.disturbance_set, Single(0, 1, false));
    ASSERT_EQ(projected.size(), 1u);
    EXPECT_TRUE(projected.front().set.IsValid());
    EXPECT_TRUE(projected.front().set.IsFeasible());
}

TEST(ControlledInvariantSetGeneratorTest, TenStateSingleInputLiftedAndProjectedSetsAreInvariant) {
    constexpr Index n = 10;
    MatrixXd A = MatrixXd::Zero(n, n);
    for (Index i = 0; i < n - 1; i++) {
        A(i, i + 1) = 1.0;
    }
    MatrixXd B = MatrixXd::Zero(n, 1);
    B(n - 1, 0) = 1.0;
    MatrixXd E = MatrixXd::Zero(n, 1);
    E(n - 1, 0) = 0.08;
    E(n - 3, 0) = 0.04;
    const HPolyhedron box = SymmetricBox(n + 1, 2.0);
    // This active facet couples the first state with the physical input.
    MatrixXd G(box.Ai().rows() + 1, n + 1);
    G.topRows(box.Ai().rows()) = box.Ai();
    G.bottomRows(1).setZero();
    G(box.Ai().rows(), 0) = 0.5;
    G(box.Ai().rows(), n) = 0.5;
    VectorXd F(box.bi().size() + 1);
    F.head(box.bi().size()) = box.bi();
    F(F.size() - 1) = 1.25;
    const HPolyhedron safe(G, F);
    VectorXd lower(1), upper(1);
    lower << -0.1;
    upper << 0.15;
    const HPolyhedron W = Box(lower, upper);

    const ControlledInvariantSetGenerator nominal(A, B);
    const ControlledInvariantSetGenerator disturbed(A, B, E);
    ASSERT_EQ(nominal.Transformation().MaxControllabilityIndex(), static_cast<std::size_t>(n));
    ASSERT_EQ(disturbed.Transformation().MaxControllabilityIndex(), static_cast<std::size_t>(n));
    const auto nominal_lifted = nominal.Compute(safe, Single(0, 1));
    const auto robust_lifted = disturbed.Compute(safe, W, Single(0, 1));
    const auto nominal_projected = nominal.Compute(safe, Single(0, 1, false));
    const auto robust_projected = disturbed.Compute(safe, W, Single(0, 1, false));
    ASSERT_EQ(nominal_lifted.size(), 1u);
    ASSERT_EQ(robust_lifted.size(), 1u);
    ASSERT_EQ(nominal_projected.size(), 1u);
    ASSERT_EQ(robust_projected.size(), 1u);

    for (const auto* component : {&nominal_lifted.front(), &robust_lifted.front()}) {
        EXPECT_EQ(component->set.Dimension(), static_cast<std::size_t>(n + 1));
        EXPECT_TRUE(component->set.Contains(VectorXd::Zero(n + 1)));
    }
    for (const auto* component : {&nominal_projected.front(), &robust_projected.front()}) {
        EXPECT_EQ(component->set.Dimension(), static_cast<std::size_t>(n));
        EXPECT_TRUE(component->set.Contains(VectorXd::Zero(n)));
        EXPECT_FALSE(component->set.Contains(VectorXd::Constant(n, 3.0)));
    }
    // The disturbance removes this nominally admissible boundary state.
    VectorXd boundary_state = VectorXd::Zero(n);
    boundary_state(n - 1) = 2.0;
    VectorXd boundary_lifted = VectorXd::Zero(n + 1);
    boundary_lifted(n - 1) = 2.0;
    EXPECT_TRUE(nominal_lifted.front().set.Contains(boundary_lifted));
    EXPECT_FALSE(robust_lifted.front().set.Contains(boundary_lifted));
    EXPECT_TRUE(nominal_projected.front().set.Contains(boundary_state));
    EXPECT_FALSE(robust_projected.front().set.Contains(boundary_state));
    VerifyLiftedInvariance(nominal_lifted.front());
    VerifyLiftedInvariance(robust_lifted.front(), W);
    VerifyProjectedControlledInvariance(nominal_projected.front().set, safe, A, B);
    VerifyProjectedControlledInvariance(robust_projected.front().set, safe, A, B, E, W);
}

TEST(ControlledInvariantSetGeneratorTest, TenStateTenSampleLiftHasExpectedPeriodicity) {
    constexpr Index n = 10;
    MatrixXd A = MatrixXd::Zero(n, n);
    for (Index i = 0; i < n - 1; i++) {
        A(i, i + 1) = 1.0;
    }
    MatrixXd B = MatrixXd::Zero(n, 1);
    B(n - 1, 0) = 1.0;
    const ControlledInvariantSetGenerator generator(A, B);
    const HPolyhedron safe = SymmetricBox(n + 1, 3.0);
    const auto results = generator.Compute(safe, Single(4, 6));
    ASSERT_EQ(results.size(), 1u);
    const CISComponent& result = results.front();
    EXPECT_EQ(result.set.Dimension(), 20u);
    EXPECT_EQ(result.set.NumInequalities(), 440u);
    EXPECT_TRUE(result.set.Contains(VectorXd::Zero(20)));

    MatrixXd at_transient = MatrixXd::Identity(20, 20);
    for (int i = 0; i < 14; i++) {
        at_transient *= result.lifted_dynamics;
    }
    MatrixXd after_period = at_transient;
    for (int i = 0; i < 6; i++) {
        after_period *= result.lifted_dynamics;
    }
    ExpectMatrixNear(after_period, at_transient);

    std::mt19937 rng(405);
    std::uniform_real_distribution<double> small_distribution(-0.5, 0.5);
    std::uniform_real_distribution<double> wide_distribution(-4.0, 4.0);
    int accepted = 0;
    int rejected = 0;
    for (int trial = 0; trial < 60; trial++) {
        VectorXd initial(20);
        for (Index i = 0; i < initial.size(); i++) {
            initial(i) = trial < 30 ? small_distribution(rng) : wide_distribution(rng);
        }
        const bool expected = FiniteTrajectoryOracle(result, safe, initial, n, 20);
        accepted += expected;
        rejected += !expected;
        EXPECT_EQ(result.set.Contains(initial, 1e-7), expected) << "trial=" << trial;
    }
    EXPECT_GT(accepted, 0);
    EXPECT_GT(rejected, 0);
}

TEST(ControlledInvariantSetGeneratorTest, RandomTransformedSystemsMatchRobustTrajectoryOracle) {
    std::mt19937 rng(1783);
    std::uniform_real_distribution<double> small(-0.08, 0.08);
    std::uniform_real_distribution<double> mixed(-0.4, 0.4);
    std::uniform_real_distribution<double> near_origin(-0.15, 0.15);
    std::uniform_real_distribution<double> wide(-2.5, 2.5);
    int accepted = 0;
    int rejected = 0;
    for (int trial = 0; trial < 12; trial++) {
        const Index n = 2 + trial % 7;
        const Index m = trial % 2 == 0 ? 1 : 2;
        const Index q = 1 + trial % 5;
        const Index tau = std::min<Index>(trial % 3, q - 1);
        const Index lambda = q - tau;

        MatrixXd Ac = MatrixXd::Zero(n, n);
        MatrixXd Bc = MatrixXd::Zero(n, m);
        Index first = 0;
        for (Index channel = 0; channel < m; channel++) {
            const Index length = channel == 0 ? n - m + 1 : 1;
            for (Index i = 0; i < length - 1; i++) {
                Ac(first + i, first + i + 1) = 1.0;
            }
            Bc(first + length - 1, channel) = 1.0;
            first += length;
        }
        MatrixXd T = MatrixXd::Identity(n, n);
        MatrixXd Bm = MatrixXd::Identity(m, m);
        MatrixXd Am(m, n);
        MatrixXd Ec(n, 1);
        for (Index i = 0; i < n; i++) {
            Ec(i, 0) = small(rng);
            for (Index j = 0; j < n; j++) {
                T(i, j) += small(rng);
            }
        }
        for (Index i = 0; i < m; i++) {
            for (Index j = 0; j < m; j++) {
                Bm(i, j) += small(rng);
            }
            for (Index j = 0; j < n; j++) {
                Am(i, j) = mixed(rng);
            }
        }
        const MatrixXd A = T.inverse() * (Ac + Bc * Am) * T;
        const MatrixXd B = T.inverse() * Bc * Bm;
        const MatrixXd E = T.inverse() * Ec;
        const HPolyhedron base_safe = SymmetricBox(n + m, 2.0);
        MatrixXd G(base_safe.Ai().rows() + 2, n + m);
        G.topRows(base_safe.Ai().rows()) = base_safe.Ai();
        for (Index row = base_safe.Ai().rows(); row < G.rows(); row++) {
            for (Index col = 0; col < G.cols(); col++) {
                G(row, col) = mixed(rng);
            }
        }
        VectorXd F(base_safe.bi().size() + 2);
        F.head(base_safe.bi().size()) = base_safe.bi();
        F.tail(2).setConstant(1.5);
        const HPolyhedron safe(G, F);
        const HPolyhedron W = SymmetricBox(1, 0.1);
        const ControlledInvariantSetGenerator generator(A, B, E);
        const auto results = generator.Compute(safe, W, Single(tau, lambda));
        ASSERT_EQ(results.size(), 1u);
        const CISComponent& result = results.front();
        ASSERT_TRUE(result.set.IsFeasible()) << "system=" << trial;
        EXPECT_EQ(result.set.Dimension(), static_cast<std::size_t>(n + m * q));
        MatrixXd expected_H = MatrixXd::Zero(m, m * q);
        for (Index channel = 0; channel < m; channel++) {
            expected_H(channel, channel * q) = 1.0;
        }
        const auto& transform = generator.Transformation();
        ExpectMatrixNear(transform.InputTransformationMatrix() * result.input_from_state +
                             transform.StateFeedbackMatrix() * transform.TransformationMatrix(),
                         MatrixXd::Zero(m, n));
        ExpectMatrixNear(transform.InputTransformationMatrix() * result.input_from_virtual,
                         expected_H);
        MatrixXd T_lift = MatrixXd::Identity(result.set.Dimension(), result.set.Dimension());
        T_lift.topLeftCorner(n, n) = transform.TransformationMatrix();
        MatrixXd companion = MatrixXd::Zero(result.set.Dimension(), result.set.Dimension());
        companion.topLeftCorner(n, n) = transform.CanonicalStateMatrix();
        companion.topRightCorner(n, m * q) = transform.CanonicalInputMatrix() * expected_H;
        companion.bottomRightCorner(m * q, m * q) =
            result.lifted_dynamics.bottomRightCorner(m * q, m * q);
        ExpectMatrixNear(T_lift * result.lifted_dynamics, companion * T_lift);
        VerifyLiftedInvariance(result, W);

        const VectorXd center = VectorXd::Zero(1);
        const VectorXd radius = VectorXd::Constant(1, 0.1);
        for (int sample = 0; sample < 24; sample++) {
            VectorXd initial(result.set.Dimension());
            for (Index i = 0; i < initial.size(); i++) {
                initial(i) = sample < 12 ? near_origin(rng) : wide(rng);
            }
            const bool expected = FiniteTrajectoryOracle(result, safe, initial, n,
                                                         n + q, center, radius);
            accepted += expected;
            rejected += !expected;
            EXPECT_EQ(result.set.Contains(initial, 2e-7), expected)
                << "system=" << trial << ", sample=" << sample;
        }
    }
    EXPECT_GT(accepted, 0);
    EXPECT_GT(rejected, 0);
}

TEST(ControlledInvariantSetGeneratorTest, RepeatedCallsDoNotRetainPreviousLassoOrDisturbanceState) {
    const ControlledInvariantSetGenerator generator(MatrixXd::Zero(1, 1), MatrixXd::Identity(1, 1));
    const HPolyhedron safe = SymmetricBox(2, 1.0);
    const auto first = generator.Compute(safe, Single(0, 1));
    const auto second = generator.Compute(safe, Single(2, 2));
    const auto again = generator.Compute(safe, Single(0, 1));
    ASSERT_EQ(first.size(), 1u);
    ASSERT_EQ(second.size(), 1u);
    ASSERT_EQ(again.size(), 1u);
    EXPECT_EQ(second.front().set.Dimension(), 5u);
    ExpectSameRepresentation(first.front().set, again.front().set);
}

}
