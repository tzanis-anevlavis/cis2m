#include "brunovskytransformation.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <gtest/gtest.h>

using cis2m::BrunovskyTransformation;
using cis2m::HPolyhedron;
using Eigen::FullPivLU;
using Eigen::Index;
using Eigen::JacobiSVD;
using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {

constexpr double kTolerance = 1e-9;

VectorXd Vec(std::initializer_list<double> values) {
    return Eigen::Map<const VectorXd>(values.begin(), static_cast<Index>(values.size()));
}

MatrixXd RandomMatrix(std::mt19937& rng, Index rows, Index cols, double scale = 1.0) {
    std::uniform_real_distribution<double> distribution(-scale, scale);
    MatrixXd result(rows, cols);
    for (Index row = 0; row < rows; row++) {
        for (Index column = 0; column < cols; column++) {
            result(row, column) = distribution(rng);
        }
    }
    return result;
}

std::pair<MatrixXd, MatrixXd> CanonicalSystem(
    const std::vector<std::size_t>& chain_lengths) {
    std::size_t state_dimension = 0;
    for (std::size_t length : chain_lengths) {
        state_dimension += length;
    }

    MatrixXd Ac = MatrixXd::Zero(
        static_cast<Index>(state_dimension), static_cast<Index>(state_dimension));
    MatrixXd Bc = MatrixXd::Zero(
        static_cast<Index>(state_dimension), static_cast<Index>(chain_lengths.size()));
    Index first = 0;
    for (Index input = 0; input < static_cast<Index>(chain_lengths.size()); input++) {
        const Index length = static_cast<Index>(chain_lengths[static_cast<std::size_t>(input)]);
        for (Index state = first; state + 1 < first + length; state++) {
            Ac(state, state + 1) = 1.0;
        }
        Bc(first + length - 1, input) = 1.0;
        first += length;
    }
    return std::make_pair(Ac, Bc);
}

Index NumericalRankOracle(const MatrixXd& matrix) {
    if (matrix.size() == 0) {
        return 0;
    }

    const JacobiSVD<MatrixXd> svd(matrix);
    const VectorXd singular_values = svd.singularValues();
    if ((singular_values.size() == 0) || (singular_values(0) == 0.0)) {
        return 0;
    }

    const double threshold = static_cast<double>(std::max(matrix.rows(), matrix.cols())) *
        std::numeric_limits<double>::epsilon() * singular_values(0);
    return (singular_values.array() > threshold).count();
}

MatrixXd ControllabilityMatrixOracle(const MatrixXd& A, const MatrixXd& B) {
    MatrixXd controllability(A.rows(), A.rows() * B.cols());
    MatrixXd block = B;
    for (Index power = 0; power < A.rows(); power++) {
        controllability.middleCols(power * B.cols(), B.cols()) = block;
        block = A * block;
    }
    return controllability;
}

std::vector<std::size_t> ControllabilityIndexMultisetOracle(
    const MatrixXd& A,
    const MatrixXd& B) {
    // The rank increase after adding A^(k-1)*B equals the number of
    // controllability indices that are at least k.
    const MatrixXd controllability = ControllabilityMatrixOracle(A, B);
    std::vector<Index> rank_increments;
    Index previous_rank = 0;
    for (Index power = 0; power < A.rows(); power++) {
        const Index columns = (power + 1) * B.cols();
        const Index current_rank = NumericalRankOracle(controllability.leftCols(columns));
        rank_increments.push_back(current_rank - previous_rank);
        previous_rank = current_rank;
    }
    if (previous_rank != A.rows()) {
        return {};
    }

    std::vector<std::size_t> indices(static_cast<std::size_t>(B.cols()), 0);
    for (Index chain = 0; chain < B.cols(); chain++) {
        for (Index increment : rank_increments) {
            if (increment > chain) {
                ++indices[static_cast<std::size_t>(chain)];
            }
        }
    }
    std::sort(indices.begin(), indices.end(), std::greater<std::size_t>());
    return indices;
}

bool HasWellSeparatedRank(const MatrixXd& matrix, Index expected_rank) {
    const JacobiSVD<MatrixXd> svd(matrix);
    const VectorXd singular_values = svd.singularValues();
    if ((expected_rank <= 0) || (singular_values.size() < expected_rank)) {
        return false;
    }
    return singular_values(expected_rank - 1) > 1e-7 * singular_values(0);
}

struct ConstructedSystem {
    MatrixXd A;
    MatrixXd B;
    MatrixXd Ac;
    MatrixXd Bc;
    MatrixXd T;
    MatrixXd Am;
    MatrixXd Bm;
};

ConstructedSystem MakeSystem(
    const std::vector<std::size_t>& chain_lengths,
    const MatrixXd& T,
    const MatrixXd& Am,
    const MatrixXd& Bm) {
    const std::pair<MatrixXd, MatrixXd> canonical = CanonicalSystem(chain_lengths);
    FullPivLU<MatrixXd> state_solve(T);

    ConstructedSystem system;
    system.Ac = canonical.first;
    system.Bc = canonical.second;
    system.T = T;
    system.Am = Am;
    system.Bm = Bm;
    system.A = state_solve.solve((system.Ac + system.Bc * Am) * T);
    system.B = state_solve.solve(system.Bc * Bm);
    return system;
}

MatrixXd WellConditionedStateTransformation(Index n, int trial) {
    MatrixXd T = MatrixXd::Identity(n, n);
    for (Index row = 0; row < n; row++) {
        T(row, row) = 1.0 + 0.1 * static_cast<double>(row + 1);
        for (Index column = 0; column < row; column++) {
            T(row, column) = 0.08 * std::sin(
                static_cast<double>((trial + 1) * (row + 1) * (column + 2)));
        }
    }
    return T;
}

std::vector<std::size_t> ChainLengths(Index n, Index m, int trial) {
    std::vector<std::size_t> lengths(static_cast<std::size_t>(m), 1);
    for (Index remaining = m; remaining < n; remaining++) {
        const Index input = (remaining + trial) % m;
        ++lengths[static_cast<std::size_t>(input)];
    }
    return lengths;
}

void ExpectMatrixNear(
    const MatrixXd& actual,
    const MatrixXd& expected,
    double tolerance = kTolerance) {
    ASSERT_EQ(actual.rows(), expected.rows());
    ASSERT_EQ(actual.cols(), expected.cols());
    const double scale = std::max(1.0, std::max(actual.norm(), expected.norm()));
    EXPECT_LE((actual - expected).norm(), tolerance * scale);
}

void ExpectCanonicalStructure(
    const MatrixXd& Ac,
    const MatrixXd& Bc,
    const std::vector<std::size_t>& indices) {
    const std::pair<MatrixXd, MatrixXd> expected = CanonicalSystem(indices);
    ExpectMatrixNear(Ac, expected.first, 0.0);
    ExpectMatrixNear(Bc, expected.second, 0.0);
}

void ExpectTransformationIdentities(
    const BrunovskyTransformation& form,
    const MatrixXd& A,
    const MatrixXd& B,
    double tolerance = kTolerance) {
    ASSERT_TRUE(form.IsValid());
    const MatrixXd& T = form.TransformationMatrix();
    const MatrixXd& Ac = form.CanonicalStateMatrix();
    const MatrixXd& Bc = form.CanonicalInputMatrix();
    const MatrixXd& Am = form.StateFeedbackMatrix();
    const MatrixXd& Bm = form.InputTransformationMatrix();

    ExpectCanonicalStructure(Ac, Bc, form.ControllabilityIndices());
    ExpectMatrixNear(T * A, (Ac + Bc * Am) * T, tolerance);
    ExpectMatrixNear(T * B, Bc * Bm, tolerance);
    EXPECT_EQ(FullPivLU<MatrixXd>(T).rank(), T.rows());
    EXPECT_EQ(FullPivLU<MatrixXd>(Bm).rank(), Bm.rows());

    std::vector<std::size_t> actual_indices = form.ControllabilityIndices();
    std::sort(actual_indices.begin(), actual_indices.end(), std::greater<std::size_t>());
    EXPECT_EQ(actual_indices, ControllabilityIndexMultisetOracle(A, B));
    EXPECT_EQ(
        std::accumulate(actual_indices.begin(), actual_indices.end(), std::size_t{0}),
        static_cast<std::size_t>(A.rows()));

    const std::size_t nu = form.MaxControllabilityIndex();
    ASSERT_GT(nu, 0u);
    MatrixXd power = MatrixXd::Identity(Ac.rows(), Ac.cols());
    for (std::size_t exponent = 0; exponent < nu; exponent++) {
        power *= Ac;
    }
    ExpectMatrixNear(power, MatrixXd::Zero(Ac.rows(), Ac.cols()), 0.0);
    if (nu > 1) {
        power = MatrixXd::Identity(Ac.rows(), Ac.cols());
        for (std::size_t exponent = 1; exponent < nu; exponent++) {
            power *= Ac;
        }
        EXPECT_GT(power.norm(), 0.5);
    }
}

TEST(BrunovskyTransformationConstruction, RejectsInvalidDimensionsAndNonfiniteData) {
    EXPECT_FALSE(BrunovskyTransformation(MatrixXd(0, 0), MatrixXd(0, 0)).IsValid());
    EXPECT_FALSE(BrunovskyTransformation(MatrixXd::Zero(2, 3), MatrixXd::Zero(2, 1)).IsValid());
    EXPECT_FALSE(BrunovskyTransformation(MatrixXd::Zero(2, 2), MatrixXd::Zero(3, 1)).IsValid());
    EXPECT_FALSE(BrunovskyTransformation(MatrixXd::Zero(2, 2), MatrixXd::Zero(2, 0)).IsValid());
    EXPECT_FALSE(BrunovskyTransformation(MatrixXd::Zero(2, 2), MatrixXd::Zero(2, 3)).IsValid());

    MatrixXd A = MatrixXd::Identity(2, 2);
    MatrixXd B = MatrixXd::Identity(2, 2);
    A(0, 0) = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(BrunovskyTransformation(A, B).IsValid());
    A = MatrixXd::Identity(2, 2);
    B(1, 1) = std::numeric_limits<double>::infinity();
    EXPECT_FALSE(BrunovskyTransformation(A, B).IsValid());
}

TEST(BrunovskyTransformationConstruction, RejectsUncontrollableAndRankDeficientInputs) {
    const MatrixXd uncontrollable_A = (MatrixXd(3, 3) <<
        0, 0, 0,
        0, 1, 0,
        0, 0, 2).finished();
    const MatrixXd uncontrollable_B = Vec({1, 1, 0});
    EXPECT_FALSE(BrunovskyTransformation(uncontrollable_A, uncontrollable_B).IsValid());

    // The pair is controllable through one effective actuator, but B has two
    // dependent columns, so no invertible two-input Bm exists.
    const MatrixXd A = (MatrixXd(2, 2) << 0, 1, 0, 0).finished();
    const MatrixXd B = (MatrixXd(2, 2) << 0, 0, 1, 2).finished();
    EXPECT_FALSE(BrunovskyTransformation(A, B).IsValid());
}

TEST(BrunovskyTransformationConstruction, InvalidObjectHasEmptyTransformationData) {
    const BrunovskyTransformation form(MatrixXd::Identity(2, 2), MatrixXd::Zero(2, 1));
    ASSERT_FALSE(form.IsValid());
    EXPECT_TRUE(form.TransformationMatrix().size() == 0);
    EXPECT_TRUE(form.ControllabilityIndices().empty());
    EXPECT_EQ(form.MaxControllabilityIndex(), 0u);
    EXPECT_EQ(form.CanonicalStateMatrix().size(), 0);
    EXPECT_EQ(form.CanonicalInputMatrix().size(), 0);
    EXPECT_EQ(form.StateFeedbackMatrix().size(), 0);
    EXPECT_EQ(form.InputTransformationMatrix().size(), 0);
}

TEST(BrunovskyTransformationCanonical, CanonicalSystemsRemainCanonical) {
    const std::vector<std::vector<std::size_t>> cases = {
        {1}, {4}, {1, 1, 1}, {3, 1}, {4, 2, 1}, {3, 3, 2, 1, 1}
    };
    for (const std::vector<std::size_t>& chains : cases) {
        SCOPED_TRACE(::testing::PrintToString(chains));
        const std::pair<MatrixXd, MatrixXd> canonical = CanonicalSystem(chains);
        const BrunovskyTransformation form(canonical.first, canonical.second);
        ASSERT_TRUE(form.IsValid());
        EXPECT_EQ(form.ControllabilityIndices(), chains);
        ExpectMatrixNear(
            form.TransformationMatrix(),
            MatrixXd::Identity(canonical.first.rows(), canonical.first.cols()),
            kTolerance);
        ExpectMatrixNear(
            form.StateFeedbackMatrix(),
            MatrixXd::Zero(canonical.second.cols(), canonical.first.cols()),
            kTolerance);
        ExpectMatrixNear(
            form.InputTransformationMatrix(),
            MatrixXd::Identity(canonical.second.cols(), canonical.second.cols()),
            kTolerance);
        ExpectTransformationIdentities(form, canonical.first, canonical.second);
    }
}

TEST(BrunovskyTransformationCanonical, FullyActuatedSystemHasChainsOfLengthOne) {
    const MatrixXd A = (MatrixXd(3, 3) <<
        2, -1, 0,
        1, 3, 2,
        0, 1, -2).finished();
    const MatrixXd B = (MatrixXd(3, 3) <<
        2, 1, 0,
        0, -1, 1,
        1, 0, 2).finished();
    const BrunovskyTransformation form(A, B);
    ASSERT_TRUE(form.IsValid());
    EXPECT_EQ(form.ControllabilityIndices(), (std::vector<std::size_t>{1, 1, 1}));
    EXPECT_EQ(form.MaxControllabilityIndex(), 1u);
    ExpectMatrixNear(form.CanonicalStateMatrix(), MatrixXd::Zero(3, 3), 0.0);
    ExpectMatrixNear(form.CanonicalInputMatrix(), MatrixXd::Identity(3, 3), 0.0);
    ExpectTransformationIdentities(form, A, B);
}

TEST(BrunovskyTransformationCanonical, DenseInputMixingPreservesUnequalChains) {
    const std::vector<std::size_t> chains = {3, 1};
    const MatrixXd T = (MatrixXd(4, 4) <<
        2, 0, 0, 0,
        1, 2, 0, 0,
        0, -1, 3, 0,
        1, 0, 2, 2).finished();
    const MatrixXd Am = (MatrixXd(2, 4) <<
        1, -2, 3, 0,
        -1, 1, 0, 2).finished();
    const MatrixXd Bm = (MatrixXd(2, 2) << 2, 1, 1, -1).finished();
    const ConstructedSystem system = MakeSystem(chains, T, Am, Bm);

    const BrunovskyTransformation form(system.A, system.B);
    ASSERT_TRUE(form.IsValid());
    EXPECT_EQ(form.ControllabilityIndices(), chains);
    ExpectTransformationIdentities(form, system.A, system.B, 2e-9);
}

TEST(BrunovskyTransformationCanonical, AcceptsControllableIllConditionedInputScaling) {
    const MatrixXd A = MatrixXd::Zero(2, 2);
    const MatrixXd B = (MatrixXd(2, 2) << 1, 0, 0, 1e-15).finished();
    const BrunovskyTransformation form(A, B);
    ASSERT_TRUE(form.IsValid());
    EXPECT_TRUE(form.TransformationMatrix().allFinite());
    EXPECT_TRUE(form.InputTransformationMatrix().allFinite());
    ExpectTransformationIdentities(form, A, B, 2e-9);
}

TEST(BrunovskyTransformationCanonical, HandlesModeratelyIllConditionedNontrivialSystem) {
    const std::vector<std::size_t> chains = {3, 2};
    MatrixXd T = MatrixXd::Identity(5, 5);
    T.diagonal() << 1e-3, 1e-2, 1.0, 1e2, 1e3;
    T(2, 0) = 2e-3;
    T(4, 1) = -3.0;
    const MatrixXd Am = (MatrixXd(2, 5) <<
        0.5, -1.0, 2.0, 0.0, 1.0,
        -0.25, 0.0, 1.5, -2.0, 0.75).finished();
    const MatrixXd Bm = (MatrixXd(2, 2) <<
        1e-3, 2e-4,
        -1e2, 1e3).finished();
    const ConstructedSystem system = MakeSystem(chains, T, Am, Bm);

    const BrunovskyTransformation form(system.A, system.B);
    ASSERT_TRUE(form.IsValid());
    EXPECT_EQ(form.ControllabilityIndices(), chains);
    ExpectTransformationIdentities(form, system.A, system.B, 2e-7);
}

TEST(BrunovskyTransformationCanonical, RandomConstructedSystemsUpToTenStates) {
    std::mt19937 rng(78421);
    for (int trial = 0; trial < 80; trial++) {
        SCOPED_TRACE(trial);
        const Index n = 1 + trial % 10;
        const Index m = 1 + (trial * 3) % std::min<Index>(n, 4);
        const std::vector<std::size_t> chains = ChainLengths(n, m, trial);
        const MatrixXd T = WellConditionedStateTransformation(n, trial);
        const MatrixXd Am = RandomMatrix(rng, m, n, 0.75);
        MatrixXd Bm = MatrixXd::Zero(m, m);
        for (Index input = 0; input < m; input++) {
            const double magnitude = 0.5 + 0.2 * static_cast<double>(input + 1);
            Bm(input, input) = ((trial + input) % 2 == 0) ? magnitude : -magnitude;
        }
        const ConstructedSystem system = MakeSystem(chains, T, Am, Bm);

        const BrunovskyTransformation form(system.A, system.B);
        ASSERT_TRUE(form.IsValid());
        EXPECT_EQ(form.ControllabilityIndices(), chains);
        ExpectTransformationIdentities(form, system.A, system.B, 2e-8);

        const VectorXd x = RandomMatrix(rng, n, 1);
        const VectorXd u = RandomMatrix(rng, m, 1);
        const MatrixXd& computed_T = form.TransformationMatrix();
        const MatrixXd& Ac = form.CanonicalStateMatrix();
        const MatrixXd& Bc = form.CanonicalInputMatrix();
        const MatrixXd& computed_Am = form.StateFeedbackMatrix();
        const MatrixXd& computed_Bm = form.InputTransformationMatrix();
        const VectorXd z = computed_T * x;
        const VectorXd r = computed_Am * z + computed_Bm * u;
        ExpectMatrixNear(
            Ac * z + Bc * r,
            computed_T * (system.A * x + system.B * u),
            2e-8);
    }
}

TEST(BrunovskyTransformationCanonical, RandomDenseInputMixingPreservesChainStructure) {
    std::mt19937 rng(31907);
    for (int trial = 0; trial < 40; trial++) {
        SCOPED_TRACE(trial);
        const Index n = 3 + trial % 8;
        const Index m = 2 + trial % std::min<Index>(n - 1, 3);
        std::vector<std::size_t> chains = ChainLengths(n, m, trial);
        std::sort(chains.begin(), chains.end(), std::greater<std::size_t>());

        const MatrixXd T = WellConditionedStateTransformation(n, trial + 100);
        const MatrixXd Am = RandomMatrix(rng, m, n, 0.5);
        MatrixXd Bm = RandomMatrix(rng, m, m, 0.15);
        for (Index input = 0; input < m; input++) {
            Bm(input, input) += 1.0 + 0.2 * static_cast<double>(input);
        }
        ASSERT_EQ(FullPivLU<MatrixXd>(Bm).rank(), m);
        const ConstructedSystem system = MakeSystem(chains, T, Am, Bm);

        const BrunovskyTransformation form(system.A, system.B);
        ASSERT_TRUE(form.IsValid());
        EXPECT_EQ(form.ControllabilityIndices(), chains);
        ExpectTransformationIdentities(form, system.A, system.B, 3e-8);
    }
}

TEST(BrunovskyTransformationCanonical, RandomIndependentControllableSystemsUpToTenStates) {
    // Sample (A, B) directly so this test does not inherit the construction
    // assumptions used by MakeSystem().
    std::mt19937 rng(64013);
    for (Index n = 2; n <= 10; n++) {
        for (int trial = 0; trial < 8; trial++) {
            SCOPED_TRACE(::testing::Message() << "n=" << n << ", trial=" << trial);
            const Index max_inputs = std::min<Index>(n, 4);
            const Index m = 1 + (n + 3 * trial) % max_inputs;
            bool tested = false;

            for (int attempt = 0; (attempt < 200) && !tested; attempt++) {
                const MatrixXd A = RandomMatrix(rng, n, n, 0.6);
                const MatrixXd B = RandomMatrix(rng, n, m, 1.0);
                const MatrixXd controllability = ControllabilityMatrixOracle(A, B);
                if (!HasWellSeparatedRank(B, m) ||
                    !HasWellSeparatedRank(controllability, n)) {
                    continue;
                }

                const BrunovskyTransformation form(A, B);
                ASSERT_TRUE(form.IsValid());
                ExpectTransformationIdentities(form, A, B, 5e-8);
                tested = true;
            }
            ASSERT_TRUE(tested);
        }
    }
}

TEST(BrunovskyTransformationTransformations, StateAndInputSetsPreserveConstraintValues) {
    const std::vector<std::size_t> chains = {2, 1};
    const MatrixXd T = (MatrixXd(3, 3) << 2, 0, 0, 1, 3, 0, -1, 2, 2).finished();
    const MatrixXd Am = (MatrixXd(2, 3) << 1, -2, 0, 0, 1, 3).finished();
    const MatrixXd Bm = (MatrixXd(2, 2) << 2, 1, 1, -1).finished();
    const ConstructedSystem system = MakeSystem(chains, T, Am, Bm);
    const BrunovskyTransformation form(system.A, system.B);
    ASSERT_TRUE(form.IsValid());

    const MatrixXd Gx = (MatrixXd(5, 3) <<
        1, 0, 0,
        -1, 0, 0,
        0, 1, -2,
        2, -1, 1,
        -1, 3, 2).finished();
    const VectorXd Fx = Vec({2, 2, 4, 5, 6});
    const MatrixXd Hex = Vec({1, -2, 3}).transpose();
    const HPolyhedron state_set(Gx, Fx, Hex, Vec({0.5}));
    const HPolyhedron transformed_state = form.TransformStateSet(state_set);
    ASSERT_TRUE(transformed_state.IsValid());
    ASSERT_EQ(transformed_state.NumEqualities(), 1u);

    const MatrixXd Gu = (MatrixXd(4, 2) << 1, 0, -1, 0, 1, 2, -2, 1).finished();
    const VectorXd Fu = Vec({2, 3, 4, 5});
    const MatrixXd Heu = Vec({1, -1}).transpose();
    const HPolyhedron input_set(Gu, Fu, Heu, Vec({0.25}));
    const HPolyhedron transformed_input = form.TransformInputSet(input_set);
    ASSERT_TRUE(transformed_input.IsValid());
    ASSERT_EQ(transformed_input.NumEqualities(), 1u);

    std::mt19937 rng(8182);
    const MatrixXd& computed_Bm = form.InputTransformationMatrix();
    for (int sample = 0; sample < 30; sample++) {
        const VectorXd x = RandomMatrix(rng, 3, 1);
        const VectorXd u = RandomMatrix(rng, 2, 1);
        const VectorXd z = form.TransformationMatrix() * x;
        const VectorXd s = computed_Bm * u;
        ExpectMatrixNear(transformed_state.Ai() * z, state_set.Ai() * x);
        ExpectMatrixNear(transformed_state.Ae() * z, state_set.Ae() * x);
        ExpectMatrixNear(transformed_input.Ai() * s, input_set.Ai() * u);
        ExpectMatrixNear(transformed_input.Ae() * s, input_set.Ae() * u);
    }
}

TEST(BrunovskyTransformationTransformations, JointSetPreservesMixedConstraintsAndEqualities) {
    const std::vector<std::size_t> chains = {3, 1};
    const MatrixXd T = WellConditionedStateTransformation(4, 7);
    const MatrixXd Am = (MatrixXd(2, 4) <<
        1, -2, 0.5, 1,
        -1, 0, 2, 3).finished();
    const MatrixXd Bm = (MatrixXd(2, 2) << 2, 1, -1, 2).finished();
    const ConstructedSystem system = MakeSystem(chains, T, Am, Bm);
    const BrunovskyTransformation form(system.A, system.B);
    ASSERT_TRUE(form.IsValid());

    const MatrixXd Gi = (MatrixXd(8, 6) <<
        1, 0, 0, 0, 0, 0,
        -1, 0, 0, 0, 0, 0,
        0, 1, -2, 0, 1, 0,
        2, 0, 1, -1, -2, 3,
        -1, 3, 0, 2, 1, -1,
        0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, -1, 0,
        1, -1, 2, -2, 3, 1).finished();
    const VectorXd Fi = Vec({2, 2, 5, 8, 7, 3, 3, 10});
    const MatrixXd Ge = (MatrixXd(2, 6) <<
        1, 0, -1, 2, 1, -2,
        0, 2, 1, -1, -3, 1).finished();
    const VectorXd Fe = Vec({0.5, -1.0});
    const HPolyhedron joint_set(Gi, Fi, Ge, Fe);
    const HPolyhedron transformed = form.TransformJointStateInputSet(joint_set);
    ASSERT_TRUE(transformed.IsValid());
    ASSERT_EQ(transformed.Dimension(), 6u);
    ASSERT_EQ(transformed.NumEqualities(), 2u);
    ExpectMatrixNear(transformed.bi(), Fi, 0.0);
    ExpectMatrixNear(transformed.be(), Fe, 0.0);

    std::mt19937 rng(17731);
    const MatrixXd& computed_Am = form.StateFeedbackMatrix();
    const MatrixXd& computed_Bm = form.InputTransformationMatrix();
    for (int sample = 0; sample < 50; sample++) {
        const VectorXd x = RandomMatrix(rng, 4, 1);
        const VectorXd u = RandomMatrix(rng, 2, 1);
        const VectorXd z = form.TransformationMatrix() * x;
        const VectorXd r = computed_Am * z + computed_Bm * u;
        VectorXd original(6);
        original << x, u;
        VectorXd canonical(6);
        canonical << z, r;
        ExpectMatrixNear(transformed.Ai() * canonical, joint_set.Ai() * original, 2e-9);
        ExpectMatrixNear(transformed.Ae() * canonical, joint_set.Ae() * original, 2e-9);
    }
}

TEST(BrunovskyTransformationTransformations, RandomHighDimensionalSetsMatchCoordinateSubstitution) {
    std::mt19937 rng(51077);
    for (int trial = 0; trial < 30; trial++) {
        SCOPED_TRACE(trial);
        const Index n = 5 + trial % 6;
        const Index m = 2 + trial % std::min<Index>(n - 1, 3);
        const std::vector<std::size_t> chains = ChainLengths(n, m, trial + 20);
        const MatrixXd seed_T = WellConditionedStateTransformation(n, trial + 200);
        const MatrixXd seed_Am = RandomMatrix(rng, m, n, 0.5);
        MatrixXd seed_Bm = RandomMatrix(rng, m, m, 0.1);
        for (Index input = 0; input < m; input++) {
            seed_Bm(input, input) += 1.0 + 0.25 * static_cast<double>(input);
        }
        const ConstructedSystem system = MakeSystem(chains, seed_T, seed_Am, seed_Bm);
        const BrunovskyTransformation form(system.A, system.B);
        ASSERT_TRUE(form.IsValid());

        const MatrixXd& T = form.TransformationMatrix();
        const MatrixXd& Am = form.StateFeedbackMatrix();
        const MatrixXd& Bm = form.InputTransformationMatrix();
        const MatrixXd T_inverse = FullPivLU<MatrixXd>(T).solve(MatrixXd::Identity(n, n));
        const MatrixXd Bm_inverse =
            FullPivLU<MatrixXd>(Bm).solve(MatrixXd::Identity(m, m));

        const Index dimension = n + m;
        const Index random_inequalities = 8;
        // Start with a bounding box and add dense mixed state-input cuts.
        MatrixXd Gi = MatrixXd::Zero(2 * dimension + random_inequalities, dimension);
        Gi.topRows(dimension) = MatrixXd::Identity(dimension, dimension);
        Gi.middleRows(dimension, dimension) = -MatrixXd::Identity(dimension, dimension);
        Gi.bottomRows(random_inequalities) =
            RandomMatrix(rng, random_inequalities, dimension, 1.0);
        const VectorXd center = RandomMatrix(rng, dimension, 1, 0.5);
        VectorXd slack = VectorXd::Ones(Gi.rows());
        for (Index row = 0; row < slack.size(); row++) {
            slack(row) += 0.05 * static_cast<double>((row + trial) % 7);
        }
        const VectorXd Fi = Gi * center + slack;
        const MatrixXd Ge = RandomMatrix(rng, 3, dimension, 1.0);
        const VectorXd Fe = Ge * center;
        const HPolyhedron joint_set(Gi, Fi, Ge, Fe);
        ASSERT_TRUE(joint_set.Contains(center, 1e-10));

        const HPolyhedron transformed = form.TransformJointStateInputSet(joint_set);
        ASSERT_TRUE(transformed.IsValid());

        // [x; u] = [T^-1, 0; -Bm^-1*Am, Bm^-1] * [z; r].
        MatrixXd canonical_to_original = MatrixXd::Zero(dimension, dimension);
        canonical_to_original.topLeftCorner(n, n) = T_inverse;
        canonical_to_original.bottomLeftCorner(m, n) = -Bm_inverse * Am;
        canonical_to_original.bottomRightCorner(m, m) = Bm_inverse;
        ExpectMatrixNear(transformed.Ai(), Gi * canonical_to_original, 5e-9);
        ExpectMatrixNear(transformed.Ae(), Ge * canonical_to_original, 5e-9);
        ExpectMatrixNear(transformed.bi(), Fi, 0.0);
        ExpectMatrixNear(transformed.be(), Fe, 0.0);

        const VectorXd x = center.head(n);
        const VectorXd u = center.tail(m);
        VectorXd canonical_point(dimension);
        canonical_point << T * x, Am * T * x + Bm * u;
        EXPECT_TRUE(transformed.Contains(canonical_point, 1e-8));

        const MatrixXd Gx = RandomMatrix(rng, 2 * n + 3, n, 1.0);
        const MatrixXd Hex = RandomMatrix(rng, 2, n, 1.0);
        const HPolyhedron state_set(Gx, Gx * x + VectorXd::Ones(Gx.rows()), Hex, Hex * x);
        const HPolyhedron transformed_state = form.TransformStateSet(state_set);
        ASSERT_TRUE(transformed_state.IsValid());
        ExpectMatrixNear(transformed_state.Ai(), Gx * T_inverse, 5e-9);
        ExpectMatrixNear(transformed_state.Ae(), Hex * T_inverse, 5e-9);

        const MatrixXd Gu = RandomMatrix(rng, 2 * m + 3, m, 1.0);
        const MatrixXd Heu = RandomMatrix(rng, 1, m, 1.0);
        const HPolyhedron input_set(Gu, Gu * u + VectorXd::Ones(Gu.rows()), Heu, Heu * u);
        const HPolyhedron transformed_input = form.TransformInputSet(input_set);
        ASSERT_TRUE(transformed_input.IsValid());
        ExpectMatrixNear(transformed_input.Ai(), Gu * Bm_inverse, 5e-9);
        ExpectMatrixNear(transformed_input.Ae(), Heu * Bm_inverse, 5e-9);
    }
}

TEST(BrunovskyTransformationTransformations, DisturbanceTransformationPreservesDynamics) {
    const std::vector<std::size_t> chains = {4, 2, 1};
    const Index n = 7;
    const Index m = 3;
    const MatrixXd T = WellConditionedStateTransformation(n, 11);
    const MatrixXd Am = (MatrixXd(3, 7) <<
        1, -2, 0, 1, 0, 1, -1,
        0, 1, 2, 0, 1, 0, 1,
        -1, 0, 1, 2, 0, 1, 3).finished();
    const MatrixXd Bm = (MatrixXd(3, 3) << 2, 0, 0, 0, -1, 0, 0, 0, 3).finished();
    const ConstructedSystem system = MakeSystem(chains, T, Am, Bm);
    const BrunovskyTransformation form(system.A, system.B);
    ASSERT_TRUE(form.IsValid());

    const MatrixXd E = (MatrixXd(7, 2) <<
        1, 0,
        0, -1,
        2, 1,
        -1, 2,
        1, 1,
        0, 2,
        -2, 1).finished();
    const MatrixXd Ec = form.TransformDisturbanceMatrix(E);
    ExpectMatrixNear(Ec, form.TransformationMatrix() * E);

    std::mt19937 rng(9021);
    const VectorXd x = RandomMatrix(rng, n, 1);
    const VectorXd u = RandomMatrix(rng, m, 1);
    const VectorXd w = RandomMatrix(rng, 2, 1);
    const MatrixXd& computed_T = form.TransformationMatrix();
    const MatrixXd& Ac = form.CanonicalStateMatrix();
    const MatrixXd& Bc = form.CanonicalInputMatrix();
    const MatrixXd& computed_Am = form.StateFeedbackMatrix();
    const MatrixXd& computed_Bm = form.InputTransformationMatrix();
    const VectorXd z = computed_T * x;
    const VectorXd r = computed_Am * z + computed_Bm * u;
    ExpectMatrixNear(
        Ac * z + Bc * r + Ec * w,
        computed_T * (system.A * x + system.B * u + E * w),
        2e-8);

    const MatrixXd no_disturbance(n, 0);
    const MatrixXd transformed_empty = form.TransformDisturbanceMatrix(no_disturbance);
    EXPECT_EQ(transformed_empty.rows(), n);
    EXPECT_EQ(transformed_empty.cols(), 0);
}

TEST(BrunovskyTransformationTransformations, RejectsInvalidTransformationInputs) {
    const std::pair<MatrixXd, MatrixXd> canonical = CanonicalSystem({2, 1});
    const BrunovskyTransformation form(canonical.first, canonical.second);
    ASSERT_TRUE(form.IsValid());

    EXPECT_FALSE(form.TransformStateSet(HPolyhedron()).IsValid());
    EXPECT_FALSE(form.TransformInputSet(HPolyhedron()).IsValid());
    EXPECT_FALSE(form.TransformJointStateInputSet(HPolyhedron()).IsValid());
    EXPECT_FALSE(form.TransformStateSet(HPolyhedron::FullSpace(2)).IsValid());
    EXPECT_FALSE(form.TransformInputSet(HPolyhedron::FullSpace(1)).IsValid());
    EXPECT_FALSE(form.TransformJointStateInputSet(HPolyhedron::FullSpace(4)).IsValid());
    EXPECT_EQ(form.TransformDisturbanceMatrix(MatrixXd::Zero(2, 1)).size(), 0);
    EXPECT_EQ(form.TransformDisturbanceMatrix(MatrixXd(0, 0)).size(), 0);

    MatrixXd nonfinite = MatrixXd::Zero(3, 1);
    nonfinite(0, 0) = std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(form.TransformDisturbanceMatrix(nonfinite).size(), 0);

    const BrunovskyTransformation invalid(MatrixXd::Identity(2, 2), MatrixXd::Zero(2, 1));
    EXPECT_FALSE(invalid.TransformStateSet(HPolyhedron::FullSpace(2)).IsValid());
    EXPECT_FALSE(invalid.TransformInputSet(HPolyhedron::FullSpace(1)).IsValid());
    EXPECT_FALSE(invalid.TransformJointStateInputSet(HPolyhedron::FullSpace(3)).IsValid());
    EXPECT_EQ(invalid.TransformDisturbanceMatrix(MatrixXd::Zero(2, 1)).size(), 0);
}

TEST(BrunovskyTransformationTransformations, FullSpaceAndEmptySetsPreserveTheirMeaning) {
    const std::pair<MatrixXd, MatrixXd> canonical = CanonicalSystem({2, 1});
    const BrunovskyTransformation form(canonical.first, canonical.second);
    ASSERT_TRUE(form.IsValid());

    EXPECT_TRUE(form.TransformStateSet(HPolyhedron::FullSpace(3)).IsFullSpace());
    EXPECT_TRUE(form.TransformInputSet(HPolyhedron::FullSpace(2)).IsFullSpace());
    EXPECT_TRUE(form.TransformJointStateInputSet(HPolyhedron::FullSpace(5)).IsFullSpace());
    EXPECT_FALSE(form.TransformStateSet(HPolyhedron::EmptySet(3)).IsFeasible());
    EXPECT_FALSE(form.TransformInputSet(HPolyhedron::EmptySet(2)).IsFeasible());
    EXPECT_FALSE(form.TransformJointStateInputSet(HPolyhedron::EmptySet(5)).IsFeasible());
}

}
