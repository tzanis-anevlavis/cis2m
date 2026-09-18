#include "brunovskytransformation.hpp"

#include <algorithm>
#include <limits>

using Eigen::ColPivHouseholderQR;
using Eigen::Index;
using Eigen::MatrixXd;
using Eigen::RowVectorXd;

namespace cis2m {
namespace {

// Use a scale-relative rank threshold close to Eigen and MATLAB defaults.
double RankThreshold(const MatrixXd& matrix) {
    const Index dimension = std::max(matrix.rows(), matrix.cols());
    return static_cast<double>(dimension) * std::numeric_limits<double>::epsilon();
}

// Compute numerical rank using the same threshold throughout the construction.
Index NumericalRank(const MatrixXd& matrix) {
    ColPivHouseholderQR<MatrixXd> qr(matrix);
    qr.setThreshold(RankThreshold(matrix));
    return qr.rank();
}

// Solve X*denominator = numerator without explicitly forming denominator^-1.
bool SolveRight(const MatrixXd& numerator, const MatrixXd& denominator, MatrixXd& result) {
    if ((denominator.rows() == 0) ||
        (denominator.rows() != denominator.cols()) ||
        (numerator.cols() != denominator.rows()) ||
        !numerator.allFinite() ||
        !denominator.allFinite()) {
        return false;
    }

    ColPivHouseholderQR<MatrixXd> qr(denominator.transpose());
    qr.setThreshold(RankThreshold(denominator));
    if (qr.rank() != denominator.rows()) {
        return false;
    }

    result = qr.solve(numerator.transpose()).transpose();
    if (!result.allFinite()) {
        return false;
    }

    const double scale = std::max(1.0, numerator.norm());
    return ((result * denominator) - numerator).norm() <= 1e-9 * scale;
}

// Compare matrices with a norm-relative tolerance.
bool MatricesApprox(const MatrixXd& lhs, const MatrixXd& rhs, double tolerance = 1e-8) {
    if ((lhs.rows() != rhs.rows()) || (lhs.cols() != rhs.cols())) {
        return false;
    }
    const double scale = std::max(1.0, std::max(lhs.norm(), rhs.norm()));
    return (lhs - rhs).norm() <= tolerance * scale;
}

}

BrunovskyTransformation::BrunovskyTransformation(const MatrixXd& A, const MatrixXd& B) {
    if ((A.rows() == 0) ||
        (A.rows() != A.cols()) ||
        (B.rows() != A.rows()) ||
        (B.cols() == 0) ||
        (B.cols() > A.rows()) ||
        !A.allFinite() ||
        !B.allFinite()) {
        return;
    }

    const Index n = A.rows();
    const Index m = B.cols();
    if (NumericalRank(B) != m) {
        return;
    }

    // Compute the controllability matrix.
    const MatrixXd C = ControllabilityMatrix(A, B);
    if (NumericalRank(C) != n) {
        return;
    }
    // Compute the controllability indices: populates `controllability_indices_`.
    if (!ComputeControllabilityIndices(C, m)) {
        return;
    }

    // Get the corresponding basis.
    const MatrixXd basis = ControllabilityBasis(A, B, controllability_indices_);
    if (NumericalRank(basis) != n) {
        controllability_indices_.clear();
        return;
    }

    // Compute the similarity transformation `T_`.
    MatrixXd basis_inverse;
    if (!SolveRight(MatrixXd::Identity(n, n), basis, basis_inverse)) {
        controllability_indices_.clear();
        return;
    }
    T_.setZero(n, n);
    Index first = 0;
    for (Index input = 0; input < m; input++) {
        const Index length = static_cast<Index>(controllability_indices_[static_cast<std::size_t>(input)]);
        const Index last = first + length - 1;
        RowVectorXd row = basis_inverse.row(last);
        for (Index state = first; state <= last; state++) {
            T_.row(state) = row;
            row = row * A;
        }
        first += length;
    }
    if ((NumericalRank(T_) != n) || !T_.allFinite()) {
        T_.resize(0, 0);
        controllability_indices_.clear();
        return;
    }

    // Build the feedback and input transformations (`Am_`, `Bm_`).
    MatrixXd transformed_A;
    if (!SolveRight(T_ * A, T_, transformed_A)) {
        T_.resize(0, 0);
        controllability_indices_.clear();
        return;
    }
    const MatrixXd transformed_B = T_ * B;

    Am_.resize(m, n);
    Bm_.resize(m, m);
    first = 0;
    for (Index input = 0; input < m; input++) {
        const Index length = static_cast<Index>(controllability_indices_[static_cast<std::size_t>(input)]);
        const Index last = first + length - 1;
        Am_.row(input) = transformed_A.row(last);
        Bm_.row(input) = transformed_B.row(last);
        first += length;
    }
    if ((NumericalRank(Bm_) != m) || !Am_.allFinite() || !Bm_.allFinite()) {
        T_.resize(0, 0);
        Am_.resize(0, 0);
        Bm_.resize(0, 0);
        controllability_indices_.clear();
        return;
    }

    // Build the canonical system matrices (`Ac_`, `Bc_`).
    BuildCanonicalSystem();

    // Check that the resulting transformations are sound.
    const MatrixXd reconstructed_A = (Ac_ + Bc_ * Am_) * T_;
    const MatrixXd reconstructed_B = Bc_ * Bm_;
    if (!MatricesApprox(T_ * A, reconstructed_A) || !MatricesApprox(T_ * B, reconstructed_B)) {
        Ac_.resize(0, 0);
        Bc_.resize(0, 0);
        T_.resize(0, 0);
        Am_.resize(0, 0);
        Bm_.resize(0, 0);
        controllability_indices_.clear();
        return;
    }
    // Extract max controllability index.
    max_controllability_index_ = *std::max_element(controllability_indices_.begin(), controllability_indices_.end());

    // If we reach here, everything computed successfully. Set validity to `true`.
    valid_ = true;
}

HPolyhedron BrunovskyTransformation::TransformStateSet(const HPolyhedron& state_set) const {
    if (!valid_ || !state_set.IsValid() ||
        (state_set.Dimension() != static_cast<std::size_t>(T_.cols()))) {
        return {};
    }

    MatrixXd Ai;
    MatrixXd Ae;
    if (!SolveRight(state_set.Ai(), T_, Ai) || !SolveRight(state_set.Ae(), T_, Ae)) {
        return {};
    }
    return HPolyhedron(Ai, state_set.bi(), Ae, state_set.be());
}

HPolyhedron BrunovskyTransformation::TransformInputSet(const HPolyhedron& input_set) const {
    if (!valid_ || !input_set.IsValid() ||
        (input_set.Dimension() != static_cast<std::size_t>(Bm_.cols()))) {
        return {};
    }

    MatrixXd Ai;
    MatrixXd Ae;
    if (!SolveRight(input_set.Ai(), Bm_, Ai) || !SolveRight(input_set.Ae(), Bm_, Ae)) {
        return {};
    }
    return HPolyhedron(Ai, input_set.bi(), Ae, input_set.be());
}

HPolyhedron BrunovskyTransformation::TransformJointStateInputSet(const HPolyhedron& joint_set) const {
    if (!valid_ || !joint_set.IsValid() ||
        (joint_set.Dimension() != static_cast<std::size_t>(T_.cols() + Bm_.cols()))) {
        return {};
    }

    const Index n = T_.cols();
    const Index m = Bm_.cols();
    const auto transform_block = [this, n, m](const MatrixXd& G, MatrixXd& transformed) {
        const MatrixXd Gx = G.leftCols(n);
        const MatrixXd Gu = G.rightCols(m);
        MatrixXd Gx_transformed;
        MatrixXd Gr;
        if (!SolveRight(Gx, T_, Gx_transformed) || !SolveRight(Gu, Bm_, Gr)) {
            return false;
        }
        transformed.resize(G.rows(), n + m);
        transformed.leftCols(n) = Gx_transformed - Gr * Am_;
        transformed.rightCols(m) = Gr;
        return transformed.allFinite();
    };

    MatrixXd Ai;
    MatrixXd Ae;
    if (!transform_block(joint_set.Ai(), Ai) || !transform_block(joint_set.Ae(), Ae)) {
        return {};
    }
    return HPolyhedron(Ai, joint_set.bi(), Ae, joint_set.be());
}

MatrixXd BrunovskyTransformation::TransformDisturbanceMatrix(const MatrixXd& E) const {
    if (!valid_ || !E.allFinite()) {
        return {};
    }
    if (E.size() == 0) {
        if (E.rows() != T_.cols()) {
            return {};
        }
        return T_ * E;
    }
    if (E.rows() != T_.cols()) {
        return {};
    }
    return T_ * E;
}

// Private functions.

MatrixXd BrunovskyTransformation::ControllabilityMatrix(const MatrixXd& A, const MatrixXd& B) const {
    // Construct [B, A B, ..., A^(n-1) B].
    const Index n = A.rows();
    const Index m = B.cols();
    MatrixXd C(n, n * m);
    MatrixXd block = B;
    for (Index power = 0; power < n; power++) {
        C.middleCols(power * m, m) = block;
        block = A * block;
    }
    return C;
}

bool BrunovskyTransformation::ComputeControllabilityIndices(const MatrixXd& C, Index m) {
    // Count the independent Krylov columns assigned to each input channel.
    const Index n = C.rows();
    MatrixXd selected = MatrixXd::Zero(n, n);
    Index rank = 0;
    controllability_indices_.assign(static_cast<std::size_t>(m), 0);

    for (Index column = 0; (column < C.cols()) && (rank < n); column++) {
        selected.col(rank) = C.col(column);
        const Index candidate_rank = NumericalRank(selected.leftCols(rank + 1));
        if (candidate_rank > rank) {
            const Index input = column % m;
            ++(controllability_indices_)[static_cast<std::size_t>(input)];
            ++rank;
        }
    }

    if (rank != n) {
        return false;
    }
    return std::all_of(
        controllability_indices_.begin(),
        controllability_indices_.end(),
        [](std::size_t index) {
            return index > 0;
        });
}

MatrixXd BrunovskyTransformation::ControllabilityBasis(
    const MatrixXd& A,
    const MatrixXd& B,
    const std::vector<std::size_t>& indices) const {
    // Reorder the selected Krylov chains by input channel.
    const Index n = A.rows();
    MatrixXd basis(n, n);
    Index column = 0;
    for (Index input = 0; input < B.cols(); input++) {
        MatrixXd power = MatrixXd::Identity(n, n);
        for (std::size_t k = 0; k < indices[static_cast<std::size_t>(input)]; k++) {
            basis.col(column++) = power * B.col(input);
            power = power * A;
        }
    }
    return basis;
}

void BrunovskyTransformation::BuildCanonicalSystem() {
    // Construct exact nilpotent Brunovsky chains from their lengths.
    std::size_t state_dimension = 0;
    for (std::size_t index : controllability_indices_) {
        state_dimension += index;
    }

    Ac_.setZero(static_cast<Index>(state_dimension), static_cast<Index>(state_dimension));
    Bc_.setZero(static_cast<Index>(state_dimension), static_cast<Index>(controllability_indices_.size()));
    Index first = 0;
    for (Index input = 0; input < static_cast<Index>(controllability_indices_.size()); input++) {
        const Index length = static_cast<Index>(controllability_indices_[static_cast<std::size_t>(input)]);
        for (Index state = first; state + 1 < first + length; state++) {
            Ac_(state, state + 1) = 1.0;
        }
        Bc_(first + length - 1, input) = 1.0;
        first += length;
    }
}

}
