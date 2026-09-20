#include "cis_generator.hpp"

#include <algorithm>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>

using Eigen::Index;
using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace cis2m {
namespace {

// Append free input coordinates to constraints on x alone.
HPolyhedron JointSafeSet(const HPolyhedron& safe_set, Index n, Index m) {
    if (safe_set.Dimension() == static_cast<std::size_t>(n + m)) {
        return safe_set;
    }

    // The validated state-only set leaves the input coordinates free.
    MatrixXd Ai = MatrixXd::Zero(safe_set.Ai().rows(), n + m);
    MatrixXd Ae = MatrixXd::Zero(safe_set.Ae().rows(), n + m);
    Ai.leftCols(n) = safe_set.Ai();
    Ae.leftCols(n) = safe_set.Ae();
    return HPolyhedron(Ai, safe_set.bi(), Ae, safe_set.be());
}

// Construct the projection matrix `H` and the eventually periodic matrix `P`
// for the sequence of virtual inputs `v = [v_1; ...; v_mq]`.
std::pair<MatrixXd, MatrixXd> LassoMatrices(Index m, Index q, Index tau) {
    MatrixXd H = MatrixXd::Zero(m, m * q);
    MatrixXd P = MatrixXd::Zero(m * q, m * q);
    for (Index channel = 0; channel < m; channel++) {
        const Index first = channel * q;
        H(channel, first) = 1.0;
        for (Index sample = 0; sample < q - 1; sample++) {
            P(first + sample, first + sample + 1) = 1.0;
        }
        P(first + q - 1, first + tau) = 1.0;
    }
    return {H, P};
}

// Substitute `r = H * v` into all joint state-input constraints.
HPolyhedron LiftJointSet(const HPolyhedron& joint_set, Index n, const MatrixXd& H) {
    // For each joint normal [Gz, Gr], substitution gives [Gz, Gr * H]
    // acting on [z; v]. Apply the same map to inequality and equality rows.
    const Index mq = H.cols();
    MatrixXd Ai(joint_set.Ai().rows(), n + mq);
    MatrixXd Ae(joint_set.Ae().rows(), n + mq);
    Ai.leftCols(n) = joint_set.Ai().leftCols(n);
    Ai.rightCols(mq) = joint_set.Ai().rightCols(H.rows()) * H;
    Ae.leftCols(n) = joint_set.Ae().leftCols(n);
    Ae.rightCols(mq) = joint_set.Ae().rightCols(H.rows()) * H;
    return HPolyhedron(Ai, joint_set.bi(), Ae, joint_set.be());
}

// Reject dimensions that cannot be represented by Eigen or by the lifted stack.
void ValidateLiftDimensions(Index n, Index m, std::size_t nu, std::size_t q, std::size_t inequalities, std::size_t equalities) {
    const auto limit = static_cast<std::size_t>(std::numeric_limits<Index>::max());
    if ((q == 0) || (q > limit) || (nu > limit - q)) {
        throw std::invalid_argument("Lasso length exceeds supported dimensions.");
    }
    if ((static_cast<std::size_t>(m) > limit / q) || (static_cast<std::size_t>(n) > limit - static_cast<std::size_t>(m) * q)) {
        throw std::invalid_argument("Lifted state dimension exceeds supported dimensions.");
    }
    const std::size_t count = nu + q;
    if ((inequalities > limit / count) || (equalities > limit / count)) {
        throw std::invalid_argument("Lifted constraint count exceeds supported dimensions.");
    }
}

}

// Constructors.
ControlledInvariantSetGenerator::ControlledInvariantSetGenerator(
    const MatrixXd& A,
    const MatrixXd& B)
    : ControlledInvariantSetGenerator(A, B, MatrixXd(A.rows(), 0)) {}

ControlledInvariantSetGenerator::ControlledInvariantSetGenerator(
    const MatrixXd& A,
    const MatrixXd& B,
    const MatrixXd& E)
    : A_(A), B_(B), E_(E), transformation_(A, B) {
    if (!transformation_.IsValid()) {
        throw std::invalid_argument("(A, B) must be finite, controllable, and have full-column-rank B.");
    }
    if ((E.rows() != A.rows()) || !E.allFinite()) {
        throw std::invalid_argument("E must have one finite row per state.");
    }
}

// Public functions.
std::vector<CISComponent> ControlledInvariantSetGenerator::Compute(
    const HPolyhedron& safe_set,
    const HPolyhedron& disturbance_set,
    const CISOptions& options) const {
    return Compute_(safe_set, disturbance_set, options);
}

std::vector<CISComponent> ControlledInvariantSetGenerator::Compute(
    const MatrixXd& Gxu,
    const VectorXd& Fxu,
    const MatrixXd& Gw,
    const VectorXd& Fw,
    const CISOptions& options) const {
    return Compute(HPolyhedron(Gxu, Fxu), HPolyhedron(Gw, Fw), options);
}

std::vector<CISComponent> ControlledInvariantSetGenerator::Compute(
    const HPolyhedron& safe_set,
    const CISOptions& options) const {
    const std::size_t disturbance_dimension = static_cast<std::size_t>(std::max<Index>(E_.cols(), 1));
    return Compute_(safe_set, HPolyhedron::EmptySet(disturbance_dimension), options);
}

std::vector<CISComponent> ControlledInvariantSetGenerator::Compute(
    const MatrixXd& Gxu,
    const VectorXd& Fxu,
    const CISOptions& options) const {
    return Compute(HPolyhedron(Gxu, Fxu), options);
}

// Private functions.
std::vector<CISComponent> ControlledInvariantSetGenerator::Compute_(
    const HPolyhedron& safe_set,
    const HPolyhedron& disturbance_set,
    const CISOptions& options) const {

    // Validate compatibility of provided `disturbance_set` and system.
    // Notice that the public functions that do not accept a disturbance set, still
    // provide a valid empty one through their implementation.
    if (!disturbance_set.IsValid()) {
        throw std::invalid_argument(
            "ControlledInvariantSetGenerator::Compute: Disturbance set is invalid.");
    }
    const bool is_accepting_disturbance = E_.cols() != 0;
    const bool is_empty_disturbance_set = !disturbance_set.IsFeasible();
    if (is_accepting_disturbance && (disturbance_set.Dimension() != static_cast<std::size_t>(E_.cols()))) {
        throw std::invalid_argument(
            "ControlledInvariantSetGenerator::Compute: Disturbance set dimension must match E.cols().");
    }
    if (!is_empty_disturbance_set && !disturbance_set.IsBounded()) {
        throw std::invalid_argument(
            "ControlledInvariantSetGenerator::Compute: Nonempty disturbance set must be bounded.");
    }
    if (!is_accepting_disturbance && !is_empty_disturbance_set) {
        std::clog
            << "ControlledInvariantSetGenerator::Compute: ignoring disturbance set because E has no columns.\n";
    }
    const bool is_disturbed = is_accepting_disturbance && !is_empty_disturbance_set;

    // Validate the `CISOptions`.
    const bool is_hierarchy_specified = options.hierarchy_level != 0;
    if (!is_hierarchy_specified && (options.lambda == 0)) {
        throw std::invalid_argument("Specify a positive lambda or hierarchy_level.");
    }
    const std::size_t limit = static_cast<std::size_t>(std::numeric_limits<Index>::max());
    if (!is_hierarchy_specified && ((options.lambda > limit) || (options.tau > limit - options.lambda))) {
        throw std::invalid_argument("tau + lambda exceeds supported dimensions.");
    }
    const std::size_t q = is_hierarchy_specified ? options.hierarchy_level : options.tau + options.lambda;

    // Validate the safe set before appending any free input coordinates.
    const Index n = A_.rows();
    const Index m = B_.cols();
    const std::size_t nu = transformation_.MaxControllabilityIndex();
    if (!safe_set.IsValid()) {
        throw std::invalid_argument("Safe set is invalid.");
    }
    ValidateLiftDimensions(n, m, nu, q, safe_set.NumInequalities(), safe_set.NumEqualities());
    if ((safe_set.Dimension() != static_cast<std::size_t>(n)) &&
        (safe_set.Dimension() != static_cast<std::size_t>(n + m))) {
        throw std::invalid_argument("Safe set must constrain x or [x; u].");
    }
    const bool is_safe_set_full_space = safe_set.IsFullSpace();
    if (is_safe_set_full_space) {
        std::clog
            << "ControlledInvariantSetGenerator::Compute: safe set is full space; returning full-space components.\n";
    }

    // Transform the joint safe set to Brunovsky coordinates.
    const HPolyhedron Sxu = JointSafeSet(safe_set, n, m);
    const HPolyhedron Sc = transformation_.TransformJointStateInputSet(Sxu);
    if (!Sc.IsValid()) {
        throw std::runtime_error("Failed to transform the safe set to Brunovsky coordinates.");
    }

    // Compute the successive safe sets.
    const MatrixXd& Ac = transformation_.CanonicalStateMatrix();
    const MatrixXd& Bc = transformation_.CanonicalInputMatrix();
    std::vector<HPolyhedron> shrunk_sets;
    shrunk_sets.reserve(is_disturbed ? nu + 1 : 1);
    shrunk_sets.push_back(Sc); // S_0 = Sc.
    if (is_disturbed && !is_safe_set_full_space) {
        const MatrixXd Ec = transformation_.TransformDisturbanceMatrix(E_);
        MatrixXd A_power = MatrixXd::Identity(n, n);
        for (std::size_t t = 1; t <= nu; t++) {
            MatrixXd disturbance_map = MatrixXd::Zero(n + m, E_.cols());
            disturbance_map.topRows(n) = A_power * Ec;
            HPolyhedron next_shrunk_set = shrunk_sets.back().PontryaginDifferenceOfLinearImage(disturbance_set, disturbance_map);
            if (!next_shrunk_set.IsValid()) {
                throw std::runtime_error("Failed to shrink the safe set by accumulated disturbance.");
            }
            shrunk_sets.push_back(std::move(next_shrunk_set));
            A_power = A_power * Ac;
        }
    }

    // Compute the lifted RCIS.
    const Index mq = m * static_cast<Index>(q);
    const Index lifted_dimension = n + mq;
    const std::size_t component_count = is_hierarchy_specified ? q : 1;
    std::vector<CISComponent> components;
    components.reserve(component_count);
    for (std::size_t i = 0; i < component_count; i++) {
        // Construct the projection matrix `H` and the eventually periodic matrix `P`.
        const std::size_t lambda = is_hierarchy_specified ? i + 1 : options.lambda;
        const std::size_t tau = q - lambda;
        const auto lasso = LassoMatrices(m, static_cast<Index>(q), static_cast<Index>(tau));
        const MatrixXd& H = lasso.first;
        const MatrixXd& P = lasso.second;

        // Construct the lifted companion dynamical system.
        MatrixXd A_lifted = MatrixXd::Zero(lifted_dimension, lifted_dimension);
        A_lifted.topLeftCorner(n, n) = Ac;
        A_lifted.topRightCorner(n, mq) = Bc * H;
        A_lifted.bottomRightCorner(mq, mq) = P;
        if (!A_lifted.allFinite()) {
            throw std::runtime_error("Lifted companion dynamics are nonfinite.");
        }

        CISComponent result = InitializeComponent_(tau, lambda, options.is_implicit, H, P, is_disturbed);

        if (is_safe_set_full_space) {
            // If the safe set is the full space, then so is the RCIS.
            result.set = HPolyhedron::FullSpace(static_cast<std::size_t>(options.is_implicit ? lifted_dimension : n));
            components.push_back(std::move(result));
            continue;
        }

        const Index N = static_cast<Index>(nu + q);
        // Time t uses S_t until nilpotency, then S_nu. Count rows first:
        // robustifying an equality can replace S_t with the canonical empty
        // set, whose row counts differ from those of S_0.
        Index ni = 0;
        Index ne = 0;
        for (Index t = 0; t < N; t++) {
            const auto idx_t =
                is_disturbed ? std::min<std::size_t>(static_cast<std::size_t>(t), nu) : 0;
            const HPolyhedron& Sc_t = shrunk_sets[idx_t];
            if ((Sc_t.Ai().rows() > std::numeric_limits<Index>::max() - ni) ||
                (Sc_t.Ae().rows() > std::numeric_limits<Index>::max() - ne)) {
                throw std::invalid_argument("Lifted constraint count exceeds supported dimensions.");
            }
            ni += Sc_t.Ai().rows();
            ne += Sc_t.Ae().rows();
        }

        // Construct the RCIS in the canonical lifted space.
        MatrixXd Ai(ni, lifted_dimension);
        MatrixXd Ae(ne, lifted_dimension);
        VectorXd bi(ni);
        VectorXd be(ne);
        Index inequality_row = 0;
        Index equality_row = 0;
        MatrixXd power = MatrixXd::Identity(lifted_dimension, lifted_dimension);
        for (Index t = 0; t < N; t++) {
            const auto idx_t =
                is_disturbed ? std::min<std::size_t>(static_cast<std::size_t>(t), nu) : 0;
            const HPolyhedron& Sc_t = shrunk_sets[idx_t];
            const HPolyhedron Sc_lifted_t = LiftJointSet(Sc_t, n, H);
            const Index irows = Sc_lifted_t.Ai().rows();
            const Index erows = Sc_lifted_t.Ae().rows();
            Ai.middleRows(inequality_row, irows) = Sc_lifted_t.Ai() * power;
            Ae.middleRows(equality_row, erows) = Sc_lifted_t.Ae() * power;
            bi.segment(inequality_row, irows) = Sc_lifted_t.bi();
            be.segment(equality_row, erows) = Sc_lifted_t.be();
            inequality_row += irows;
            equality_row += erows;
            power = power * A_lifted;
        }
        const HPolyhedron rcis_canonical_lifted(Ai, bi, Ae, be);

        // Transform the RCIS from the canonical lifted space to the original state coordinates lifted space.
        MatrixXd T_lift = MatrixXd::Identity(lifted_dimension, lifted_dimension);
        T_lift.topLeftCorner(n, n) = transformation_.TransformationMatrix();
        HPolyhedron rcis_lifted = rcis_canonical_lifted.Preimage(T_lift);
        if (!rcis_lifted.IsValid()) {
            throw std::runtime_error("Failed to express the lifted set in original-state coordinates.");
        }

        // If necessary, project to the state coordinate space.
        if (options.is_implicit) {
            result.set = std::move(rcis_lifted);
        } else {
            std::vector<std::size_t> state_coordinates(static_cast<std::size_t>(n));
            for (Index state = 0; state < n; state++) {
                state_coordinates[static_cast<std::size_t>(state)] = static_cast<std::size_t>(state);
            }
            result.set = rcis_lifted.ProjectOnto(state_coordinates);
            if (!result.set.IsValid()) {
                throw std::runtime_error("Projection of the implicit RCIS failed.");
            }
        }
        components.push_back(std::move(result));
    }
    return components;
}

CISComponent ControlledInvariantSetGenerator::InitializeComponent_(
    std::size_t tau,
    std::size_t lambda,
    bool is_implicit,
    const MatrixXd& H,
    const MatrixXd& P,
    bool is_disturbed) const {

    const Index n = A_.rows();
    const Index mq = H.cols();
    const Index lifted_dimension = n + mq;

    // The physical input satisfies r = Am*T*x + Bm*u = H*v.
    CISComponent result;
    result.tau = tau;
    result.lambda = lambda;
    result.is_implicit = is_implicit;
    const auto Bm_qr = transformation_.InputTransformationMatrix().colPivHouseholderQr();
    result.input_from_state = -Bm_qr.solve(
        transformation_.StateFeedbackMatrix() * transformation_.TransformationMatrix());
    result.input_from_virtual = Bm_qr.solve(H);
    result.lifted_dynamics = MatrixXd::Zero(lifted_dimension, lifted_dimension);
    result.lifted_dynamics.topLeftCorner(n, n) = A_ + B_ * result.input_from_state;
    result.lifted_dynamics.topRightCorner(n, mq) = B_ * result.input_from_virtual;
    result.lifted_dynamics.bottomRightCorner(mq, mq) = P;
    result.lifted_disturbance = MatrixXd::Zero(lifted_dimension, is_disturbed ? E_.cols() : 0);
    if (is_disturbed) {
        result.lifted_disturbance.topRows(n) = E_;
    }
    if (!result.input_from_state.allFinite() ||
        !result.input_from_virtual.allFinite() ||
        !result.lifted_dynamics.allFinite()) {
        throw std::runtime_error("Lifted coordinate transformation is nonfinite.");
    }
    return result;
}

}
