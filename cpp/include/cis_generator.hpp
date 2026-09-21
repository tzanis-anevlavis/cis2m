#pragma once

#include <cstddef>
#include <vector>

#include <Eigen/Dense>

#include "brunovskytransformation.hpp"
#include "hpolyhedron.hpp"

namespace cis2m {

/**
 * \brief Parameters of a (tau, lambda)-lasso component.
 *
 *   - `lambda` must be positive and `tau` nonnegative.
 *   - A positive `hierarchy_level` computes all `lambda = 1, ..., hierarchy_level`,
 *     with `tau = hierarchy_level - lambda`;
 *
 * If hierarchy level is specified then any given tau and lambda are ignored.
 */
struct CISOptions {
    std::size_t tau = 0;
    std::size_t lambda = 0;
    std::size_t hierarchy_level = 0;
    bool is_implicit = true;
};

/**
 * \brief One component of the controlled invariant set construction.
 *
 * When `is_implicit` is true, set is in [x; v] coordinates, where `v` has `m * q`
 * entries grouped by input channel and `q = tau + lambda`.
 * Otherwise, set is the projection onto x.
 * The lifted dynamics and input maps always use [x; v] coordinates:
 *
 *     [x+; v+] = lifted_dynamics * [x; v] + lifted_disturbance * w,
 *            u = input_from_state * x + input_from_virtual * v.
 * In nominal mode, lifted_disturbance has zero columns, even if E was
 * supplied when constructing the generator.
 */
struct CISComponent {
    HPolyhedron set;
    Eigen::MatrixXd lifted_dynamics;
    Eigen::MatrixXd lifted_disturbance;
    Eigen::MatrixXd input_from_state;
    Eigen::MatrixXd input_from_virtual;
    std::size_t tau = 0;
    std::size_t lambda = 0;
    bool is_implicit = true;
};

/**
 * \brief Closed-form implicit RCIS generator for a controllable linear system.
 *
 * The system is `x+ = A x + B u + E w`. Omit `E` for a nominal system.
 * Construction rejects invalid, uncontrollable systems or rank-deficient `B`.
 * Compute rejects invalid constraints and numerical failures with exceptions.
 */
class ControlledInvariantSetGenerator {
public:
    ControlledInvariantSetGenerator(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B);
    ControlledInvariantSetGenerator(
        const Eigen::MatrixXd& A,
        const Eigen::MatrixXd& B,
        const Eigen::MatrixXd& E);

    /**
     * \brief Computes one RCIS or all RCISs of a hierarchy.
     *
     * \param[in] safe_set Constraints on [x; u], or on x alone (free u)
     * \param[in] disturbance_set Constraints on w. With E, a valid empty set
     * selects nominal mode; otherwise it must be bounded and match E.cols().
     * Without E, a bounded nonempty set is ignored with a warning.
     * \param[in] options Lasso parameters and output representation
     * \return The resulting RCIS(s) as CISComponent object(s) ordered by increasing lambda
     */
    std::vector<CISComponent> Compute(
        const HPolyhedron& safe_set,
        const HPolyhedron& disturbance_set,
        const CISOptions& options) const;

    /** \brief Matrix-inequality overload with disturbance constraints. */
    std::vector<CISComponent> Compute(
        const Eigen::MatrixXd& Gxu,
        const Eigen::VectorXd& Fxu,
        const Eigen::MatrixXd& Gw,
        const Eigen::VectorXd& Fw,
        const CISOptions& options) const;

    /**
     * \brief Computes one CIS or all CISs of a hierarchy.
     * Disturbance set omitted, the computation is nominal even if E was supplied.
     *
     * \param[in] safe_set Constraints on [x; u], or on x alone (free u)
     * \param[in] options Lasso parameters and output representation
     * \return The resulting CIS(s) as CISComponent object(s) ordered by increasing lambda
     */
    std::vector<CISComponent> Compute(
        const HPolyhedron& safe_set,
        const CISOptions& options) const;

    /** \brief Matrix-inequality overload for a nominal system. */
    std::vector<CISComponent> Compute(
        const Eigen::MatrixXd& Gxu,
        const Eigen::VectorXd& Fxu,
        const CISOptions& options) const;

    /** \brief Returns the Brunovsky transformation object. */
    const BrunovskyTransformation& Transformation() const {
        return transformation_;
    }

private:
    /**
     * \brief Computes one RCIS or all RCISs of a hierarchy.
     * Core method. All the public functions route to this.
     *
     * \param[in] safe_set Constraints on [x; u], or on x alone (free u)
     * \param[in] disturbance_set Valid bounded or empty disturbance set.
     * With E it must match E.cols(); without E it is ignored.
     * \param[in] options Lasso parameters and output representation
     * \return The resulting RCIS(s) as CISComponent object(s) ordered by increasing lambda
     */
    std::vector<CISComponent> Compute_(
        const HPolyhedron& safe_set,
        const HPolyhedron& disturbance_set,
        const CISOptions& options) const;

    /** \brief Initializes component metadata, input maps, and lifted dynamics. */
    CISComponent InitializeComponent_(
        std::size_t tau,
        std::size_t lambda,
        bool is_implicit,
        const Eigen::MatrixXd& H,
        const Eigen::MatrixXd& P,
        bool is_disturbed) const;

    Eigen::MatrixXd A_;
    Eigen::MatrixXd B_;
    Eigen::MatrixXd E_;
    BrunovskyTransformation transformation_;
};

}
