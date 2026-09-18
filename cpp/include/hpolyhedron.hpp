#pragma once

#include <Eigen/Dense>
#include <cstddef>
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace cis2m {

/**
 * @brief H-representation of a polyhedron {x | Ai*x <= bi, Ae*x = be},
 * in a positive ambient dimension.
 * Validity checks finite data and dimensions, not feasibility. The default
 * object is invalid; use EmptySet(n) for an empty set and FullSpace(n) for R^n.
 * Failed set operations return an invalid object. Numerical LP failures in
 * IsFeasible throw; ComputeSupport reports NaN for infeasibility or failure.
 */
class HPolyhedron {
public:
    /**
     * \brief Constructor of polyhedron, default invalid object
     */
    HPolyhedron();

    /**
     * \brief Constructor of polyhedron providing inequality constraints Ax <= b
     * \param[in] A Matrix A
     * \param[in] b Matrix b
     */
    HPolyhedron(const MatrixXd& A, const VectorXd& b);

    /**
     * \brief Constructor of a polyhedron providing inequality and equality
     * constraints:
     *  Ai x <= bi
     *  Ae x = be
     * \param[in] Ai Inequality constraints matrix Ai
     * \param[in] bi Inequality constraints vector bi
     * \param[in] Ae Equality constraints matrix Ae
     * \param[in] be Equality constraints vector be
     */
    HPolyhedron(
        const MatrixXd& Ai,
        const VectorXd& bi,
        const MatrixXd& Ae,
        const VectorXd& be);

    /* Copy constructor */
    HPolyhedron(const HPolyhedron&) = default;
    /* Copy assignment operator */
    HPolyhedron& operator=(const HPolyhedron&) = default;

    /* Move constructor */
    HPolyhedron(HPolyhedron&& other) noexcept;
    /* Move assignment operator */
    HPolyhedron& operator=(HPolyhedron&& other) noexcept;

    /**
     * \brief Constructor for full-space polyhedron R^n.
     * The resulting object is valid and has no inequality/equality constraints.
     * \param[in] n The dimension of the polyhedron
     */
    static HPolyhedron FullSpace(size_t n);

    /**
     * \brief Constructor for an empty polyhedron in R^n.
     * The resulting object is valid and has the contradictory constraint 0 <= -1.
     * \param[in] n The dimension of the polyhedron
     */
    static HPolyhedron EmptySet(size_t n);

    /**
     * \brief Prints information about the HPolyhedron object
     */
    void Show() const;

    /**
     * \brief Normalizes the inequality and equality constraints.
     * Divides rows with |b| > tol by |b|, giving RHS +1 or -1.
     * Zero/small RHS rows remain unchanged; inequality signs are preserved.
     * \param[in] tol Numerical tolerance for the normalization operation
     */
    void NormalizeConstraints(double tol = 1e-6);

    /**
     * \brief Sequential redundancy constraint removal, preserving equalities.
     * tol is an absolute allowance after normalizing each constraint normal.
     *
     * \param[in] tol Absolute allowance after normalizing each constraint normal
     */
    void RemoveRedundantConstraints(double tol = 1e-6);

    /**
     * \brief Checks if a point belongs to the Polyhedron
     * \param[in] point Point for which to check containment
     * \param[in] tol Numerical absolute constraint-residual tolerance
     * \return Boolean indicating containment (true) or not (false)
     */
    bool Contains(const VectorXd& point, double tol = 1e-6) const;

    /**
     * \brief Check whether the polyhedron constraints are feasible
     * \return Boolean indicating whether the constraints are feasible (true) or not (false)
     */
    bool IsFeasible() const;

    /**
     * \brief Check whether the polyhedron is bounded.
     * The empty set is bounded. LP failures in feasibility checks propagate.
     * \return Boolean indicating whether is is bounded (true) or not (false)
     */
    bool IsBounded() const;

    /**
     * \brief Checks if this polyhedron is subset of another
     * \param[in] other A polyhedron object
     * \param[in] tol Numerical absolute constraint-residual tolerance
     * \return Boolean indicating containment (true) or not (false)
     */
    bool IsSubsetOf(const HPolyhedron& other, double tol = 1e-6) const;

    /**
     * \brief Compute the support of the polyhedron along the directions provided in A_other.
     * Returns +infinity for unbounded support and NaN for invalid input, infeasibility, or solver failure.
     * h_P(0) = 0 if P is nonempty. The constraint model is built once for the whole batch.
     * \param[in] directions Matrix containing the directions as rows
     * \return Vector with the support along the required directions
     */
    VectorXd ComputeSupport(const MatrixXd& directions) const;

    /**
     * \brief Computes the intersection between this polyhedron and another
     * \param[in] other A polyhedron object
     * \return A polyhedron object that is the intersection of the two polyhedra
     */
    HPolyhedron Intersection(const HPolyhedron& other) const;

    /**
     * \brief Computes the inverse image of this polyhedron through a linear map T.
     * Preimage supports rectangular/singular maps without projection.
     * \param[in] T The linear map
     * \return A polyhedron object {x | T*x belongs to this polyhedron}
     */
    HPolyhedron Preimage(const MatrixXd& T) const;

    /**
     * \brief Projection of the polyhedron on the space with a component removed.
     * It removes one zero-based coordinate using equality substitution or Fourier-Motzkin elimination.
     * Nonzero coefficients are never discarded solely for being small.
     * Zero-dimensional output and more than 100000 generated rows are rejected.
     * \param[in] var_index Index of the element to be removed
     * \param[in] base_tol Numerical tolerance for redundancy removal
     * \return The resulting polyhedron object after the projection
     */
    HPolyhedron Projection(size_t var_index, double base_tol = 1e-6) const;

    /**
     * \brief Projection of the polyhedron on a set of coordinates.
     * Keeps exactly these coordinates, in the requested order.
     * \param[in] coordinates Coordinates to project onto
     * \param[in] base_tol Numerical tolerance for redundancy removal
     * \return The resulting polyhedron object after the projection
     */
    HPolyhedron ProjectOnto(const std::vector<size_t>& coordinates, double base_tol = 1e-6) const;

    /**
     * \brief Affine transformation of a polyhedron through a linear map:
     * {T*x | x belongs to this set};
     * QR handles square invertible maps and falls back to SVD otherwise
     * \param[in] T The linear map
     * \param[in] condition_threshold Threshold for the matrix condition number in order
     * to select the method for the affine transformation function
     * \return The transformed polyhedron
     */
    HPolyhedron AffineTransform(const MatrixXd& T, double condition_threshold = 1e8) const;

    /**
     * \brief Affine transformation of a polyhedron through a linear map using SVD
     * to compute the pseudo-inverse; SVD enforces the image subspace and projects
     * out kernel coordinates
     *
     * 1. Compute T⁺ (pseudo-inverse) using SVD: T = U*Σ*V^T
     * 2. Find nullspace basis N of T
     * 3. Express original variables as x = T⁺*y + N*z
     * 4. Substitute into constraints and project out auxiliary variables z
     *
     * \param[in] T The linear map
     * \return The transformed polyhedron
     */
    HPolyhedron AffineTransform_SVD(const MatrixXd& T) const;

    /**
     * \brief Affine transformation of a polyhedron through a linear map using QR
     * decomposition with column pivoting for numerical stability and handling of
     * rank-deficient transformations
     * \param[in] T The linear map
     * \return The transformed polyhedron
     */
    HPolyhedron AffineTransform_QR(const MatrixXd& T) const;

    /**
     * \brief Computes P minus T*W, using h_W(G*T) directly, without forming T*W.
     * W must be valid and nonempty. T may be rectangular, singular or zero.
     * Inequality normals are retained; equalities are robustified too. Equality-image
     * widths up to 1e-9 times the normal's infinity norm are treated as numerical zero.
     *
     * \param[in] W A polyhedron object
     * \param[in] T A linear map
     * \return A polyhedron object that is the P - T * W
     */
    HPolyhedron PontryaginDifferenceOfLinearImage(const HPolyhedron& W, const MatrixXd& T) const;

    /**
     * \brief Exact Minkowski sum between this polyhedron and another.
     * Operands must be valid and nonempty (including canonical FullSpace).
     * No feasibility probes are performed to enforce this precondition.
     * Invalid/mismatched operands invalidate the result.
     * \param[in] other A polyhedron object
     * \return A polyhedron object corresponding to the Minkowski sum
     */
    HPolyhedron& operator+=(const HPolyhedron& other);

    /**
     * \brief Exact Pontryagin difference between this polyhedron and another.
     * Operands must be valid and nonempty (including canonical FullSpace).
     * No feasibility probes are performed to enforce this precondition.
     * Invalid/mismatched operands invalidate the result.
     * \param[in] other A polyhedron object
     * \return A polyhedron object corresponding to the Pontryagin difference
     */
    HPolyhedron& operator-=(const HPolyhedron& other);

    /* Getters */
    const MatrixXd& Ai() const { return Ai_; }
    const VectorXd& bi() const { return bi_; }
    const MatrixXd& Ae() const { return Ae_; }
    const VectorXd& be() const { return be_; }
    size_t NumInequalities() const { return static_cast<size_t>(Ai_.rows()); }
    size_t NumEqualities() const { return static_cast<size_t>(Ae_.rows()); }
    size_t Dimension() const { return static_cast<size_t>(Ai_.cols()); }
    bool IsValid() const { return valid_; }

    /**
     * \brief Check if the polyhedron represents the full space R^n.
     * Representation check only: redundant descriptions of R^n need not
     * satisfy this predicate until RemoveRedundantConstraints is called.
     * \return Boolean indicating whether the polyhedron represents the full space (true) or not (false)
     */
    inline bool IsFullSpace() const {
        return valid_ && Ai_.rows() == 0 && Ae_.rows() == 0;
    }

private:
    bool valid_ = false;
    MatrixXd Ai_;
    VectorXd bi_;
    MatrixXd Ae_;
    VectorXd be_;
};

// Copy only the modified operand.
inline HPolyhedron operator+(HPolyhedron lhs, const HPolyhedron& rhs) {
    lhs += rhs;
    return lhs;
}

inline HPolyhedron operator-(HPolyhedron lhs, const HPolyhedron& rhs) {
    lhs -= rhs;
    return lhs;
}

}
