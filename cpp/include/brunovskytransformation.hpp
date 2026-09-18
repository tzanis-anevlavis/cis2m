#pragma once

#include <cstddef>
#include <vector>

#include <Eigen/Dense>

#include "hpolyhedron.hpp"

namespace cis2m {

/**
 * \brief Brunovsky normal-form transformation for a controllable system.
 *
 * For x+ = A x + B u, the class computes
 *
 *     z = T x,
 *     r = Am z + Bm u,
 *     z+ = Ac z + Bc r,
 *
 * where (Ac, Bc) consists of consecutive Brunovsky chains. The input matrix
 * B must have full column rank so that Bm is invertible.
 */
class BrunovskyTransformation {
public:
    /**
     * \brief Construct the normal-form transformation for (A, B).
     * Invalid dimensions, nonfinite data, rank-deficient B, and
     * uncontrollable systems produce an invalid object.
     * \param[in] A State matrix
     * \param[in] B Input matrix
     */
    BrunovskyTransformation(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B);

    /**
     * \brief Transform a state set from x coordinates to z = T x.
     * \param[in] state_set Polyhedron in the original state coordinates
     * \return The same set represented in Brunovsky state coordinates
     */
    HPolyhedron TransformStateSet(const HPolyhedron& state_set) const;

    /**
     * \brief Transform an input set from u coordinates to s = Bm u.
     * This transforms only the linear input-coordinate part. The complete
     * virtual input is r = Am z + s.
     * \param[in] input_set Polyhedron in the physical input coordinates
     * \return The same set represented in s coordinates
     */
    HPolyhedron TransformInputSet(const HPolyhedron& input_set) const;

    /**
     * \brief Transform a joint state-input set from (x, u) to (z, r).
     * Both inequality and equality constraints are preserved under
     * x = T^-1 z and u = Bm^-1 (r - Am z).
     * \param[in] joint_set Polyhedron in joint (x, u) coordinates
     * \return The same set represented in joint (z, r) coordinates
     */
    HPolyhedron TransformJointStateInputSet(const HPolyhedron& joint_set) const;

    /**
     * \brief Transform the disturbance matrix as Ec = T E.
     * \param[in] E Disturbance matrix in the original state coordinates
     * \return Disturbance matrix in Brunovsky state coordinates
     */
    Eigen::MatrixXd TransformDisturbanceMatrix(const Eigen::MatrixXd& E) const;

    /**
     * \brief Get the canonical state matrix Ac.
     */
    const Eigen::MatrixXd& CanonicalStateMatrix() const {
        return Ac_;
    }

    /**
     * \brief Get the canonical input matrix Bc.
     */
    const Eigen::MatrixXd& CanonicalInputMatrix() const {
        return Bc_;
    }

    /**
     * \brief Get the controllability index of each input channel.
     */
    const std::vector<std::size_t>& ControllabilityIndices() const {
        return controllability_indices_;
    }

    /**
     * \brief Get the nilpotency index of Ac.
     */
    std::size_t MaxControllabilityIndex() const {
        return max_controllability_index_;
    }

    /**
     * \brief Get the state transformation T in z = T x.
     */
    const Eigen::MatrixXd& TransformationMatrix() const {
        return T_;
    }

    /**
     * \brief Get the state-feedback matrix Am in r = Am*z + Bm*u.
     */
    const Eigen::MatrixXd& StateFeedbackMatrix() const {
        return Am_;
    }

    /**
     * \brief Get the input transformation matrix Bm in r = Am*z + Bm*u.
     */
    const Eigen::MatrixXd& InputTransformationMatrix() const {
        return Bm_;
    }

    /**
     * \brief Check whether the transformation was constructed successfully.
     */
    bool IsValid() const {
        return valid_;
    }

private:
    /**
     * \brief Constructs the Controllability matrix of a linear system described by the
     * pair of (A, B) matrices.
     * \param[in] A State matrix
     * \param[in] B Input matrix
     * \return Controllability matrix of the linear system (A, B)
     */
    Eigen::MatrixXd ControllabilityMatrix(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B) const;

    /**
     * \brief Computes the controllability indices given a Controllability matrix and the
     * dimension of the input. It populates the internal `controllability_indices_`.
     * \param[in] C The controllability matrix
     * \param[in] m Dimension of the input
     * \return Boolean flag indicating whether the operation was successful (true) or not (false)
     */
    bool ComputeControllabilityIndices(const Eigen::MatrixXd& C, Eigen::Index m);

    /**
     * \brief Computes a controllability basis based on the state matrix A, input matrix B, and the
     * controllability indices.
     * \param[in] A State matrix
     * \param[in] B Input matrix
     * \param[in] indices The controllability indices
     * \return The controllability basis as columns of a matrix
     */
    Eigen::MatrixXd ControllabilityBasis(
        const Eigen::MatrixXd& A,
        const Eigen::MatrixXd& B,
        const std::vector<std::size_t>& indices) const;

    /**
     * \brief Builds the internal Ac_, Bc_ matrices.
     */
    void BuildCanonicalSystem();

    bool valid_ = false;
    Eigen::MatrixXd Ac_;
    Eigen::MatrixXd Bc_;
    Eigen::MatrixXd T_;
    Eigen::MatrixXd Am_;
    Eigen::MatrixXd Bm_;
    std::vector<std::size_t> controllability_indices_;
    std::size_t max_controllability_index_ = 0;
};

}
