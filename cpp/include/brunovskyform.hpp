#pragma once

#include "hpolyhedron.hpp"
#include <Eigen/Dense>
#include <vector>

namespace cis2m {
class BrunovskyForm {
  public:
    /// Constructor
    BrunovskyForm(const Eigen::MatrixXd &Ad, const Eigen::MatrixXd &Bd);

    BrunovskyForm(const Eigen::MatrixXd &Ad, const Eigen::MatrixXd &Bd,
                  const Eigen::MatrixXd &Ed);

    /// Destructor
    ~BrunovskyForm();

    /**
     * Get the Brunovsky Form A and B matrices
     */
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> GetBrunovskySystem();

    /**
     * Get the original A and B matrices
     */
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> GetOriginalSystem();

    /**
     * Transform the dynamics constraints expressed with
     * the Controllability Form basis
     */
    HPolyhedron GetStateConstraints(const HPolyhedron &StateCnstr);

    /**
     * Transform the input constraints considering the transformation
     * used to get the Brunovksy form.
     */
    HPolyhedron GetInputConstraints(const HPolyhedron &InputCnstr);

    /**
     * Return the disturbance matrix of the system in Brunovsky space
     */
    Eigen::MatrixXd GetDisturbanceMatrixBrunovsky();

    /**
     * Return the disturbance matrix of the system in original space
     */
    Eigen::MatrixXd GetDisturbanceMatrixOriginal();

    /**
     * Get the Am, Bm matrices used in the transformation from Controller Form
     * to Brunovksy Form. A_brunovsky = A_controllable + B_brunovsky * Am
     * B_brunovksy = B_controllable * Bm
     */
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> GetIntermediateMatrices();

    /**
     * Return the controllability indexes
     */
    std::vector<int> GetControllabilityIndexes();

    /**
     *  Get the Transformation Matrix (Change of basis to get the
     *  Controllability Form)
     */
    Eigen::MatrixXd GetTransformationMatrix();

    /**
     * Return the max controllability index
     */
    int GetMaxControllabilityIndex();

    /**
     * Return whether the system has disturbance
     */
    bool hasDisturbance();

  private:
    /// System in original space
    Eigen::MatrixXd Ad_;
    Eigen::MatrixXd Bd_;
    Eigen::MatrixXd Ed_;

    /// Dynamics A matrix in controllability form
    Eigen::MatrixXd A_cntrl_form_;

    /// Dynamics B matrix in controllability form
    Eigen::MatrixXd B_cntrl_form_;

    /// System in Brunovsky form
    Eigen::MatrixXd A_brunovsky_form_;
    Eigen::MatrixXd B_brunovsky_form_;
    Eigen::MatrixXd E_brunovsky_form_;

    /// Controllability Indexes
    std::vector<int> controllability_indexes_;

    /// Max Controllability Index
    int max_controllability_index_;

    /// Transformation Matrix
    Eigen::MatrixXd TransformationMatrix_;

    /// Brunovsky Form Bm Matrix
    Eigen::MatrixXd Bm_;

    /// Brunovsky Form Am Matrix
    Eigen::MatrixXd Am_;

    /// Indicator if the system has disturbance
    bool hasDisturbance_;
};
} // namespace cis2m
