#pragma once

#include "brunovskyform.hpp"
#include "hpolyhedron.hpp"
#include <Eigen/Dense>

namespace cis2m {
class CISGenerator {
public:
  /**
   * \brief Constructor with only dynamics information
   * x_new = Ad * x_old + Bd * u
   *
   * \param[in]	Ad	Dynamics matrix of the system
   * \param[in]	Bd	Input matrix of the system
   *
   */
  CISGenerator(const int L, const int T, const Eigen::MatrixXd &Ad,
               const Eigen::MatrixXd &Bd);

  /**
   * \brief Constructor with only dynamics information
   * x_new = Ad * x_old + Bd * u + Ed * w
   *
   * \param[in]	Ad	Dynamics matrix of the system
   * \param[in]	Bd	Input matrix of the system
   * \param[in]	Ed	Disturbance matrix
   *
   */
  CISGenerator(const int L, const int T, const Eigen::MatrixXd &Ad,
               const Eigen::MatrixXd &Bd, const Eigen::MatrixXd &Ed);

  /// Destructor
  ~CISGenerator();

  /**
   * \brief Add disturbance information
   *
   * \param[in]	DisturbanceSet	Disturbance set described as a HPolyhedron
   *
   * \return 	void
   *
   */
  void AddDisturbanceSet(const HPolyhedron &DisturbanceSet);

  /**
   * \brief Add input constraints set
   *
   * \param[in]	ics	Input Constraints described as a HPolyhedron
   *
   * \return 	void
   *
   */
  void AddInputConstraintsSet(const HPolyhedron &ics);

  /**
   * \brief Compute the Control Invariant Set
   *
   * \param[in]	SafeSet 	Safe set described as a HPolyhedron
   * \param[in]	L		Period of control policy
   * \param[in]	T		Transient of control policy
   *
   * \return	Control Invariant Set as a HPolyhedron
   *
   */
  void computeCIS(const HPolyhedron &SafeSet, int L, int T);

  /**
   * \brief Fetch the state part of the Constraint Coefficients
   *
   * \return	Control Invariant Set as a HPolyhedron
   *
   */
  HPolyhedron Fetch_CIS();

  /**
   * \brief Fetch the state part of the Constraint Coefficients
   *
   * \return	Eigen::MatrixXd
   *
   */
  Eigen::MatrixXd Fetch_A_State();

  /**
   * \brief Fetch the input part of the Constraint Coefficients
   *
   * \return	Eigen::MatrixXd
   *
   */
  Eigen::MatrixXd Fetch_A_Input();

  /**
   * \brief Fetch the virtual input part of the Constraint Coefficients
   *
   * \return	Eigen::MatrixXd
   *
   */
  Eigen::MatrixXd Fetch_A_Virtual();

  /**
   * \brief Fetch the matrix A of the lifted companion system.
   *
   * \return	Eigen::MatrixXd
   *
   */
  Eigen::MatrixXd Fetch_A_lifted();

  /**
   * \brief Compute the input in the Brunovksy coordinates
   *
   * \param[in]	u	Input in the original coordinates
   * \param[in]	x	State in the original coordinates
   *
   * \return	Eigen::VectorXd	Input in Brunovksy coordinates
   */
  Eigen::VectorXd TransformU2B(const Eigen::VectorXd &u,
                               const Eigen::VectorXd &x);

  /**
   * \brief Compute the input in the Original coordinates
   *
   * \param[in]	u	Input in the Brunovksy coordinates
   * \param[in]	x	State in the original coordinates
   *
   * \return	Eigen::VectorXd	Input in Original coordinates
   */
  Eigen::VectorXd TransformU2O(const Eigen::VectorXd &u,
                               const Eigen::VectorXd &x);

  /**
   * \brief Get the size of the virtual input
   */
  int GetExtendedDim() const;

  /**
   * \brief Get the size of the state
   */
  int GetStateDim() const;

  int GetLevel() const;

  /**
   * \brief Get the Brunovsky Form Transformation Class
   */
  BrunovskyForm *getBrunovskyForm();

private:
  int StateDim_;
  int NumberInputs_;
  int DisturbanceDim_;

  int ExtendedInputDim_;

  int Loop_;
  int Transient_;

  bool ThereAreInputConstraints_;

  HPolyhedron InputCnstrSet_;
  HPolyhedron DisturbanceSet_;

  Eigen::MatrixXd A_lifted_;

  HPolyhedron CIS_;

  Eigen::MatrixXd Ed_BF_;

  Eigen::MatrixXd ExtendedU2U_;

  /**
   * \brief Reference to the Brunovsky Form Transformation Class
   */
  BrunovskyForm *brunovsky_form_;

  void Reset();

  void ComputeLiftedSystem(int L, int T);

  void NewComputeLiftedSystem(int L, int T,
                              std::vector<Eigen::MatrixXd> &SysMatrices);

  void ComputeExtended(HPolyhedron &SafeSet_BF,
                       std::vector<Eigen::MatrixXd> &SysMatrices_BF);

  HPolyhedron NewMinkDiff(HPolyhedron &P, HPolyhedron &S);

  std::vector<HPolyhedron>
  ComputeSafeSetsSequence(const HPolyhedron &SafeSet,
                          std::vector<Eigen::MatrixXd> &SysMatrices);

  void GenerateBrunovksyForm(const Eigen::MatrixXd &A,
                             const Eigen::MatrixXd &B);

  void GenerateBrunovksyForm(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B,
                             const Eigen::MatrixXd &E);

  bool cis_computed_;
};
} // namespace cis2m
