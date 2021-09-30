#include "brunovskyform.hpp"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace cis2m {
// Helper Functions

// Compute the sigmas
std::vector<int>
ComputeSigmas(const std::vector<int> &controllability_indexes) {
  int N = controllability_indexes.size();
  std::vector<int> sigma_v(N, 0);

  if (N == 0)
    return sigma_v;

  sigma_v[0] = controllability_indexes[0];
  for (int i = 1; i < controllability_indexes.size(); i++) {
    sigma_v[i] = sigma_v[i - 1] + controllability_indexes[i];
  }
  return sigma_v;
}

// Compute Controllability Indexes
std::pair<MatrixXd, std::vector<int>>
ComputeControllabilityIndexes(const MatrixXd &ControllabilityMatrix,
                              int ControlNum) {

  int StateSize = ControllabilityMatrix.rows();
  int maxCols = ControllabilityMatrix.cols();
  int indexC = 0;

  // Allocate the temp matrix
  MatrixXd S(MatrixXd::Zero(StateSize, StateSize));
  int indexS = 0;
  int rank = 0;

  std::vector<std::vector<int>> indexes(ControlNum, std::vector<int>{});

  int watchdog = 0;
  // Selecting the first linearly independent StateSize columns from the
  // controllabity matrix
  while (rank < StateSize && watchdog < maxCols) {
    S.col(indexS) = ControllabilityMatrix.col(indexC);
    // Compute the rank
    Eigen::ColPivHouseholderQR<MatrixXd> lu_decomp(S);
    lu_decomp.setThreshold(1e-10);
    int nrank = lu_decomp.rank();
    if (nrank > rank) {
      rank = nrank;
      indexes[indexC % ControlNum].push_back(indexS);
      indexS++;
    }
    indexC++;
  }

  if (rank < StateSize) {
    // This should not happen!
    std::cerr << "The system is not controllable!" << std::endl;
  }

  // Reordering the column of S to make the Cbar matrix [ref. Antsklis Linear
  // Systems pg. 283]
  MatrixXd Cbar(MatrixXd::Zero(StateSize, StateSize));
  std::vector<int> ctrl_indexes;
  int counter = 0;
  for (auto &vel : indexes) {
    int ctrl_ind_i = 0;
    for (auto &el : vel) {
      Cbar.col(counter++) = ControllabilityMatrix.col(el);
      ctrl_ind_i++;
    }
    ctrl_indexes.push_back(ctrl_ind_i);
  }

  // Check the condition number of the matrix before attempting an inversion
  Eigen::JacobiSVD<MatrixXd> svd(Cbar);
  double cond = svd.singularValues()(0) /
                svd.singularValues()(svd.singularValues().size() - 1);

  if (cond > 1e10) {
    std::cerr << "Matrix singularity, unsafe to invert!" << std::endl;
  }

  std::pair<MatrixXd, std::vector<int>> output(Cbar, ctrl_indexes);

  return output;
}

// Compute the change of basis tranformation
MatrixXd ComputeTransformation(const MatrixXd &A, const MatrixXd &B,
                               const MatrixXd &Cbar, std::vector<int> &sigmas) {
  int StateSize = A.cols();

  // Invert the C
  MatrixXd invCbar = Cbar.inverse();

  MatrixXd P(MatrixXd::Zero(StateSize, StateSize));
  int counter = 0;
  for (int i = 0; i < sigmas.size(); i++) {
    P.row(counter++) = invCbar.row(sigmas[i] - 1);
    int mu = i > 0 ? sigmas[i] - sigmas[i - 1] : sigmas[i];
    for (int j = 1; j < mu; j++) {
      P.row(counter) = P.row(counter - 1) * A;
      counter++;
    }
  }

  return P;
}

bool isControllable(const MatrixXd &C, int fullrank, double tol) {
  bool out = false;
  Eigen::ColPivHouseholderQR<MatrixXd> lu_decomp(C);
  lu_decomp.setThreshold(tol);
  if (lu_decomp.rank() != fullrank) {
    std::cerr << "[" << __FILE__ << " @ " << __LINE__
              << "] System is not Controllable!" << std::endl;
  } else {
    out = true;
  }
  return out;
}

MatrixXd ComputeCTRL(const MatrixXd &A, const MatrixXd &B) {
  const int num_states = B.rows();
  const int num_inputs = B.cols();

  MatrixXd CTRL(MatrixXd::Zero(num_states, num_states * num_inputs));
  CTRL.leftCols(num_inputs) = B;
  for (int i = 1; i < num_states; i++) {
    CTRL.middleCols(num_inputs * i, num_inputs) =
        A * CTRL.middleCols(num_inputs * (i - 1), num_inputs);
  }

  if (!isControllable(CTRL, num_states, 0.001)) {
    throw;
  }

  return CTRL;
}

// ============================================================
// CLASS
BrunovskyForm::BrunovskyForm(const MatrixXd &Ad, const MatrixXd &Bd)
    : A_cntrl_form_(MatrixXd::Zero(Ad.rows(), Ad.cols())),
      B_cntrl_form_(MatrixXd::Zero(Bd.rows(), Bd.cols())),
      TransformationMatrix_(MatrixXd::Zero(Ad.rows(), Ad.rows())),
      Am_(MatrixXd::Zero(Bd.cols(), Ad.cols())),
      Bm_(MatrixXd::Zero(Bd.cols(), Bd.cols())) {

  // Compute the controllability matrix
  MatrixXd CTRL = ComputeCTRL(Ad, Bd);

  // Compute the controllability indexes and the intermediate Cbar matrix
  std::pair<MatrixXd, std::vector<int>> data =
      ComputeControllabilityIndexes(CTRL, Bd.cols());
  MatrixXd Cbar = data.first;
  controllability_indexes_ = data.second;
  max_controllability_index_ = *std::max_element(
      controllability_indexes_.begin(), controllability_indexes_.end());
  std::vector<int> sigmas = ComputeSigmas(controllability_indexes_);

  TransformationMatrix_ = ComputeTransformation(Ad, Bd, Cbar, sigmas);

  // Compute the controller forms
  A_cntrl_form_ = TransformationMatrix_ * Ad * TransformationMatrix_.inverse();
  B_cntrl_form_ = TransformationMatrix_ * Bd;

  for (int i = 0; i < sigmas.size(); i++) {
    Bm_.row(i) = B_cntrl_form_.row(sigmas[i] - 1);
    Am_.row(i) = A_cntrl_form_.row(sigmas[i] - 1);
  }

  // Generate the Brunovksy Forms after Feedback
  // A_brunovsky_form_ = MatrixXd::Zero(Ad.rows(), Ad.cols());
  // B_brunovsky_form_ = MatrixXd::Zero(Bd.rows(), Bd.cols());
  // int start_index = 0;
  // for (int i = 0; i < controllability_indexes_.size(); ++i) {
  //   int mu_i = controllability_indexes_[i];
  //   A_brunovsky_form_.block(start_index, start_index + 1, mu_i - 1, mu_i - 1)
  //   =
  //       Eigen::MatrixXd::Identity(mu_i - 1, mu_i - 1);
  //   B_brunovsky_form_(start_index + mu_i - 1, i) = 1.0;
  //   start_index += mu_i;
  // }
  //
  // Use the forms:
  B_brunovsky_form_ = TransformationMatrix_ * Bd * Bm_.inverse();
  A_brunovsky_form_ =
      TransformationMatrix_ * Ad * TransformationMatrix_.inverse() -
      B_brunovsky_form_ * Am_;
}

BrunovskyForm::~BrunovskyForm(){};

std::pair<MatrixXd, MatrixXd> BrunovskyForm::GetDynSystem() {
  std::pair<MatrixXd, MatrixXd> output(A_brunovsky_form_, B_brunovsky_form_);
  return output;
}

HPolyhedron BrunovskyForm::GetDynConstraints(const HPolyhedron &DynCnstr) {
  // Return the HPolyhedron [F, b] F x < b
  MatrixXd F_state_cnstr = DynCnstr.Ai() * TransformationMatrix_.inverse();
  VectorXd b_state_cnstr = DynCnstr.bi();

  return HPolyhedron(F_state_cnstr, b_state_cnstr);
}

HPolyhedron BrunovskyForm::GetInputConstraints(const HPolyhedron &InputCnstr) {
  MatrixXd F_input_cnstr = InputCnstr.Ai() * Bm_.inverse();
  VectorXd b_input_cnstr = InputCnstr.bi();

  return HPolyhedron(F_input_cnstr, b_input_cnstr);
}

MatrixXd
BrunovskyForm::GetDisturbanceMatrix(const MatrixXd &DisturbanceMatrix) {
  MatrixXd E_cntrl_form = TransformationMatrix_ * DisturbanceMatrix;
  return E_cntrl_form;
}

std::pair<MatrixXd, MatrixXd> BrunovskyForm::GetIntermediateMatrixes() {
  std::pair<MatrixXd, MatrixXd> output(Am_, Bm_);
  return output;
}

MatrixXd BrunovskyForm::GetTransformationMatrix() {
  return TransformationMatrix_;
}

std::vector<int> BrunovskyForm::GetControllabilityIndexes() {
  return controllability_indexes_;
}

int BrunovskyForm::GetMaxControllabilityIndex() {
  return max_controllability_index_;
}

} // namespace cis2m
