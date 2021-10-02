#include "cis_generator.hpp"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace cis2m {
// Helper Functions
MatrixXd blkdiag(const MatrixXd &A, int count) {
  MatrixXd bdm = MatrixXd::Zero(A.rows() * count, A.cols() * count);
  for (int i = 0; i < count; ++i) {
    bdm.block(i * A.rows(), i * A.cols(), A.rows(), A.cols()) = A;
  }
  return bdm;
}

// ========================================================================
// CLASS
CISGenerator::CISGenerator(const MatrixXd &Ad, const MatrixXd &Bd) {
  Reset();

  StateDim_ = Ad.rows();
  NumberInputs_ = Bd.cols();
  GenerateBrunovksyForm(Ad, Bd);
}

CISGenerator::CISGenerator(const MatrixXd &Ad, const MatrixXd &Bd,
                           const MatrixXd &Ed) {
  Reset();

  StateDim_ = Ad.rows();
  NumberInputs_ = Bd.cols();
  DisturbanceDim_ = Ed.cols();

  GenerateBrunovksyForm(Ad, Bd);
  MatrixXd Transform = brunovsky_form_->GetTransformationMatrix();
  Ed_BR_ = Transform * Ed;
}

CISGenerator::~CISGenerator(){};

void CISGenerator::Reset() {
  DisturbanceDim_ = 0;

  Level_ = -1;
  Transient_ = -1;
  ExtendedInputDim_ = -1;

  ThereAreInputConstraints_ = false;

  cis_computed_ = false;
}

void CISGenerator::GenerateBrunovksyForm(const MatrixXd &A, const MatrixXd &B) {
  brunovsky_form_ = new BrunovskyForm(A, B);
}

void CISGenerator::AddDisturbanceSet(const HPolyhedron &ds) {
  DisturbanceSet_ = ds;
}

void CISGenerator::AddInputConstraintsSet(const HPolyhedron &ics) {
  ThereAreInputConstraints_ = true;
  InputCnstrSet_ = ics;
}

std::vector<HPolyhedron>
CISGenerator::ComputeShrinkedSafeSetsSequence(const HPolyhedron &ss) {
  int nmax = brunovsky_form_->GetMaxControllabilityIndex();

  std::pair<MatrixXd, MatrixXd> pairAB = brunovsky_form_->GetDynSystem();
  MatrixXd Ad_BF = pairAB.first;
  HPolyhedron SafeSet_BF = brunovsky_form_->GetDynConstraints(ss);

  std::vector<HPolyhedron> SafeSet_seq;
  SafeSet_seq.push_back(SafeSet_BF);
  if (Ed_BR_.size() > 0) {
    MatrixXd DynMat = MatrixXd::Identity(StateDim_, StateDim_);
    for (int i = 1; i < nmax; i++) {
      HPolyhedron Sub(DisturbanceSet_.affineT(DynMat * Ed_BR_));
      SafeSet_seq.push_back(SafeSet_BF - Sub);
      DynMat *= Ad_BF;
    }
  }

#ifdef CIS2M_DEBUG
  std::cout << __FILE__ << std::endl;
  std::cout << "Computation of Safe set sequence" << std::endl;
  int i = 0;
  for (auto &el : SafeSet_seq) {
    std::cout << "S_k " << el.Ai().rows() << " x " << el.Ai().cols()
              << std::endl;
    std::cout << "SafeSet [" << i << "] Base A: " << std::endl
              << el.Ai() << std::endl;
    std::cout << "SafeSet [" << i++ << "]Base B: " << std::endl
              << el.bi() << std::endl;
    std::cout << std::endl;
  }
#endif
  return SafeSet_seq;
}

void CISGenerator::ComputeLiftedSystem(int L, int T) {
  // Update the parameters
  Level_ = L;
  Transient_ = T;
  int length = Transient_ + Level_;
  ExtendedInputDim_ = length * NumberInputs_;

  // Construct the High-Dimensional system
  MatrixXd Ki(MatrixXd::Zero(1, length));
  Ki(0) = 1.0;
  MatrixXd Pi(MatrixXd::Zero(length, length));
  Pi.block(0, 1, length - 1, length - 1) =
      MatrixXd::Identity(length - 1, length - 1);
  Pi(length - 1, Transient_) = 1.0;
  ExtendedU2U_ = blkdiag(Ki, NumberInputs_);
  MatrixXd P = blkdiag(Pi, NumberInputs_);

  std::pair<MatrixXd, MatrixXd> pairAB = brunovsky_form_->GetDynSystem();
  MatrixXd Ad_BF = pairAB.first;
  MatrixXd Bd_BF = pairAB.second;

  int Nrow_hd = StateDim_ + ExtendedInputDim_;
  int Ncol_hd = StateDim_ + P.cols();
  MatrixXd Ahd(MatrixXd::Zero(Nrow_hd, Ncol_hd));
  Ahd.block(0, 0, StateDim_, StateDim_) = Ad_BF;
  Ahd.block(0, StateDim_, StateDim_, ExtendedU2U_.cols()) =
      Bd_BF * ExtendedU2U_;
  Ahd.block(StateDim_, StateDim_, P.rows(), P.cols()) = P;

  A_lifted_ = Ahd;
#ifdef CIS2M_DEBUG
  std::cout << __FILE__ << std::endl;
  std::cout << "Computation of Lifted System" << std::endl;
  std::cout << "A_BF: " << std::endl << Ad_BF << std::endl;
  std::cout << "B_BF: " << std::endl << Bd_BF << std::endl;
  std::cout << "A_lifted: " << std::endl << Ahd << std::endl;
  std::cout << "A_lifted Size: " << Ahd.rows() << " x " << Ahd.cols()
            << std::endl;
  std::cout << std::endl;
#endif
}

void CISGenerator::computeCIS(const HPolyhedron &SafeSet, int L, int T) {
  if (Level_ != L || Transient_ != T)
    ComputeLiftedSystem(L, T);

  // Compute the shrinked Free Space
  std::vector<HPolyhedron> seq = ComputeShrinkedSafeSetsSequence(SafeSet);

  HPolyhedron polyhedron_Gb = brunovsky_form_->GetDynConstraints(SafeSet);

  int NDynconstr = polyhedron_Gb.Ai().rows();
  int length = Level_ + Transient_;
  int mu_max = brunovsky_form_->GetMaxControllabilityIndex();

  // Compute the mcisA; mcisB matrixes
  int mcisA_rows = NDynconstr * (mu_max + length);
  int mcisA_cols = StateDim_ + ExtendedInputDim_;
  MatrixXd mcisA(MatrixXd::Zero(mcisA_rows, mcisA_cols));
  VectorXd mcisb(VectorXd::Zero(mcisA_rows));

#ifdef CIS2M_DEBUG
  std::cout << " mcisA: " << mcisA_rows << " x " << mcisA_cols << std::endl;
#endif
  mcisA.block(0, 0, NDynconstr, StateDim_) = polyhedron_Gb.Ai();
  mcisb.head(NDynconstr) = polyhedron_Gb.bi();

  MatrixXd Acurr = A_lifted_;
  MatrixXd TempA(MatrixXd::Zero(NDynconstr, mcisA_cols));

  MatrixXd mcisA_ctrl;
  VectorXd mcisb_ctrl;
  MatrixXd AA;
  MatrixXd BB;
  HPolyhedron FGu;

  if (ThereAreInputConstraints_) {
    FGu = brunovsky_form_->GetInputConstraints(InputCnstrSet_);
    mcisA_ctrl = MatrixXd::Zero(ExtendedInputDim_, mcisA_cols);
    mcisb_ctrl = VectorXd::Zero(ExtendedInputDim_);

    std::pair<MatrixXd, MatrixXd> AmBm =
        brunovsky_form_->GetIntermediateMatrixes();

    AA = -FGu.Ai() * AmBm.first;
    BB = FGu.Ai() * ExtendedU2U_;
    mcisA_ctrl.block(0, 0, NumberInputs_, StateDim_) = AA;
    mcisA_ctrl.block(0, StateDim_, NumberInputs_, mcisA_cols - StateDim_) = BB;
    mcisb_ctrl.head(NumberInputs_) = FGu.bi();
  }

  for (int t = 1; t < mu_max + length; t++) {
    int tbar = t < seq.size() ? t : seq.size() - 1;
    TempA.leftCols(StateDim_) = seq[tbar].Ai();

    mcisA.block(NDynconstr * t, 0, NDynconstr, mcisA_cols) = TempA * Acurr;
    mcisb.segment(NDynconstr * t, NDynconstr) = seq[tbar].bi();

    if (ThereAreInputConstraints_) {
      if (t < length) {
        mcisA_ctrl.block(t * NumberInputs_, 0, NumberInputs_, StateDim_) = AA;
        mcisA_ctrl.block(t * NumberInputs_, StateDim_, NumberInputs_,
                         mcisA_cols - StateDim_) =
            BB * Acurr.block(StateDim_, StateDim_, ExtendedInputDim_,
                             ExtendedInputDim_);
        mcisb_ctrl.segment(t * NumberInputs_, NumberInputs_) = FGu.bi();
      }
    }
    Acurr *= A_lifted_;
  }
  std::cout << std::endl << A_lifted_ << std::endl;

  MatrixXd mcisA_tot(mcisA.rows() + mcisA_ctrl.rows(), mcisA.cols());
  VectorXd mcisb_tot(mcisA_tot.rows());

  if (ThereAreInputConstraints_) {
    mcisA_tot << mcisA, mcisA_ctrl;
    mcisb_tot << mcisb, mcisb_ctrl;
  } else {
    mcisA_tot << mcisA;
    mcisb_tot << mcisb;
  }
  // Get the transformation
  //   MatrixXd Transform = brunovsky_form_->GetTransformationMatrix();
  //   mcisA_tot.block(0, 0, mcisA_tot.rows(), StateDim_) =
  //       mcisA_tot.block(0, 0, mcisA_tot.rows(), StateDim_) * Transform;
  //   // Note: The CIS is expressed in the original basis
  CIS_ = HPolyhedron(mcisA_tot, mcisb_tot);
  cis_computed_ = true;

#ifdef CIS2M_DEBUG
  std::cout << __FILE__ << std::endl;
  std::cout << "Computation of CIS" << std::endl;
  std::cout << "Transf: " << std::endl << Transform << std::endl;
  std::cout << "CIS Size: " << mcisA.rows() << " x " << mcisA.cols()
            << std::endl;
  std::cout << "CIS_A_BF: " << std::endl << mcisA << std::endl;
  std::cout << "CIS_b_BF: " << std::endl << mcisb.transpose() << std::endl;
  std::cout << std::endl;
#endif
}

HPolyhedron CISGenerator::Fetch_CIS() {
  if (cis_computed_) {
    return CIS_;
  } else {
    return {};
  }
}

MatrixXd CISGenerator::Fetch_A_State() {
  if (cis_computed_) {
    return CIS_.Ai().leftCols(StateDim_);
  } else {
    return {};
  }
}

MatrixXd CISGenerator::Fetch_A_Input() {
  if (cis_computed_) {
    int NumRows = CIS_.Ai().rows();
    return CIS_.Ai().block(0, StateDim_, NumRows, NumberInputs_);
  } else {
    return {};
  }
}

MatrixXd CISGenerator::Fetch_A_Virtual() {
  if (cis_computed_) {
    int NumCols = CIS_.Ai().cols();
    return CIS_.Ai().rightCols(NumCols - (NumberInputs_ + StateDim_));
  } else {
    return {};
  }
}

MatrixXd CISGenerator::Fetch_A_lifted() { return A_lifted_; }

VectorXd CISGenerator::TransformU2B(const VectorXd &u, const VectorXd &x) {
  MatrixXd T = brunovsky_form_->GetTransformationMatrix();
  std::pair<MatrixXd, MatrixXd> AmBm =
      brunovsky_form_->GetIntermediateMatrixes();
  return AmBm.second * u + AmBm.first * T * x;
}

VectorXd CISGenerator::TransformU2O(const VectorXd &u, const VectorXd &x) {
  MatrixXd T = brunovsky_form_->GetTransformationMatrix();
  std::pair<MatrixXd, MatrixXd> AmBm =
      brunovsky_form_->GetIntermediateMatrixes();
  return AmBm.second.inverse() * (u - AmBm.first * T * x);
}

int CISGenerator::GetExtendedDim() const { return NumberInputs_ * Level_; }

int CISGenerator::GetStateDim() const { return StateDim_; }

int CISGenerator::GetLevel() const { return Level_; }

BrunovskyForm *CISGenerator::getBrunovskyForm() { return brunovsky_form_; }
} // namespace cis2m
