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
CISGenerator::CISGenerator(const int L, const int T, const MatrixXd &Ad,
                           const MatrixXd &Bd) {
  Reset();

  Loop_ = L;
  Transient_ = T;
  StateDim_ = Ad.rows();
  NumberInputs_ = Bd.cols();
  GenerateBrunovksyForm(Ad, Bd);
}

CISGenerator::CISGenerator(const int L, const int T, const MatrixXd &Ad,
                           const MatrixXd &Bd, const MatrixXd &Ed) {
  Reset();

  Loop_ = L;
  Transient_ = T;
  StateDim_ = Ad.rows();
  NumberInputs_ = Bd.cols();
  DisturbanceDim_ = Ed.cols();

  GenerateBrunovksyForm(Ad, Bd, Ed);
}

CISGenerator::~CISGenerator(){};

void CISGenerator::Reset() {
  DisturbanceDim_ = 0;

  Loop_ = -1;
  Transient_ = -1;
  ExtendedInputDim_ = -1;

  ThereAreInputConstraints_ = false;

  cis_computed_ = false;
}

void CISGenerator::GenerateBrunovksyForm(const MatrixXd &A, const MatrixXd &B) {
  brunovsky_form_ = new BrunovskyForm(A, B);
}

void CISGenerator::GenerateBrunovksyForm(const MatrixXd &A, const MatrixXd &B,
                                         const MatrixXd &E) {
  brunovsky_form_ = new BrunovskyForm(A, B, E);
}

void CISGenerator::AddDisturbanceSet(const HPolyhedron &DisturbanceSet) {
  DisturbanceSet_ = DisturbanceSet;
}

void CISGenerator::AddInputConstraintsSet(const HPolyhedron &ics) {
  ThereAreInputConstraints_ = true;
  InputCnstrSet_ = ics;
}

HPolyhedron CISGenerator::NewMinkDiff(HPolyhedron &P, HPolyhedron &S) {
  // Returns P - S

  if (P.GetSpaceDim() != S.GetSpaceDim()) {
    std::cerr << "The polyhedra should have the same dimension" << std::endl;
    return P;
  }

  VectorXd supp = S.ComputeSupport(P.Ai());

  MatrixXd tempA(P.Ai().rows() + P.Ai().rows(), P.Ai().cols());
  tempA << P.Ai(), P.Ai();

  VectorXd tempb(P.bi().rows() + P.bi().rows());
  tempb << P.bi(), (P.bi() - supp);

  HPolyhedron MinkDiff(tempA, tempb);

  return MinkDiff;
}

std::vector<HPolyhedron>
CISGenerator::ComputeSafeSetsSequence(const HPolyhedron &SafeSet,
                                      std::vector<MatrixXd> &SysMatrices) {
  int nmax = brunovsky_form_->GetMaxControllabilityIndex();

  // std::pair<MatrixXd, MatrixXd> pairAB = brunovsky_form_->GetSystem();
  //   MatrixXd Ad_BF = brunovsky_form_->GetSystem().first;
  //   HPolyhedron SafeSet_BF = brunovsky_form_->GetStateConstraints(ss);

  std::vector<HPolyhedron> SafeSet_seq;
  SafeSet_seq.push_back(SafeSet);
  if (brunovsky_form_->hasDisturbance()) {
    MatrixXd A_curr = MatrixXd::Identity(StateDim_, StateDim_);
    for (int i = 1; i < nmax + Transient_ + Loop_; i++) {
      if (i <= nmax) {
        HPolyhedron Sub(DisturbanceSet_.affineT(A_curr * SysMatrices[2]));
        // SafeSet_seq.push_back(SafeSet_seq[i - 1] - Sub);
        HPolyhedron Diff = SafeSet_seq[i - 1] - Sub;
        // HPolyhedron Diff = NewMinkDiff(SafeSet_seq[i - 1], Sub);
        SafeSet_seq.push_back(Diff);
        A_curr *= SysMatrices[0];

        std::cout << i << " - A: " << std::endl << Sub.Ai() << std::endl;
        std::cout << "b: " << std::endl << Sub.bi() << std::endl;
      } else {
        SafeSet_seq.push_back(SafeSet_seq[nmax]);
      }
    }
  } else {
    for (int i = 1; i < nmax + Transient_ + Loop_; i++) {
      SafeSet_seq.push_back(SafeSet);
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
  Loop_ = L;
  Transient_ = T;
  int length = Transient_ + Loop_;
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

  std::pair<MatrixXd, MatrixXd> pairAB = brunovsky_form_->GetSystem();
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

void CISGenerator::NewComputeLiftedSystem(int L, int T,
                                          std::vector<MatrixXd> &SysMatrices) {
  // Update the parameters
  Loop_ = L;
  Transient_ = T;
  int length = Transient_ + Loop_;
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

  //   std::pair<MatrixXd, MatrixXd> pairAB = brunovsky_form_->GetSystem();
  //   MatrixXd Ad_BF = pairAB.first;
  //   MatrixXd Bd_BF = pairAB.second;
  MatrixXd Ad_BF = SysMatrices[0];
  MatrixXd Bd_BF = SysMatrices[1];

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

void CISGenerator::ComputeExtended(HPolyhedron &SafeSet_BF,
                                   std::vector<MatrixXd> &SysMatrices_BF) {

  // Compute the extended system
  MatrixXd A_BF_extended(SysMatrices_BF[0].cols() + SysMatrices_BF[1].cols(),
                         SysMatrices_BF[0].cols() + SysMatrices_BF[1].cols());
  A_BF_extended << SysMatrices_BF[0], SysMatrices_BF[1],
      MatrixXd::Zero(SysMatrices_BF[1].cols(),
                     SysMatrices_BF[0].cols() + SysMatrices_BF[1].cols());
  MatrixXd B_BF_extended(SysMatrices_BF[0].cols() + SysMatrices_BF[1].cols(),
                         SysMatrices_BF[1].cols());
  B_BF_extended << MatrixXd::Zero(SysMatrices_BF[0].cols(),
                                  SysMatrices_BF[1].cols()),
      MatrixXd::Identity(SysMatrices_BF[1].cols(), SysMatrices_BF[1].cols());

  std::cout << std::endl << "Gotcha" << std::endl;

  if (SysMatrices_BF.size() == 3) {
    MatrixXd E_BF_extended(SysMatrices_BF[2].rows() + SysMatrices_BF[1].cols(),
                           SysMatrices_BF[2].cols());
    E_BF_extended << SysMatrices_BF[2],
        MatrixXd::Zero(SysMatrices_BF[1].cols(), SysMatrices_BF[2].cols()),

        SysMatrices_BF[2] = E_BF_extended;
  }
  SysMatrices_BF[1] = B_BF_extended;
  SysMatrices_BF[0] = A_BF_extended;

  // Compute the extended safe set
  MatrixXd Gxu(SafeSet_BF.Ai().rows() + InputCnstrSet_.Ai().rows(),
               SafeSet_BF.Ai().cols() + InputCnstrSet_.Ai().cols());
  Gxu << SafeSet_BF.Ai(),
      MatrixXd::Zero(SafeSet_BF.Ai().rows(), InputCnstrSet_.Ai().cols()),
      InputCnstrSet_.Ai();
  //   Ge = [Gc sparse(size(Gc, 1), size(Gu, 2)); Gu * alpha_e];
  MatrixXd Fxu(SafeSet_BF.Ai().rows() + InputCnstrSet_.Ai().rows(), 1);
  Fxu << SafeSet_BF.bi(), InputCnstrSet_.bi();

  HPolyhedron SafeSet_BF_extended(Gxu, Fxu);
  SafeSet_BF = SafeSet_BF_extended;
}

void CISGenerator::computeCIS(const HPolyhedron &SafeSet, int L, int T) {
  //   if (Loop_ != L || Transient_ != T)
  //     ComputeLiftedSystem(L, T);

  // Constructor has generated the Brunovsky form.
  std::vector<MatrixXd> SysMatrices_BF;
  SysMatrices_BF.push_back(brunovsky_form_->GetSystem().first);
  SysMatrices_BF.push_back(brunovsky_form_->GetSystem().second);
  if (brunovsky_form_->hasDisturbance())
    SysMatrices_BF.push_back(brunovsky_form_->GetDisturbanceMatrix());

  HPolyhedron SafeSet_BF = brunovsky_form_->GetStateConstraints(SafeSet);

  // Extend the system if there are input constraints:
  if (ThereAreInputConstraints_)
    ComputeExtended(SafeSet_BF, SysMatrices_BF);

  // Compute the shrinked Free Space
  std::vector<HPolyhedron> SafeSetSeq =
      ComputeSafeSetsSequence(SafeSet_BF, SysMatrices_BF);

  // Compute the lifted system
  NewComputeLiftedSystem(L, T, SysMatrices_BF);

  //   // // OLD CODE:
  //   int NStateConstr = SafeSet_BF.Ai().rows();
  //   int length = Loop_ + Transient_;
  //   int mu_max = brunovsky_form_->GetMaxControllabilityIndex();
  //   // Compute the mcisA; mcisB matrixes
  //   int mcisA_rows = NStateConstr * (mu_max + length);
  //   int mcisA_cols = StateDim_ + ExtendedInputDim_;
  //   MatrixXd mcisA(MatrixXd::Zero(mcisA_rows, mcisA_cols));
  //   VectorXd mcisb(VectorXd::Zero(mcisA_rows));
  // #ifdef CIS2M_DEBUG
  //   std::cout << " mcisA: " << mcisA_rows << " x " << mcisA_cols <<
  //   std::endl;
  // #endif
  //   mcisA.block(0, 0, NStateConstr, StateDim_) = SafeSet_BF.Ai();
  //   mcisb.head(NStateConstr) = SafeSet_BF.bi();
  //   MatrixXd Acurr = A_lifted_;
  //   MatrixXd TempA(MatrixXd::Zero(NStateConstr, mcisA_cols));
  //   MatrixXd mcisA_ctrl;
  //   VectorXd mcisb_ctrl;
  //   MatrixXd AA;
  //   MatrixXd BB;
  //   HPolyhedron FGu;
  //   if (ThereAreInputConstraints_) {
  //     FGu = brunovsky_form_->GetInputConstraints(InputCnstrSet_);
  //     mcisA_ctrl = MatrixXd::Zero(ExtendedInputDim_, mcisA_cols);
  //     mcisb_ctrl = VectorXd::Zero(ExtendedInputDim_);
  //     std::pair<MatrixXd, MatrixXd> AmBm =
  //         brunovsky_form_->GetIntermediateMatrixes();
  //     AA = -FGu.Ai() * AmBm.first;
  //     BB = FGu.Ai() * ExtendedU2U_;
  //     mcisA_ctrl.block(0, 0, NumberInputs_, StateDim_) = AA;
  //     mcisA_ctrl.block(0, StateDim_, NumberInputs_, mcisA_cols - StateDim_) =
  //     BB; mcisb_ctrl.head(NumberInputs_) = FGu.bi();
  //   }
  //   for (int t = 1; t < mu_max + length; t++) {
  //     int tbar = t < SafeSetSeq.size() ? t : SafeSetSeq.size() - 1;
  //     TempA.leftCols(StateDim_) = SafeSetSeq[tbar].Ai();
  //     mcisA.block(NStateConstr * t, 0, NStateConstr, mcisA_cols) = TempA *
  //     Acurr; mcisb.segment(NStateConstr * t, NStateConstr) =
  //     SafeSetSeq[tbar].bi(); if (ThereAreInputConstraints_) {
  //       if (t < length) {
  //         mcisA_ctrl.block(t * NumberInputs_, 0, NumberInputs_, StateDim_) =
  //         AA; mcisA_ctrl.block(t * NumberInputs_, StateDim_, NumberInputs_,
  //                          mcisA_cols - StateDim_) =
  //             BB * Acurr.block(StateDim_, StateDim_, ExtendedInputDim_,
  //                              ExtendedInputDim_);
  //         mcisb_ctrl.segment(t * NumberInputs_, NumberInputs_) = FGu.bi();
  //       }
  //     }
  //     Acurr *= A_lifted_;
  //   }
  //   MatrixXd mcisA_tot(mcisA.rows() + mcisA_ctrl.rows(), mcisA.cols());
  //   VectorXd mcisb_tot(mcisA_tot.rows());
  //   if (ThereAreInputConstraints_) {
  //     mcisA_tot << mcisA, mcisA_ctrl;
  //     mcisb_tot << mcisb, mcisb_ctrl;
  //   } else {
  //     mcisA_tot << mcisA;
  //     mcisb_tot << mcisb;
  //   }
  // Get the transformation
  //   MatrixXd Transform = brunovsky_form_->GetTransformationMatrix();
  //   mcisA_tot.block(0, 0, mcisA_tot.rows(), StateDim_) =
  //       mcisA_tot.block(0, 0, mcisA_tot.rows(), StateDim_) * Transform;
  //   // Note: The CIS is expressed in the original basis

  // // NEW CODE:
  MatrixXd A_curr = A_lifted_;
  int nmax = brunovsky_form_->GetMaxControllabilityIndex();
  int space_dim = A_lifted_.cols();
  int totalConstraintNum = 0;
  for (int i = 0; i < SafeSetSeq.size(); i++) {
    HPolyhedron set = SafeSetSeq[i];
    totalConstraintNum += set.GetNumInequalities();
  }
  MatrixXd mcisA_tot(totalConstraintNum, space_dim);
  VectorXd mcisb_tot(totalConstraintNum, 1);
  for (int i = 0; i < nmax + Transient_ + Loop_; i++) {
    int runningIdx = 0;
    MatrixXd tempAi(SafeSetSeq[i].GetNumInequalities(), space_dim);
    tempAi << SafeSetSeq[i].Ai(),
        MatrixXd::Zero(SafeSetSeq[i].GetNumInequalities(),
                       space_dim - SafeSetSeq[0].GetSpaceDim());
    if (i == 0) {
      mcisA_tot.block(runningIdx, 0, SafeSetSeq[0].GetNumInequalities(),
                      space_dim)
          << tempAi;
      mcisb_tot.block(runningIdx, 0, SafeSetSeq[0].GetNumInequalities(), 1)
          << SafeSetSeq[0].bi();
    } else {
      runningIdx += SafeSetSeq[i - 1].GetNumInequalities();
      mcisA_tot.block(runningIdx, 0, SafeSetSeq[i].GetNumInequalities(),
                      space_dim)
          << tempAi * A_curr;
      mcisb_tot.block(runningIdx, 0, SafeSetSeq[i].GetNumInequalities(), 1)
          << SafeSetSeq[i].bi();
    }
    A_curr *= A_lifted_;
  }

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

int CISGenerator::GetExtendedDim() const { return NumberInputs_ * Loop_; }

int CISGenerator::GetStateDim() const { return StateDim_; }

int CISGenerator::GetLevel() const { return Loop_; }

BrunovskyForm *CISGenerator::getBrunovskyForm() { return brunovsky_form_; }
} // namespace cis2m
