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
                           const MatrixXd &Bd, const bool &original_space) {
    Reset();

    Loop_ = L;
    Transient_ = T;
    StateDim_ = Ad.rows();
    InputDim_ = Bd.cols();
    original_space_ = original_space;

    GenerateBrunovksyForm(Ad, Bd);
}

CISGenerator::CISGenerator(const int L, const int T, const MatrixXd &Ad,
                           const MatrixXd &Bd, const MatrixXd &Ed,
                           const bool &original_space) {
    Reset();

    Loop_ = L;
    Transient_ = T;
    StateDim_ = Ad.rows();
    InputDim_ = Bd.cols();
    DisturbanceDim_ = Ed.cols();
    original_space_ = original_space;

    GenerateBrunovksyForm(Ad, Bd, Ed);
}

CISGenerator::~CISGenerator(){};

void CISGenerator::Reset() {
    DisturbanceDim_ = 0;

    Loop_ = -1;
    Transient_ = -1;
    VirtualInputDim_ = -1;

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

std::vector<HPolyhedron>
CISGenerator::ComputeSafeSetsSequence(const HPolyhedron &SafeSet,
                                      std::vector<MatrixXd> &SysMatrices) {
    int nmax = brunovsky_form_->GetMaxControllabilityIndex();
    int StateDim = StateDim_;
    // Now we extend the state when we have input constraints.
    if (ThereAreInputConstraints_) {
        StateDim = StateDim_ + InputDim_;
        nmax += 1;
    } else {
        StateDim = StateDim_;
    }

    std::vector<HPolyhedron> SafeSet_seq;
    SafeSet_seq.push_back(SafeSet);
    if (brunovsky_form_->hasDisturbance()) {
        MatrixXd A_curr = MatrixXd::Identity(StateDim, StateDim);
        for (int i = 1; i < nmax + Transient_ + Loop_; i++) {
            if (i < nmax) {
                HPolyhedron Sub(
                    DisturbanceSet_.affineT(A_curr * SysMatrices[2]));
                SafeSet_seq.push_back(SafeSet_seq.back() - Sub);
                A_curr *= SysMatrices[0];

            } else {
                SafeSet_seq.push_back(SafeSet_seq.back());
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

void CISGenerator::ComputeLiftedSystem(std::vector<MatrixXd> &SysMatrices) {
    int Level = Transient_ + Loop_;
    VirtualInputDim_ = Level * InputDim_;

    int StateDim = StateDim_;
    // Now we extend the state when we have input constraints.
    if (ThereAreInputConstraints_) {
        StateDim = StateDim_ + InputDim_;
    } else {
        StateDim = StateDim_;
    }

    // Construct the High-Dimensional system:
    // A_lifted_ =
    //      [A    B * K]
    //      [0      P  ]
    MatrixXd Ki(MatrixXd::Zero(1, Level));
    Ki(0) = 1.0;
    MatrixXd Pi(MatrixXd::Zero(Level, Level));
    Pi.block(0, 1, Level - 1, Level - 1) =
        MatrixXd::Identity(Level - 1, Level - 1);
    Pi(Level - 1, Transient_) = 1.0;
    ExtendedU2U_ = blkdiag(Ki, InputDim_);
    MatrixXd P = blkdiag(Pi, InputDim_);

    int Nrow_hd = StateDim + VirtualInputDim_;
    int Ncol_hd = StateDim + P.cols();
    MatrixXd Ahd(MatrixXd::Zero(Nrow_hd, Ncol_hd));
    Ahd.block(0, 0, StateDim, StateDim) = SysMatrices[0];
    Ahd.block(0, StateDim, StateDim, ExtendedU2U_.cols()) =
        SysMatrices[1] * ExtendedU2U_;
    Ahd.block(StateDim, StateDim, P.rows(), P.cols()) = P;

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

    if (brunovsky_form_->hasDisturbance()) {
        MatrixXd E_BF_extended(SysMatrices_BF[2].rows() +
                                   SysMatrices_BF[1].cols(),
                               SysMatrices_BF[2].cols());
        E_BF_extended << SysMatrices_BF[2],
            MatrixXd::Zero(SysMatrices_BF[1].cols(), SysMatrices_BF[2].cols());

        SysMatrices_BF[2] = E_BF_extended;
    }
    SysMatrices_BF[1] = B_BF_extended;
    SysMatrices_BF[0] = A_BF_extended;

    // Compute the extended safe set
    MatrixXd Gxu(SafeSet_BF.Ai().rows() + InputCnstrSet_.Ai().rows(),
                 SafeSet_BF.Ai().cols() + InputCnstrSet_.Ai().cols());

    std::pair<MatrixXd, MatrixXd> AmBm =
        brunovsky_form_->GetIntermediateMatrices();
    MatrixXd FeedbackConstraints(InputCnstrSet_.Ai().cols(),
                                 AmBm.first.cols() + AmBm.second.cols());
    FeedbackConstraints << -AmBm.second.inverse() * AmBm.first, AmBm.second;

    Gxu << SafeSet_BF.Ai(),
        MatrixXd::Zero(SafeSet_BF.Ai().rows(), InputCnstrSet_.Ai().cols()),
        InputCnstrSet_.Ai() * FeedbackConstraints;
    MatrixXd Fxu(SafeSet_BF.Ai().rows() + InputCnstrSet_.Ai().rows(), 1);
    Fxu << SafeSet_BF.bi(), InputCnstrSet_.bi();

    HPolyhedron SafeSet_BF_extended(Gxu, Fxu);
    SafeSet_BF = SafeSet_BF_extended;
}

void CISGenerator::TransformToOriginal() {
    MatrixXd Ai_OG(CIS_.GetNumInequalities(), CIS_.GetSpaceDim());
    if (ThereAreInputConstraints_) {
        // Fetch matrices:
        MatrixXd T = brunovsky_form_->GetTransformationMatrix();
        std::pair<MatrixXd, MatrixXd> AmBm =
            brunovsky_form_->GetIntermediateMatrices();
        // Transform CIS_ to original space:
        Ai_OG << Fetch_A_State() * T + Fetch_A_Input() * AmBm.first * T,
            Fetch_A_Input() * AmBm.second, Fetch_A_Virtual();

        // Construct A_lifted_ in original space:
        // A_lifted_ =
        //      [A                     B 0_{StateDim_,VirtualInputDim_}]
        //      [-inv(Bm)*Am*T*A  -inv(Bm)*Am*T*B           inv(Bm)*K ]
        //      [0_{VirtualInputDim_,StateDim_+InputDim_}   P ]

        int Level = Transient_ + Loop_;
        // Construct K (ExtendedU2U_):
        MatrixXd Ki(MatrixXd::Zero(1, Level));
        Ki(0) = 1.0;
        ExtendedU2U_ = blkdiag(Ki, InputDim_);
        // Construct P:
        MatrixXd Pi(MatrixXd::Zero(Level, Level));
        Pi.block(0, 1, Level - 1, Level - 1) =
            MatrixXd::Identity(Level - 1, Level - 1);
        Pi(Level - 1, Transient_) = 1.0;
        MatrixXd P = blkdiag(Pi, InputDim_);

        int LiftedDim = StateDim_ + InputDim_ + VirtualInputDim_;
        MatrixXd Ahd(MatrixXd::Zero(LiftedDim, LiftedDim));
        std::pair<MatrixXd, MatrixXd> AB = brunovsky_form_->GetOriginalSystem();
        // x+ = Ax + Bu:
        Ahd.block(0, 0, StateDim_, StateDim_) = AB.first;
        Ahd.block(0, StateDim_, StateDim_, InputDim_) = AB.second;
        // u+ = -inv(Bm)*Am*T*A x -inv(Bm)*Am*T*B u + inv(Bm)*K {virtual}
        Ahd.block(StateDim_, 0, InputDim_, StateDim_) =
            -AmBm.second.inverse() * AmBm.first * T * AB.first;
        Ahd.block(StateDim_, StateDim_, InputDim_, InputDim_) =
            -AmBm.second.inverse() * AmBm.first * T * AB.second;
        Ahd.block(StateDim_, StateDim_ + InputDim_, InputDim_,
                  VirtualInputDim_) = AmBm.second.inverse() * ExtendedU2U_;
        // {virtual}+ = P {virtual}
        Ahd.block(StateDim_ + InputDim_, StateDim_ + InputDim_, P.rows(),
                  P.cols()) = P;

        A_lifted_ = Ahd;
    } else {
        // Transform CIS_ to original space:
        Ai_OG << Fetch_A_State() * brunovsky_form_->GetTransformationMatrix(),
            Fetch_A_Input(), Fetch_A_Virtual();

        // Construct A_lifted_ in original space:
        // A_lifted_ =
        //      [A    B * K]
        //      [0      P  ]
        int Level = Transient_ + Loop_;
        MatrixXd Ki(MatrixXd::Zero(1, Level));
        Ki(0) = 1.0;
        MatrixXd Pi(MatrixXd::Zero(Level, Level));
        Pi.block(0, 1, Level - 1, Level - 1) =
            MatrixXd::Identity(Level - 1, Level - 1);
        Pi(Level - 1, Transient_) = 1.0;
        ExtendedU2U_ = blkdiag(Ki, InputDim_);
        MatrixXd P = blkdiag(Pi, InputDim_);

        int Nrow_hd = StateDim_ + VirtualInputDim_;
        int Ncol_hd = StateDim_ + P.cols();
        MatrixXd Ahd(MatrixXd::Zero(Nrow_hd, Ncol_hd));
        Ahd.block(0, 0, StateDim_, StateDim_) =
            brunovsky_form_->GetOriginalSystem().first;
        Ahd.block(0, StateDim_, StateDim_, ExtendedU2U_.cols()) =
            brunovsky_form_->GetOriginalSystem().second * ExtendedU2U_;
        Ahd.block(StateDim_, StateDim_, P.rows(), P.cols()) = P;

        A_lifted_ = Ahd;
    }

    CIS_ = HPolyhedron(Ai_OG, CIS_.bi());
}

void CISGenerator::computeImplicitCIS(const HPolyhedron &SafeSet, int L,
                                      int T) {

    // Constructor has generated the Brunovsky Form (BF).
    // Get the system matrices in BF:
    std::vector<MatrixXd> SysMatrices_BF;
    SysMatrices_BF.push_back(brunovsky_form_->GetBrunovskySystem().first);
    SysMatrices_BF.push_back(brunovsky_form_->GetBrunovskySystem().second);
    if (brunovsky_form_->hasDisturbance())
        SysMatrices_BF.push_back(
            brunovsky_form_->GetDisturbanceMatrixBrunovsky());
    // and the safe set in BF:
    HPolyhedron SafeSet_BF = brunovsky_form_->GetStateConstraints(SafeSet);

    // Extend the system if there are input constraints:
    if (ThereAreInputConstraints_)
        ComputeExtended(SafeSet_BF, SysMatrices_BF);

    // Compute the safe set sequence:
    std::vector<HPolyhedron> SafeSetSeq =
        ComputeSafeSetsSequence(SafeSet_BF, SysMatrices_BF);

    // Compute the lifted system:
    ComputeLiftedSystem(SysMatrices_BF);

    // // Below Luigi was considering all sets in SafeSetSeq to have the same
    // NumIneqs, but with my new MinkDiff it's not the case.
    // // So I changed the below.
    // //  // New MinkDiff
    // int mcisA_rows = 0;
    // for (auto set : SafeSetSeq){
    //     mcisA_rows += set.GetNumInequalities();
    // }
    // int mcisA_cols = SysMatrices_BF[0].cols() + VirtualInputDim_;

    int NStateConstr =
        SafeSet_BF.Ai()
            .rows(); // Assumes that all sets in sequence have the same ineq
                     // num. || Will change it in the loop later.
    int Level = Loop_ + Transient_;
    int nmax = brunovsky_form_->GetMaxControllabilityIndex();
    if (ThereAreInputConstraints_) {
        nmax += 1;
    }
    // Compute the mcisA; mcisB matrixes
    int mcisA_rows = NStateConstr * (nmax + Level);
    int mcisA_cols = SysMatrices_BF[0].cols() + VirtualInputDim_;
    MatrixXd mcisA(MatrixXd::Zero(mcisA_rows, mcisA_cols));
    VectorXd mcisb(VectorXd::Zero(mcisA_rows));
#ifdef CIS2M_DEBUG
    std::cout << " mcisA: " << mcisA_rows << " x " << mcisA_cols << std::endl;
#endif
    mcisA.block(0, 0, NStateConstr, SysMatrices_BF[0].cols()) = SafeSet_BF.Ai();
    mcisb.head(NStateConstr) = SafeSet_BF.bi();
    MatrixXd Acurr = A_lifted_;
    MatrixXd TempA(MatrixXd::Zero(NStateConstr, mcisA_cols));
    //   MatrixXd mcisA_ctrl;
    //   VectorXd mcisb_ctrl;
    MatrixXd AA;
    MatrixXd BB;
    HPolyhedron FGu;

    //  // New MinkDiff
    //   for (int t = 1; t < nmax + Level; t++) {
    //     MatrixXd TempAi(MatrixXd::Zero(SafeSetSeq[t].Ai().rows(),
    //     mcisA_cols)); TempAi.leftCols(SafeSetSeq[t].Ai().cols()) =
    //     SafeSetSeq[t].Ai(); mcisA.block(NStateConstr, 0,
    //     SafeSetSeq[t].Ai().rows(), mcisA_cols) = TempAi * Acurr;
    //     mcisb.segment(NStateConstr, SafeSetSeq[t].Ai().rows()) =
    //     SafeSetSeq[t].bi(); Acurr *= A_lifted_; NStateConstr +=
    //     SafeSetSeq[t].Ai().rows();
    //   }

    //   if (ThereAreInputConstraints_) {
    //     FGu = brunovsky_form_->GetInputConstraints(InputCnstrSet_);
    //     mcisA_ctrl = MatrixXd::Zero(VirtualInputDim_, mcisA_cols);
    //     mcisb_ctrl = VectorXd::Zero(VirtualInputDim_);
    //     std::pair<MatrixXd, MatrixXd> AmBm =
    //     brunovsky_form_->GetIntermediateMatrices(); AA = -FGu.Ai() *
    //     AmBm.first; BB = FGu.Ai() * ExtendedU2U_; mcisA_ctrl.block(0, 0,
    //     InputDim_, StateDim_) = AA; mcisA_ctrl.block(0, StateDim_, InputDim_,
    //     mcisA_cols - StateDim_) = BB; mcisb_ctrl.head(InputDim_) = FGu.bi();
    //   }

    for (int t = 1; t < nmax + Level; t++) {
        int tbar = t < SafeSetSeq.size() ? t : SafeSetSeq.size() - 1;
        TempA.leftCols(SysMatrices_BF[0].cols()) = SafeSetSeq[tbar].Ai();
        mcisA.block(NStateConstr * t, 0, NStateConstr, mcisA_cols) =
            TempA * Acurr;
        mcisb.segment(NStateConstr * t, NStateConstr) = SafeSetSeq[tbar].bi();
        // if (ThereAreInputConstraints_) {
        //   if (t < Level) {
        //     mcisA_ctrl.block(t * InputDim_, 0, InputDim_, StateDim_) =
        //     AA; mcisA_ctrl.block(t * InputDim_, StateDim_, InputDim_,
        //                      mcisA_cols - StateDim_) =
        //         BB * Acurr.block(StateDim_, StateDim_, VirtualInputDim_,
        //                          VirtualInputDim_);
        //     mcisb_ctrl.segment(t * InputDim_, InputDim_) = FGu.bi();
        //   }
        // }
        Acurr *= A_lifted_;
    }
    //   MatrixXd mcisA_tot(mcisA.rows() + mcisA_ctrl.rows(), mcisA.cols());
    //   VectorXd mcisb_tot(mcisA_tot.rows());
    //   if (ThereAreInputConstraints_) {
    //     mcisA_tot << mcisA, mcisA_ctrl;
    //     mcisb_tot << mcisb, mcisb_ctrl;
    //   } else {
    //     mcisA_tot << mcisA;
    //     mcisb_tot << mcisb;
    //   }

    // // Get the transformation
    // //   MatrixXd Transform = brunovsky_form_->GetTransformationMatrix();
    // //   mcisA_tot.block(0, 0, mcisA_tot.rows(), StateDim_) =
    // //       mcisA_tot.block(0, 0, mcisA_tot.rows(), StateDim_) * Transform;
    //   // Note: The CIS is expressed in the original basis

    //   CIS_ = HPolyhedron(mcisA_tot, mcisb_tot);
    CIS_ = HPolyhedron(mcisA, mcisb);
    cis_computed_ = true;

    // Transform CIS_ and A_lifted_ from BF to original space:
    if (original_space_)
        TransformToOriginal();

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
        return CIS_.Ai().block(0, StateDim_, NumRows, InputDim_);
    } else {
        return {};
    }
}

MatrixXd CISGenerator::Fetch_A_Virtual() {
    if (cis_computed_) {
        // int NumCols = CIS_.Ai().cols();
        // return CIS_.Ai().rightCols(NumCols - (InputDim_ + StateDim_));
        // return CIS_.Ai().rightCols(NumCols - VirtualInputDim_);
        return CIS_.Ai().rightCols(VirtualInputDim_);
    } else {
        return {};
    }
}

MatrixXd CISGenerator::Fetch_A_lifted() { return A_lifted_; }

VectorXd CISGenerator::TransformU2B(const VectorXd &u, const VectorXd &x) {
    // The transformation from original coordinates (x,u) to Brunovsky
    // coordinates (z,v) is: x = inv(T) * z, u = -inv(Bm) * Am * T * x + inv(Bm)
    // * v then: v = Am * T * x + Bm * u.
    MatrixXd T = brunovsky_form_->GetTransformationMatrix();
    std::pair<MatrixXd, MatrixXd> AmBm =
        brunovsky_form_->GetIntermediateMatrices();
    return AmBm.second * u + AmBm.first * T * x;
}

VectorXd CISGenerator::TransformU2O(const VectorXd &u, const VectorXd &x) {
    // The transformation from original coordinates (x,u) to Brunovsky
    // coordinates (z,v) is: x = inv(T) * z, u = -inv(Bm) * Am * T * x + inv(Bm)
    // * v
    MatrixXd T = brunovsky_form_->GetTransformationMatrix();
    std::pair<MatrixXd, MatrixXd> AmBm =
        brunovsky_form_->GetIntermediateMatrices();
    return AmBm.second.inverse() * (u - AmBm.first * T * x);
}

int CISGenerator::GetExtendedDim() const { return VirtualInputDim_; }

int CISGenerator::GetStateDim() const { return StateDim_; }

int CISGenerator::GetInputDim() const { return InputDim_; }

int CISGenerator::GetLevel() const { return Loop_; }

HPolyhedron CISGenerator::GetDisturbanceSet() { return DisturbanceSet_; }

HPolyhedron CISGenerator::GetInputConstraints() { return InputCnstrSet_; }

bool CISGenerator::isInOriginalSpace() { return original_space_; }

BrunovskyForm *CISGenerator::getBrunovskyForm() { return brunovsky_form_; }
} // namespace cis2m
