#include "cis_generator.hpp"
#include <iostream>


namespace cis2m {
	// Helper Functions
	MatrixXd blkdiag(const MatrixXd& A, int count) {
		MatrixXd bdm = MatrixXd::Zero(A.rows() * count, A.cols() * count);
		for (int i = 0; i < count; ++i) {
			bdm.block(i * A.rows(), i * A.cols(), A.rows(), A.cols()) = A;
		}
		return bdm;
	}


	// ========================================================================
	// CLASS
	CISGenerator::CISGenerator(const MatrixXd& Ad,const MatrixXd& Bd) {
		StateDim_ = Ad.rows();
		NumberInputs_ = Bd.cols();
		cis_computed_ = false;

		brunovsky_form_ = new BrunovskyForm(Ad, Bd);  
	}

	CISGenerator::~CISGenerator() {};

	HPolyhedron CISGenerator::computeCIS(const HPolyhedron& SafeSet, int L, int T) {
		// Get the A, B pair in the Brunovsky form
		std::pair<MatrixXd, MatrixXd> pairAB = brunovsky_form_->GetDynSystem();

		int length = T + L - 1;
		// Construct the High-Dimensional system
		MatrixXd Ki(MatrixXd::Zero(1, length));
		Ki(0) = 1.0;
		MatrixXd Pi(MatrixXd::Zero(T + L, T + L));
		Pi.block(0, 1, length, length) = MatrixXd::Identity(length, length);
		MatrixXd K = blkdiag(Ki, NumberInputs_);
		MatrixXd P = blkdiag(Pi, NumberInputs_); 

		int Nrow_hd = StateDim_ + (T + L)  * NumberInputs_;
		int Ncol_hd = StateDim_ + P.cols();
		MatrixXd Ahd(MatrixXd::Zero(Nrow_hd, Ncol_hd));

		Ahd.block(0, 0, StateDim_, StateDim_) = pairAB.first;
		Ahd.block(0, StateDim_, StateDim_, K.cols()) = pairAB.second * K;
		Ahd.block(StateDim_, StateDim_, P.rows(), P.cols()) = P;

		HPolyhedron polyhedron_Gb = brunovsky_form_->GetDynConstraints(SafeSet);
		int NDynconstr = polyhedron_Gb.Ai().rows();

		int mu_max = brunovsky_form_->GetMaxControllabilityIndex();
		int mcisA_rows =  NDynconstr * (mu_max + T + L); 
		int mcisA_cols = StateDim_ + NumberInputs_ * (T + L);
		MatrixXd mcisA(MatrixXd::Zero(mcisA_rows, mcisA_cols));
		MatrixXd mcisb(MatrixXd::Zero(mcisA_rows, 1)); 

		mcisA.block(0, 0, NDynconstr, StateDim_) = polyhedron_Gb.Ai(); 
		mcisb.block(0, 0, NDynconstr, 1) = polyhedron_Gb.bi();

		for (int t = 1; t < mu_max + (T + L - 1); t++) {
			mcisA.block(NDynconstr * t, 0, NDynconstr, mcisA_cols) =
				mcisA.block(NDynconstr * (t - 1), 0, NDynconstr, mcisA_cols) * Ahd; 

			mcisb.block(NDynconstr * t, 0, NDynconstr, 1) =
				polyhedron_Gb.bi();
		}

		// Get the transformation
		MatrixXd Transform = brunovsky_form_->GetTransformationMatrix();

		std::cout << "Transf: "  << std::endl << Transform << std::endl;

		mcisA.block(0, 0, mcisA.rows(), StateDim_) = mcisA.block(0, 0, mcisA.rows(), StateDim_) * 
			Transform;

		// Note: The CIS is expressed in the original basis
		CIS_ = new HPolyhedron(mcisA, mcisb);
		cis_computed_ = true;
		return *CIS_;
	}

	std::optional<HPolyhedron> CISGenerator::FetchCIS() {
		if (cis_computed_) {
			return *CIS_;
		}
		else{
			return {};	
		}
	}

}
