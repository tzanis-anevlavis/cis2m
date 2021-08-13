#include "hpolyhedron.hpp"
#include <ortools/linear_solver/linear_solver.h>

using operations_research::MPSolver;
using operations_research::MPConstraint;
using operations_research::MPVariable;
using operations_research::MPObjective;

using Eigen::VectorXd;
using Eigen::MatrixXd;

namespace cis2m {

	// CLASS ===================================================================
	HPolyhedron::HPolyhedron() {
		SpaceDim_ = 0;
		NumIneqs_ = 0;
	}

	HPolyhedron::~HPolyhedron() {}

	HPolyhedron::HPolyhedron(const Eigen::MatrixXd& A, const Eigen::MatrixXd& b) : Ai_(A), bi_(b) {
		SpaceDim_ = Ai_.cols();
		NumIneqs_ = Ai_.rows();
	}

	HPolyhedron::HPolyhedron(
			const Eigen::MatrixXd& Ai,
			const Eigen::MatrixXd& bi,
			const Eigen::MatrixXd& Ae,
			const Eigen::MatrixXd& be) : Ai_(Ai), bi_(bi), Ae_(Ae), be_(be)  {
		SpaceDim_ = Ai_.cols();
		NumIneqs_ = Ai_.rows();
	}

	// Support
	VectorXd HPolyhedron::ComputeSupport(const MatrixXd& A_other) const {
		// Number of hyperplanes
		int NHyperPlanes = A_other.rows(); 

		VectorXd supp(NHyperPlanes);

		std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("SCIP"));
		if (!solver) {
			std::cerr << "SCIP Solver not available" << std::endl;
			return supp;
		}
		const double infinity = solver->infinity();

		for (int hp = 0; hp < NHyperPlanes; hp++) {
			solver->Clear();

			// Create the X variable
			int NumVars = SpaceDim_;
			std::vector<MPVariable*> x(NumVars);
			for (int i = 0; i < NumVars; i++) {
				x[i] = solver->MakeNumVar(-infinity, infinity, "");
			}

			// Create the inequalities constraints
			for (int i = 0; i < NumIneqs_; i++) {
				MPConstraint* con = solver->MakeRowConstraint(-infinity, bi_(i), "");
				for (int j = 0; j < NumVars; j++) {
					con->SetCoefficient(x[j],  Ai_(i, j));
				}
			}

			// Define the objective
			MPObjective* const obj = solver->MutableObjective();
			obj->SetMaximization();

			VectorXd cost_coeff = VectorXd(A_other.row(hp));
			for (int j = 0; j < NumVars; j++) {
				obj->SetCoefficient(x[j], cost_coeff(j));
			}

			const MPSolver::ResultStatus result_status = solver->Solve();

			if (result_status != MPSolver::OPTIMAL) {
				std::cerr << "Did not get the optimal solution!" << std::endl;
				std::cerr << "Return Result: " << result_status << std::endl;
			}
			supp(hp) = obj->Value();
		}
		return supp;
	}

	HPolyhedron HPolyhedron::Projection(int el_index, double tol) {
		MatrixXd el_col = Ai_.col(el_index);

		std::vector<int> posIndex;
		std::vector<int> negIndex;
		std::vector<int> zeroIndex;

		// Normalize the rows using the elements on the column
		// we want to eliminate.
		for (int i = 0; i < el_col.rows(); i++) {
			double val = el_col(i);

			if (val < -tol) {
				negIndex.push_back(i);
				Ai_.row(i) = Ai_.row(i) / (-val);
				bi_.row(i) = bi_.row(i) / (-val);
			} else {

				if (val > tol) {
					posIndex.push_back(i);
					Ai_.row(i) = Ai_.row(i) / val;
					bi_.row(i) = bi_.row(i) / val;
				} else {
					zeroIndex.push_back(i);
				}
			}
		}

		// Create a new system composing the rows of the old one
		int num_rows = zeroIndex.size() + (posIndex.size() * negIndex.size());
		int num_cols = Ai_.cols() - 1;
		MatrixXd A_new(num_rows, num_cols);
		MatrixXd B_new(num_rows, 1);

		int left_n_elements = el_index;
		int right_n_elements = Ai_.cols() - el_index - 1;

		int cum_counter = 0;
		for (auto elp : posIndex) {
			for (auto eln : negIndex) {
				A_new.block(cum_counter, 0, 1, left_n_elements) = Ai_.block(elp, 0, 1, left_n_elements) + Ai_.block(eln, 0, 1, left_n_elements);
				A_new.block(cum_counter, el_index, 1, right_n_elements) = Ai_.block(elp, el_index + 1, 1, right_n_elements) + Ai_.block(eln, el_index + 1, 1, right_n_elements);
				B_new.block(cum_counter++, 0, 1, 1) = bi_.block(elp, 0, 1, 1) + bi_.block(eln, 0, 1, 1);
			}
		}

		for (int i : zeroIndex) {
			A_new.block(cum_counter, 0, 1, left_n_elements) = Ai_.block(i, 0, 1, left_n_elements);
			A_new.block(cum_counter, el_index, 1, right_n_elements) = Ai_.block(i, el_index + 1, 1, right_n_elements);
			B_new.block(cum_counter++, 0, 1, 1) = bi_.block(i, 0, 1, 1);
		}

		HPolyhedron output(A_new, B_new);
		return output;
	}


	HPolyhedron HPolyhedron::affineT(const MatrixXd& T) {
		int Nc = T.cols();
		int Nr = T.rows();

		// Full LU decomposition
		Eigen::FullPivLU<MatrixXd> lu_decomposition(T);
		int rank = lu_decomposition.rank();

		MatrixXd LLUU  = lu_decomposition.matrixLU();
		int NLLUU_r = LLUU.rows();
		int NLLUU_c = LLUU.cols();

		// Extract L
		MatrixXd ll(MatrixXd::Identity(NLLUU_r, NLLUU_c));
		ll.block(0, 0, NLLUU_r, NLLUU_c).triangularView<Eigen::StrictlyLower>() = LLUU;
		MatrixXd l(ll.block(0, 0, Nr, rank));

		// Extract U
		MatrixXd uu = lu_decomposition.matrixLU().triangularView<Eigen::Upper>();
		MatrixXd u(uu.block(0,  0, rank, Nc));

		// Compose LU
		MatrixXd LU = l * u;

		MatrixXd LU_11 = LU.block(0,0, rank, rank);
		MatrixXd LU_12 = LU.block(0, rank, rank, Nc - rank); 
		MatrixXd LU_21 = LU.block(rank, 0, Nr - rank, rank);

		// Build the generic inverse for the system
		MatrixXd LU_11_inv = LU_11.inverse();
		MatrixXd invT = MatrixXd::Identity(Nc, Nc);
		invT.block(0,0, rank, rank) = LU_11_inv;
		invT.block(0, rank, rank, Nc - rank) = -LU_11_inv * LU_12; 

		// Inequality matrix in the permutated coordinates
		MatrixXd Q = lu_decomposition.permutationQ();
		MatrixXd A_prime = Ai_ * Q * invT; 
		HPolyhedron P_prime = HPolyhedron(A_prime, bi_);

		// Project the polyhedron on the first 'rank' coordinates
		for (int i = 0; i < Nc - rank; i++) {
			int deleted_col = P_prime.Ai().cols() - 1;
			P_prime = P_prime.Projection(deleted_col);
		}

		// Build the full system in the new permuted coordinates
		MatrixXd Ai_final(MatrixXd::Zero(P_prime.Ai().rows(), Nr));
		Ai_final.block(0,0, P_prime.Ai().rows(), rank) = P_prime.Ai();
		MatrixXd Ae_final(MatrixXd::Zero(Nr - rank, Nr));
		Ae_final.block(0, 0, Nr - rank, rank) = LU_21 * LU_11_inv;
		Ae_final.block(0, rank, Nr - rank, Nr - rank) = -MatrixXd::Identity(Nr - rank, Nr - rank);

		// Remove the permutation
		MatrixXd P = lu_decomposition.permutationP();
		Ai_final = Ai_final * P;
		Ae_final = Ae_final * P;
		MatrixXd bi_final = P_prime.bi();
		MatrixXd be_final = MatrixXd::Zero(Nr - rank, 1);

		HPolyhedron P_final(Ai_final, bi_final, Ae_final, be_final);

#ifdef CIS2M_DEBUG 
		std::cout << __FILE__ << std::endl;
		std::cout << "LU Decomposition" << std::endl;
		std::cout << "Original: " << std::endl << T << std::endl;
		std::cout << "P: " << std::endl << P << std::endl;
		std::cout << "Q: " << std::endl << Q << std::endl;
		std::cout << "LLUU: " << std::endl << LLUU << std::endl;
		std::cout << "L: " << std::endl << l << std::endl;
		std::cout << "U: " << std::endl << u << std::endl;
		std::cout << "LU: " << std::endl << LU << std::endl;
		std::cout << "invT: " << std::endl << invT << std::endl;
		std::cout << "A_prime: " << std::endl << A_prime << std::endl;
		std::cout << "A_prime_proj: " << std::endl << P_prime.Ai() << std::endl; 
		std::cout << "A_final: " << std::endl << Ai_final << std::endl;
		std::cout << "Reconstruction: " << std::endl << P.inverse() * l * u * Q.inverse() << std::endl;
		std::cout << std::endl;
#endif

		return P_final;
	}

	// Operators
	HPolyhedron& HPolyhedron::operator=(const HPolyhedron& other) {
		Ai_ = other.Ai();
		bi_ = other.bi();

		Ae_ = other.Ae();
		be_ = other.be();

		return *this;
	}


	HPolyhedron& HPolyhedron::operator+=(const HPolyhedron& rhs) {
		// Substract the rhs.Support(this.A) from the this.b vector
		// max <a_i, x>
		// subj. to x \in rhs
		// where a_i is the vector representing a hyperplane of *this polyhedron.
		VectorXd supp = rhs.ComputeSupport(Ai_);
		bi_ = bi_ + supp;

		return *this;
	}


	HPolyhedron& HPolyhedron::operator-=(const HPolyhedron& rhs) {
		// Substract the rhs.Support(this.A) from the this.b vector
		// max <a_i, x>
		// subj. to x \in rhs
		// where a_i is the vector representing a hyperplane of *this polyhedron.
		VectorXd supp = rhs.ComputeSupport(Ai_);
		bi_ = bi_ - supp;

		return *this;
	}


	// Getters
	Eigen::MatrixXd HPolyhedron::Ai() const {
		return Ai_;
	}

	Eigen::VectorXd HPolyhedron::bi() const {
		return bi_;
	}

	Eigen::MatrixXd HPolyhedron::Ae() const {
		return Ae_;
	}

	Eigen::VectorXd HPolyhedron::be() const {
		return be_;
	}



	// Operators
	HPolyhedron operator+(HPolyhedron Pl, HPolyhedron Pr) {
		return Pl += Pr;
	}

	HPolyhedron operator-(HPolyhedron Pl, HPolyhedron Pr) {
		return Pl -= Pr;
	}
}
