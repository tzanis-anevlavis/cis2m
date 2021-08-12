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
	HPolyhedron::~HPolyhedron() {}

	HPolyhedron::HPolyhedron(const Eigen::MatrixXd& A, const Eigen::MatrixXd& b) : Ai_(A), bi_(b) {
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
			}
			supp(hp) = obj->Value();
		}
		return supp;
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



	HPolyhedron operator+(HPolyhedron Pl, HPolyhedron Pr) {
		return Pl += Pr;
	}

	HPolyhedron operator-(HPolyhedron Pl, HPolyhedron Pr) {
		return Pl -= Pr;
	}
}
