#include "hpolyhedron.hpp"
#include <cfloat>
#include <iostream>
#include <ortools/linear_solver/linear_solver.h>
#include <string>

using operations_research::MPConstraint;
using operations_research::MPObjective;
using operations_research::MPSolver;
using operations_research::MPVariable;

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace cis2m {

struct optimization_result {
    bool solved;
    float objective_value;
    VectorXd solution;
    // Constructor:
    optimization_result(bool b = false, float val = 0.0,
                        VectorXd vec = VectorXd::Zero(0))
        : solved(b), objective_value(val), solution(vec){};
};

// Helper function: Solve a linear program.
optimization_result solveLP(const VectorXd &c, const MatrixXd &A,
                            const VectorXd &b, const bool &maximize = false,
                            const std::string &solver_name = "GUROBI") {

    optimization_result res;

    // Select solver:
    std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver(solver_name));
    if (!solver) {
        LOG(WARNING) << "Solver unavailable.";
        return res;
    }
    solver->Clear();
    // Set silent:
    solver->SuppressOutput();

    // Check correctness of dimensions:
    if (A.rows() != b.rows()) {
        std::cerr
            << "Matrix 'A' and vector 'b' must have the same number of rows."
            << std::endl;
        return {};
    }
    if (A.cols() != c.rows()) {
        std::cerr
            << "Matrix 'A' must have the same number of columns as the size "
               "of vector 'c'."
            << std::endl;
        return {};
    }

    // Create the program variable:
    const double infinity = solver->infinity();
    int NumVars = A.cols();
    std::vector<MPVariable *> x(NumVars);
    for (int i = 0; i < NumVars; i++) {
        x[i] = solver->MakeNumVar(-infinity, infinity, "");
    }

    // Create the inequalities constraints
    int NumIneqs = A.rows();
    for (int i = 0; i < NumIneqs; i++) {
        MPConstraint *con = solver->MakeRowConstraint(-infinity, b(i));
        for (int j = 0; j < NumVars; j++) {
            con->SetCoefficient(x[j], A(i, j));
        }
    }

    // Define the objective
    MPObjective *const obj = solver->MutableObjective();
    for (int j = 0; j < NumVars; j++) {
        obj->SetCoefficient(x[j], c(j));
    }
    if (maximize) {
        obj->SetMaximization();
    } else {
        obj->SetMinimization();
    }

    const MPSolver::ResultStatus result_status = solver->Solve();
    if (result_status == MPSolver::INFEASIBLE) {
        res.solved = false;
        res.objective_value = DBL_MAX;
        res.solution = {};
    } else if (result_status == MPSolver::OPTIMAL) {
        res.solved = true;
        res.objective_value = obj->Value();
        res.solution = VectorXd::Zero(NumVars);
        for (int j = 0; j < NumVars; j++) {
            res.solution(j) = x[j]->solution_value();
        }
    } else {
        // Catch if we get something unbounded or abnormal.
        LOG(INFO) << "Solution status: " << result_status;
        LOG(INFO) << "Optimal objective value = " << obj->Value();
        LOG(INFO) << "Solution:";
        for (int j = 0; j < NumVars; j++) {
            LOG(INFO) << x[j]->solution_value();
        }
    }
    return res;
}

// CLASS ===================================================================
HPolyhedron::HPolyhedron() {
    SpaceDim_ = 0;
    NumIneqs_ = 0;
    NumEqs_ = 0;
    initialized_ = false;
}

HPolyhedron::~HPolyhedron() {}

HPolyhedron::HPolyhedron(const Eigen::MatrixXd &A, const Eigen::VectorXd &b)
    : Ai_(A), bi_(b) {
    // Run checks
    if (Ai_.rows() != bi_.rows()) {
        std::cerr
            << "Matrices 'Ai_' and 'bi_' must have the same number of rows."
            << std::endl;
    }
    SpaceDim_ = Ai_.cols();
    NumIneqs_ = Ai_.rows();
    NumEqs_ = 0;
    initialized_ = true;
}

HPolyhedron::HPolyhedron(const Eigen::MatrixXd &Ai, const Eigen::VectorXd &bi,
                         const Eigen::MatrixXd &Ae, const Eigen::VectorXd &be)
    : Ai_(Ai), bi_(bi), Ae_(Ae), be_(be) {
    // Run checks
    if (Ai_.rows() != bi_.rows()) {
        std::cerr
            << "Matrices 'Ai_' and 'bi_' must have the same number of rows."
            << std::endl;
    }
    if (Ae_.rows() != be_.rows()) {
        std::cerr
            << "Matrices 'Ae_' and 'be_' must have the same number of rows."
            << std::endl;
    }
    if (Ai_.cols() != Ae_.cols()) {
        std::cerr
            << "Matrices 'Ai_' and 'Ae_' must have the same number of columns."
            << std::endl;
    }
    SpaceDim_ = Ai_.cols();
    NumIneqs_ = Ai_.rows();
    NumEqs_ = Ae_.rows();
    initialized_ = true;
}

// Check emptiness of HPolyhedron.
bool HPolyhedron::isEmpty() const {
    // Argument checks
    if (!initialized_) {
        std::cerr << "Calling a method on an uninitialized HPolyhedron!"
                  << std::endl;
        return true;
    }
    return !solveLP(VectorXd::Zero(SpaceDim_), Ai_, bi_).solved;
}

// Check if HPolyhedron contains another HPolyhedron.
bool HPolyhedron::Contains(const HPolyhedron &P) const {
    // Argument checks
    if (!initialized_) {
        std::cerr << "Calling a method on an uninitialized HPolyhedron!"
                  << std::endl;
        return false;
    }
    if (SpaceDim_ != P.SpaceDim_) {
        std::cerr << "The HPolyhedron objects have different dimensions!"
                  << std::endl;
        return false;
    }
    if (P.isEmpty()) {
        std::cout << "HPolyhedron is empty and, hence, trivially contained."
                  << std::endl;
        return true;
    }
    // If the HPolyhedron contains P, then maximizing the left-hand side of
    // all inequalities of the HPolyhedron, when constrained in P, still
    // satisfies the right-hand side of the inequalities.

    // Select solver
    std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("GUROBI"));
    if (!solver) {
        LOG(WARNING) << "Solver unavailable.";
        return true;
    }
    solver->Clear();
    // Set silent
    solver->SuppressOutput();

    const double infinity = solver->infinity();
    // Create the X variable
    int NumVars = SpaceDim_;
    std::vector<MPVariable *> x(NumVars);
    for (int i = 0; i < NumVars; i++) {
        x[i] = solver->MakeNumVar(-infinity, infinity, "");
    }
    // Create the inequality constraints: x \in P
    Eigen::MatrixXd Y_Ai = P.Ai();
    Eigen::VectorXd Y_bi = P.bi();
    for (int i = 0; i < NumIneqs_; i++) {
        MPConstraint *con = solver->MakeRowConstraint(-infinity, Y_bi(i));
        for (int j = 0; j < NumVars; j++) {
            con->SetCoefficient(x[j], Y_Ai(i, j));
        }
    }
    // Define the objective: each inequality in the HPolyhedron
    MPObjective *const obj = solver->MutableObjective();
    bool violated = false;
    for (int i = 0; i < NumIneqs_; i++) {
        for (int j = 0; j < NumVars; j++) {
            obj->SetCoefficient(x[j], Ai_(i, j));
        }
        obj->SetMaximization();

        const MPSolver::ResultStatus result_status = solver->Solve();
        if (result_status == MPSolver::INFEASIBLE) {
            LOG(INFO) << "Inequality: " << i
                      << " - Solution status: " << result_status;
            return false;
        }
        // If the maximization value exceed the right-hand side value of the
        // inequality, then return false.
        if (obj->Value() > bi_(i) + 1e-6) {
            LOG(INFO) << "Inequality: " << i
                      << " - Solution status: " << result_status;
            std::cout << "obj->Value: " << obj->Value()
                      << " > bi_(i): " << bi_(i) << std::endl;
            violated = true;
        }
    }
    // Return here so that we see all the violated inequalities.
    return !violated;
}

// Check if HPolyhedron contains a point.
bool HPolyhedron::Contains(const Eigen::VectorXd &point) const {
    bool output = true;
    // Argument checks
    if (!initialized_) {
        std::cerr << "Calling a method on an uninitialized HPolyhedron!"
                  << std::endl;
        output = false;
        return output;
    }
    if (point.rows() != Ai_.cols()) {
        std::cerr << "The vector dimensions is not compatbile!" << std::endl;
        std::cerr << "Expected: " << Ai_.cols() << " | "
                  << "Provided: " << point.rows() << std::endl;
        output = false;
        return output;
    }
    // Check if point satisfies all inequalities in HPolyhedron
    VectorXd ineq = Ai_ * point - bi_;
    for (int i = 0; i < ineq.size(); i++) {
        if (ineq(i) > 0) {
            output = false;
            break;
        }
    }
    return output;
}

// Check if HPolyhedron is (Robustly) Positively Invariant.
bool HPolyhedron::isPositivelyInvariant(const Eigen::MatrixXd &A) const {
    // Argument checks
    if (!initialized_) {
        std::cerr << "Calling a method on an uninitialized HPolyhedron!"
                  << std::endl;
        return false;
    }
    if (A.cols() != A.rows()) {
        std::cerr << "Matrix 'A' must be square." << std::endl;
        return false;
    }
    if (A.cols() != SpaceDim_) {
        std::cerr << "Matrix 'A' cols must be the same as the HPolyhedron's "
                     "dimension."
                  << std::endl;
        return false;
    }

    Eigen::MatrixXd PreA = Ai_ * A;
    HPolyhedron Pre(PreA, bi_);
    HPolyhedron This(Ai_, bi_);
    return Pre.Contains(This);
}

bool HPolyhedron::isPositivelyInvariant(const Eigen::MatrixXd &A,
                                        const Eigen::MatrixXd &E,
                                        HPolyhedron W) const {
    // Argument checks
    if (!initialized_) {
        std::cerr << "Calling a method on an uninitialized HPolyhedron!"
                  << std::endl;
        return false;
    }
    if (A.cols() != A.rows()) {
        std::cerr << "Matrix 'A' must be square." << std::endl;
        return false;
    }
    if (A.cols() != SpaceDim_) {
        std::cerr << "Matrix 'A' cols must be the same as the HPolyhedron's "
                     "dimension."
                  << std::endl;
        return false;
    }
    if (A.rows() != E.rows()) {
        std::cerr << "Matrices 'A' and 'E' must have the same number of rows."
                  << std::endl;
        return false;
    }
    if (E.cols() != W.GetSpaceDim()) {
        std::cerr
            << "Matrix 'E' cols must be the same as the 'W' HPolyhedron's "
               "dimension."
            << std::endl;
        return false;
    }
    HPolyhedron This(Ai_, bi_);

    return (This - W.affineT(E)).isPositivelyInvariant(A);
}

// Compute support of HPolyhedron along each row of A.
VectorXd HPolyhedron::ComputeSupport(const MatrixXd &A) const {
    std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("GUROBI"));
    const double infinity = solver->infinity();

    // Number of hyperplanes
    int NHyperPlanes = A.rows();
    VectorXd supp(NHyperPlanes);
    for (int i = 0; i < NHyperPlanes; i++) {
        supp(i) = infinity;
    }

    // Argument checks
    if (!initialized_) {
        std::cerr << "Calling a method on an uninitialized HPolyhedron!"
                  << std::endl;
        return supp;
    }

#ifdef CIS2M_DEBUG
    std::cout << "S: " << std::endl << Ai_ << std::endl;
    std::cout << "Hyperplanes: " << std::endl << A << std::endl;
#endif

    for (int hp = 0; hp < NHyperPlanes; hp++) {
        solver->Clear();
        // Set silent
        solver->SuppressOutput();
        // Create the X variable
        int NumVars = SpaceDim_;
        std::vector<MPVariable *> x(NumVars);
        for (int i = 0; i < NumVars; i++) {
            x[i] = solver->MakeNumVar(-infinity, infinity, "");
        }
        // Create the inequalities constraints
        for (int i = 0; i < NumIneqs_; i++) {
            MPConstraint *con =
                solver->MakeRowConstraint(-infinity, bi_(i), "");
            for (int j = 0; j < NumVars; j++) {
                con->SetCoefficient(x[j], Ai_(i, j));
            }
        }
        // Add equality constraints, if any
        for (int i = 0; i < NumEqs_; i++) {
            MPConstraint *con = solver->MakeRowConstraint(be_(i), be_(i), "");
            for (int j = 0; j < NumVars; j++) {
                con->SetCoefficient(x[j], Ae_(i, j));
            }
        }
        // Define the objective
        MPObjective *const obj = solver->MutableObjective();
        VectorXd cost_coeff = VectorXd(A.row(hp));
        for (int j = 0; j < NumVars; j++) {
            obj->SetCoefficient(x[j], cost_coeff(j));
        }
        obj->SetMaximization();

        const MPSolver::ResultStatus result_status = solver->Solve();

        if (result_status != MPSolver::OPTIMAL &&
            result_status != MPSolver::UNBOUNDED) {
            std::cerr << "ERROR!" << std::endl;
            std::cerr << "Return Result: " << result_status << std::endl;
            std::cerr << "Hyperplane: " << std::endl
                      << cost_coeff.transpose() << std::endl;
            std::cerr << "Polyhedron: " << std::endl << Ai_ << std::endl;
            std::cerr << std::endl;
        } else if (result_status == MPSolver::UNBOUNDED) {
            supp(hp) = 0;
            std::cout << "Always Good!" << std::endl;
        } else {
            supp(hp) = obj->Value();
            /*
            std::cout << "Optimal Value = " << supp(hp) << std::endl;
            std::cout << "x: ";
            for (auto xi : x) {
                    std::cout << xi->solution_value() << " ";
            }
            std::cout << std::endl;
            */
        }
    }
    return supp;
}

// Eliminate the el_index dimension of the HPolyhedron.
HPolyhedron HPolyhedron::Projection(int el_index, double tol) {
    // Argument checks
    if (!initialized_) {
        std::cerr << "Calling a method on an uninitialized HPolyhedron!"
                  << std::endl;
        return *this;
    }
    if (el_index >= Ai_.cols()) {
        std::cerr << "Elimination index out of bounds!" << std::endl;
        return *this;
    }

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
            A_new.block(cum_counter, 0, 1, left_n_elements) =
                Ai_.block(elp, 0, 1, left_n_elements) +
                Ai_.block(eln, 0, 1, left_n_elements);
            A_new.block(cum_counter, el_index, 1, right_n_elements) =
                Ai_.block(elp, el_index + 1, 1, right_n_elements) +
                Ai_.block(eln, el_index + 1, 1, right_n_elements);
            B_new.block(cum_counter++, 0, 1, 1) =
                bi_.block(elp, 0, 1, 1) + bi_.block(eln, 0, 1, 1);
        }
    }

    for (int i : zeroIndex) {
        A_new.block(cum_counter, 0, 1, left_n_elements) =
            Ai_.block(i, 0, 1, left_n_elements);
        A_new.block(cum_counter, el_index, 1, right_n_elements) =
            Ai_.block(i, el_index + 1, 1, right_n_elements);
        B_new.block(cum_counter++, 0, 1, 1) = bi_.block(i, 0, 1, 1);
    }

    HPolyhedron output(A_new, B_new);
    return output;
}

// Compute affine transformation T * HPolyhedron.
HPolyhedron HPolyhedron::affineT(const MatrixXd &T) {
    // Argument checks
    if (!initialized_) {
        std::cerr << "Calling a method on an uninitialized HPolyhedron!"
                  << std::endl;
        return *this;
    }
    if (T.cols() != SpaceDim_) {
        std::cerr
            << "Columns of T must be the same as the dimension of HPolyhedron!"
            << std::endl;
        return *this;
    }

    int Nc = T.cols();
    int Nr = T.rows();

    // Full LU decomposition
    Eigen::FullPivLU<MatrixXd> lu_decomposition(T);
    int rank = lu_decomposition.rank();

    MatrixXd LLUU = lu_decomposition.matrixLU();
    int NLLUU_r = LLUU.rows();
    int NLLUU_c = LLUU.cols();

    // Extract L
    MatrixXd ll(MatrixXd::Identity(NLLUU_r, NLLUU_c));
    ll.block(0, 0, NLLUU_r, NLLUU_c).triangularView<Eigen::StrictlyLower>() =
        LLUU;
    MatrixXd l(ll.block(0, 0, Nr, rank));

    // Extract U
    MatrixXd uu = lu_decomposition.matrixLU().triangularView<Eigen::Upper>();
    MatrixXd u(uu.block(0, 0, rank, Nc));

    // Compose LU
    MatrixXd LU = l * u;

    MatrixXd LU_11 = LU.block(0, 0, rank, rank);
    MatrixXd LU_12 = LU.block(0, rank, rank, Nc - rank);
    MatrixXd LU_21 = LU.block(rank, 0, Nr - rank, rank);

    // Build the generic inverse for the system
    MatrixXd LU_11_inv = LU_11.inverse();
    MatrixXd invT = MatrixXd::Identity(Nc, Nc);
    invT.block(0, 0, rank, rank) = LU_11_inv;
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
    Ai_final.block(0, 0, P_prime.Ai().rows(), rank) = P_prime.Ai();
    MatrixXd Ae_final(MatrixXd::Zero(Nr - rank, Nr));
    Ae_final.block(0, 0, Nr - rank, rank) = LU_21 * LU_11_inv;
    Ae_final.block(0, rank, Nr - rank, Nr - rank) =
        -MatrixXd::Identity(Nr - rank, Nr - rank);

    // Remove the permutation
    MatrixXd P = lu_decomposition.permutationP();
    Ai_final = Ai_final * P;
    Ae_final = Ae_final * P;
    MatrixXd bi_final = P_prime.bi();
    MatrixXd be_final = MatrixXd::Zero(Nr - rank, 1);

    HPolyhedron P_final(Ai_final, bi_final, Ae_final, be_final);

#ifdef CIS2M_DEBUG
    /*
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
    std::cout << "Reconstruction: " << std::endl << P.inverse() * l * u *
    Q.inverse() << std::endl; std::cout << std::endl;
    */
#endif

    return P_final;
}

// Eliminate redundant inequalities.
void HPolyhedron::RemoveRedundantIneqs() {
    // Argument checks
    if (!initialized_) {
        std::cerr << "Calling a method on an uninitialized HPolyhedron!"
                  << std::endl;
        return;
    }

    std::vector<int> indicesToKeep;
    // For each row of Ai_ ..
    for (int i = 0; i < NumIneqs_; i++) {
        MatrixXd TempA(Ai_.rows() + 1, Ai_.cols());
        TempA << Ai_, Ai_.row(i);
        VectorXd Tempb(bi_.rows() + 1);
        Tempb << bi_, bi_(i) + 1.0;
        // Maximize Ai_(i) * x and see if it is less/equal to bi_(i) ..
        if (solveLP(Ai_.row(i), TempA, Tempb, true).objective_value > bi_(i)) {
            // if not, keep the inequality.
            indicesToKeep.push_back(i);
        }
    }

    // Keep only irredundant rows
    Eigen::VectorXi indicesToKeepVector =
        Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(indicesToKeep.data(),
                                                      indicesToKeep.size());
    MatrixXd NewAi = Ai_;
    NewAi = NewAi(indicesToKeepVector, Eigen::all);
    VectorXd Newbi = bi_;
    Newbi = Newbi(indicesToKeepVector, Eigen::all);
    HPolyhedron IrredunantPoly(NewAi, Newbi);
    // Replace HPolyhedron with the IrredunantPoly
    *this = IrredunantPoly;
}

// Operators.
HPolyhedron &HPolyhedron::operator=(const HPolyhedron &P) {
    Ai_ = P.Ai();
    bi_ = P.bi();

    Ae_ = P.Ae();
    be_ = P.be();

    SpaceDim_ = P.SpaceDim_;
    NumIneqs_ = P.NumIneqs_;
    NumEqs_ = P.NumEqs_;

    initialized_ = P.initialized_;

    return *this;
}

HPolyhedron &HPolyhedron::operator+=(const HPolyhedron &P) {
    // Substract the P.Support(this.A) from the this.b vector
    // max <a_i, x>
    // subj. to x \in P
    // where a_i is the vector representing a hyperplane of *this polyhedron.
    if (SpaceDim_ != P.SpaceDim_) {
        std::cerr << "The HPolyhedron objects should have the same dimension!"
                  << std::endl;
        return *this;
    }

    VectorXd supp = P.ComputeSupport(Ai_);
    bi_ = bi_ + supp;

    return *this;
}

HPolyhedron &HPolyhedron::operator-=(const HPolyhedron &P) {
    // Substract the P.Support(this.A) from the this.b vector
    // max <a_i, x>
    // subj. to x \in P
    // where a_i is the vector representing a hyperplane of *this polyhedron.
    if (SpaceDim_ != P.SpaceDim_) {
        std::cerr << "The HPolyhedron objects should have the same dimension!"
                  << std::endl;
        return *this;
    }

    VectorXd supp = P.ComputeSupport(Ai_);
    bi_ = bi_ - supp;

    //   MatrixXd Adiff(2 * Ai_.rows(), Ai_.cols());
    //   Adiff.block(0,0,Ai_.rows(), Ai_.cols()) = Ai_;
    //   Adiff.block(Ai_.rows(),0,Ai_.rows(), Ai_.cols()) = Ai_;
    //
    //   VectorXd bdiff(2 * bi_.rows(), 1);
    //   bdiff << bi_, (bi_ - supp);
    //
    //   HPolyhedron This(Adiff, bdiff);
    //   *this = This;
    return *this;
}

// Getters
Eigen::MatrixXd HPolyhedron::Ai() const { return Ai_; }

Eigen::VectorXd HPolyhedron::bi() const { return bi_; }

Eigen::MatrixXd HPolyhedron::Ae() const { return Ae_; }

Eigen::VectorXd HPolyhedron::be() const { return be_; }

int HPolyhedron::GetSpaceDim() { return SpaceDim_; }

int HPolyhedron::GetNumInequalities() { return NumIneqs_; }

int HPolyhedron::GetNumEqualities() { return NumEqs_; }

bool HPolyhedron::isInitialized() const { return initialized_; }

bool HPolyhedron::hasEqualityConstraints() const { return (NumEqs_ > 0); }

// Operators
HPolyhedron operator+(HPolyhedron Pl, HPolyhedron Pr) { return Pl += Pr; }

HPolyhedron operator-(HPolyhedron Pl, HPolyhedron Pr) { return Pl -= Pr; }
} // namespace cis2m
