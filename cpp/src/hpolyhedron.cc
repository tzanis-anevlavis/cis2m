#include "hpolyhedron.hpp"

#include <ortools/linear_solver/linear_solver.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>

namespace cis2m {
namespace {

using Eigen::Index;
using operations_research::MPConstraint;
using operations_research::MPSolver;
using operations_research::MPVariable;

constexpr size_t kMaxProjectionRows = 100000;

// Helper functions used by the HPolyhedron class.

bool ValidTolerance(double tol) {
    // Reject negative and nonfinite tolerances before numerical comparisons.
    return std::isfinite(tol) && tol >= 0.0;
}

double RowInfinityNorm(const MatrixXd& A, Index row) {
    // Infinity norm of a constraint normal; zero identifies a constant row.
    assert((row >= 0) && (row < A.rows()));
    assert(A.cols() > 0);
    return A.row(row).cwiseAbs().maxCoeff();
}

void NormalizeRowsByNormal(MatrixXd& A, VectorXd& b) {
    // Divide each nonzero normal and its RHS by the same positive row norm.
    // This preserves the set and makes residual tolerances independent of the
    // original row scale. Constant rows remain for feasibility handling.
    assert(A.rows() == b.size());
    for (Index i = 0; i < A.rows(); i++) {
        const double scale = RowInfinityNorm(A, i);
        if (scale > 0.0) {
            A.row(i) /= scale;
            b(i) /= scale;
        }
    }
}

void NormalizeRowsByRightHandSide(MatrixXd& A, VectorXd& b, double tol) {
    // Divide rows with a sufficiently large RHS by its absolute value, making
    // the normalized RHS +1 or -1 while preserving the constraint direction.
    assert(A.rows() == b.size());
    assert(ValidTolerance(tol));
    for (Index i = 0; i < A.rows(); i++) {
        const double scale = std::abs(b(i));
        if (scale > tol) {
            A.row(i) /= scale;
            b(i) /= scale;
        }
    }
}

MatrixXd CopyWithoutColumn(const MatrixXd& A, Index column) {
    // Copy all but one zero-based column, preserving the order of the others.
    assert(A.cols() > 1);
    assert((column >= 0) && (column < A.cols()));
    MatrixXd result(A.rows(), A.cols() - 1);
    result.leftCols(column) = A.leftCols(column);
    result.rightCols(A.cols() - column - 1) = A.rightCols(A.cols() - column - 1);
    return result;
}

class LinearProgram {
    // Internal LP model shared by feasibility, support, and redundancy checks.
    // It is built from an original polyhedron `A * x {<=,=} b`.
    //
    // First, it performs coordinate/data normalization: `x = D * z`.
    // This aims to make the magnitudes of each coordinate comperable and is performed
    // given the data of all constraints holistically.
    // Second, it normalizes each constraint `A_i * D * z {<=,=} b_i` independently.
    // The two steps are complementary conditioning steps to improve numerical robustness;
    // both preserve the represented polyhedron.
    //
    // Further the model performs objective normalization to further improve numerical
    // robustness.
public:
    explicit LinearProgram(const HPolyhedron& P) : solver_(MPSolver::CreateSolver("SCIP")) {
        if (!solver_) {
            return;
        }

        // Perform coordinate/data normalization: `x = D * z`.
        // Extract the scales for each coordinate:
        //      coordinate_norms_j = max_i |A_ij| / s_i, where:
        //                      s_i = |b_i| if non-zero, else |A_i|_infty
        // Scales are stored in variable_scales_ = 1 / coordinate_norms.
        VectorXd coordinate_norms = VectorXd::Zero(P.Dimension());
        for (const auto& block : {std::make_pair(&P.Ai(), &P.bi()), std::make_pair(&P.Ae(), &P.be())}) {
            const auto& lhs = block.first;
            const auto& rhs = block.second;
            for (Index i = 0; i < lhs->rows(); i++) {
                // For each row: scale is the RHS if non-zero, else the largest LHS coefficient.
                const double scale =
                    ((*rhs)(i) != 0.0) ? std::abs((*rhs)(i)) : RowInfinityNorm(*lhs, i);
                if (scale == 0.0) {
                    continue;
                }
                // Store the current max over rows all parsed rows.
                coordinate_norms = coordinate_norms.cwiseMax(lhs->row(i).cwiseAbs().transpose() / scale);
            }
        }
        variable_scales_ = VectorXd::Ones(P.Dimension());
        for (Index j = 0; j < coordinate_norms.size(); j++) {
            if (coordinate_norms(j) > 0.0) {
                variable_scales_(j) = 1.0 / coordinate_norms(j);
            }
        }

        if (!coordinate_norms.allFinite() || !variable_scales_.allFinite() || (variable_scales_.array() == 0.0).any()) {
            // Set the LP invalid if there is a numerical issue.
            model_valid_ = false;
            return;
        }

        // Instantiate the LP.
        const double inf = solver_->infinity();
        for (size_t j = 0; j < P.Dimension(); j++) {
            // Initialize variables.
            variables_.push_back(solver_->MakeNumVar(-inf, inf, ""));
        }
        for (Index i = 0; i < P.Ai().rows(); i++) {
            // Instantiate inequality constraints.
            // Save inequality bounds so redundancy checks can temporarily
            // disable individual rows without changing P or rebuilding the model.
            inequalities_.push_back(AddRow(P.Ai(), P.bi(), i, false));
            inequality_bounds_.push_back(inequalities_.back()->ub());
        }
        for (Index i = 0; i < P.Ae().rows(); i++) {
            // Instantiate equality constraints.
            AddRow(P.Ae(), P.be(), i, true);
        }
    }

    MPSolver::ResultStatus Maximize(const Eigen::RowVectorXd& cost_coeffs) {
        // Replace the objective with `cost_coeffs * x`, expressed as `cost_coeffs * D * z`.
        // Normalize it for the solver and retain its scale for Support.
        // A zero `cost_coeffs` is a feasibility query; unsafe scaling is NOT_SOLVED.
        if (!Ready()) {
            return MPSolver::NOT_SOLVED;
        }

        // Rebuild the extracted SCIP model once after any temporary constraint
        // bound changes. This avoids incorrect certificates observed when the
        // incrementally modified model is reused across redundancy checks.
        if (solver_reset_required_) {
            solver_->Reset();
            solver_reset_required_ = false;
        }

        auto* objective = solver_->MutableObjective();
        objective->Clear();
        // cost = cost_coeffs * D
        const Eigen::RowVectorXd cost = cost_coeffs.array() * variable_scales_.transpose().array();
        if (!cost.allFinite() || ((cost.array() == 0.0) && (cost_coeffs.array() != 0.0)).any()) {
            return MPSolver::NOT_SOLVED;
        }

        // Normalize cost
        cost_scale_ = cost.cwiseAbs().maxCoeff();
        if (cost_scale_ > 0.0) {
            for (Index j = 0; j < cost_coeffs.size(); j++) {
                objective->SetCoefficient(variables_[j], cost(j) / cost_scale_);
            }
        }

        objective->SetMaximization();
        return solver_->Solve();
    }

    double Support(const Eigen::RowVectorXd& direction) {
        // Return the supremum in the original coordinates. Only an optimal
        // solve certifies a finite value; unboundedness gives +infinity, while
        // infeasibility or any unresolved solver status gives NaN.
        switch (Maximize(direction)) {
            case MPSolver::OPTIMAL:
                return solver_->Objective().Value() * cost_scale_;
            case MPSolver::UNBOUNDED:
                return std::numeric_limits<double>::infinity();
            default:
                return std::numeric_limits<double>::quiet_NaN();
        }
    }

    bool Ready() const {
        // Model construction succeeded; this does not assert feasibility.
        return solver_ && model_valid_;
    }

    void Disable(Index row) {
        // Relax one inequality in place, retaining its coefficients and handle.
        // row indexes P's inequality block, not the combined constraint list.
        inequalities_[row]->SetBounds(-solver_->infinity(), solver_->infinity());
        solver_reset_required_ = true;
    }

    void Restore(Index row) {
        // Reinstate the saved scaled bound when the disabled row is essential.
        inequalities_[row]->SetBounds(-solver_->infinity(), inequality_bounds_[row]);
        solver_reset_required_ = true;
    }

private:
    MPConstraint* AddRow(const MatrixXd& A, const VectorXd& b, Index row, bool is_equality) {
        // Substitute `x = D * z` and normalize the row by a positive coefficient
        // norm. Inequalities get only an upper bound; equalities get identical
        // lower and upper bounds. Unsafe arithmetic invalidates the whole LP
        // instead of allowing the solver to certify an altered constraint.
        const Eigen::RowVectorXd coefficients = A.row(row).array() * variable_scales_.transpose().array();
        double scale = coefficients.cwiseAbs().maxCoeff();
        // A constant contradiction must not disappear into solver tolerance.
        if (scale == 0.0) {
            scale = (b(row) == 0.0) ? 1.0 : std::abs(b(row));
        }
        const double bound = b(row) / scale;

        const bool is_row_data_finite = coefficients.allFinite() && std::isfinite(bound);
        if (!is_row_data_finite || ((coefficients.array() == 0.0) && (A.row(row).array() != 0.0)).any()) {
            // Keep the row handle usable, but never solve this invalid model.
            model_valid_ = false;
            return solver_->MakeRowConstraint(-solver_->infinity(), solver_->infinity(), "");
        }

        MPConstraint* constraint = solver_->MakeRowConstraint(is_equality ? bound : -solver_->infinity(), bound, "");
        for (Index j = 0; j < A.cols(); j++) {
            constraint->SetCoefficient(variables_[j], coefficients(j) / scale);
        }
        return constraint;
    }

    std::unique_ptr<MPSolver> solver_;          // OR-Tools LP solver instance.
    std::vector<MPVariable*> variables_;        // Solver variables representing scaled coordinates.
    std::vector<MPConstraint*> inequalities_;   // Inequality handles used for temporary disabling.
    std::vector<double> inequality_bounds_;     // Original scaled upper bounds for restoration.
    VectorXd variable_scales_;                  // Diagonal coordinate scaling: x = D*z.
    bool model_valid_ = true;                   // False if scaling or model construction failed.
    bool solver_reset_required_ = false;        // True after constraint bounds change between solves.
    double cost_scale_ = 1.0;                   // Factor used to recover the original objective value.
};

}

// Constructors.

HPolyhedron::HPolyhedron() = default;

// Move constructor.
HPolyhedron::HPolyhedron(HPolyhedron&& other) noexcept
    : valid_(std::exchange(other.valid_, false)), // Move and invalidate `other`
      Ai_(std::move(other.Ai_)),
      bi_(std::move(other.bi_)),
      Ae_(std::move(other.Ae_)),
      be_(std::move(other.be_)) {}

// Move assignment operator.
HPolyhedron& HPolyhedron::operator=(HPolyhedron&& other) noexcept {
    if (this != &other) {
        Ai_ = std::move(other.Ai_);
        bi_ = std::move(other.bi_);
        Ae_ = std::move(other.Ae_);
        be_ = std::move(other.be_);
        valid_ = std::exchange(other.valid_, false); // Move and invalidate `other`
    }
    return *this;
}

// Specialized constructors.
HPolyhedron::HPolyhedron(const MatrixXd& A, const VectorXd& b)
    : HPolyhedron(A, b, MatrixXd(0, A.cols()), VectorXd(0)) {}

HPolyhedron::HPolyhedron(const MatrixXd& Ai, const VectorXd& bi, const MatrixXd& Ae, const VectorXd& be)
    : Ai_(Ai), bi_(bi), Ae_(Ae), be_(be) {
    const Index n = std::max(Ai.cols(), Ae.cols());
    // Accept an omitted 0x0 block, but not a supplied block of the wrong width.
    if (Ai_.rows() == 0 && Ai_.cols() == 0) {
        Ai_.resize(0, n);
    }
    if (Ae_.rows() == 0 && Ae_.cols() == 0) {
        Ae_.resize(0, n);
    }
    valid_ =
        (n > 0) &&
        (Ai_.cols() == Ae_.cols()) &&
        (Ai_.rows() == bi_.size()) &&
        (Ae_.rows() == be_.size()) &&
        Ai_.allFinite() &&
        bi_.allFinite() &&
        Ae_.allFinite() &&
        be_.allFinite();
}

HPolyhedron HPolyhedron::FullSpace(size_t n) {
    if (n == 0 || n > static_cast<size_t>(std::numeric_limits<Index>::max())) {
        return {};
    }
    return HPolyhedron(MatrixXd(0, static_cast<Index>(n)), VectorXd(0));
}

HPolyhedron HPolyhedron::EmptySet(size_t n) {
    if (n == 0 || n > static_cast<size_t>(std::numeric_limits<Index>::max())) {
        return {};
    }
    return HPolyhedron(MatrixXd::Zero(1, static_cast<Index>(n)), VectorXd::Constant(1, -1.0));
}

// Public functions.

void HPolyhedron::Show() const {
    std::cout << "HPolyhedron: valid=" << IsValid() << ", dimension=" << Dimension()
              << ", inequalities=" << NumInequalities() << ", equalities=" << NumEqualities()
              << "\nAi:\n" << Ai_ << "\nbi: " << bi_.transpose()
              << "\nAe:\n" << Ae_ << "\nbe: " << be_.transpose() << '\n';
}

void HPolyhedron::NormalizeConstraints(double tol) {
    if (!valid_ || !ValidTolerance(tol)) {
        return;
    }
    MatrixXd Ai = Ai_;
    VectorXd bi = bi_;
    MatrixXd Ae = Ae_;
    VectorXd be = be_;
    NormalizeRowsByRightHandSide(Ai, bi, tol);
    NormalizeRowsByRightHandSide(Ae, be, tol);
    // Update object.
    *this = HPolyhedron(Ai, bi, Ae, be);
}

bool HPolyhedron::Contains(const VectorXd& point, double tol) const {
    if (!valid_ || point.size() != Ai_.cols() || !point.allFinite() || !ValidTolerance(tol)) {
        return false;
    }
    const VectorXd inequality_residual = Ai_ * point - bi_;
    const VectorXd equality_residual = Ae_ * point - be_;
    return inequality_residual.allFinite() &&
        equality_residual.allFinite() &&
        (inequality_residual.array() <= tol).all() &&
        (equality_residual.array().abs() <= tol).all();
}

bool HPolyhedron::IsFeasible() const {
    if (!valid_) {
        return false;
    }
    if (IsFullSpace()) {
        return true;
    }
    // Feasibility is an LP with zero objective.
    LinearProgram lp(*this);
    const auto status = lp.Maximize(Eigen::RowVectorXd::Zero(Ai_.cols()));
    if (status == MPSolver::OPTIMAL || status == MPSolver::FEASIBLE) {
        return true;
    }
    if (status == MPSolver::INFEASIBLE) {
        return false;
    }
    throw std::runtime_error("HPolyhedron::IsFeasible: LP did not determine feasibility.");
}

VectorXd HPolyhedron::ComputeSupport(const MatrixXd& directions) const {
    // Initialize output to NaN.
    VectorXd result = VectorXd::Constant(directions.rows(), std::numeric_limits<double>::quiet_NaN());
    if (!valid_ || (directions.cols() != Ai_.cols()) || !directions.allFinite() || (directions.rows() == 0)) {
        return result;
    }
    if (IsFullSpace()) {
        for (Index i = 0; i < directions.rows(); i++) {
            result(i) = directions.row(i).isZero(0.0) ? 0.0 : std::numeric_limits<double>::infinity();
        }
        return result;
    }
    // Solve an LP over the polyhedron maximizing along the `directions`.
    LinearProgram lp(*this);
    for (Index i = 0; i < directions.rows(); i++) {
        result(i) = lp.Support(directions.row(i));
    }
    return result;
}

bool HPolyhedron::IsBounded() const {
    if (!valid_) {
        return false;
    }
    if (!IsFeasible()) {
        return true; // The empty set is bounded.
    }
    // Finite support in both directions of every coordinate axis gives a
    // containing box, hence boundedness, including lower-dimensional sets.
    MatrixXd directions(2 * Ai_.cols(), Ai_.cols());
    directions.topRows(Ai_.cols()).setIdentity();
    directions.bottomRows(Ai_.cols()) = -MatrixXd::Identity(Ai_.cols(), Ai_.cols());
    return ComputeSupport(directions).allFinite();
}

bool HPolyhedron::IsSubsetOf(const HPolyhedron& other, double tol) const {
    if (!valid_ || !other.valid_ || (Dimension() != other.Dimension()) || !ValidTolerance(tol)) {
        return false;
    }
    if (other.IsFullSpace() || !IsFeasible()) {
        return true;
    }
    // Both signs of each equality are needed to test its whole affine hull.
    MatrixXd directions(other.Ai_.rows() + 2 * other.Ae_.rows(), Ai_.cols());
    directions << other.Ai_, other.Ae_, -other.Ae_;
    VectorXd bounds(directions.rows());
    bounds << other.bi_, other.be_, -other.be_;
    NormalizeRowsByNormal(directions, bounds);
    const VectorXd support = ComputeSupport(directions);
    return support.allFinite() && ((support - bounds).array() <= tol).all();
}

HPolyhedron HPolyhedron::Intersection(const HPolyhedron& other) const {
    if (!valid_ || !other.valid_ || (Dimension() != other.Dimension())) {
        return {};
    }
    MatrixXd Ai(Ai_.rows() + other.Ai_.rows(), Ai_.cols());
    MatrixXd Ae(Ae_.rows() + other.Ae_.rows(), Ae_.cols());
    VectorXd bi(bi_.size() + other.bi_.size()), be(be_.size() + other.be_.size());
    Ai << Ai_, other.Ai_;
    bi << bi_, other.bi_;
    Ae << Ae_, other.Ae_;
    be << be_, other.be_;
    return HPolyhedron(Ai, bi, Ae, be);
}

HPolyhedron HPolyhedron::Preimage(const MatrixXd& T) const {
    if (!valid_ || (T.rows() != Ai_.cols()) || (T.cols() == 0) || !T.allFinite()) {
        return {};
    }
    return HPolyhedron(Ai_ * T, bi_, Ae_ * T, be_);
}

void HPolyhedron::RemoveRedundantConstraints(double tol) {
    if (!valid_ || !ValidTolerance(tol) || (Ai_.rows() == 0)) {
        return;
    }

    // Get a copy of the polyhedron and normalize its rows.
    HPolyhedron normalized = *this;
    NormalizeRowsByNormal(normalized.Ai_, normalized.bi_);
    if (!normalized.Ai_.allFinite() || !normalized.bi_.allFinite()) {
        return;
    }

    // Test each inequality against the remaining constraints. Removed rows
    // stay excluded from later tests, so duplicate rows cannot all certify
    // one another as redundant and disappear together.
    std::vector<bool> removed(Ai_.rows(), false);
    LinearProgram lp(normalized);
    if (!lp.Ready()) {
        return;
    }
    for (Index i = 0; i < Ai_.rows(); i++) {
        if (normalized.Ai_.row(i).isZero(0.0)) {
            // Check if constraint is trivially satisfied and mark it removed.
            if (bi_(i) < 0.0) {
                *this = EmptySet(Dimension());
                return;
            }
            removed[i] = true;
            lp.Disable(i);
            continue;
        }

        // Compute the support over the directions of the i-th constraint,
        // while disabling said constraint.
        lp.Disable(i);
        const double support = lp.Support(normalized.Ai_.row(i));
        // Keep the constraint if the LP is unbounded, infeasible, or inconclusive;
        // only a finite support value can certify redundancy here.
        removed[i] = std::isfinite(support) && (support <= (normalized.bi_(i) + tol));
        if (!removed[i]) {
            lp.Restore(i);
        }
    }

    const Index kept = static_cast<Index>(std::count(removed.begin(), removed.end(), false));
    MatrixXd A(kept, Ai_.cols());
    VectorXd b(kept);
    Index row = 0;
    for (Index i = 0; i < Ai_.rows(); i++) {
        if (!removed[i]) {
            A.row(row) = Ai_.row(i);
            b(row) = bi_(i);
            row++;
        }
    }
    Ai_ = std::move(A);
    bi_ = std::move(b);
}

HPolyhedron HPolyhedron::Projection(size_t var_index, double base_tol) const {
    if (!valid_ || (Dimension() <= 1) || (var_index >= Dimension()) || !ValidTolerance(base_tol)) {
        return {};
    }
    const Index column = static_cast<Index>(var_index);
    MatrixXd Ai = Ai_, Ae = Ae_;
    VectorXd bi = bi_, be = be_;
    NormalizeRowsByNormal(Ai, bi);
    NormalizeRowsByNormal(Ae, be);
    if (!Ai.allFinite() || !bi.allFinite() || !Ae.allFinite() || !be.allFinite()) {
        return {};
    }

    // An equality can eliminate the coordinate directly. Prefer the largest
    // coefficient after row normalization to reduce division by small pivots.
    Index pivot = -1;
    double largest = 0.0;
    for (Index i = 0; i < Ae.rows(); i++) {
        if (std::abs(Ae(i, column)) > largest) {
            pivot = i;
            largest = std::abs(Ae(i, column));
        }
    }
    if (pivot >= 0) {
        // Substitute once, explicitly zeroing the eliminated coordinate in
        // every remaining row. It must not undergo Fourier-Motzkin again.
        const Eigen::RowVectorXd equation = Ae.row(pivot) / Ae(pivot, column);
        const double rhs = be(pivot) / Ae(pivot, column);
        for (Index i = 0; i < Ai.rows(); i++) {
            const double coefficient = Ai(i, column);
            Ai.row(i) -= coefficient * equation;
            bi(i) -= coefficient * rhs;
            Ai(i, column) = 0.0;
        }
        MatrixXd remaining_Ae(Ae.rows() - 1, Ae.cols());
        VectorXd remaining_be(Ae.rows() - 1);
        Index row = 0;
        for (Index i = 0; i < Ae.rows(); i++) {
            if (i == pivot) {
                continue;
            }
            remaining_Ae.row(row) = Ae.row(i) - Ae(i, column) * equation;
            remaining_be(row) = be(i) - Ae(i, column) * rhs;
            remaining_Ae(row++, column) = 0.0;
        }
        HPolyhedron result(
            CopyWithoutColumn(Ai, column), bi,
            CopyWithoutColumn(remaining_Ae, column), remaining_be);
        result.RemoveRedundantConstraints(base_tol);
        return result;
    }

    // Without an equality pivot, separate upper and lower bounds on the
    // eliminated coordinate. Normalize its coefficients to +1 and -1 so
    // adding each opposite-sign pair cancels that coordinate exactly.
    std::vector<Index> positive, negative, independent;
    for (Index i = 0; i < Ai.rows(); i++) {
        const double coefficient = Ai(i, column);
        if (coefficient > 0.0) {
            positive.push_back(i);
        }
        else if (coefficient < 0.0) {
            negative.push_back(i);
        }
        else {
            independent.push_back(i);
        }
        if (coefficient != 0.0) {
            Ai.row(i) /= std::abs(coefficient);
            bi(i) /= std::abs(coefficient);
        }
    }
    // Check the total before multiplication/allocation, including independent rows.
    if ((independent.size() > kMaxProjectionRows) ||
        (!negative.empty() && (positive.size() > (kMaxProjectionRows - independent.size()) / negative.size()))) {
        return {};
    }
    const size_t count = independent.size() + positive.size() * negative.size();
    MatrixXd A(static_cast<Index>(count), Ai.cols());
    VectorXd b(static_cast<Index>(count));
    Index row = 0;
    for (Index p : positive) {
        for (Index n : negative) {
            A.row(row) = Ai.row(p) + Ai.row(n);
            b(row++) = bi(p) + bi(n);
        }
    }
    // A one-sided family imposes no restriction on the retained coordinates.
    for (Index i : independent) {
        A.row(row) = Ai.row(i);
        b(row++) = bi(i);
    }
    HPolyhedron result(CopyWithoutColumn(A, column), b, CopyWithoutColumn(Ae, column), be);
    result.RemoveRedundantConstraints(base_tol);
    return result;
}

HPolyhedron HPolyhedron::ProjectOnto(const std::vector<size_t>& coordinates, double base_tol) const {
    if (!valid_ || coordinates.empty() || !ValidTolerance(base_tol)) {
        return {};
    }
    std::vector<bool> keep(Dimension(), false);
    for (size_t coordinate : coordinates) {
        if (coordinate >= Dimension() || keep[coordinate]) {
            return {};
        }
        keep[coordinate] = true;
    }
    HPolyhedron result = *this;
    // Eliminate in descending order so earlier coordinate indices stay valid.
    for (Index j = static_cast<Index>(Dimension()) - 1; j >= 0; j--) {
        const size_t coordinate = static_cast<size_t>(j);
        if (!keep[coordinate]) {
            result = result.Projection(coordinate, base_tol);
            if (!result.valid_) {
                return result;
            }
        }
    }
    MatrixXd Ai(result.Ai_.rows(), coordinates.size());
    MatrixXd Ae(result.Ae_.rows(), coordinates.size());
    // Surviving columns are in original index order; reorder them to match
    // the caller's requested coordinates.
    for (size_t j = 0; j < coordinates.size(); j++) {
        const Index source = std::count(keep.begin(), keep.begin() + coordinates[j], true);
        Ai.col(j) = result.Ai_.col(source);
        Ae.col(j) = result.Ae_.col(source);
    }
    return HPolyhedron(Ai, result.bi_, Ae, result.be_);
}

HPolyhedron HPolyhedron::AffineTransform(const MatrixXd& T, double condition_threshold) const {
    if (!valid_ || (T.cols() != Ai_.cols()) || (T.rows() == 0) ||
        !T.allFinite() || !std::isfinite(condition_threshold) || (condition_threshold < 1.0)) {
        return {};
    }
    if (T.rows() == T.cols()) {
        Eigen::JacobiSVD<MatrixXd> svd(T);
        const auto& values = svd.singularValues();
        if ((values(values.size() - 1) > 0.0) && (values(0) / values(values.size() - 1) <= condition_threshold)) {
            return AffineTransform_QR(T);
        }
    }
    return AffineTransform_SVD(T);
}

HPolyhedron HPolyhedron::AffineTransform_QR(const MatrixXd& T) const {
    if (!valid_ || (T.cols() != Ai_.cols()) || (T.rows() == 0)|| !T.allFinite()) {
        return {};
    }
    if (T.rows() != T.cols()) {
        return AffineTransform_SVD(T);
    }
    Eigen::ColPivHouseholderQR<MatrixXd> qr(T.transpose());
    if (qr.rank() != T.cols()) {
        return AffineTransform_SVD(T);
    }
    // Solve T' * G' = A' without explicitly constructing an inverse.
    const MatrixXd Ai = qr.solve(Ai_.transpose()).transpose();
    const MatrixXd Ae = qr.solve(Ae_.transpose()).transpose();
    return HPolyhedron(Ai, bi_, Ae, be_);
}

HPolyhedron HPolyhedron::AffineTransform_SVD(const MatrixXd& T) const {
    if (!valid_ || (T.cols() != Ai_.cols()) || (T.rows() == 0) || !T.allFinite()) {
        return {};
    }
    const Index n = T.cols(), m = T.rows();
    Eigen::JacobiSVD<MatrixXd> svd(T, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Index rank = svd.rank();
    const MatrixXd inverse =
        svd.matrixV().leftCols(rank) *
        svd.singularValues().head(rank).cwiseInverse().asDiagonal() *
        svd.matrixU().leftCols(rank).transpose();
    const MatrixXd kernel = svd.matrixV().rightCols(n - rank);
    const Index lifted_dim = m + n - rank;

    // x = T^+ y + kernel*s is valid only for y in Im(T). Enforce the
    // orthogonal complement of Im(T) for every rank, including the zero map.
    MatrixXd Ai(Ai_.rows(), lifted_dim);
    Ai << Ai_ * inverse, Ai_ * kernel;
    MatrixXd Ae = MatrixXd::Zero(Ae_.rows() + m - rank, lifted_dim);
    Ae.topRows(Ae_.rows()).leftCols(m) = Ae_ * inverse;
    Ae.topRows(Ae_.rows()).rightCols(n - rank) = Ae_ * kernel;
    Ae.bottomRows(m - rank).leftCols(m) = svd.matrixU().rightCols(m - rank).transpose();
    VectorXd be(Ae.rows());
    be << be_, VectorXd::Zero(m - rank);
    HPolyhedron result(Ai, bi_, Ae, be);
    for (Index j = lifted_dim; j > m && result.valid_; --j) {
        result = result.Projection(static_cast<size_t>(j - 1), 1e-10);
    }
    return result;
}

HPolyhedron& HPolyhedron::operator+=(const HPolyhedron& other) {
    if (!valid_ || !other.valid_ || (Dimension() != other.Dimension())) {
        *this = HPolyhedron();
        return *this;
    }
    if (IsFullSpace() || other.IsFullSpace()) {
        *this = FullSpace(Dimension());
        return *this;
    }

    // z = x + y: P(z-y), Q(y). Eliminating y gives the exact sum, including
    // facets absent from P. Updating only b with support_Q(A) is not exact.
    const Index n = Ai_.cols();
    MatrixXd Ai = MatrixXd::Zero(Ai_.rows() + other.Ai_.rows(), 2 * n);
    Ai.topRows(Ai_.rows()).leftCols(n) = Ai_;
    Ai.topRows(Ai_.rows()).rightCols(n) = -Ai_;
    Ai.bottomRows(other.Ai_.rows()).rightCols(n) = other.Ai_;
    VectorXd bi(bi_.size() + other.bi_.size());
    bi << bi_, other.bi_;
    MatrixXd Ae = MatrixXd::Zero(Ae_.rows() + other.Ae_.rows(), 2 * n);
    Ae.topRows(Ae_.rows()).leftCols(n) = Ae_;
    Ae.topRows(Ae_.rows()).rightCols(n) = -Ae_;
    Ae.bottomRows(other.Ae_.rows()).rightCols(n) = other.Ae_;
    VectorXd be(be_.size() + other.be_.size());
    be << be_, other.be_;
    HPolyhedron result(Ai, bi, Ae, be);
    for (Index j = 2 * n; j > n && result.valid_; --j) {
        result = result.Projection(static_cast<size_t>(j - 1), 1e-10);
    }
    *this = std::move(result);
    return *this;
}

HPolyhedron HPolyhedron::PontryaginDifferenceOfLinearImage(const HPolyhedron& W, const MatrixXd& T) const {
    if (!valid_ || !W.valid_ || (T.rows() != Ai_.cols()) || (T.cols() != W.Ai_.cols()) || !T.allFinite()) {
        return {};
    }
    if (IsFullSpace()) {
        return *this;
    }
    // For every inequality g*x <= b, require g*x <= b - h_W(g*T).
    // Evaluate both signs of equality normals in the same support batch to
    // determine whether the disturbance is constant in each such direction.
    MatrixXd normals(Ai_.rows() + 2 * Ae_.rows(), Ai_.cols());
    normals << Ai_, Ae_, -Ae_;
    const VectorXd support = W.ComputeSupport(normals * T);
    if (support.array().isNaN().any()) {
        return {};
    }
    // Infinite support along a constrained direction makes the robustified
    // set empty; a numerical failure above must remain a distinct invalid set.
    if (!support.allFinite()) {
        return EmptySet(Dimension());
    }

    VectorXd bi = bi_ - support.head(Ai_.rows());
    VectorXd be = be_;
    for (Index i = 0; i < Ae_.rows(); i++) {
        const double upper = support(Ai_.rows() + i);
        const double lower = -support(Ai_.rows() + Ae_.rows() + i);
        // Both directions must have the same support location. Scaling by
        // the normal avoids accepting a varying image of a tiny equality.
        const double scale = RowInfinityNorm(Ae_, i);
        const double tolerance = 1e-9 * scale;
        if (std::abs(upper - lower) > tolerance) {
            return EmptySet(Dimension());
        }
        be(i) -= 0.5 * upper + 0.5 * lower;
    }
    return HPolyhedron(Ai_, bi, Ae_, be);
}

HPolyhedron& HPolyhedron::operator-=(const HPolyhedron& other) {
    if (!valid_ || !other.valid_ || (Dimension() != other.Dimension())) {
        *this = HPolyhedron();
        return *this;
    }
    *this = PontryaginDifferenceOfLinearImage(other, MatrixXd::Identity(Ai_.cols(), Ai_.cols()));
    return *this;
}

}
