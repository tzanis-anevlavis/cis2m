#pragma once

#include <Eigen/Dense>

namespace cis2m {
class HPolyhedron {
public:
  /**
   * \brief Constructor of empty polyhedron
   */
  HPolyhedron();

  /**
   * \brief Constructor of polyhedron providing the inequality constraints Ax <=
   * b \param[in]	A	Matrix A \param[in]	b	Matrix b
   */
  HPolyhedron(const Eigen::MatrixXd &A, const Eigen::VectorXd &b);

  /**
   * \brief Constructor of a polyhedron providing the inequality and equality
   * constraints:
   * 	Ai x <= bi
   * 	Ae x = be
   * \param[in] Ai	Matrix Ai
   * \param[in] bi	Vector bi
   * \param[in] Ae	Matrix Ae
   * \param[in] be	Vector be
   */
  HPolyhedron(const Eigen::MatrixXd &Ai, const Eigen::VectorXd &bi,
              const Eigen::MatrixXd &Ae, const Eigen::VectorXd &be);
  ~HPolyhedron();

  /**
   * \brief Check if a point belongs to the Polyhedron
   *
   * \param[in]	point	Point to be tested
   *
   * return	Boolean
   */
  bool Contains(const Eigen::VectorXd &point) const;

  /**
   * \brief Check if the Polyhedron is empty
   *
   * return	Boolean
   */
  bool isEmpty() const;

  /**
   * \brief Compute the support of the polyhedron along
   * 	the directions provided in A_other
   * \param[in]	A_other	Matrix containing the directions as rows
   *
   * \return	Vector with the support along the required directions
   */
  Eigen::VectorXd ComputeSupport(const Eigen::MatrixXd &A_other) const;

  /**
   * \brief Projection of the polyhedron on the space with a component removed
   *
   * \param[in]	el_index	Index of the element to be removed
   * \param[in]	tol		Tollerance
   *
   * return 	HPolyhedron
   */
  HPolyhedron Projection(int el_index, double tol = 1e-6);

  /**
   * \brief Affine transformation of a HPolyhedron
   *
   * \param[in]	T	Affine transformation
   */
  HPolyhedron affineT(const Eigen::MatrixXd &T);

  /**
   * \brief Operator=
   */
  HPolyhedron &operator=(const HPolyhedron &other);

  /**
   * \brief Operator+=
   */
  HPolyhedron &operator+=(const HPolyhedron &P);

  /**
   * \brief Operator-=
   */
  HPolyhedron &operator-=(const HPolyhedron &P);

  // ===========================================
  // Getters

  /**
   * \brief Get the Ai matrix
   */
  Eigen::MatrixXd Ai() const;

  /**
   * \brief Get the bi vector
   */
  Eigen::VectorXd bi() const;

  /**
   * \brief Get the Ae matrix
   */
  Eigen::MatrixXd Ae() const;

  /**
   * \brief Get the be vector
   */
  Eigen::VectorXd be() const;

  /**
   * \brief Get the number of inequalities composing the system
   */
  int GetNumInequalities();

  /**
   * \brief Get the size of the state
   */
  int GetSpaceDim();

  /**
   * \brief Check if the HPolyedron is defined
   */
  bool isValid() const;

private:
  bool valid_;

  Eigen::MatrixXd Ai_;
  Eigen::VectorXd bi_;

  Eigen::MatrixXd Ae_;
  Eigen::VectorXd be_;

  int SpaceDim_;
  int NumIneqs_;
};

// Operators  outside che class
HPolyhedron operator+(HPolyhedron Pl, HPolyhedron Pr);
HPolyhedron operator-(HPolyhedron Pl, HPolyhedron Pr);
} // namespace cis2m
