#pragma once

#include <Eigen/Dense>
#include "brunovskyform.hpp"
#include "hpolyhedron.hpp"

namespace cis2m {
class CISGenerator {
	public:
		/**
		 * \brief Constructor with only dynamics information
		 * x_new = Ad * x_old + Bd * u
		 *
		 * \param[in]	Ad	Dynamics matrix of the system		
		 * \param[in]	Bd	Input matrix of the system		
		 *
		 */
		CISGenerator(const Eigen::MatrixXd& Ad,const Eigen::MatrixXd& Bd);


		/**
		 * \brief Constructor with only dynamics information
		 * x_new = Ad * x_old + Bd * u + Ed * w
		 *
		 * \param[in]	Ad	Dynamics matrix of the system		
		 * \param[in]	Bd	Input matrix of the system		
		 * \param[in]	Ed	Disturbance matrix
		 *
		 */
		CISGenerator(const Eigen::MatrixXd& Ad,const Eigen::MatrixXd& Bd, const Eigen::MatrixXd& Ed);


		/// Destructor
		~CISGenerator();


		/**
		 * \brief Add disturbance information 
		 *
		 * \param[in]	ds	Disturbance set described as a HPolyhedron		
		 *
		 * \return 	void	
		 *
		 */
		void AddDisturbanceSet(const HPolyhedron& ds);


		/**
		 * \brief Add input constraints set
		 *
		 * \param[in]	ics	Input Constraints described as a HPolyhedron		
		 *
		 * \return 	void	
		 *
		 */
		void AddInputConstraintsSet(const HPolyhedron& ics);


		/**
		 * \brief Compute the Control Invariant Set
		 *
		 * \param[in]	SafeSet 	Safe set described as a HPolyhedron	
		 * \param[in]	L		Level of Hierarchy 
		 * \param[in]	T		Transient before L 
		 *
		 * \return	Control Invariant Set as a HPolyhedron
		 *
		 */
		void computeCIS(const HPolyhedron& SafeSet, int L, int T);


		/**
		 * \brief Fetch the state part of the Constraint Coefficients
		 *
		 * \return	Control Invariant Set as a HPolyhedron
		 *
		 */
		HPolyhedron Fetch_CIS();


		/**
		 * \brief Fetch the state part of the Constraint Coefficients
		 *
		 * \return	Eigen::MatrixXd 
		 *
		 */
		Eigen::MatrixXd Fetch_A_State();


		/**
		 * \brief Fetch the input part of the Constraint Coefficients
		 *
		 * \return	Eigen::MatrixXd 
		 *
		 */
		Eigen::MatrixXd Fetch_A_Input();

	
		/**
		 * \brief Fetch the virtual input part of the Constraint Coefficients
		 *
		 * \return	Eigen::MatrixXd 
		 *
		 */
		Eigen::MatrixXd Fetch_A_Virtual();


		/**
		 * \brief Compute the input in the Brunovksy coordinates
		 * 
		 * \param[in]	u	Input in the original coordinates	
		 * \param[in]	x	State in the original coordinates	
		 *
		 * \return	Eigen::VectorXd	Input in Brunovksy coordinates
		 */
		Eigen::VectorXd TransformU2B(const Eigen::VectorXd& u, const Eigen::VectorXd& x);


		/**
		 * \brief Compute the input in the Original coordinates
		 * 
		 * \param[in]	u	Input in the Brunovksy coordinates	
		 * \param[in]	x	State in the original coordinates	
		 *
		 * \return	Eigen::VectorXd	Input in Original coordinates
		 */
		Eigen::VectorXd TransformU2O(const Eigen::VectorXd& u, const Eigen::VectorXd& x);


		
		/**
		 * \brief Get the size of the virtual input 
		 */
		int GetExtendedDim() const;

		/**
		 * \brief Get the size of the state
		 */
		int GetStateDim() const;

		int GetLevel() const;

	private:

		int StateDim_;
		int NumberInputs_;
		int DisturbanceDim_;

		int ExtendedInputDim_;
		
		int Level_;
		int Transient_;

		bool ThereAreInputConstraints_;

		HPolyhedron InputCnstrSet_;
		HPolyhedron DisturbanceSet_;

		Eigen::MatrixXd A_lifted_;

		HPolyhedron CIS_;

		Eigen::MatrixXd Ed_BR_;

		Eigen::MatrixXd ExtendedU2U_;

		/**
		 * \brief Reference to the Brunovsky Form Transformation Class
		 */
		BrunovskyForm* brunovsky_form_;

		void Reset();

		void ComputeLiftedSystem(int L, int T);
	
		std::vector<HPolyhedron> ComputeShrinkedSafeSetsSequence(const HPolyhedron& ss);

		void GenerateBrunovksyForm(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B);

		bool cis_computed_;
};
}
