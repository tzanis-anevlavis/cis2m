#include <Eigen/Dense>
#include "brunovskyform.hpp"
#include "hpolyhedron.hpp"

using Eigen::MatrixXd;

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
		CISGenerator(const MatrixXd& Ad,const MatrixXd& Bd);


		/**
		 * \brief Constructor with only dynamics information
		 * x_new = Ad * x_old + Bd * u + Ed * w
		 *
		 * \param[in]	Ad	Dynamics matrix of the system		
		 * \param[in]	Bd	Input matrix of the system		
		 * \param[in]	Ed	Disturbance matrix
		 *
		 */
		CISGenerator(const MatrixXd& Ad,const MatrixXd& Bd, const MatrixXd& Ed);


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
		 * \return	MatrixXd 
		 *
		 */
		MatrixXd Fetch_A_State();


		/**
		 * \brief Fetch the input part of the Constraint Coefficients
		 *
		 * \return	MatrixXd 
		 *
		 */
		MatrixXd Fetch_A_Input();

	
		/**
		 * \brief Fetch the virtual input part of the Constraint Coefficients
		 *
		 * \return	MatrixXd 
		 *
		 */
		MatrixXd Fetch_A_Virtual();

	private:

		int StateDim_;
		int NumberInputs_;
		int DisturbanceDim_;
		
		int Level_;
		int Transient_;

		HPolyhedron InputCnstrSet_;
		HPolyhedron DisturbanceSet_;

		MatrixXd A_lifted_;

		HPolyhedron CIS_;

		MatrixXd Ed_;


		/**
		 * \brief Reference to the Brunovsky Form Transformation Class
		 */
		BrunovskyForm* brunovsky_form_;



		void ComputeLiftedSystem(int L, int T);
	
		std::vector<HPolyhedron> ComputeShrinkedSafeSetsSequence(const HPolyhedron& ss);

		void GenerateBrunovksyForm(const MatrixXd& A, const MatrixXd& B);

		bool cis_computed_;
};
}
