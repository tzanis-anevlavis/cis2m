#include <Eigen/Dense>
#include "brunovskyform.hpp"
#include "hpolyhedron.hpp"

using Eigen::MatrixXd;

namespace cis2m {
class CISGenerator {
	public:
		/// Constructor
		CISGenerator(const MatrixXd& Ad,const MatrixXd& Bd);
		/// Destructor
		~CISGenerator();

		/**
		 * \brief Compute the Control Invariant Set
		 *
		 * \param[in]	SafeSet 	Safe set described as a HPolyhedron	
		 * \param[in]	InputCnstr	Input constraints described as a HPolyhedron		
		 * \param[in]	Disturbance	Disturbance set described as a HPolyhedron		
		 * \param[in]	L		Level of Hierarchy 
		 * \param[in]	T		Transient before L 
		 *
		 * \return	Control Invariant Set as a HPolyhedron
		 *
		 */
		HPolyhedron computeCIS(const HPolyhedron& SafeSet, int L, int T);
	private:

		int StateDim_;
		int NumberInputs_;

		HPolyhedron* CIS_;

		BrunovskyForm* brunovsky_form_;
		std::optional<HPolyhedron> FetchCIS();

		bool cis_computed_;
};
}
