/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_PHAMILTONIAN__
#define __INCLUDE_PHAMILTONIAN__
#include "pUtils.h"
REAL8 EOBHamiltonian(const REAL8 eta,
                     REAL8Vector * x,
                     REAL8Vector * p,
                     REAL8Vector * s1Vec,
                     REAL8Vector * s2Vec,
                     REAL8Vector * sigmaKerr,
                     REAL8Vector * sigmaStar,
                     INT    tortoise,
                     SpinEOBHCoeffs *coeffs);
void CalculateSpinEOBHSACoeffs(REAL8 m1, REAL8 m2, REAL8 s1z, REAL8 s2z, SpinEOBHSACoeffs *coeffs);
REAL8 EOBSAHamiltonian(REAL8 r, REAL8 prT, REAL8 pphi, SpinEOBHSACoeffs *coeffs, REAL8 *invcsi);
int XLALSpinAlignedHcapDerivative_SAConserve(
                  double t,          /**< UNUSED */
                  const REAL8   values[],   /**< dynamical varables */
                  REAL8         dvalues[],  /**< time derivative of dynamical variables */
                  void         *funcParams  /**< EOB parameters */
                  );

INT EOBCalculateSpinEOBHamCoeffs(SpinEOBHCoeffs *coeffs,
                                 const REAL8 eta,
                                 const REAL8 a,
                                 HyperParams *params);
INT XLALSimIMRCalculateSpinPrecEOBHCoeffs_v2(
    SpinEOBHCoeffs *coeffs,           /**<< OUTPUT, EOB parameters including pre-computed coefficients */
    const REAL8 eta,                  /**<< symmetric mass ratio */
    REAL8 a,                          /**<< Normalized deformed Kerr spin */
    REAL8 chi,                        /**<< The augmented spin, with correct aligned-spin limit */
    const UINT SpinAlignedEOBversion, /**<< 4 for SEOBNRv4P; Possible to extend this later later */
    HyperParams *params
);

double GSLSpinHamiltonianWrapper( double x, void *params );
double GSLSpinAlignedHamiltonianWrapper_SA (double x, void *params);

REAL8 XLALSpinHcapNumDerivWRTParam(
                 const INT paramIdx,      /**<< Index of the parameters */
                 const REAL8 values[],     /**<< Dynamical variables */
                 SpinEOBParams *funcParams /**<< EOB Parameters */
);

REAL8
XLALSimIMRSpinEOBHamiltonianDeltaT (SpinEOBHCoeffs * coeffs,
				/**<< Pre-computed coefficients which appear in the function */
				    const REAL8 r,
				/**<< Current orbital radius (in units of total mass) */
				    const REAL8 eta,
				/**<< Symmetric mass ratio */
				    const REAL8 a
				/**<< Normalized deformed Kerr spin */
  );

REAL8
XLALSimIMRSpinEOBHamiltonianDeltaR (SpinEOBHCoeffs * coeffs,
				/**<< Pre-computed coefficients which appear in the function */
				    const REAL8 r,
				/**<< Current orbital radius (in units of total mass) */
				    const REAL8 eta,
				/**<< Symmetric mass ratio */
				    const REAL8 a
				/**<< Normalized deformed Kerr spin */
  );

REAL8
XLALSimIMRSpinAlignedEOBNonKeplerCoeff (const REAL8 values[],
							/**<< Dynamical variables */
					SpinEOBParams * funcParams
							/**<< EOB parameters */
);

REAL8
XLALSimIMRSpinAlignedEOBCalcOmega (const REAL8 values[],/**<< Dynamical variables */
				   SpinEOBParams * funcParams,
							/**<< EOB parameters */
                    REAL8 STEP_SIZE /**<< Step size for numerical derivation of H */
  );

int XLALSpinAlignedHcapDerivative(
                  double t,          /**< UNUSED */
                  const REAL8   values[],   /**< dynamical varables */
                  REAL8         dvalues[],  /**< time derivative of dynamical variables */
                  void         *funcParams  /**< EOB parameters */
);

int XLALSpinAlignedHcapDerivative_SA(
                  double t,          /**< UNUSED */
                  const REAL8   values[],   /**< dynamical varables */
                  REAL8         dvalues[],  /**< time derivative of dynamical variables */
                  void         *funcParams  /**< EOB parameters */
                  );

int XLALSpinAlignedHcapDerivative_Conserve(
                  double t,          /**< UNUSED */
                  const REAL8   values[],   /**< dynamical varables */
                  REAL8         dvalues[],  /**< time derivative of dynamical variables */
                  void         *funcParams  /**< EOB parameters */
                  );

REAL8 XLALSimIMRSpinPrecEOBHamiltonian(
               const REAL8    eta,                  /**<< Symmetric mass ratio */
               REAL8Vector    * x,         /**<< Position vector */
               REAL8Vector    * p,	    /**<< Momentum vector (tortoise radial component pr*) */
               REAL8Vector    * s1Vec,     /**<< Spin vector 1 */
               REAL8Vector    * s2Vec,     /**<< Spin vector 2 */
               REAL8Vector    * sigmaKerr, /**<< Spin vector sigma_kerr */
               REAL8Vector    * sigmaStar, /**<< Spin vector sigma_star */
               INT4                      tortoise,  /**<< flag to state whether the momentum is the tortoise co-ord */
	            SpinEOBHCoeffs *coeffs,               /**<< Structure containing various coefficients */
                HyperParams *hParams
               );

REAL8
GSLSpinPrecHamiltonianWrapper(double x, void *params);

REAL8
XLALSpinPrecHcapNumDerivWRTParam(
			     const INT paramIdx,	/**<< Index of the parameters */
			     const REAL8 values[],	/**<< Dynamical variables */
			     SpinEOBParams * funcParams	/**<< EOB Parameters */
);

int	XLALSpinPrecHcapNumericalDerivative(
                    double	t,	/**<< UNUSED */
                    const	REAL8	values[],	/**<< Dynamical variables */
                    REAL8	dvalues[],	/**<< Time derivatives of variables (returned) */
                    void    *funcParams	/**<< EOB parameters */
);

int	XLALSpinPrecHcapNumericalDerivative_Conserve(
                    double	t,	/**<< UNUSED */
                    const	REAL8	values[],	/**<< Dynamical variables */
                    REAL8	dvalues[],	/**<< Time derivatives of variables (returned) */
                    void    *funcParams	/**<< EOB parameters */
);
int XLALSpinAlignedHcapDerivative_invConserve(
                  double t,          /**< UNUSED */
                  const REAL8   values[],   /**< dynamical varables */
                  REAL8         dvalues[],  /**< time derivative of dynamical variables */
                  void         *funcParams  /**< EOB parameters */
                  );

REAL8 XLALSimIMRSpinPrecEOBHamiltonianDeltaT(
        SpinEOBHCoeffs *coeffs, /**<< Pre-computed coefficients which appear in the function */
        const REAL8    r,       /**<< Current orbital radius (in units of total mass) */
        const REAL8    eta,     /**<< Symmetric mass ratio */
        const REAL8    a        /**<< Normalized deformed Kerr spin */
        );

REAL8 XLALSimIMRSpinPrecEOBHamiltonianDeltaR(
        SpinEOBHCoeffs *coeffs, /**<< Pre-computed coefficients which appear in the function */
        const REAL8    r,       /**<< Current orbital radius (in units of total mass) */
        const REAL8    eta,     /**<< Symmetric mass ratio */
        const REAL8    a        /**<< Normalized deformed Kerr spin */
        );

REAL8 XLALSimIMRSpinPrecEOBCalcOmega(
                      const REAL8           values[],   /**<< Dynamical variables */
                      SpinEOBParams         *funcParams /**<< EOB parameters */
                      );

int XLALSpinPrecHcapRvecDerivative(
            double     t,         /**<< UNUSED */
            const  REAL8      values[],  /**<< Dynamical variables */
            REAL8             dvalues[], /**<< Time derivatives of variables (returned) */
            void             *funcParams /**<< EOB parameters */
                               );

REAL8 XLALCalculateSphHamiltonianDeriv2(
                 const int      idx1,     /**<< Derivative w.r.t. index 1 */
                 const int      idx2,     /**<< Derivative w.r.t. index 2 */
                 const REAL8    values[], /**<< Dynamical variables in spherical coordinates */
                 SpinEOBParams *params    /**<< Spin EOB Parameters */
                 );

#endif

