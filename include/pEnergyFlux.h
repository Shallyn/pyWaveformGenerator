/**
 * Writer: Xiaolin.liu
 * xiaolin.liu@mail.bnu.edu.cn
 *
 * This module contains basic functions for  calculation.
 * Functions list:
 * Kernel:
 * 20xx.xx.xx, LOC
 **/

#ifndef __INCLUDE_PENERGYFLUX__
#define __INCLUDE_PENERGYFLUX__
#include "pUtils.h"
REAL8
XLALInspiralSpinFactorizedFlux(REAL8Vector *polvalues, /**< polar dynamical variables */
                               REAL8Vector *values,    /**< cart dynamical variables */
                               EOBNonQCCoeffs *nqcCoeffs,
                               /**< pre-computed NQC coefficients */
                               const REAL8 omega,                                   /**< orbital frequency */
                               const REAL8 dr, const REAL8 ncrv, SpinEOBParams *ak, /**< physical parameters */
                               const REAL8 H,                                       /**< real Hamiltonian */
                               const INT lMax /**< upper limit of the summation over l */
);

REAL8
XLALInspiralSpinFactorizedFlux_SA(REAL8Vector *polvalues, /**< polar dynamical variables */
                                  EOBNonQCCoeffs *nqcCoeffs,
                                  /**< pre-computed NQC coefficients */
                                  const REAL8 omega,                                   /**< orbital frequency */
                                  const REAL8 dr, const REAL8 ncrv, SpinEOBParams *ak, /**< physical parameters */
                                  const REAL8 H,                                       /**< real Hamiltonian */
                                  const INT lMax /**< upper limit of the summation over l */
);

REAL8
XLALInspiralPrecSpinFactorizedFlux(REAL8Vector *polvalues,    /**< \f$(r,\phi,p_r,p_\phi)\f$ */
                                   REAL8Vector *values,       /**< dynamical variables */
                                   EOBNonQCCoeffs *nqcCoeffs, /**< pre-computed NQC coefficients */
                                   const REAL8 omega,         /**< orbital frequency */
                                   const REAL8 dr, const REAL8 ncrv, SpinEOBParams *ak, /**< physical parameters */
                                   const REAL8 H,                                       /**< real Hamiltonian */
                                   const INT4 lMax,                 /**< upper limit of the summation over l */
                                   const UINT SpinAlignedEOBversion /**< 1 for SEOBNRv1, 2 for SEOBNRv2 */
);

INT CalculateRRForceCoeffs(RRForceCoeffs *coeffsFf, RRForceCoeffs *coeffsFr, SpinEOBParams *params);

INT CalculateRRForce(RRForceCoeffs *coeffsFf, RRForceCoeffs *coeffsFr, REAL8 *rrFf, REAL8 *rrFr, REAL8 r, REAL8 pr,
                     REAL8 p2, REAL8 pf);

INT CalculateRRForceSpinCoeffs(RRForceCoeffs *coeffsFf, RRForceCoeffs *coeffsFr, REAL8 m1, REAL8 m2, REAL8 chi1,
                               REAL8 chi2);
void CalculateEccCorrectionCoeffs(REAL8 eta, REAL8 chi1, REAL8 chi2, EccCorrectionCoeffs *coeffs);

INT CalculateEccCorrectionToFlux(REAL8 eta, REAL8 chi1, REAL8 chi2, REAL8 r, REAL8 pr, REAL8 prDot, REAL8 *cFr,
                                 REAL8 *cFf);

INT CalculateEccCorrectionToFluxV2(REAL8 eta, REAL8 chi1, REAL8 chi2, REAL8 r, REAL8 prt, REAL8 prDot, REAL8 *cFr,
                                   REAL8 *cFf);

INT CalculateEccCorrectionToFluxV3(REAL8 eta, REAL8 chi1, REAL8 chi2, REAL8 r, REAL8 vr, REAL8 prDot, REAL8 *cFr,
                                   REAL8 *cFf, REAL8 e0);

INT CalculateEccCorrectionToFluxV3X(REAL8 r, REAL8 vr, REAL8 prDot, REAL8 *cFr, REAL8 *cFf, REAL8 e0,
                                    EccCorrectionCoeffs *coeffs);

INT CalculateEccCorrectionToFluxV4(REAL8 eta, REAL8 chi1, REAL8 chi2, REAL8 r, REAL8 prt, REAL8 prDot, REAL8 *cFr,
                                   REAL8 *cFf);

#endif
