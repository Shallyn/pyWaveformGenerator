/**
 * Writer: Xiaolin.liu
 * xiaolin.liu@mail.bnu.edu.cn
 *
 * This module contains basic functions for  calculation.
 * Functions list:
 * Kernel:
 * 20xx.xx.xx, LOC
 **/

#ifndef __INCLUDE_PFACTORIZEDWAVEFORM__
#define __INCLUDE_PFACTORIZEDWAVEFORM__
#include "pUtils.h"

int XLALSimIMRSpinEOBCalculateNewtonianMultipole(
    COMPLEX16 *multipole, /**<< OUTPUT, Newtonian multipole */
    REAL8 x,              /**<< Dimensionless parameter \f$\equiv v^2\f$ */
    REAL8 r,              /**<< Orbital separation (units of total mass M */
    REAL8 phi,            /**<< Orbital phase (in radians) */
    UINT l,               /**<< Mode l */
    INT m,                /**<< Mode m */
    SpinEOBParams *params /**<< Pre-computed coefficients, parameters, etc. */
);

int XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole(
    COMPLEX16 *multipole, /**<< OUTPUT, Newtonian multipole */
    REAL8 x,              /**<< Dimensionless parameter \f$\equiv v^2\f$ */
    REAL8 r,              /**<< Orbital separation (units of total mass M */
    REAL8 phi,            /**<< Orbital phase (in radians) */
    UINT l,               /**<< Mode l */
    INT m,                /**<< Mode m */
    SpinEOBParams *params /**<< Pre-computed coefficients, parameters, etc. */
);

INT XLALSimIMREOBComputeNewtonMultipolePrefixes(
    NewtonMultipolePrefixes *prefix, /**<< OUTPUT Structure containing the coeffs */
    const REAL8 m1,                  /**<< Mass of first component */
    const REAL8 m2                   /**<< Nass of second component */
);

INT XLALSimIMREOBCalcSpinFacWaveformCoefficients(FacWaveformCoeffs *const coeffs, SpinEOBParams *params, REAL8 a,
                                                 const REAL8 chiS, const REAL8 chiA);

INT XLALSimIMRSpinEOBGetSpinFactorizedWaveform(COMPLEX16 *ret,
                                               /**< OUTPUT, hlm waveforms */
                                               REAL8Vector *values,
                                               /**< dyanmical variables */
                                               const REAL8 v,
                                               /**< velocity */
                                               const REAL8 Hreal,
                                               /**< real Hamiltonian */
                                               const INT4 l,
                                               /**< l mode index */
                                               const INT4 m,
                                               /**< m mode index */
                                               SpinEOBParams *params,
                                               /**< Spin EOB parameters */
                                               INT ret_type
                                               /**< ret_type: 0:full, 1:amp, 2:amp_res */
);

int XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(
    FacWaveformCoeffs *const coeffs, /**< OUTPUT, pre-computed waveform coefficients */
    const REAL8 m1,                  /**< mass 1 */
    const REAL8 m2,                  /**< mass 2 */
    const REAL8 eta,                 /**< symmetric mass ratio */
    const REAL8 tmpa,                /**< Kerr spin parameter for test-particle terms */
    const REAL8 chiS,                /**< (chi1+chi2)/2 */
    const REAL8 chiA,                /**< (chi1-chi2)/2 */
    UINT SpinAlignedEOBversion       /**< 1 for SEOBNRv1; 2 for SEOBNRv2; 4 for the
      coefficients in the flux of v4P and v4Pwave for the coefficients in the
      waveform of v4P */
);

INT XLALSimIMRSpinEOBFluxGetPrecSpinFactorizedWaveform(
    COMPLEX16 *hlmTab,       /**< OUTPUT, hlm waveforms */
    REAL8Vector *values,     /**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
    REAL8Vector *cartvalues, /**< dyanmical variables */
    const REAL8 v,           /**< velocity */
    const REAL8 Hreal,       /**< real Hamiltonian */
    const INT4 lMax,         /**< maximum l mode to compute, compute 0 < m <= lMax */
    SpinEOBParams *params    /**< Spin EOB parameters */
);

INT XLALSimIMRSpinEOBFluxGetSASpinFactorizedWaveform(
    COMPLEX16 *hlmTab,    /**< OUTPUT, hlm waveforms */
    REAL8Vector *values,  /**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
    const REAL8 v,        /**< velocity */
    const REAL8 Hreal,    /**< real Hamiltonian */
    const INT4 lMax,      /**< maximum l mode to compute, compute 0 < m <= lMax */
    SpinEOBParams *params /**< Spin EOB parameters */
);

int XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(
    FacWaveformCoeffs *const coeffs, /**< OUTPUT, pre-computed waveform coefficients */
    const REAL8 m1,                  /**< mass 1 */
    const REAL8 m2,                  /**< mass 2 */
    const REAL8 eta,                 /**< symmetric mass ratio */
    const REAL8 tmpa,                /**< Kerr spin parameter for test-particle terms */
    const REAL8 chiS,                /**< (chi1+chi2)/2 */
    const REAL8 chiA,                /**< (chi1-chi2)/2 */
    UINT SpinAlignedEOBversion       /**< 1 for SEOBNRv1; 2 for SEOBNRv2; 4 for the
      coefficients in the flux of v4P and v4Pwave for the coefficients in the
      waveform of v4P */
);

INT XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform(
    COMPLEX16 *hlm,          /**< OUTPUT, hlm waveforms */
    REAL8Vector *values,     /**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
    REAL8Vector *cartvalues, /**< dyanmical variables */
    const REAL8 v,           /**< velocity */
    const REAL8 Hreal,       /**< real Hamiltonian */
    const INT l,             /**< l mode index */
    const INT m,             /**< m mode index */
    SpinEOBParams *params    /**< Spin EOB parameters */
);

INT XLALSimIMRSpinEOBGetSASpinFactorizedWaveform(
    COMPLEX16 *hlm,       /**< OUTPUT, hlm waveforms */
    REAL8Vector *values,  /**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
    const REAL8 v,        /**< velocity */
    const REAL8 Hreal,    /**< real Hamiltonian */
    const INT l,          /**< l mode index */
    const INT m,          /**< m mode index */
    SpinEOBParams *params /**< Spin EOB parameters */
);

INT XLALSimIMRSpinEOBGetAmplitudeResidualPrec(COMPLEX16 *rholmpwrl, const REAL8 v, const REAL8 Hreal, const INT modeL,
                                              const INT modeM, SpinEOBParams *params);

REAL8 XLALSimIMRSpinPrecEOBNonKeplerCoeff(const REAL8 values[],     /**<< Dynamical variables */
                                          SpinEOBParams *funcParams /**<< EOB parameters */
);

#endif
