/**
 * Writer: Xiaolin.liu
 * xiaolin.liu@mail.bnu.edu.cn
 *
 * This module contains basic functions for  calculation.
 * Functions list:
 * Kernel:
 * 20xx.xx.xx, LOC
 **/

#ifndef __INCLUDE_PNQCCORRECTION__
#define __INCLUDE_PNQCCORRECTION__
#include "pUtils.h"

int XLALSimIMREOBNonQCCorrection(COMPLEX16 *nqc, /**<< OUTPUT, The NQC correction */
                                 REAL8Vector *values,
                                 /**<< Dynamics r, phi, pr, pphi */
                                 const REAL8 omega, /**<< Angular frequency */
                                 EOBNonQCCoeffs *coeffs
                                 /**<< NQC coefficients */
);

int XLALSimIMRSpinEOBNonQCCorrection(COMPLEX16 *nqc,
                                     /**<< OUTPUT, The NQC correction */
                                     REAL8Vector *values,
                                     /**<< Dynamics r, phi, pr, pphi */
                                     const REAL8 omega, /**<< Angular frequency */
                                     const REAL8 t, const REAL8 tWind, const REAL8 wWind, EOBNonQCCoeffs *coeffs
                                     /**<< NQC coefficients */
);

int XLALSimIMRSpinEOBSANonQCCorrection(COMPLEX16 *nqc,
                                       /**<< OUTPUT, The NQC correction */
                                       REAL8Vector *values,
                                       /**<< Dynamics r, phi, pr, pphi */
                                       const REAL8 omega, /**<< Angular frequency */
                                       const REAL8 t, const REAL8 tWind, const REAL8 wWind, EOBNonQCCoeffs *coeffs
                                       /**<< NQC coefficients */
);

REAL8
XLALSimIMREOBGetNRSpinPeakDeltaTv4(INT4 l,     /**<< Mode l */
                                   INT m,      /**<< Mode m */
                                   REAL8 m1,   /**<< mass 1 */
                                   REAL8 m2,   /**<< mass 2 */
                                   REAL8 chi1, /**<< Dimensionless spin1 */
                                   REAL8 chi2, /**<< Dimensionless spin2 */
                                   HyperParams *hparams);

REAL8
XLALSimIMREOBGetNRSpinPeakAmplitudeV4(INT4 modeL, INT4 modeM, REAL8 m1, REAL8 m2, REAL8 chiS, REAL8 chiA);
REAL8
XLALSimIMREOBGetNRSpinPeakADotV4(INT4 modeL, INT4 modeM, REAL8 m1, REAL8 m2, REAL8 chiS, REAL8 chiA);
REAL8
XLALSimIMREOBGetNRSpinPeakADDotV4(INT4 modeL, INT4 modeM, REAL8 m1, REAL8 m2, REAL8 chiS, REAL8 chiA);

REAL8 XLALSimIMREOBGetNRSpinPeakOmegaV4(INT4 modeL, INT4 modeM, REAL8 eta, REAL8 a);
REAL8
XLALSimIMREOBGetNRSpinPeakOmegaDotV4(INT4 modeL, INT4 modeM, REAL8 eta, REAL8 a);

int XLALSimIMRSpinEOBCalculateNQCCoefficientsV4(REAL8Vector *amplitude,   /**<< Waveform amplitude, func of time */
                                                REAL8Vector *phase,       /**<< Waveform phase(rad), func of time */
                                                REAL8Vector *rVec,        /**<< Position-vector, function of time */
                                                REAL8Vector *prVec,       /**<< Momentum vector, function of time */
                                                REAL8Vector *orbOmegaVec, /**<< Orbital frequency, func of time */
                                                INT4 modeL,               /**<< Mode index l */
                                                INT4 modeM,               /**<< Mode index m */
                                                REAL8 timePeak,           /**<< Time of peak orbital frequency */
                                                REAL8 timeStart,          /**<< Start time */
                                                REAL8 deltaT,             /**<< Sampling interval */
                                                REAL8 m1,                 /**<< Component mass 1 */
                                                REAL8 m2,                 /**<< Component mass 2 */
                                                REAL8 chiA,  /**<< Assymmetric dimensionless spin combination */
                                                REAL8 chiS,  /**<< Symmetric dimensionless spin combination */
                                                REAL8 tWind, /**<< Location of the NQC window */
                                                REAL8 wWind, /**<< Width of the NQC window */
                                                EOBNonQCCoeffs *coeffs, /**<< OUTPUT, NQC coefficients */
                                                SpinEOBParams *ak);

#endif
