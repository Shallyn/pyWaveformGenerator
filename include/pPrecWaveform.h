/**
 * Writer: Xiaolin.liu
 * xiaolin.liu@mail.bnu.edu.cn
 *
 * This module contains basic functions for  calculation.
 * Functions list:
 * Kernel:
 * 20xx.xx.xx, LOC
 **/

#ifndef __INCLUDE_PPRECWAVEFORM__
#define __INCLUDE_PPRECWAVEFORM__
#include "pFactorizedWaveform.h"
#include "pPrecUtils.h"

void prec_CalculateSEOBPrecWaveformVariables(SEOBPrecWaveformVariables *vars, REAL8 nchia, REAL8 nchis, REAL8 lchia,
                                             REAL8 lchis, REAL8 echia, REAL8 echis, REAL8 chi1chi1, REAL8 chi1chi2,
                                             REAL8 chi2chi2, REAL8 Jn, REAL8 Jl, REAL8 Je, REAL8 r, REAL8 prT,
                                             REAL8 prTDot);

INT prec_EOBGetPrecEccSpinFactorizedWaveform_v1(
    COMPLEX16 *hlm,                                        /**< OUTPUT, hlm waveforms */
    REAL8Vector *values,                                   /**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
    REAL8Vector *cartvalues,                               /**< dyanmical variables */
    const REAL8 v,                                         /**< velocity */
    const REAL8 Hreal,                                     /**< real Hamiltonian */
    const INT l,                                           /**< l mode index */
    const INT m,                                           /**< m mode index */
    SEOBPrecWaveformVariables *vars, SpinEOBParams *params /**< Spin EOB parameters */
);

INT prec_EOBGetPrecEccSpinFactorizedWaveform_v2(
    COMPLEX16 *hlm,                                        /**< OUTPUT, hlm waveforms */
    REAL8Vector *values,                                   /**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
    REAL8Vector *cartvalues,                               /**< dyanmical variables */
    const REAL8 v,                                         /**< velocity */
    const REAL8 Hreal,                                     /**< real Hamiltonian */
    const INT l,                                           /**< l mode index */
    const INT m,                                           /**< m mode index */
    SEOBPrecWaveformVariables *vars, SpinEOBParams *params /**< Spin EOB parameters */
);

#endif
