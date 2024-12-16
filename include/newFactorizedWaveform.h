/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_NEWFACTORIZEDWAVEFORM__
#define __INCLUDE_NEWFACTORIZEDWAVEFORM__
#include "pUtils.h"

// INT
// XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2(
//             COMPLEX16 * hlm,	/**< OUTPUT, hlm waveforms */
//             REAL8Vector * values,	/**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
//             REAL8Vector * cartvalues,	/**< dyanmical variables */
//             const REAL8 v,	/**< velocity */
//             const REAL8 dr, /**< radial velocity */
//             const REAL8 ncrv, /**< angular velocity */
//             const REAL8 Hreal,	/**< real Hamiltonian */
//             const INT l,	/**< l mode index */
//             const INT m,	/**< m mode index */
//             SpinEOBParams * params	/**< Spin EOB parameters */
// );

INT
XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2(
            COMPLEX16 * hlm,	/**< OUTPUT, hlm waveforms */
            REAL8Vector * values,	/**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
            REAL8Vector * cartvalues,	/**< dyanmical variables */
            const REAL8 v,	/**< velocity */
            const REAL8 dr, /**< radial velocity */
            const REAL8 ncrv, /**< angular velocity */
            const REAL8 prDot,
            const REAL8 Hreal,	/**< real Hamiltonian */
            const INT l,	/**< l mode index */
            const INT m,	/**< m mode index */
            SpinEOBParams * params	/**< Spin EOB parameters */
);

INT
XLALSimIMRSpinEOBGetSASpinFactorizedWaveformV2(
            COMPLEX16 * hlm,	/**< OUTPUT, hlm waveforms */
            REAL8Vector * values,	/**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
            const REAL8 v,	/**< velocity */
            const REAL8 dr, /**< radial velocity */
            const REAL8 ncrv, /**< angular velocity */
            const REAL8 prDot,
            const REAL8 Hreal,	/**< real Hamiltonian */
            const INT l,	/**< l mode index */
            const INT m,	/**< m mode index */
            SpinEOBParams * params	/**< Spin EOB parameters */
);


INT CalculateFacWaveformAmpResV2(FacWaveformCoeffs *const hCoeffs,
                                 const REAL8 eta,
                                 const INT l,
                                 const INT m,
                                 const REAL8 v,
                                 const REAL8 vPhi,
                                 const REAL8 vh3,
                                 const REAL8 eulerlogxabs,
                                 const REAL8 x0,
                                 const REAL8 x1,
                                 const REAL8 x2,
                                 COMPLEX16 *facAmpRes);

INT XLALSimIMRSpinEOBGetAmplitudeResidualPrecV2(COMPLEX16 *rholmpwrl, 
                REAL8Vector *  values,	/**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
                REAL8Vector *  cartvalues,	/**< dyanmical variables */
                const REAL8 v, 
                const REAL8 dr,
                const REAL8 ncrv,
                const REAL8 Hreal, 
                const INT modeL, const INT modeM, 
                SpinEOBParams *params);

INT XLALSimIMRSpinEOBSAGetAmplitudeResidualPrecV2(COMPLEX16 *rholmpwrl, 
                REAL8Vector *  values,	/**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
                const REAL8 v, 
                const REAL8 dr,
                const REAL8 ncrv,
                const REAL8 Hreal, 
                const INT modeL, const INT modeM, 
                SpinEOBParams *params);

#endif

