/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_NEWFACTORIZEDWAVEFORMPREC__
#define __INCLUDE_NEWFACTORIZEDWAVEFORMPREC__
#include "pUtils.h"

INT
XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV3(
            COMPLEX16 * hlm,	/**< OUTPUT, hlm waveforms */
            REAL8Vector * values,	/**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
            REAL8Vector * cartvalues,	/**< dyanmical variables */
            const REAL8 v,	/**< velocity */
            const REAL8 dr, /**< radial velocity */
            const REAL8 ncrv, /**< angular velocity */
            const REAL8 Hreal,	/**< real Hamiltonian */
            const INT l,	/**< l mode index */
            const INT m,	/**< m mode index */
            SpinEOBParams * params	/**< Spin EOB parameters */
);

INT CalculateFacWaveformAmpResV3(FacWaveformCoeffs *const hCoeffs,
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

INT XLALSimIMRSpinEOBGetAmplitudeResidualPrecV3(COMPLEX16 *rholmpwrl, 
                REAL8Vector *  values,	/**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
                REAL8Vector *  cartvalues,	/**< dyanmical variables */
                const REAL8 v, 
                const REAL8 dr,
                const REAL8 ncrv,
                const REAL8 Hreal, 
                const INT modeL, const INT modeM, 
                SpinEOBParams *params);

int
EccPrec_CalcSpinPrecFacWaveformCoefficients(
                        FacWaveformCoeffs * const coeffs,	/**< OUTPUT, pre-computed waveform coefficients */
                        const REAL8 m1,	/**< mass 1 */
                        const REAL8 m2,	/**< mass 2 */
                        const REAL8 eta,	/**< symmetric mass ratio */
                        const REAL8 tmpa,	/**< Kerr spin parameter for test-particle terms */
                        const REAL8 chiS,
                        const REAL8 chiA,
                        const REAL8 csx,	/**< (chi1+chi2)/2 */
                        const REAL8 csy,
                        const REAL8 csz,
                        const REAL8 cax,	/**< (chi1-chi2)/2 */
                        const REAL8 cay,
                        const REAL8 caz,
                        UINT SpinAlignedEOBversion	/**< 1 for SEOBNRv1; 2 for SEOBNRv2; 4 for the coefficients
                        in the flux of v4P and v4Pwave for the coefficients in the waveform of v4P */
                        );

#endif

