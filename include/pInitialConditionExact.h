/**
 * Writer: Xiaolin.liu
 * xiaolin.liu@mail.bnu.edu.cn
 *
 * This module contains basic functions for  calculation.
 * Functions list:
 * Kernel:
 * 20xx.xx.xx, LOC
 **/

#ifndef __INCLUDE_PINITIALCONDITIONEXACT__
#define __INCLUDE_PINITIALCONDITIONEXACT__
#include "pCore.h"

INT SEOBComputeExactEquatorialInitialCondition(REAL8Vector *ICvalues, REAL8 MfMin, REAL8 ecc,
                                               SpinEOBParams *seobParams);

INT EstimateEquatorialEccentricity(REAL8 r, REAL8 pphi, REAL8 *MfOrb, REAL8 *ecc, REAL8 inputfMin, REAL8 inputecc,
                                   SpinEOBParams *core);

INT EOBInitialConditionsPrec_epi(REAL8Vector *initConds,         /**<< OUTPUT, Initial dynamical variables */
                                 const REAL8 mass1,              /**<< mass 1 */
                                 const REAL8 mass2,              /**<< mass 2 */
                                 const REAL8 p0,                 /**<< Initial semi-latus */
                                 const REAL8 e0, const REAL8 x0, /**<< Inclination */
                                 const REAL8 spin1[],            /**<< Initial spin vector 1 */
                                 const REAL8 spin2[],            /**<< Initial spin vector 2 */
                                 SpinEOBParams *params           /**<< Spin EOB parameters */
);

#endif
