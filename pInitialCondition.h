/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_PINITIALCONDITION__
#define __INCLUDE_PINITIALCONDITION__
#include "pUtils.h"
INT EOBInitialConditions(REAL8Vector    *initConds,
                         const REAL8    mass1,
                         const REAL8    mass2,
                         const REAL8    omega0,
                         const REAL8    ecc,
                         const REAL8    inc,
                         const REAL8    spin1[],
                         const REAL8    spin2[],
                         SpinEOBParams  *params);

INT
EOBInitialConditionsPrec(
                    REAL8Vector * initConds,	/**<< OUTPUT, Initial dynamical variables */
                    const REAL8 mass1,	/**<< mass 1 */
                    const REAL8 mass2,	/**<< mass 2 */
                    const REAL8 fMin,	/**<< Initial frequency (given) */
                    const REAL8 e0,
                    const REAL8 inc,	/**<< Inclination */
                    const REAL8 spin1[],	/**<< Initial spin vector 1 */
                    const REAL8 spin2[],	/**<< Initial spin vector 2 */
                    SpinEOBParams * params	/**<< Spin EOB parameters */
);

INT
EOBInitialConditionsPrec_Conserve(
                    REAL8Vector * initConds,	/**<< OUTPUT, Initial dynamical variables */
                    const REAL8 mass1,	/**<< mass 1 */
                    const REAL8 mass2,	/**<< mass 2 */
                    const REAL8 fMin,	/**<< Initial frequency (given) */
                    const REAL8 e0,
                    const REAL8 inc,	/**<< Inclination */
                    const REAL8 spin1[],	/**<< Initial spin vector 1 */
                    const REAL8 spin2[],	/**<< Initial spin vector 2 */
                    SpinEOBParams * params	/**<< Spin EOB parameters */
);

INT EOBInitialConditionsSA_egw(REAL8Vector    *initConds,
                         const REAL8    mass1,
                         const REAL8    mass2,
                         const REAL8    fMin,
                         const REAL8    ecc,
                         const REAL8    inc,
                         const REAL8    spin1[],
                         const REAL8    spin2[],
                         SpinEOBParams  *params);

#endif

