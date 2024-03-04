/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "myLog.h"
#include "pPrecUtils.h"
#include "pPrec.h"
#include "pPrecHam.h"
#include "pPrecWaveform.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_spline.h>


#define GSL_START \
          gsl_error_handler_t *saveGSLErrorHandler_; \
          saveGSLErrorHandler_ = gsl_set_error_handler_off();

#define GSL_END \
          gsl_set_error_handler( saveGSLErrorHandler_ );

/* We need to encapsulate the data for the GSL derivative function */
typedef struct tagPrecEulerAnglesIntegration
{
   gsl_spline *alpha_spline;
   gsl_spline *beta_spline;
   gsl_interp_accel *alpha_acc;
   gsl_interp_accel *beta_acc;
}
PrecEulerAnglesIntegration;

INT prec_CalculateFactorizedWaveformCorrection(
    COMPLEX16 * rholm,
    COMPLEX16 * flm,
    REAL8Vector *values,
    REAL8Vector *cartvalues,
    const REAL8 v,
    const REAL8 Hreal,
    const INT l,
    const INT m,
    SEOBPrecWaveformVariables *vars,
    SpinEOBParams *params
);

/**
 * Computes RHS of ODE for gamma. Eq. 10 of PRD 89, 084006 (2014)
 */
static double f_alphadotcosi( double x, void * inparams )
{
	PrecEulerAnglesIntegration* params = (PrecEulerAnglesIntegration*) inparams;

	REAL8 alphadot = gsl_spline_eval_deriv( params->alpha_spline, x, params->alpha_acc );
	REAL8 beta = gsl_spline_eval( params->beta_spline, x, params->beta_acc );

	return -1. * alphadot * cos(beta);

}

/**
 * Stopping conditions for dynamics integration for SEOBNRv4P
 */
static int PrecStopConditionBasedOnPR(double t, 
                                      const double values[],
                                      double dvalues[],
                                      void *funcParams) 
{
    int debugPK = 0;
    int debugPKverbose = 0;
    INT i;
    SpinEOBParams *params = (SpinEOBParams *)funcParams;

    REAL8 r2, pDotr = 0;
    REAL8 p[3], r[3], pdotVec[3], rdotVec[3];
    REAL8 omega, omega_xyz[3], L[3], dLdt1[3], dLdt2[3];

    memcpy(r, values, 3 * sizeof(REAL8));
    memcpy(p, values + 3, 3 * sizeof(REAL8));
    memcpy(rdotVec, dvalues, 3 * sizeof(REAL8));
    memcpy(pdotVec, dvalues + 3, 3 * sizeof(REAL8));

    r2 = inner_product3d(r, r);
    cross_product3d(values, dvalues, omega_xyz);
    omega = sqrt(inner_product3d(omega_xyz, omega_xyz)) / r2;
    pDotr = inner_product3d(p, r) / sqrt(r2);
    // if (debugPK) {
    //     XLAL_PRINT_INFO("PrecStopConditionBasedOnPR:: r = %e %e\n",
    //                     sqrt(r2), omega);
    // }
    // if (debugPK) {
    //     XLAL_PRINT_INFO(
    //         "PrecStopConditionBasedOnPR:: values = %e %e %e %e %e %e\n",
    //         values[6], values[7], values[8], values[9], values[10], values[11]);
    // }
    // if (debugPK) {
    //     XLAL_PRINT_INFO(
    //         "PrecStopConditionBasedOnPR:: dvalues = %e %e %e %e %e %e\n",
    //         dvalues[6], dvalues[7], dvalues[8], dvalues[9], dvalues[10],
    //         dvalues[11]);
    // }
    REAL8 rdot;
    // this is d(r)/dt obtained by differentiating r2 (see above)
    rdot = inner_product3d(rdotVec, r) / sqrt(r2);
    // This is d/dt(pDotr) see pDotr above.
    double prDot = -inner_product3d(p, r) * rdot / r2 +
                    inner_product3d(pdotVec, r) / sqrt(r2) +
                    inner_product3d(rdotVec, p) / sqrt(r2);

    cross_product3d(r, pdotVec, dLdt1);
    cross_product3d(rdotVec, p, dLdt2);
    cross_product3d(r, p, L);

    /* ********************************************************** */
    /* *******  Different termination conditions Follow  ******** */
    /* ********************************************************** */

    /* Table of termination conditions

        Value                   Reason
        -1                   Any of the derivatives are Nan
        0                    r < 8 and pDotr >= 0 (outspiraling)
        1                    r < 8 and rdot >= 0 (outspiraling)
        2                    r < 2 and prDot > 0 (dp_r/dt is growing)
        3                    r < 8 and |p_vec| > 10 (the momentum vector is large)
        4                    r < 8 and |p_vec| < 1e-10 (momentum vector is small)
        5                    r < 2 and omega has a another peak
        6                    r < 8 and omega < 0.04 or (r < 2. and  omega < 0.14 and
        omega has a peak) 7                    r < 8 and omega > 1 (unphysical
        omega) 8                    r < 5 and any of  |dp_i/dt| > 10 9 r < 8 and
        pphi > 10
        10                   r < 3 and rdot increases
    */

    /* Terminate if any derivative is Nan */
    for (i = 0; i < 12; i++) {
        if (isnan(dvalues[i]) || isnan(values[i])) {
        // if (debugPK) {
        //     XLAL_PRINT_INFO("\n  isnan reached. r2 = %f\n", r2);
        //     fflush(NULL);
        // }
            PRINT_LOG_INFO(LOG_CRITICAL, "nan reached at r2 = %f \n", r2);
            params->termination_reason = -1;
            return 1;
        }
    }

    /* ********************************************************** */
    /* *******  Unphysical orbital conditions  ******** */
    /* ********************************************************** */

    /* Terminate if p_r points outwards */
    if (r2 < 16 && pDotr >= 0) 
    {
        // if (debugPK) {
        // XLAL_PRINT_INFO(
        //     "\n Integration stopping, p_r pointing outwards -- out-spiraling!\n");
        // fflush(NULL);
        // }
        params->termination_reason = 0;
        return 1;
    }

    /* Terminate if rdot is >0 (OUTspiraling) for separation <4M */
    if (r2 < 16 && rdot >= 0) 
    {
        // if (debugPK) {
        // XLAL_PRINT_INFO("\n Integration stopping, dr/dt>0 -- out-spiraling!\n");
        // fflush(NULL);
        // }
        params->termination_reason = 1;
        return 1;
    }

    /* Terminate if dp_R/dt > 0, i.e. radial momentum is increasing for separation
    * <2M */
    if (r2 < 4. && prDot > 0.) 
    {
        // if (debugPK) {
        // XLAL_PRINT_INFO("\n Integration stopping as prDot = %lf at r = %lf\n",
        //                 prDot, sqrt(r2));
        // fflush(NULL);
        // }
        params->termination_reason = 2;
        return 1;
    }

    if (r2 < 16. && (sqrt(values[3] * values[3] + values[4] * values[4] +
                            values[5] * values[5]) > 10.)) 
    {
        // if (debugPK)
        // XLAL_PRINT_INFO("\n Integration stopping |pvec|> 10\n");
        fflush(NULL);
        params->termination_reason = 3;
        return 1;
    }

    if (r2 < 16. && (sqrt(values[3] * values[3] + values[4] * values[4] +
                            values[5] * values[5]) < 1.e-10)) 
    {
        // if (debugPK)
        // XLAL_PRINT_INFO("\n Integration stopping |pvec|<1e-10\n");
        fflush(NULL);
        params->termination_reason = 4;
        return 1;
    }

    /* **************************************************************** */
    /*                         Omega related                            */
    /* **************************************************************** */
    /* Terminate when omega reaches peak, and separation is < 4M */
    if (r2 < 16. && omega < params->omega)
        params->omegaPeaked = 1;

    /* If omega has gone through a second extremum, break */
    if (r2 < 4. && params->omegaPeaked == 1 &&
        omega > params->omega) 
    {
        // if (debugPK) {
        // XLAL_PRINT_INFO(
        //     "\n Integration stopping, omega reached second extremum\n");
        // fflush(NULL);
        // }
        params->termination_reason = 5;

        return 1;
    }

    /* If Momega did not evolve above 0.01 even though r < 4 or omega<0.14 for
    * r<2, break */
    if ((r2 < 16. && omega < 0.04) ||
        (r2 < 4. && omega < 0.14 && params->omegaPeaked == 1)) 
    {
        // if (debugPK) {
        // XLAL_PRINT_INFO("\n Integration stopping for omega below threshold, "
        //                 "omega=%f at r = %f\n",
        //                 omega, sqrt(r2));
        // fflush(NULL);
        // }
        params->termination_reason = 6;

        return 1;
    }

    if (r2 < 16. && omega > 1.) 
    {
        // if (debugPK) {
        // XLAL_PRINT_INFO("\n Integration stopping, omega>1 at r = %f\n", sqrt(r2));
        // fflush(NULL);
        // }
        params->termination_reason = 7;

        return 1;
    }
    params->omega = omega;

    /* **************************************************************** */
    /*              related to Numerical values of x/p/derivatives      */
    /* **************************************************************** */

    /* If momentum derivatives are too large numerically, break */
    if (r2 < 25 && (fabs(dvalues[3]) > 10 || fabs(dvalues[4]) > 10 ||
                    fabs(dvalues[5]) > 10)) 
    {
        // if (debugPK) {
        // XLAL_PRINT_INFO("\n Integration stopping, dpdt > 10 -- too large!\n");
        // fflush(NULL);
        // }
        params->termination_reason = 8;
        return 1;
    }

    /* If p_\Phi is too large numerically, break */
    if (r2 < 16. && values[5] > 10) 
    {
        // if (debugPK) {
        // XLAL_PRINT_INFO("Integration stopping, Pphi > 10 now\n\n");
        // fflush(NULL);
        // }
        params->termination_reason = 9;
        return 1;
    }
    /* If rdot inclreases, break */
    if (r2 < 9. && rdot > params->prev_dr) 
    {
        // if (debugPK) {
        // XLAL_PRINT_INFO("\n Integration stopping, dr/dt increasing!\n");
        // fflush(NULL);
        // }
        params->prev_dr = rdot;
        params->termination_reason = 10;

        return 1;
    }
    params->prev_dr = rdot;

    /* **************************************************************** */
    /*              Last resort conditions                              */
    /* **************************************************************** */

    /* Very verbose output */
    // if (debugPKverbose && r2 < 16.) {
    //     XLAL_PRINT_INFO("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", t,
    //                     values[0], values[1], values[2], values[3], values[4],
    //                     values[5], values[6], values[7], values[8], values[9],
    //                     values[10], values[11], values[12], values[13], omega);
    // }
    return GSL_SUCCESS;
}

static int PrecStopConditionBasedOnPR_withfmax(double t, 
                                      const double values[],
                                      double dvalues[],
                                      void *funcParams) 
{
    int debugPK = 0;
    int debugPKverbose = 0;
    INT i;
    SpinEOBParams *params = (SpinEOBParams *)funcParams;

    REAL8 r2, pDotr = 0;
    REAL8 p[3], r[3], pdotVec[3], rdotVec[3];
    REAL8 omega, omega_xyz[3], L[3], dLdt1[3], dLdt2[3];

    memcpy(r, values, 3 * sizeof(REAL8));
    memcpy(p, values + 3, 3 * sizeof(REAL8));
    memcpy(rdotVec, dvalues, 3 * sizeof(REAL8));
    memcpy(pdotVec, dvalues + 3, 3 * sizeof(REAL8));

    r2 = inner_product3d(r, r);
    cross_product3d(values, dvalues, omega_xyz);
    omega = sqrt(inner_product3d(omega_xyz, omega_xyz)) / r2;
    pDotr = inner_product3d(p, r) / sqrt(r2);
    REAL8 omegaMax = params->hParams->Mf_max*CST_PI;
    if (omega > omegaMax)
        return 1;
    // if (debugPK) {
    //     XLAL_PRINT_INFO("PrecStopConditionBasedOnPR:: r = %e %e\n",
    //                     sqrt(r2), omega);
    // }
    // if (debugPK) {
    //     XLAL_PRINT_INFO(
    //         "PrecStopConditionBasedOnPR:: values = %e %e %e %e %e %e\n",
    //         values[6], values[7], values[8], values[9], values[10], values[11]);
    // }
    // if (debugPK) {
    //     XLAL_PRINT_INFO(
    //         "PrecStopConditionBasedOnPR:: dvalues = %e %e %e %e %e %e\n",
    //         dvalues[6], dvalues[7], dvalues[8], dvalues[9], dvalues[10],
    //         dvalues[11]);
    // }
    REAL8 rdot;
    // this is d(r)/dt obtained by differentiating r2 (see above)
    rdot = inner_product3d(rdotVec, r) / sqrt(r2);
    // This is d/dt(pDotr) see pDotr above.
    double prDot = -inner_product3d(p, r) * rdot / r2 +
                    inner_product3d(pdotVec, r) / sqrt(r2) +
                    inner_product3d(rdotVec, p) / sqrt(r2);

    cross_product3d(r, pdotVec, dLdt1);
    cross_product3d(rdotVec, p, dLdt2);
    cross_product3d(r, p, L);

    /* ********************************************************** */
    /* *******  Different termination conditions Follow  ******** */
    /* ********************************************************** */

    /* Table of termination conditions

        Value                   Reason
        -1                   Any of the derivatives are Nan
        0                    r < 8 and pDotr >= 0 (outspiraling)
        1                    r < 8 and rdot >= 0 (outspiraling)
        2                    r < 2 and prDot > 0 (dp_r/dt is growing)
        3                    r < 8 and |p_vec| > 10 (the momentum vector is large)
        4                    r < 8 and |p_vec| < 1e-10 (momentum vector is small)
        5                    r < 2 and omega has a another peak
        6                    r < 8 and omega < 0.04 or (r < 2. and  omega < 0.14 and
        omega has a peak) 7                    r < 8 and omega > 1 (unphysical
        omega) 8                    r < 5 and any of  |dp_i/dt| > 10 9 r < 8 and
        pphi > 10
        10                   r < 3 and rdot increases
    */

    /* Terminate if any derivative is Nan */
    for (i = 0; i < 12; i++) {
        if (isnan(dvalues[i]) || isnan(values[i])) {
        // if (debugPK) {
        //     XLAL_PRINT_INFO("\n  isnan reached. r2 = %f\n", r2);
        //     fflush(NULL);
        // }
            PRINT_LOG_INFO(LOG_CRITICAL, "nan reached at r2 = %f \n", r2);
            params->termination_reason = -1;
            return 1;
        }
    }

    /* ********************************************************** */
    /* *******  Unphysical orbital conditions  ******** */
    /* ********************************************************** */

    /* Terminate if p_r points outwards */
    if (r2 < 16 && pDotr >= 0) 
    {
        // if (debugPK) {
        // XLAL_PRINT_INFO(
        //     "\n Integration stopping, p_r pointing outwards -- out-spiraling!\n");
        // fflush(NULL);
        // }
        params->termination_reason = 0;
        return 1;
    }

    /* Terminate if rdot is >0 (OUTspiraling) for separation <4M */
    if (r2 < 16 && rdot >= 0) 
    {
        // if (debugPK) {
        // XLAL_PRINT_INFO("\n Integration stopping, dr/dt>0 -- out-spiraling!\n");
        // fflush(NULL);
        // }
        params->termination_reason = 1;
        return 1;
    }

    /* Terminate if dp_R/dt > 0, i.e. radial momentum is increasing for separation
    * <2M */
    if (r2 < 4. && prDot > 0.) 
    {
        // if (debugPK) {
        // XLAL_PRINT_INFO("\n Integration stopping as prDot = %lf at r = %lf\n",
        //                 prDot, sqrt(r2));
        // fflush(NULL);
        // }
        params->termination_reason = 2;
        return 1;
    }

    if (r2 < 16. && (sqrt(values[3] * values[3] + values[4] * values[4] +
                            values[5] * values[5]) > 10.)) 
    {
        // if (debugPK)
        // XLAL_PRINT_INFO("\n Integration stopping |pvec|> 10\n");
        fflush(NULL);
        params->termination_reason = 3;
        return 1;
    }

    if (r2 < 16. && (sqrt(values[3] * values[3] + values[4] * values[4] +
                            values[5] * values[5]) < 1.e-10)) 
    {
        // if (debugPK)
        // XLAL_PRINT_INFO("\n Integration stopping |pvec|<1e-10\n");
        fflush(NULL);
        params->termination_reason = 4;
        return 1;
    }

    /* **************************************************************** */
    /*                         Omega related                            */
    /* **************************************************************** */
    /* Terminate when omega reaches peak, and separation is < 4M */
    if (r2 < 16. && omega < params->omega)
        params->omegaPeaked = 1;

    /* If omega has gone through a second extremum, break */
    if (r2 < 4. && params->omegaPeaked == 1 &&
        omega > params->omega) 
    {
        // if (debugPK) {
        // XLAL_PRINT_INFO(
        //     "\n Integration stopping, omega reached second extremum\n");
        // fflush(NULL);
        // }
        params->termination_reason = 5;

        return 1;
    }

    /* If Momega did not evolve above 0.01 even though r < 4 or omega<0.14 for
    * r<2, break */
    if ((r2 < 16. && omega < 0.04) ||
        (r2 < 4. && omega < 0.14 && params->omegaPeaked == 1)) 
    {
        // if (debugPK) {
        // XLAL_PRINT_INFO("\n Integration stopping for omega below threshold, "
        //                 "omega=%f at r = %f\n",
        //                 omega, sqrt(r2));
        // fflush(NULL);
        // }
        params->termination_reason = 6;

        return 1;
    }

    if (r2 < 16. && omega > 1.) 
    {
        // if (debugPK) {
        // XLAL_PRINT_INFO("\n Integration stopping, omega>1 at r = %f\n", sqrt(r2));
        // fflush(NULL);
        // }
        params->termination_reason = 7;

        return 1;
    }
    params->omega = omega;

    /* **************************************************************** */
    /*              related to Numerical values of x/p/derivatives      */
    /* **************************************************************** */

    /* If momentum derivatives are too large numerically, break */
    if (r2 < 25 && (fabs(dvalues[3]) > 10 || fabs(dvalues[4]) > 10 ||
                    fabs(dvalues[5]) > 10)) 
    {
        // if (debugPK) {
        // XLAL_PRINT_INFO("\n Integration stopping, dpdt > 10 -- too large!\n");
        // fflush(NULL);
        // }
        params->termination_reason = 8;
        return 1;
    }

    /* If p_\Phi is too large numerically, break */
    if (r2 < 16. && values[5] > 10) 
    {
        // if (debugPK) {
        // XLAL_PRINT_INFO("Integration stopping, Pphi > 10 now\n\n");
        // fflush(NULL);
        // }
        params->termination_reason = 9;
        return 1;
    }
    /* If rdot inclreases, break */
    if (r2 < 9. && rdot > params->prev_dr) 
    {
        // if (debugPK) {
        // XLAL_PRINT_INFO("\n Integration stopping, dr/dt increasing!\n");
        // fflush(NULL);
        // }
        params->prev_dr = rdot;
        params->termination_reason = 10;

        return 1;
    }
    params->prev_dr = rdot;

    /* **************************************************************** */
    /*              Last resort conditions                              */
    /* **************************************************************** */

    /* Very verbose output */
    // if (debugPKverbose && r2 < 16.) {
    //     XLAL_PRINT_INFO("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", t,
    //                     values[0], values[1], values[2], values[3], values[4],
    //                     values[5], values[6], values[7], values[8], values[9],
    //                     values[10], values[11], values[12], values[13], omega);
    // }
    return GSL_SUCCESS;
}

static int PrecStopConditionBasedOnPR_inverse(double t, 
                                      const double values[],
                                      double dvalues[],
                                      void *funcParams) 
{
    int debugPK = 0;
    int debugPKverbose = 0;
    INT i;
    SpinEOBParams *params = (SpinEOBParams *)funcParams;

    REAL8 r2, pDotr = 0;
    REAL8 p[3], r[3], pdotVec[3], rdotVec[3];
    REAL8 omega, omega_xyz[3], L[3], dLdt1[3], dLdt2[3];

    memcpy(r, values, 3 * sizeof(REAL8));
    memcpy(p, values + 3, 3 * sizeof(REAL8));
    memcpy(rdotVec, dvalues, 3 * sizeof(REAL8));
    memcpy(pdotVec, dvalues + 3, 3 * sizeof(REAL8));

    r2 = inner_product3d(r, r);
    cross_product3d(values, dvalues, omega_xyz);
    omega = sqrt(inner_product3d(omega_xyz, omega_xyz)) / r2;
    REAL8 omega0 = CST_PI * params->hParams->Mf_min;
    // print_debug("omega, omega0 = %.16e, %.16e\n", omega, omega0);
    if (omega < omega0)
        return 1;
    return GSL_SUCCESS;
}


void SortSEOBDynamicsArrayPrec(REAL8Array **ret_dyn, INT length, REAL8Array *dyn_inv)
{
    INT i,j;
    REAL8Array *dyn = CreateREAL8Array(2, 21, length);
    for (i=0; i<length; i++)
    {
        dyn->data[i] = -dyn_inv->data[length - 1 - i];
        for (j=1; j<21; j++)
            dyn->data[i + j*length] = dyn_inv->data[length - 1 - i + j*length];
    }
    //     FILE *out = fopen("debug_inverse_Cons1.dat" ,"w");
    // for (i=0; i<length; i++)
    // {
    //     fprintf(out, "%.16e\t%.16e\t%.16e\t%.16e\n", 
    //         dyn->data[i], dyn->data[i + length], dyn->data[i + 2*length],
    //         dyn->data[i + 3*length]);
    // }
    //     fclose(out);
    *ret_dyn = dyn;
    return;
}

void SEOBConcactInverseDynToAdaSDynPrec(REAL8Array **dyn_out, REAL8Array *dyn_inv, 
        INT *retLen_out, INT retLen_inv)
{
    REAL8Array *dyn_adas = *dyn_out;
    INT retLenAdaS = (*retLen_out);
    INT i, retLen = retLenAdaS + retLen_inv-1;
    REAL8Array *dyn_conc = CreateREAL8Array(2, 21, retLen);
    *retLen_out = retLen;
    for (i=0; i<21; i++)
    {
        memcpy(dyn_conc->data + i*retLen, dyn_inv->data + i*retLen_inv, (retLen_inv-1)*sizeof(REAL8));
        memcpy(dyn_conc->data+i*retLen+retLen_inv-1, dyn_adas->data + i*retLenAdaS, retLenAdaS*sizeof(REAL8));
    }
    // FILE *out = fopen("debug_inverse_Cons.dat" ,"w");
    REAL8 t0 = dyn_conc->data[0];
    for (i=0; i<retLen; i++)
    {
        dyn_conc->data[i] = dyn_conc->data[i] - t0;
        // fprintf(out, "%.16e\t%.16e\t%.16e\t%.16e\n", 
        //     dyn_conc->data[i], dyn_conc->data[i+retLen], 
        //     dyn_conc->data[i+2*retLen], dyn_conc->data[i+3*retLen]);
    }
    // fclose(out);
    STRUCTFREE(dyn_adas, REAL8Array);
    *dyn_out = dyn_conc;
    return;
}

INT SEOBIntegrateDynamics_prec_inverse(REAL8Array **dynamics,
                          INT *retLenOut,
                          REAL8Vector *ICvalues,
                          REAL8 EPS_ABS,
                          REAL8 EPS_REL,
                          REAL8 deltaT,
                          REAL8 deltaT_min,
                          REAL8 tstart,
                          REAL8 tend ,
                          SpinEOBParams *seobParams,
                          INT flagConstantSampling)
{
    PRINT_LOG_INFO(LOG_INFO, "Setting Integrator");
    INT retLen;
    UINT i;
    REAL8Array *dynamics_spinaligned = NULL;
    REAL8Array *dynamics_inverse = NULL;

    INT status, failed = 0;
    /* Dimensions of vectors of dynamical variables to be integrated */
    UINT nb_Hamiltonian_variables = 14;

    REAL8Vector *values = CreateREAL8Vector(nb_Hamiltonian_variables);
    if (!values) {failed = 1; goto QUIT;}
    memcpy(values->data, ICvalues->data, values->length * sizeof(REAL8));

    ARKIntegrator *integrator = NULL;
    integrator = XLALAdaptiveRungeKutta4Init(
        nb_Hamiltonian_variables, PrecHcapNumericalDerivative_inverse,
        PrecStopConditionBasedOnPR_inverse, EPS_ABS, EPS_REL);

    if (!integrator) {failed = 1; goto QUIT;}

    /* Ensure that integration stops ONLY when the stopping condition is True */
    integrator->stopontestonly = 1;
    /* When this option is set to 0, the integration can be exceedingly slow for
    * spin-aligned systems */
    integrator->retries = 1;
    /* Computing the dynamical evolution of the system */
    // Prec
    retLen = XLALAdaptiveRungeKutta4NoInterpolateWithDerivPrec(integrator, seobParams, 
        values->data, 0., tend-tstart, 
        deltaT, deltaT_min, &dynamics_inverse);
    if (retLen < 0)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Integration Failed!!!");
        failed = 1;
        goto QUIT;
    }
    PRINT_LOG_INFO(LOG_INFO, "Integration End");
    // sort dynamics
    SortSEOBDynamicsArrayPrec(dynamics, retLen, dynamics_inverse);
    STRUCTFREE(dynamics_inverse, REAL8Array);

    // NOTE: functions like XLALAdaptiveRungeKutta4 would give nans if the times
    // do not start at 0 -- we have to adjust the starting time after integration
    /* Adjust starting time */
    for ( i = 0; i < retLen; i++)
        (*dynamics)->data[i] += tstart;

QUIT:
    STRUCTFREE(values, REAL8Vector);
    STRUCTFREE(integrator, ARKIntegrator);
    *retLenOut = retLen;
    if (failed)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}

INT SEOBIntegrateDynamics_prec_withFMax(REAL8Array **dynamics,
                          INT *retLenOut,
                          REAL8Vector *ICvalues,
                          REAL8 EPS_ABS,
                          REAL8 EPS_REL,
                          REAL8 deltaT,
                          REAL8 deltaT_min,
                          REAL8 tstart,
                          REAL8 tend ,
                          SpinEOBParams *seobParams,
                          INT flagConstantSampling)
{
    INT retLen;
    UINT i;
    REAL8Array *dynamics_spinaligned = NULL;
    INT status, failed = 0;
    /* Dimensions of vectors of dynamical variables to be integrated */
    UINT nb_Hamiltonian_variables = 14;

    REAL8Vector *values = CreateREAL8Vector(nb_Hamiltonian_variables);
    if (!values) {failed = 1; goto QUIT;}
    memcpy(values->data, ICvalues->data, values->length * sizeof(REAL8));

    ARKIntegrator *integrator = NULL;
    // Prec
    if (tstart > 0) 
    {
        integrator = XLALAdaptiveRungeKutta4Init(
            nb_Hamiltonian_variables, PrecHcapNumericalDerivative,
            PrecStopConditionBasedOnPR_withfmax, EPS_ABS, EPS_REL);
    } else {
        integrator = XLALAdaptiveRungeKutta4Init(
            nb_Hamiltonian_variables, PrecHcapNumericalDerivative,
            PrecStopConditionBasedOnPR_withfmax, EPS_ABS, EPS_REL);
    }

    if (!integrator) {failed = 1; goto QUIT;}

    /* Ensure that integration stops ONLY when the stopping condition is True */
    integrator->stopontestonly = 1;
    /* When this option is set to 0, the integration can be exceedingly slow for
    * spin-aligned systems */
    integrator->retries = 1;
    /* Computing the dynamical evolution of the system */
    // Prec
    if (!flagConstantSampling)
        retLen = XLALAdaptiveRungeKutta4NoInterpolateWithDerivPrec(integrator, seobParams, 
            values->data, 0., tend-tstart, 
            deltaT, deltaT_min, dynamics);
    else
        retLen = XLALAdaptiveRungeKutta4WithDerivPrec(integrator, seobParams, 
            values->data, 0., tend-tstart, 
            deltaT, deltaT_min, dynamics);
    if (retLen < 0)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Integration Failed!!!");
        failed = 1;
        goto QUIT;
    }
    PRINT_LOG_INFO(LOG_INFO, "Integration End");
    // NOTE: functions like XLALAdaptiveRungeKutta4 would give nans if the times
    // do not start at 0 -- we have to adjust the starting time after integration
    /* Adjust starting time */
    for ( i = 0; i < retLen; i++)
        (*dynamics)->data[i] += tstart;

QUIT:
    STRUCTFREE(values, REAL8Vector);
    STRUCTFREE(integrator, ARKIntegrator);
    *retLenOut = retLen;
    if (failed)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}

INT SEOBIntegrateDynamics_prec(REAL8Array **dynamics,
                          INT *retLenOut,
                          REAL8Vector *ICvalues,
                          REAL8 EPS_ABS,
                          REAL8 EPS_REL,
                          REAL8 deltaT,
                          REAL8 deltaT_min,
                          REAL8 tstart,
                          REAL8 tend ,
                          SpinEOBParams *seobParams,
                          INT flagConstantSampling)
{
    INT retLen;
    UINT i;
    REAL8Array *dynamics_spinaligned = NULL;
    INT status, failed = 0;
    /* Dimensions of vectors of dynamical variables to be integrated */
    UINT nb_Hamiltonian_variables = 14;

    REAL8Vector *values = CreateREAL8Vector(nb_Hamiltonian_variables);
    if (!values) {failed = 1; goto QUIT;}
    memcpy(values->data, ICvalues->data, values->length * sizeof(REAL8));

    ARKIntegrator *integrator = NULL;
    // Prec
    if (tstart > 0) 
    {
        integrator = XLALAdaptiveRungeKutta4Init(
            nb_Hamiltonian_variables, PrecHcapNumericalDerivative,
            PrecStopConditionBasedOnPR, EPS_ABS, EPS_REL);
    } else {
        integrator = XLALAdaptiveRungeKutta4Init(
            nb_Hamiltonian_variables, PrecHcapNumericalDerivative,
            PrecStopConditionBasedOnPR, EPS_ABS, EPS_REL);
    }

    if (!integrator) {failed = 1; goto QUIT;}

    /* Ensure that integration stops ONLY when the stopping condition is True */
    integrator->stopontestonly = 1;
    /* When this option is set to 0, the integration can be exceedingly slow for
    * spin-aligned systems */
    integrator->retries = 1;
    /* Computing the dynamical evolution of the system */
    // Prec
    if (!flagConstantSampling)
        retLen = XLALAdaptiveRungeKutta4NoInterpolateWithDerivPrec(integrator, seobParams, 
            values->data, 0., tend-tstart, 
            deltaT, deltaT_min, dynamics);
    else
        retLen = XLALAdaptiveRungeKutta4WithDerivPrec(integrator, seobParams, 
            values->data, 0., tend-tstart, 
            deltaT, deltaT_min, dynamics);
    if (retLen < 0)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Integration Failed!!!");
        failed = 1;
        goto QUIT;
    }
    PRINT_LOG_INFO(LOG_INFO, "Integration End");
    // NOTE: functions like XLALAdaptiveRungeKutta4 would give nans if the times
    // do not start at 0 -- we have to adjust the starting time after integration
    /* Adjust starting time */
    for ( i = 0; i < retLen; i++)
        (*dynamics)->data[i] += tstart;

QUIT:
    STRUCTFREE(values, REAL8Vector);
    STRUCTFREE(integrator, ARKIntegrator);
    *retLenOut = retLen;
    if (failed)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}

INT CutSEOBPrecdynamics(SEOBPrecdynamics **eobdyn, REAL8 MfMin)
{
    INT i, len_old, len_new;
    REAL8 omega0 = CST_PI*MfMin;
    SEOBPrecdynamics *dyn_old = *eobdyn;
    len_old = dyn_old->length;
    for(i=0; i<len_old; i++)
    {
        // print_debug("omega = %.16e\n", dyn_old->omegaVec[i]);
        if (dyn_old->omegaVec[i] >= omega0)
            break;
    }
    if (i==0) return CEV_SUCCESS;
    if(i>0) i--;
    len_new = len_old - i;
    if (len_new == 0)
        return CEV_FAILURE;
    SEOBPrecdynamics *dyn_new = CreateSEOBPrecdynamics(len_new);
// print_debug("here: (%d, %d)\n", dyn_new->array->dimLength->data[0], dyn_new->array->dimLength->data[1]);
    for (int j=0; j<dyn_new->array->dimLength->data[0]; j++)
        memcpy(dyn_new->array->data + j*len_new, dyn_old->array->data + i + j*len_old, len_new*(sizeof(REAL8)));
    REAL8 t0 = dyn_new->tVec[0];
    for (int j=0; j<len_new; j++)
        dyn_new->tVec[j] = dyn_new->tVec[j] - t0;
    STRUCTFREE(dyn_old, SEOBPrecdynamics);
    *eobdyn = dyn_new;
    return CEV_SUCCESS;
}

INT SEOBComputeExtendedSEOBPrecdynamics(SEOBPrecdynamics **seobdynamics,
                                    REAL8Array *dynamics,
                                    INT retLen,
                                    SpinEOBParams *seobParams)
{
    SEOBPrecdynamics *seobdyn = NULL;
    INT status, failed = 0;
    UINT i, j;
    seobdyn = CreateSEOBPrecdynamics(retLen);
    // print_debug("retLen = %d, size = %zu, col = %zu\n", retLen, seobdyn->array->size, seobdyn->array->size/retLen);
    // print_debug("sizeof dyn = (%zu, %zu)\n", dynamics->dimLength->data[0], dynamics->dimLength->data[1]);
    memcpy(seobdyn->array->data, dynamics->data, v4PrecEvolvedynamicsVariables * retLen * sizeof(REAL8));
    REAL8 rVec[3] = {0, 0, 0};
    REAL8 pTVec[3] = {0, 0, 0};
    REAL8 vVec[3] = {0, 0, 0};
    REAL8 rcrossv[3] = {0, 0, 0};
    REAL8 rcrossp[3] = {0, 0, 0};
    // REAL8 LNhat[3] = {0, 0, 0};
    REAL8 chi1Vec[3] = {0, 0, 0};
    REAL8 chi2Vec[3] = {0, 0, 0};
    REAL8 chiaVec[3] = {0, 0, 0};
    REAL8 chisVec[3] = {0, 0, 0};
    REAL8 ehat[3] = {0, 0, 0};
    REAL8 nhat[3] = {0, 0, 0};
    REAL8 lhat[3] = {0, 0, 0};
    REAL8 eta = seobParams->eta;
    REAL8 dm, s1NormFac, s2NormFac;
    dm = 1. - 4.*eta;
    if (dm < 0.0)
        dm = 0.0;
    else
        dm = sqrt(dm);
    s1NormFac = 0.25*(1. + dm)*(1. + dm);
    s2NormFac = 0.25*(1. - dm)*(1. - dm);
    for (i=0; i<retLen; i++)
    {
        rVec[0] = seobdyn->posVecx[i];
        rVec[1] = seobdyn->posVecy[i];
        rVec[2] = seobdyn->posVecz[i];

        pTVec[0] = seobdyn->momTVecx[i];
        pTVec[1] = seobdyn->momTVecy[i];
        pTVec[2] = seobdyn->momTVecz[i];

        vVec[0] = seobdyn->velVecx[i];
        vVec[1] = seobdyn->velVecy[i];
        vVec[2] = seobdyn->velVecz[i];

        chi1Vec[0] = seobdyn->s1Vecx[i] / s1NormFac;
        chi1Vec[1] = seobdyn->s1Vecy[i] / s1NormFac;
        chi1Vec[2] = seobdyn->s1Vecz[i] / s1NormFac;

        chi2Vec[0] = seobdyn->s2Vecx[i] / s2NormFac;
        chi2Vec[1] = seobdyn->s2Vecy[i] / s2NormFac;
        chi2Vec[2] = seobdyn->s2Vecz[i] / s2NormFac;

        REAL8 polarr;
        seobdyn->polarrVec[i] = polarr = sqrt(inner_product3d(rVec, rVec));
        cross_product3d(rVec, pTVec, rcrossp);
        cross_product3d(rVec, vVec, rcrossv);
        REAL8 rcrossvNorm = sqrt(inner_product3d(rcrossv, rcrossv));
        REAL8 magL = sqrt(inner_product3d(rcrossp, rcrossp));
        for (j = 0; j < 3; j++)
        {
            if (PREC_FLAG == 2)
                ehat[j] = rcrossp[j] / magL;
            else
                ehat[j] = rcrossv[j] / rcrossvNorm;
            nhat[j] = rVec[j] / polarr;
            chisVec[j] = 0.5*(chi1Vec[j] + chi2Vec[j]);
            chiaVec[j] = 0.5*(chi1Vec[j] - chi2Vec[j]);
        }
        // cross_product3d(LNhat, nhat, lhat);
        cross_product3d(ehat, nhat, lhat);
        seobdyn->polarphiVec[i] = seobdyn->phiDMod[i] + seobdyn->phiMod[i];
        seobdyn->polarprTVec[i] = inner_product3d(rVec, pTVec) / seobdyn->polarrVec[i];
        seobdyn->polarpphiVec[i] = magL;
        seobdyn->omegaVec[i] = rcrossvNorm / (polarr * polarr);

        seobdyn->nchiaVec[i] = inner_product3d(nhat, chiaVec);
        seobdyn->nchisVec[i] = inner_product3d(nhat, chisVec);
        seobdyn->lchiaVec[i] = inner_product3d(lhat, chiaVec);
        seobdyn->lchisVec[i] = inner_product3d(lhat, chisVec);
        seobdyn->echiaVec[i] = inner_product3d(ehat, chiaVec);
        seobdyn->echisVec[i] = inner_product3d(ehat, chisVec);

        seobdyn->chi1chi1[i] = inner_product3d(chi1Vec, chi1Vec);
        seobdyn->chi1chi2[i] = inner_product3d(chi1Vec, chi2Vec);
        seobdyn->chi2chi2[i] = inner_product3d(chi2Vec, chi2Vec);

        seobdyn->JnVec[i] = inner_product3d(nhat, seobParams->J0Vec->data);
        seobdyn->JlVec[i] = inner_product3d(lhat, seobParams->J0Vec->data);
        seobdyn->JeVec[i] = inner_product3d(ehat, seobParams->J0Vec->data);

        seobdyn->s1dotZVec[i] = (rcrossp[0]*seobdyn->s1Vecx[i] + rcrossp[1]*seobdyn->s1Vecy[i] + rcrossp[2]*seobdyn->s1Vecz[i]) / magL;
        seobdyn->s2dotZVec[i] = (rcrossp[0]*seobdyn->s2Vecx[i] + rcrossp[1]*seobdyn->s2Vecy[i] + rcrossp[2]*seobdyn->s2Vecz[i]) / magL;
        // Debug:
        // print_debug("chi2 dot chi2 = %.16e\n", inner_product3d(chi2Vec, chi2Vec));
    }
    *seobdynamics = seobdyn;
    return CEV_SUCCESS;
}

void OrbitalPhaseReducePrec(SEOBPrecdynamics *dyn, REAL8 phiD, REAL8 phiM, REAL8 phi)
{
    INT i, length = dyn->length;
    for (i=0; i<length; i++)
    {
        dyn->phiDMod[i] -= phiD;
        dyn->phiMod[i] -= phiM;
        dyn->polarphiVec[i] -= phi;
    }
    return;
}

INT SetZeroPhaseAtTimePrec(SEOBPrecdynamics *dyn, REAL8 t, REAL8 *ret_dphiD, REAL8 *ret_dphiM, REAL8 *ret_dphi)
{
    REAL8 phiD, phiM, phi;
    REAL8Vector *dynVec = NULL;
    INT status;
    status = SEOBInterpolatePrecDynamicsAtTime(&dynVec, t, dyn);
    if (status != CEV_SUCCESS)
        return CEV_FAILURE;
    *ret_dphiD = phiD = dynVec->data[13];
    *ret_dphiM = phiM = dynVec->data[14];
    *ret_dphi = phi = dynVec->data[22];
    OrbitalPhaseReducePrec(dyn, phiD, phiM, phi);
    STRUCTFREE(dynVec, REAL8Vector);
    return CEV_SUCCESS;
}

int SEOBInterpolatePrecDynamicsAtTime(
    REAL8Vector **seobdynamics_values, /**<< Output: pointer to vector for
                                          seobdynamics interpolated values */
    REAL8 t,                           /**<< Input: time at which to evaluate */
    SEOBPrecdynamics *seobdynamics         /**<< Input: SEOB dynamics */
) 
{
    UINT j;
    int is_failed = 0;
    /* Create output vector */
    if (!((*seobdynamics_values) = CreateREAL8Vector(v4PrecdynamicsVariables))) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failed to allocate REAL8Vector ");
        return CEV_FAILURE;
    }
    memset((*seobdynamics_values)->data, 0, ((*seobdynamics_values)->length) * sizeof(REAL8));

    /* Check that the time asked for is in range */

    UINT retLen = seobdynamics->length;
    REAL8 *tVec = seobdynamics->tVec;
    if ((t < tVec[0]) || (t > tVec[retLen - 1])) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "time %f for interpolation is out of range of (%f, %f)", t, tVec[0], tVec[retLen - 1]);
        return CEV_FAILURE;;
    }

    /* Get the start and end indices that we will use to interpolate */
    /* indext max index such that tVec[indext] <= t */
    UINT indext = 0;
    while ((indext < retLen - 1) && (tVec[indext + 1] <= t))
        indext++;
    INT4 indexstart = indext - 20 > 0 ? indext - 20 : 0;
    INT4 indexend = indext + 20 < retLen - 1 ? indext + 20 : retLen - 1;
    INT4 interp_length = indexend - indexstart + 1;
    if (interp_length <= 0) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "not finding a strictly positive number of ");
        return CEV_FAILURE;;
    }
    GSL_START;
    /* Spline allocation */
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, interp_length);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    int status;
    /* Interpolate all quantities */
    (*seobdynamics_values)->data[0] = t;
    for ( j = 1; j < v4PrecdynamicsVariables; j++) {
        status = gsl_spline_init(spline, &(tVec[indexstart]),
                        &(seobdynamics->array->data[j * retLen + indexstart]),
                        interp_length);
        if (status != GSL_SUCCESS)
        {
            is_failed = 1;
            goto QUIT;
        }
        (*seobdynamics_values)->data[j] = gsl_spline_eval(spline, t, acc);
    }

QUIT:
    /* Cleanup */
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    GSL_END;
    if (is_failed)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}


int SEOBPrecLocateTimePeakOmega(
    REAL8 *tPeakOmega, /**<< Output: time of peak of Omega if found (see inside
        XLALSimLocateOmegaTime for what is returned otherwise) */
    INT *foundPeakOmega, /**<< Output: flag indicating wether tPeakOmega has been found */
    REAL8Array *dynamics,      /**<< Input: array for dynamics */
    SEOBPrecdynamics *seobdynamics,       /**<< Input: SEOB dynamics object */
    UINT retLen,                     /**<< Input: length of dynamics */
    SpinEOBParams *seobParams /**<< SEOB params */
)
{ 
    REAL8Vector tVec;
    tVec.length = retLen;
    tVec.data = seobdynamics->tVec;
    REAL8Vector omegaVec;
    omegaVec.length = retLen;
    omegaVec.data = seobdynamics->omegaVec;
    XLALEOBFindRobustPeak(tPeakOmega, &tVec, &omegaVec, 3);
    *foundPeakOmega = 1;
    return CEV_SUCCESS;
}


INT SEOBPrecCalculateNQCWindowFactorsFromDyn(SEOBPrecdynamics *dyn,
                                         REAL8 thPeak,
                                         REAL8 tr6M,
                                         REAL8 tHiStart,
                                         REAL8 tThresh,
                                         INT is_first,
                                         SpinEOBParams *seobParams)
{
    REAL8 tp1, tp2;
    REAL8 pr_n, pr_o;
    INT i, ip1, ip2, count = 0;
    ip1 = ip2 = 0;
    REAL8 dr_n, dr_o;
    // find_dy_indpr0(dy, 6, ip1, ip2)
    for (i = dyn->length - 2; i>=0; i--)
    {
        if (dyn->polarrVec[i] < tThresh) continue;
        if (is_first)
        {
            dr_n = (dyn->posVecx[i]*dyn->velVecx[i] + 
                dyn->posVecy[i]*dyn->velVecy[i] + 
                dyn->posVecz[i]*dyn->velVecz[i]) / dyn->polarrVec[i];
            dr_o = (dyn->posVecx[i+1]*dyn->velVecx[i+1] + 
                dyn->posVecy[i+1]*dyn->velVecy[i+1] + 
                dyn->posVecz[i+1]*dyn->velVecz[i+1]) / dyn->polarrVec[i+1];
            if (dr_n * dr_o < 0)
            {
                ip1 = i;
                ip2 = i;
                break;
            }
            continue;
        }
        pr_n = dyn->polarprTVec[i];
        pr_o = dyn->polarprTVec[i+1];
        if (pr_n * pr_o < 0)
        {
            if (!count)
            {
                ip1 = i;
                count++;
                continue;
            }
            else
            {
                ip2 = i;
                break;
            }
        }
    }
    tp1 = dyn->tVec[ip1];
    tp2 = dyn->tVec[ip2];
    if (is_first)
        tp1 = (tp1 + thPeak)/2.;
    seobParams->tWind = tp1;
    seobParams->wWind = log(2.)*50 / GET_MAX(fabs(thPeak - tHiStart), fabs(tr6M - tp1));
    // print_debug("ip1 = %d, ip2 = %d\n", ip1, ip2);
    // print_debug("tp1 = %g, tp2 = %g, thPeak = %g, tHi0 = %g, tr6 = %g, rp = %g\n",
    //     tp1-thPeak, tp2-thPeak, thPeak, tHiStart, tr6M, dyn->polarrVec[ip1]);
    PRINT_LOG_INFO(LOG_DEBUG, "tWind = %.16e, wWind = %g\n", seobParams->tWind, seobParams->wWind);
    return CEV_SUCCESS;
}

static INT XLALSimIMREOBPrecCalcCalibCoefficientHigherModesPrec (
               SpinEOBParams * params, /**Output **/
               const UINT modeL, /*<< Mode index L */
               const UINT modeM, /*<< Mode index M */
               SEOBPrecdynamics *seobdynamics, /*<< Dynamics vector */
               const REAL8 timeorb, /*<< Time of the peak of the orbital frequency */
               const REAL8 m1, /**Component mass 1 */
               const REAL8 m2, /**Component mass 2 */
               const REAL8 deltaT  /**<< Sampling interval */
)
{
    // PRINT_LOG_INFO(LOG_DEBUG, "timeorb = %g\n", timeorb);
    /* Status of function calls */
    UINT i, j;
    UINT debugRC = 0;
    INT failed = 0;
    UINT SEOBWaveformVersion = 451;
    /** Physical quantities */
    REAL8Vector *timeVec = NULL, *hLMdivrholmVec = NULL;
    REAL8Vector polarDynamics, values;
    REAL8 valuesdata[14] = {0.};
    REAL8 polarDynamicsdata[4] = {0.};
    polarDynamics.length = 4;
                    values.length = 14;
    polarDynamics.data = polarDynamicsdata;
                    values.data = valuesdata;
    REAL8 omegaAttachmentPoint, hLMdivrholmAttachmentPoint, rholmNRAttachmentPoint, rholmBeforeCalibAttachmentPoint;
    REAL8 rAttachmentPoint, hLMdivrholmEAttachmentPoint, rholmENRAttachmentPoint, rholmERealBeforeCalAP, rholmEImagBeforeCalAP;
    REAL8 v;
    COMPLEX16Vector *hLMVec = NULL, *rholmpwrlVec = NULL;
    COMPLEX16 hLME, rholmpwrlE;
    REAL8Vector *rholmpwrlVecReal = NULL;
    REAL8Vector *rholmpwrlEVecReal = NULL;
    REAL8Vector *rholmpwrlEVecImag = NULL;
    REAL8Vector *hLMdivrholmEVec = NULL;
    REAL8Vector orbOmegaVec, HamVec, rVec;
    UINT retLen = seobdynamics->length;

    REAL8 hLMrealAP, hLMimagAP, hLMrealAPDot, hLMrealAPDDot, hLMimagAPDot, hLMimagAPDDot, x0AP, x0DotAP, x0DDotAP;
    REAL8Vector *hLMrealVec = NULL;
    REAL8Vector *hLMimagVec = NULL;
    REAL8Vector *drVec = NULL;
    REAL8Vector *ncrvVec = NULL;
    /** Find the vaulues of the final spins */
    REAL8 mtot = m1+m2;
    REAL8 eta = params->eta;
    REAL8 s1dotZ = 0, s2dotZ = 0, chi1dotZ = 0, chi2dotZ = 0;
    REAL8 chiS = 0;
    REAL8 chiA = 0;
    REAL8 tplspin = 0;
    REAL8 dr, ncrv; // New
    REAL8 spin1z_omegaPeak = params->spin1z_omegaPeak;
    REAL8 spin2z_omegaPeak = params->spin2z_omegaPeak;
    REAL8 chiS_omegaPeak = 0.5*(spin1z_omegaPeak+ spin2z_omegaPeak);
    REAL8 chiA_omegaPeak = 0.5*(spin1z_omegaPeak-spin2z_omegaPeak);
    SEOBPrecWaveformVariables vars;
/* Create dynamical arrays */

    orbOmegaVec.length = HamVec.length = rVec.length = retLen;

    hLMVec = CreateCOMPLEX16Vector (retLen);
    rholmpwrlVec = CreateCOMPLEX16Vector (retLen);
    timeVec = CreateREAL8Vector (retLen);
    hLMdivrholmVec = CreateREAL8Vector (retLen);
    rholmpwrlVecReal = CreateREAL8Vector (retLen);
    hLMdivrholmEVec = CreateREAL8Vector (retLen);
    rholmpwrlEVecReal = CreateREAL8Vector (retLen);
    rholmpwrlEVecImag = CreateREAL8Vector (retLen);
    drVec = CreateREAL8Vector(retLen);
    ncrvVec = CreateREAL8Vector(retLen);
    if (!hLMVec || !rholmpwrlVec|| !timeVec|| !rholmpwrlVec || !rholmpwrlVecReal || !hLMdivrholmEVec || !rholmpwrlEVecReal)
    {failed = 1; goto QUIT;}

    orbOmegaVec.data = seobdynamics->omegaVec;
    HamVec.data = seobdynamics->HamVec;
    rVec.data = seobdynamics->polarrVec;


    /* Stuff for interpolating function */
    GSL_START;
    gsl_spline *spline = NULL;
    gsl_interp_accel *acc = NULL;

    /* The calibration parameter is only used for 21 and 55 modes, if you try to use this function for other modes, you get an error */

    if (!((modeL == 2 && modeM == 1) || (modeL == 5 && modeM == 5)))
    {
        PRINT_LOG_INFO (LOG_CRITICAL, "Mode %d,%d is not supported by this function.", modeL, modeM);
        failed = 1; 
        goto QUIT;
    }
    /* Populate time vector as necessary */
    for (i = 0; i < timeVec->length; i++)
    {
        timeVec->data[i] = i * deltaT;
    }
    /**Initializing stuff for interpolation */
    spline = gsl_spline_alloc (gsl_interp_cspline, orbOmegaVec.length);
    acc = gsl_interp_accel_alloc ();
    /* Calculation of the frequency at the attachment point */
    REAL8 timewavePeak = timeorb-XLALSimIMREOBGetNRSpinPeakDeltaTv4 (modeL, modeM, m1, m2, spin1z_omegaPeak, spin2z_omegaPeak, params->hParams);
    // print_debug("timewavePeak = %g\n", timewavePeak);
    INT status;
    status = gsl_spline_init (spline, timeVec->data, orbOmegaVec.data, orbOmegaVec.length);
    if (status != GSL_SUCCESS)
    {
        PRINT_LOG_INFO (LOG_CRITICAL, "gsl_spline_init failed.");
        failed = 1; 
        goto QUIT;
    }
    gsl_interp_accel_reset (acc);
    omegaAttachmentPoint = gsl_spline_eval (spline, timewavePeak, acc);
    // REAL8 domAP, ddomAP;
    // domAP = gsl_spline_eval_deriv (spline, timewavePeak, acc);
    // ddomAP = gsl_spline_eval_deriv2 (spline, timewavePeak, acc);
    /** Calculate x0 at the attachment point */
    REAL8 drAP, ddrAP;
    status = gsl_spline_init (spline, timeVec->data, rVec.data, rVec.length);
    if (status != GSL_SUCCESS)
    {
        PRINT_LOG_INFO (LOG_CRITICAL, "gsl_spline_init failed.");
        failed = 1; 
        goto QUIT;
    }
    gsl_interp_accel_reset (acc);
    rAttachmentPoint = gsl_spline_eval (spline, timewavePeak, acc);
    // print_debug("rAttachPoint = %g\n", rAttachmentPoint);
    drAP = gsl_spline_eval_deriv (spline, timewavePeak, acc);
    ddrAP = gsl_spline_eval_deriv2 (spline, timewavePeak, acc);
    x0AP = 1./sqrt(rAttachmentPoint);
    x0DotAP = -0.5*drAP / pow(rAttachmentPoint,1.5);
    x0DDotAP = 3*drAP*drAP / (4.*pow(rAttachmentPoint, 2.5)) - 0.5*ddrAP / pow(rAttachmentPoint, 1.5);
    hLMrealVec = CreateREAL8Vector (retLen);
    hLMimagVec = CreateREAL8Vector (retLen);

// print_debug("timewave = %g, Final Time = %g\n", timewavePeak, timeVec->data[orbOmegaVec.length-1]);
    /** Calculation and interpolation at the matching point of rho_lm^l + f_lm */
    for(i=0; i<orbOmegaVec.length; i++)
    {
        for ( j=0; j<14; j++) 
        {
            values.data[j] = seobdynamics->array->data[i + (j+1)*retLen];
        }
        s1dotZ = seobdynamics->s1dotZVec[i];
        s2dotZ = seobdynamics->s2dotZVec[i];
        polarDynamics.data[0] = seobdynamics->polarrVec[i];
        polarDynamics.data[1] = seobdynamics->polarphiVec[i];
        polarDynamics.data[2] = seobdynamics->polarprTVec[i];
        polarDynamics.data[3] = seobdynamics->polarpphiVec[i];
        chi1dotZ = s1dotZ * mtot*mtot / (m1*m1);
        chi2dotZ = s2dotZ * mtot*mtot / (m2*m2);
        chiS = 0.5*(chi1dotZ+chi2dotZ);
        chiA = 0.5*(chi1dotZ-chi2dotZ);
        tplspin = (1.-2.*eta) * chiS + (m1 - m2)/(m1 + m2) * chiA;
        v = cbrt (orbOmegaVec.data[i]);
        if ( XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients( params->hCoeffs, m1, m2, eta, tplspin, chiS, chiA, SEOBWaveformVersion ) == CEV_FAILURE )
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients at step %d of the loop.", i);
            failed = 1; 
            goto QUIT;
        }
        // if (CODE_VERSION == 0)
        // {
        if (PREC_FLAG == 2)
        {
            prec_CalculateSEOBPrecWaveformVariables(&vars, 
                seobdynamics->nchiaVec[i], seobdynamics->nchisVec[i],
                seobdynamics->lchiaVec[i], seobdynamics->lchisVec[i],
                seobdynamics->echiaVec[i], seobdynamics->echisVec[i],
                seobdynamics->chi1chi1[i], seobdynamics->chi1chi2[i], seobdynamics->chi2chi2[i],
                seobdynamics->JnVec[i], seobdynamics->JlVec[i], seobdynamics->JeVec[i],
                seobdynamics->polarrVec[i], seobdynamics->polarprTVec[i], seobdynamics->prTDotVec[i]);
            if (prec_EOBGetPrecEccSpinFactorizedWaveform_v1(&(hLMVec->data[i]), &polarDynamics, &values, v, HamVec.data[i], modeL, modeM, &vars, params) != CEV_SUCCESS)
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in prec_EOBGetPrecEccSpinFactorizedWaveform at step %d of the loop.", i);
                return CEV_FAILURE;
            }
        } else if ( PREC_FLAG > 2) {
            // print_debug("here\n");
            prec_CalculateSEOBPrecWaveformVariables(&vars, 
                seobdynamics->nchiaVec[i], seobdynamics->nchisVec[i],
                seobdynamics->lchiaVec[i], seobdynamics->lchisVec[i],
                seobdynamics->echiaVec[i], seobdynamics->echisVec[i],
                seobdynamics->chi1chi1[i], seobdynamics->chi1chi2[i], seobdynamics->chi2chi2[i],
                seobdynamics->JnVec[i], seobdynamics->JlVec[i], seobdynamics->JeVec[i],
                seobdynamics->polarrVec[i], seobdynamics->polarprTVec[i], seobdynamics->prTDotVec[i]);
            if (prec_EOBGetPrecEccSpinFactorizedWaveform_v2(&(hLMVec->data[i]), &polarDynamics, &values, v, HamVec.data[i], modeL, modeM, &vars, params) != CEV_SUCCESS)
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in prec_EOBGetPrecEccSpinFactorizedWaveform at step %d of the loop.", i);
                return CEV_FAILURE;
            }
        } else if ( XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform( &(hLMVec->data[i]), &polarDynamics, &values, v, HamVec.data[i], modeL, modeM, params ) == CEV_FAILURE )
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform at step %d of the loop.", i);
            failed = 1; 
            goto QUIT;
        }
        if (XLALSimIMRSpinEOBGetAmplitudeResidualPrec (&(rholmpwrlVec->data[i]), v, HamVec.data[i], modeL, modeM, params) == CEV_FAILURE)
        //RC: for the 21 and 55 mode rholmpwrlVec is always real. This is not true for the 33 mode. For this reason, in order to make this function general, we use a complex variable for it.
        {
            /* TODO: Clean-up */
            failed = 1; 
            goto QUIT;
        }

        rholmpwrlVecReal->data[i] = (REAL8)creal(rholmpwrlVec->data[i]);
        hLMdivrholmVec->data[i] = ((REAL8)cabs(hLMVec->data[i]))/fabs(rholmpwrlVecReal->data[i]);

    }
    status = gsl_spline_init (spline, timeVec->data, hLMdivrholmVec->data, hLMdivrholmVec->length);
    if (status != GSL_SUCCESS)
    {
        PRINT_LOG_INFO (LOG_CRITICAL, "gsl_spline_init failed.");
        failed = 1; 
        goto QUIT;
    }
    gsl_interp_accel_reset (acc);
    hLMdivrholmAttachmentPoint = gsl_spline_eval (spline, timewavePeak, acc);


    REAL8 nra = 0, nraDot = 0, nraDDot = 0;
    REAL8 nrOmega, nrOmegaDot;
    nra = XLALSimIMREOBGetNRSpinPeakAmplitudeV4 (modeL, modeM, m1, m2,chiS_omegaPeak, chiA_omegaPeak);
    nraDot = XLALSimIMREOBGetNRSpinPeakADotV4(modeL, modeM, m1, m2, chiS_omegaPeak, chiA_omegaPeak);
    nraDDot = XLALSimIMREOBGetNRSpinPeakADDotV4(modeL, modeM, m1, m2, chiS_omegaPeak, chiA_omegaPeak);
    nrOmega = XLALSimIMREOBGetNRSpinPeakOmegaV4 (modeL, modeM, eta, chiS_omegaPeak + chiA_omegaPeak * (m1 - m2) / (m1 + m2) / (1. - 2. * eta));
    nrOmegaDot = XLALSimIMREOBGetNRSpinPeakOmegaDotV4(modeL, modeM, eta, chiS_omegaPeak + chiA_omegaPeak * (m1 - m2) / (m1 + m2) / (1. - 2. * eta));
    PRINT_LOG_INFO(LOG_DEBUG, "hNR_%d%d = %g", modeL, modeM, nra);
    if((fabs(nra/eta)< 3e-2) && ((modeL == 2) && (modeM == 1)))
    {
        //R.C.: safeguard to avoid the 21 mode to go to 0
        nra = GSL_SIGN(nra)*eta*3e-2;
    }
    if((fabs(nra/eta)< 1e-4) && ((modeL == 5) && (modeM == 5)))
    {
        //R.C.: safeguard to avoid the 55 mode to go to 0
        nra = GSL_SIGN(nra)*eta*1e-4;
    }
    rholmNRAttachmentPoint = nra/hLMdivrholmAttachmentPoint;
    status = gsl_spline_init (spline, timeVec->data, rholmpwrlVecReal->data, rholmpwrlVecReal->length);
    if (status != GSL_SUCCESS)
    {
        PRINT_LOG_INFO (LOG_CRITICAL, "gsl_spline_init failed.");
        failed = 1; 
        goto QUIT;
    }
    gsl_interp_accel_reset (acc);
    /* rho_lm^l + f_lm before the calibration parameter is set */
    rholmBeforeCalibAttachmentPoint = gsl_spline_eval (spline, timewavePeak, acc);

    if ((modeL == 2) && (modeM == 1))
    {
        /* Here we compute ((rho_lm^l + f_lm + CalPar*omega^7/3)_NR - (rho_lm^l + f_lm)_EOB)/omega^7/3 to get CalPar.
        The factor rholmpwrlVecReal->data[0])/cabs(rholmpwrlVecReal->data[0]) is used to know the sign of the function (rho_lm^l + f_lm + CalPar*omega^7/3)_NR which is computed as absolute value */
        params->cal21 = (rholmNRAttachmentPoint - rholmBeforeCalibAttachmentPoint)/(pow(omegaAttachmentPoint,7./3.));
    }
    if ((modeL == 5) && (modeM == 5))
    {
        /* Here we compute ((rho_lm^l + f_lm + CalPar*omega^7/3)_NR - (rho_lm^l + f_lm)_EOB)/omega^7/3 to get CalPar.
        The factor rholmpwrlVecReal->data[0])/cabs(rholmpwrlVecReal->data[0]) is used to know the sign of the function (rho_lm^l + f_lm + CalPar*omega^7/3)_NR which is computed as absolute value */
        params->cal55 = (rholmNRAttachmentPoint - rholmBeforeCalibAttachmentPoint)/(pow(omegaAttachmentPoint,5./3.));
                            //printf("params->cal55 = %.16f\n",params->cal55);
    }
    if (isnan(creal(params->cal21)) || isnan(cimag(params->cal21)) ||
        isnan(creal(params->cal55)) || isnan(cimag(params->cal55)))
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "(l,m) = (%d, %d) calibration get nan", modeL, modeM);
        failed = 1;
    }
QUIT:
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
    GSL_END;
    STRUCTFREE(hLMdivrholmVec, REAL8Vector);
    STRUCTFREE(hLMdivrholmEVec, REAL8Vector);
    STRUCTFREE(hLMVec, COMPLEX16Vector);
    STRUCTFREE(timeVec, REAL8Vector);
    STRUCTFREE(rholmpwrlVec, COMPLEX16Vector);
    STRUCTFREE(rholmpwrlVecReal, REAL8Vector);
    STRUCTFREE(rholmpwrlEVecReal, REAL8Vector);
    STRUCTFREE(rholmpwrlEVecImag, REAL8Vector);
    STRUCTFREE(hLMrealVec, REAL8Vector);
    STRUCTFREE(hLMimagVec, REAL8Vector);
    STRUCTFREE(drVec, REAL8Vector);
    STRUCTFREE(ncrvVec, REAL8Vector);
    if (failed)
        return CEV_FAILURE;

    return CEV_SUCCESS;
}


INT prec_calculateSEOBFactorizedWaveformCorrectionFromDynVectors(
    SpinEOBParams * params,
    const INT l,
    const INT m,
    REAL8Vector *xVec,
    REAL8Vector *yVec,
    REAL8Vector *zVec,
    REAL8Vector *vxVec,
    REAL8Vector *vyVec,
    REAL8Vector *vzVec,
    REAL8Vector *pTxVec,
    REAL8Vector *pTyVec,
    REAL8Vector *pTzVec,
    REAL8Vector *s1xVec,
    REAL8Vector *s1yVec,
    REAL8Vector *s1zVec,
    REAL8Vector *s2xVec,
    REAL8Vector *s2yVec,
    REAL8Vector *s2zVec,
    REAL8Vector *hamVec,
    REAL8Vector *prTDotVec,
    REAL8Vector **rholm_real,
    REAL8Vector **rholm_imag,
    REAL8Vector **flm_real,
    REAL8Vector **flm_imag
)
{
    INT i, j, length;
    length = xVec->length;
    REAL8Vector *rholm_val_re = CreateREAL8Vector(length);
    REAL8Vector *rholm_val_im = CreateREAL8Vector(length);
    REAL8Vector *flm_val_re = CreateREAL8Vector(length);
    REAL8Vector *flm_val_im = CreateREAL8Vector(length);
    REAL8 rVec[3] = {0, 0, 0}, vVec[3] = {0, 0, 0};
    REAL8 ham, pTVec[3] = {0, 0, 0};
    REAL8 s1Vec[3] = {0, 0, 0}, s2Vec[3] = {0, 0, 0};
    REAL8 omega, s1dotZ, s2dotZ;
    REAL8 eta, m1, m2, mtot;
    m1 = params->m1;
    m2 = params->m2;
    mtot = m1+m2;
    eta = params->eta;
    REAL8 m1sq = m1*m1/(m1+m2)/(m1+m2);
    REAL8 m2sq = m2*m2/(m1+m2)/(m1+m2);
    REAL8 rcrossrdot[3] = {0, 0, 0};
    REAL8 rcrossp[3] = {0, 0, 0};
    REAL8 LNhat[3] = {0, 0, 0};
    REAL8 Lhat[3] = {0, 0, 0};

    REAL8 nhat[3] = {0,0,0}, chi1Vec[3] = {0,0,0}, chiSVec[3] = {0,0,0};
    REAL8 lhat[3] = {0,0,0}, chi2Vec[3] = {0,0,0}, chiAVec[3] = {0,0,0};
    REAL8Vector *values = NULL;
    REAL8Vector *polarDynamics = NULL;
    // REAL8Vector *dvalues = NULL;
    values = CreateREAL8Vector(14);
    polarDynamics = CreateREAL8Vector(4);
    // dvalues = CreateREAL8Vector(14);
    memset(values->data, 0, (values->length) * sizeof(REAL8));
    // memset(dvalues->data, 0, (dvalues->length) * sizeof(REAL8));
    SEOBPrecWaveformVariables vars;
    for (i = 0; i < length; i++) 
    {
        rVec[0] = xVec->data[i];
        rVec[1] = yVec->data[i];
        rVec[2] = zVec->data[i];
        vVec[0] = vxVec->data[i];
        vVec[1] = vyVec->data[i];
        vVec[2] = vzVec->data[i];
        pTVec[0] = pTxVec->data[i];
        pTVec[1] = pTyVec->data[i];
        pTVec[2] = pTzVec->data[i];
        s1Vec[0] = s1xVec->data[i];
        s1Vec[1] = s1yVec->data[i];
        s1Vec[2] = s1zVec->data[i];
        s2Vec[0] = s2xVec->data[i];
        s2Vec[1] = s2yVec->data[i];
        s2Vec[2] = s2zVec->data[i];
        REAL8 prTDot = prTDotVec->data[i];
        ham = hamVec->data[i];
        // print_debug("s1Vec = (%.5e, %.5e, %.5e), s2Vec = (%.5e, %.5e, %.5e)\n", 
        //     s1Vec[0], s1Vec[1], s1Vec[2],
        //     s2Vec[0], s2Vec[1], s2Vec[2]);
        // print_debug("prTDot = %.5e, ham = %.5e\n", prTDot, ham);
#if 1
        cross_product3d(rVec, pTVec, rcrossp);
        cross_product3d(rVec, vVec, rcrossrdot);
        REAL8 rcrossrdotNorm = sqrt(inner_product3d(rcrossrdot, rcrossrdot));
        for (j = 0; j < 3; j++) 
            LNhat[j] = rcrossrdot[j] / rcrossrdotNorm;
        REAL8 magL = sqrt(inner_product3d(rcrossp, rcrossp));
        for (j = 0; j < 3; j++)
            Lhat[j] = rcrossp[j] / magL;
        REAL8 polarr = sqrt(inner_product3d(rVec, rVec));
        REAL8 polarpr = inner_product3d(rVec, pTVec) / polarr;
        /* Compute waveform coefficients */
        REAL8 omega = rcrossrdotNorm / (polarr * polarr);
        s1dotZ = inner_product3d(s1Vec, Lhat);
        s2dotZ = inner_product3d(s2Vec, Lhat);
        REAL8 chi1dotZ = s1dotZ * mtot * mtot / (m1 * m1);
        REAL8 chi2dotZ = s2dotZ * mtot * mtot / (m2 * m2);
        REAL8 chiS = SEOBCalculateChiS(chi1dotZ, chi2dotZ);
        REAL8 chiA = SEOBCalculateChiA(chi1dotZ, chi2dotZ);
        for (j=0; j<3; j++)
        {
            chi1Vec[j] = s1Vec[j] * mtot * mtot / (m1 * m1);
            chi2Vec[j] = s2Vec[j] * mtot * mtot / (m2 * m2);
            chiSVec[j] = 0.5*(chi1Vec[j] + chi2Vec[j]);
            chiAVec[j] = 0.5*(chi1Vec[j] - chi2Vec[j]);
            nhat[j] = rVec[j] / polarr;
        }
        cross_product3d(LNhat, nhat, lhat);
        REAL8 tplspin = SEOBCalculatetplspin(m1, m2, eta, s1dotZ, s2dotZ);

        if (XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(
                params->hCoeffs, m1, m2, eta, tplspin, chiS, chiA,
                451) == CEV_FAILURE) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients at step %d of the loop.", i);
            return CEV_FAILURE;
        }
        params->hCoeffs->f21v7c = 0.0;
        params->hCoeffs->f21v7cEff00 = 0.0;
        params->hCoeffs->f21v7cEff10 = 0.0;
        params->hCoeffs->f21v7cEff11 = 0.0;
        params->hCoeffs->f21v7cEff01 = 0.0;
        params->hCoeffs->f21v7cEff02 = 0.0;
        params->hCoeffs->f55v5c = 0.0;
        // print_debug("f21v7c = %.16f\n", seobParams->hCoeffs->f21v7c);
        /* Dynamics, polar dynamics, omega */
        memcpy( values->data,   rVec,  3*sizeof(REAL8) );
        memcpy( values->data+3, pTVec, 3*sizeof(REAL8) );
        memcpy( values->data+6, s1Vec, 3*sizeof(REAL8) );
        memcpy( values->data+9, s2Vec, 3*sizeof(REAL8) );
        // for (j = 0; j < 14; j++)
        //     values.data[j] = seobdynamics->array->data[i + (j + 1) * retLen];
        polarDynamics->data[0] = polarr;
        polarDynamics->data[1] = 0.0;
        polarDynamics->data[2] = polarpr;
        polarDynamics->data[3] = magL;
        REAL8 nchia, nchis;
        REAL8 lchia, lchis;
        REAL8 echia, echis;
        REAL8 chi1chi1, chi1chi2, chi2chi2;
        nchia = inner_product3d(chiAVec, nhat);
        nchis = inner_product3d(chiSVec, nhat);
        lchia = inner_product3d(chiAVec, lhat);
        lchis = inner_product3d(chiSVec, lhat);
        echia = inner_product3d(chiAVec, LNhat);
        echis = inner_product3d(chiSVec, LNhat);
        chi1chi1 = inner_product3d(chi1Vec, chi1Vec);
        chi1chi2 = inner_product3d(chi1Vec, chi2Vec);
        chi2chi2 = inner_product3d(chi2Vec, chi2Vec);
        REAL8 v = cbrt(omega);
        COMPLEX16 rholm = 0.;
        COMPLEX16 flm = 0.;
        prec_CalculateSEOBPrecWaveformVariables(&vars, 
            nchia, nchis,
            lchia, lchis,
            echia, echis,
            chi1chi1, chi1chi2, chi2chi2,
            0, 0, 0,
            polarr, polarpr, prTDot);
        if (prec_CalculateFactorizedWaveformCorrection(&rholm, &flm, polarDynamics, values, v, ham, l, m, &vars, params) != CEV_SUCCESS)
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "failure in prec_EOBGetPrecEccSpinFactorizedWaveform at step %d of the loop.", i);
            return CEV_FAILURE;
        }
        rholm_val_re->data[i] = creal(rholm);
        rholm_val_im->data[i] = cimag(rholm);
        flm_val_re->data[i] = creal(flm);
        flm_val_im->data[i] = cimag(flm);
#endif
    }
    STRUCTFREE(values, REAL8Vector);
    STRUCTFREE(polarDynamics, REAL8Vector);
    *rholm_real = rholm_val_re;
    *rholm_imag = rholm_val_im;
    *flm_real = flm_val_re;
    *flm_imag = flm_val_im;
    return CEV_SUCCESS;
}

/**
 * This function generates a waveform mode for a given SEOB dynamics.
 */
// NOTE: as is written here, the step
// XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients in the loop will be repeated
// across modes -- would be more efficient to loop on modes inside the loop on
// times
static int SEOBPrecCalculatehlmAmpPhase(
    CAmpPhaseSequence *
        *hlm, /**<< Output: hlm in complex amplitude / phase form */
    INT4 l,   /**<< Input: mode index l */
    INT4 m,   /**<< Input: mode index m */
    SEOBPrecdynamics *seobdynamics, /**<< Input: SEOB dynamics */
    EOBNonQCCoeffs *nqcCoeffs,  /**<< Input: NQC coeffs */
    SpinEOBParams *seobParams,  /**<< SEOB params */
    // UINT4 SpinsAlmostAligned, /**<< flag to decide wether to fall back to
    // aligned spins  */
    UINT includeNQC /**<< flag to choose wether or not to include NQC */
) {
    /* Check that the input double pointer are not NULL */
    if (!hlm) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "pointer to CAmpPhaseSequence hlm is NULL.");
        return CEV_FAILURE;
    }

    /* Masses */
    REAL8 m1 = seobParams->m1;
    REAL8 m2 = seobParams->m2;
    REAL8 deltaT = seobdynamics->tVec[1] - seobdynamics->tVec[0];
    REAL8 eta = seobParams->eta;
    REAL8 mtot = m1 + m2;
    REAL8 dr, ncrv;
    SEOBPrecWaveformVariables vars;
    // UINT SpinAlignedEOBversion = seobParams->SpinAlignedEOBversion;
    UINT SpinAlignedEOBversion = 4;
    UINT SpinAlignedEOBversionWaveform; // RC: I use this different variable
                                        // because the PN terms in the waveform
                                        // are different from those in the flux

    /* Length of dynamics data and sampling step */
    UINT retLen = seobdynamics->length;

    /* Allocate structure for complex amplitude and phase */
    if (CAmpPhaseSequence_Init(hlm, retLen) == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in CAmpPhaseSequence_Init.");
        return CEV_FAILURE;
    }

    /* Workspace vectors */
    REAL8Vector values, polarDynamics;
    REAL8 valuesdata[14] = {0.};
    REAL8 polarDynamicsdata[4] = {0.};
    values.length = 14;
    polarDynamics.length = 4;
    values.data = valuesdata;
    polarDynamics.data = polarDynamicsdata;
    REAL8 tPeakOmega = seobParams->tPeakOmega;

    /* Calibration parameter */
    if (includeNQC == 0 && seobParams->hParams->Mf_max < seobParams->hParams->Mf_min)
    {
        if (((l == 2) && (m == 1)) || ((l == 5) && (m == 5))) 
        {
            if (XLALSimIMREOBPrecCalcCalibCoefficientHigherModesPrec(
                    seobParams, l, m, seobdynamics,
                    tPeakOmega - seobdynamics->tVec[0], m1, m2,
                    deltaT) == CEV_FAILURE) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBPrecCalcCalibCoefficientHigherModesPrec.");
                return CEV_FAILURE;
            }
        }
    }

    /* Loop to compute compute amplitude and phase of the hlm mode */
    REAL8 s1dotZ, s2dotZ, chiS, chiA, tplspin, t, omega, ham, v;
    UINT i, j;
    for (i = 0; i < retLen; i++) 
    {
        /* Compute waveform coefficients */
        t = seobdynamics->tVec[i];
        omega = seobdynamics->omegaVec[i];
        ham = seobdynamics->HamVec[i];
        s1dotZ = seobdynamics->s1dotZVec[i];
        s2dotZ = seobdynamics->s2dotZVec[i];
        REAL8 chi1dotZ = s1dotZ * mtot * mtot / (m1 * m1);
        REAL8 chi2dotZ = s2dotZ * mtot * mtot / (m2 * m2);
        chiS = SEOBCalculateChiS(chi1dotZ, chi2dotZ);
        chiA = SEOBCalculateChiA(chi1dotZ, chi2dotZ);

        tplspin = SEOBCalculatetplspin(m1, m2, eta, s1dotZ, s2dotZ);
        if (SpinAlignedEOBversion == 4) 
        {
            SpinAlignedEOBversionWaveform = 451;
        } else 
        {
            SpinAlignedEOBversionWaveform = SpinAlignedEOBversion;
        }

        if (XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(
                seobParams->hCoeffs, m1, m2, eta, tplspin, chiS, chiA,
                SpinAlignedEOBversionWaveform) == CEV_FAILURE) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients at step %d of the loop.", i);
            return CEV_FAILURE;
        }
        seobParams->hCoeffs->f21v7c = seobParams->cal21;
        seobParams->hCoeffs->f21v7cEff00 = seobParams->cal21E;
        seobParams->hCoeffs->f21v7cEff10 = seobParams->cal21E1;
        seobParams->hCoeffs->f21v7cEff11 = seobParams->cal21E2;
        seobParams->hCoeffs->f21v7cEff01 = seobParams->cal21E3;
        seobParams->hCoeffs->f21v7cEff02 = seobParams->cal21E4;
        seobParams->hCoeffs->f55v5c = seobParams->cal55;
        // print_debug("f21v7c = %.16f\n", seobParams->hCoeffs->f21v7c);
        /* Dynamics, polar dynamics, omega */
        for (j = 0; j < 14; j++) 
            values.data[j] = seobdynamics->array->data[i + (j + 1) * retLen];
        polarDynamics.data[0] = seobdynamics->polarrVec[i];
        polarDynamics.data[1] = seobdynamics->polarphiVec[i];
        polarDynamics.data[2] = seobdynamics->polarprTVec[i];
        polarDynamics.data[3] = seobdynamics->polarpphiVec[i];
        v = cbrt(omega);
        COMPLEX16 hlm_val = 0.;
        if (PREC_FLAG == 2)
        {
            prec_CalculateSEOBPrecWaveformVariables(&vars, 
                seobdynamics->nchiaVec[i], seobdynamics->nchisVec[i],
                seobdynamics->lchiaVec[i], seobdynamics->lchisVec[i],
                seobdynamics->echiaVec[i], seobdynamics->echisVec[i],
                seobdynamics->chi1chi1[i], seobdynamics->chi1chi2[i], seobdynamics->chi2chi2[i],
                seobdynamics->JnVec[i], seobdynamics->JlVec[i], seobdynamics->JeVec[i],
                seobdynamics->polarrVec[i], seobdynamics->polarprTVec[i], seobdynamics->prTDotVec[i]);
            if (prec_EOBGetPrecEccSpinFactorizedWaveform_v1(&hlm_val, &polarDynamics, &values, v, ham, l, m, &vars, seobParams) != CEV_SUCCESS)
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in prec_EOBGetPrecEccSpinFactorizedWaveform at step %d of the loop.", i);
                return CEV_FAILURE;
            }
        } else if (PREC_FLAG > 2) {
            prec_CalculateSEOBPrecWaveformVariables(&vars, 
                seobdynamics->nchiaVec[i], seobdynamics->nchisVec[i],
                seobdynamics->lchiaVec[i], seobdynamics->lchisVec[i],
                seobdynamics->echiaVec[i], seobdynamics->echisVec[i],
                seobdynamics->chi1chi1[i], seobdynamics->chi1chi2[i], seobdynamics->chi2chi2[i],
                seobdynamics->JnVec[i], seobdynamics->JlVec[i], seobdynamics->JeVec[i],
                seobdynamics->polarrVec[i], seobdynamics->polarprTVec[i], seobdynamics->prTDotVec[i]);
            if (prec_EOBGetPrecEccSpinFactorizedWaveform_v2(&hlm_val, &polarDynamics, &values, v, ham, l, m, &vars, seobParams) != CEV_SUCCESS)
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in prec_EOBGetPrecEccSpinFactorizedWaveform at step %d of the loop.", i);
                return CEV_FAILURE;
            }
        } else if (XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform(
                &hlm_val, &polarDynamics, &values, v, ham, l, m, seobParams) != CEV_SUCCESS) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform at step %d of the loop.", i);
            return CEV_FAILURE;
        }
        /* NQC correction */
        COMPLEX16 factor_nqc = 1.;
        if (includeNQC) 
        {
            if (XLALSimIMRSpinEOBNonQCCorrection(&factor_nqc, &values, omega, t, 
                                    seobParams->tWind, seobParams->wWind,
                                    nqcCoeffs) == CEV_FAILURE) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL,
                    "failure in XLALSimIMRSpinEOBNonQCCorrection at step %d of the loop.", i);
                return CEV_FAILURE;
            }
        } else
        factor_nqc = 1.;
        /* Result and output */
        COMPLEX16 hlmNQC = hlm_val * factor_nqc;
        (*hlm)->xdata->data[i] = t; /* Copy times */
        (*hlm)->camp_real->data[i] = cabs(hlmNQC);
        (*hlm)->camp_imag->data[i] = 0.; /* We use only real amplitudes */
        (*hlm)->phase->data[i] = carg(hlmNQC);
        // if (IS_DEBUG)
        //     print_out("%.16e\t%.16e\t%.16e\n", t, creal(hlmNQC), cimag(hlmNQC));
    }

    /* Unwrap the phase vector, in place */
    XLALREAL8VectorUnwrapAngle((*hlm)->phase, (*hlm)->phase);

    return CEV_SUCCESS;
}


int SEOBPrecCalculateSphHarmListNQCCoefficientsV4(
    SphHarmListEOBNonQCCoeffs *
        *nqcCoeffsList, /**<< Output: non-quasi-circular coefficients as a list
                           for each mode */
    INT modes[][2],    /**<< Input: array of modes (l,m) */
    UINT nmodes,       /**<< Input: number of modes (l,m) */
    REAL8 tPeakOmega,   /**<< Input: time of peak of Omega */
    SEOBPrecdynamics *seobdynamics,  /**<< Input: SEOB dynamics */
    SpinEOBParams *seobParams,   /**<< Input: SEOB params */
    REAL8Vector *chi1_omegaPeak, /**<< Input: dimensionless spin 1 at peak of
                                    omega in L_N frame */
    REAL8Vector *chi2_omegaPeak /**<< Input: dimensionless spin 2 at peak of
                                    omega in L_N frame */
)
{
    INT failed = 0;
    /* Masses */
    REAL8 m1 = seobParams->m1;
    REAL8 m2 = seobParams->m2;
    seobParams->tPeakOmega = tPeakOmega;
    seobParams->spin1z_omegaPeak = chi1_omegaPeak->data[2];
    seobParams->spin2z_omegaPeak = chi2_omegaPeak->data[2];
    // REAL8 mtot = m1 + m2;
    // UINT4 SpinAlignedEOBversion =
    // seobParams->seobCoeffs->SpinAlignedEOBversion;

    /* Length of dynamics data and sampling step */
    UINT retLen = seobdynamics->length;
    REAL8 deltaT = seobdynamics->tVec[1] - seobdynamics->tVec[0];

    /* Create vectors from dynamics */
    REAL8Vector r, pr, omega;
    r.length = pr.length = omega.length = retLen;
    r.data = seobdynamics->polarrVec;
    pr.data = seobdynamics->polarprTVec;
    omega.data = seobdynamics->omegaVec;

    /* Workspace vectors */
    /*
    REAL8Vector s1Vec, s2Vec, sigmaKerr;
    REAL8 s1Vecdata[3] = {0.};
    REAL8 s2Vecdata[3] = {0.};
    REAL8 sigmaKerrdata[3] = {0.};
    s1Vec.length = s2Vec.length = sigmaKerr.length = 3;
    s1Vec.data = s1Vecdata;
    s2Vec.data = s2Vecdata;
    sigmaKerr.data = sigmaKerrdata;
    */
    /* Final values for a and for the spins projected onto Z */
    /*
    REAL8 s1dotZfinal = seobdynamics->s1dotZVec[retLen-1];
    REAL8 s2dotZfinal = seobdynamics->s2dotZVec[retLen-1];
    REAL8 chi1dotZfinal = s1dotZfinal * mtot*mtot / (m1*m1);
    REAL8 chi2dotZfinal = s2dotZfinal * mtot*mtot / (m2*m2);
    REAL8 chiSfinal = SEOBCalculateChiS( chi1dotZfinal, chi2dotZfinal );
    REAL8 chiAfinal = SEOBCalculateChiA( chi1dotZfinal, chi2dotZfinal );
    s1Vec.data[0] = seobdynamics->s1Vecx[retLen-1];
    s1Vec.data[1] = seobdynamics->s1Vecy[retLen-1];
    s1Vec.data[2] = seobdynamics->s1Vecz[retLen-1];
    s2Vec.data[0] = seobdynamics->s2Vecx[retLen-1];
    s2Vec.data[1] = seobdynamics->s2Vecy[retLen-1];
    s2Vec.data[2] = seobdynamics->s2Vecz[retLen-1];
    SEOBCalculateSigmaKerr( &sigmaKerr, &s1Vec, &s2Vec );
    */
    // REAL8 afinal = sqrt( inner_product( sigmaKerr.data, sigmaKerr.data ) );

    REAL8 chi1dotZfinal = chi1_omegaPeak->data[2];
    REAL8 chi2dotZfinal = chi2_omegaPeak->data[2];
    REAL8 chiSfinal = SEOBCalculateChiS(chi1dotZfinal, chi2dotZfinal);
    REAL8 chiAfinal = SEOBCalculateChiA(chi1dotZfinal, chi2dotZfinal);
    REAL8 q = m1/m2;
    //printf("chiA = %.16f\n",chiAfinal);
    /* Time elapsed from the start of the dynamics to tPeakOmega */
    REAL8 tPeakOmegaFromStartDyn = tPeakOmega - seobdynamics->tVec[0];
// print_debug("tPeak = %f, t0 = %f\n", tPeakOmega, seobdynamics->tVec[0]);
    /* Compute NQC coefficients - output is nqcCoeffs */
    /* NOTE: internally, XLALSimIMRSpinEOBCalculateNQCCoefficientsV4 builds a time
    * vector as i*deltaT, starting at 0 - thus tPeakOmega has to be measured from
    * the start of the dynamics */
    CAmpPhaseSequence *hlm = NULL;
    REAL8Vector *hlm_amp = NULL;
    REAL8Vector *hlm_phase = NULL;

    /* Modes are to be generated without NQC */
    UINT includeNQC = 0, nmode;
    /* Loop over modes */
    for ( nmode = 0; nmode < nmodes; nmode++) 
    {

        INT l = modes[nmode][0];
        INT m = modes[nmode][1];

        EOBNonQCCoeffs *nqcCoeffs = MYMalloc(sizeof(EOBNonQCCoeffs));
        memset(nqcCoeffs, 0, sizeof(EOBNonQCCoeffs));
        /* In the equal mass equal spins case the odd-m modes are 0, so we set the NQCs to 0 */
        if (q<1.005 && (m % 2 != 0) && (fabs(chiAfinal) < 0.15)) 
        { /* In this case, set NQC coeffs to 0 for odd m */
            nqcCoeffs->a1 = 0.;
            nqcCoeffs->a2 = 0.;
            nqcCoeffs->a3 = 0.;
            nqcCoeffs->a3S = 0.;
            nqcCoeffs->a4 = 0.;
            nqcCoeffs->a5 = 0.;
            nqcCoeffs->b1 = 0.;
            nqcCoeffs->b2 = 0.;
            nqcCoeffs->b3 = 0.;
            nqcCoeffs->b4 = 0.;
        } 
        else 
        {

        /* Compute amplitude and phase of the mode hlm pre-NQC */
        /* Mode hlm is to be generated without NQC */
        includeNQC = 0;
        hlm = NULL;
        // PRINT_LOG_INFO(LOG_DEBUG, "here, calc Amp Phase: (l,m) = (%d,%d)", l, m);
        // if (l==2 && m==1)
        //     DEBUG_START;
        SEOBPrecCalculatehlmAmpPhase(&hlm, l, m, seobdynamics, nqcCoeffs, seobParams,
                                includeNQC);
        // if (l==2 && m==1)
        //     DEBUG_END;
        /* Cast to amp/phase vectors */
        hlm_amp = NULL;
        if (!(hlm_amp = CreateREAL8Vector(hlm->xdata->length))) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "failed to allocate REAL8Vector hlm_amp.");
            MYFree(nqcCoeffs);
            failed = 1;
            goto QUIT;
        }
        hlm_phase = NULL;
        if (!(hlm_phase = CreateREAL8Vector(hlm->xdata->length))) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "failed to allocate REAL8Vector hlm_phase.");
            MYFree(nqcCoeffs);
            failed = 1;
            goto QUIT;
        }
        memcpy(hlm_amp->data, hlm->camp_real->data,
                hlm->xdata->length * sizeof(REAL8));
        memcpy(hlm_phase->data, hlm->phase->data,
                hlm->xdata->length * sizeof(REAL8));

        /* Compute NQC */
        // PRINT_LOG_INFO(LOG_DEBUG, "here, calc NQC start: (l,m) = (%d,%d)", l, m);
        if (XLALSimIMRSpinEOBCalculateNQCCoefficientsV4(
                hlm_amp, hlm_phase, &r, &pr, &omega, l, m, tPeakOmegaFromStartDyn, seobdynamics->tVec[0],
                deltaT, m1, m2, chiAfinal, chiSfinal, seobParams->tWind, seobParams->wWind,
                nqcCoeffs, seobParams) == CEV_FAILURE) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "failure for the mode (l,m) = (%d, %d)", l, m);
            MYFree(nqcCoeffs);
            failed = 1;
            goto QUIT;
        }
        // PRINT_LOG_INFO(LOG_DEBUG, "here, calc NQC end: (l,m) = (%d,%d)", l, m);
        /* Cleanup */
        STRUCTFREE(hlm, CAmpPhaseSequence);
        STRUCTFREE(hlm_amp, REAL8Vector);
        STRUCTFREE(hlm_phase, REAL8Vector);
        }

        /* Add computed NQCs to the list */
        // PRINT_LOG_INFO(LOG_DEBUG, "here, Add computed NQCs to the list: (l,m) = (%d,%d)", l, m);
        SphHarmListEOBNonQCCoeffs_AddMode(nqcCoeffsList, nqcCoeffs, l, m);
    }
QUIT:
    STRUCTFREE(hlm, CAmpPhaseSequence);
    STRUCTFREE(hlm_amp, REAL8Vector);
    STRUCTFREE(hlm_phase, REAL8Vector);

    if (failed)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}




int SEOBPrecCalculateSphHarmListhlmAmpPhase(
    SphHarmListCAmpPhaseSequence **listhlm,               /**<< Output: list of modes for hlm */
    INT modes[][2],            /**<< Input: array of modes (l,m) */
    UINT nmodes,               /**<< Input: number of modes (l,m) */
    SEOBPrecdynamics *seobdynamics, /**<< Input: SEOB dynamics */
    SphHarmListEOBNonQCCoeffs *listnqcCoeffs, /**<< Input: list of NQCs */
    SpinEOBParams *seobParams,                /**<< SEOB params */
    UINT flagNQC /**<< flag to choose wether or not to include NQC */
) 
{
    /* Read version of SEOB to be used */
    // UINT SpinAlignedEOBversion = seobParams->seobCoeffs->SpinAlignedEOBversion;
    UINT SpinAlignedEOBversion = 4;

    /* Flag for inclusion of NQC, useful for higher modes that have no NQC
    * implemented for some SpinAlignedEOBversion */
    UINT includeNQC = 0, nmode;

    /* Loop over modes */
    for ( nmode = 0; nmode < nmodes; nmode++) 
    {

        INT l = modes[nmode][0];
        INT m = modes[nmode][1];

        if ((!(l == 2 && m == 2)) && (SpinAlignedEOBversion == 3)) 
        {
            includeNQC = 0; /* For HM beyond 22, no NQC available for
                                SpinAlignedEOBversion==3 */
        } else
            includeNQC = flagNQC;

        EOBNonQCCoeffs *nqcCoeffs = NULL;
        if (includeNQC)
            nqcCoeffs = SphHarmListEOBNonQCCoeffs_GetMode(listnqcCoeffs, l, m)->nqcCoeffs;

        CAmpPhaseSequence *hlm = NULL;
        SEOBPrecCalculatehlmAmpPhase(&hlm, l, m, seobdynamics, nqcCoeffs, seobParams,
                                    includeNQC);

        SphHarmListCAmpPhaseSequence_AddMode(listhlm, hlm, l, m);
    }

    return CEV_SUCCESS;
}


INT SEOBPrecAttachRDToSphHarmListhPlm(
    SphHarmListCAmpPhaseSequence **listhPlm_RDattached, /**<< Output: list of extended modes hlm with RD
                                 attached */
    COMPLEX16Vector **sigmaQNMlm0, /**<< Output: list of QNM complex frequency for modes lm,
                         0th overtone (dimensionless) */
    INT4 modes[][2],  /**<< Input: array of modes (l,m) */
    UINT nmodes,     /**<< Input: number of modes (l,m) */
    REAL8 finalMass,  /**<< Input: final mass computed from fit (scaled by M) */
    REAL8
        finalSpin, /**<< Input: final spin computed from fit (dimensionless) */
    SphHarmListCAmpPhaseSequence *listhPlm, /**<< Input: list of modes hlm */
    REAL8 deltaT,                           /**<< Input: time step */
    UINT retLen,        /**<< Input: length of the input modes and dynamics */
    UINT retLenRDPatch, /**<< Input: length of the ringdown patch */
    REAL8 tAttach,       /**<< Input: time of attachment */
    // REAL8 tStart, /**<< Input: starting time (of the HiS part) */
    REAL8Vector *seobvalues, /**<< Input: vector for dynamics values at time of
                                peak of omega */
    SEOBPrecdynamics *seobdynamics,            /**<< Input: SEOB dynamics */
    SpinEOBParams *seobParams             /**<< SEOB params */
) 
{
    /* Check that the input list pointer and vector pointer are NULL */
    if (!(*listhPlm_RDattached == NULL)) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "output pointer for the list hPlm_RDattached is not NULL.");
        return CEV_FAILURE;
    }
    if (!(*sigmaQNMlm0 == NULL)) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "output pointer for the vector sigmaQNMlm0 is not NULL.");
        return CEV_FAILURE;
    }

    /* Create list of output modes */
    UINT retLen_RDattached = retLen + retLenRDPatch, nmode, i;
    for ( nmode = 0; nmode < nmodes; nmode++) 
    {
        INT4 l = modes[nmode][0];
        INT4 m = modes[nmode][1];
        CAmpPhaseSequence *hPlm_RDattached = NULL;
        if (CAmpPhaseSequence_Init(&hPlm_RDattached, retLen_RDattached) == CEV_FAILURE) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "failed to allocate CAmpPhaseSequence hlm for mode (l,m) = (%d,%d).", l, m);
            return CEV_FAILURE;
        }

        SphHarmListCAmpPhaseSequence_AddMode(listhPlm_RDattached, hPlm_RDattached, l, m);
    }

    /* Masses */
    INT failed = 0;
    REAL8 m1 = seobParams->m1;
    REAL8 m2 = seobParams->m2;
    // REAL8 eta = seobParams->eobParams->eta;
    REAL8 mTotal = m1 + m2;
    REAL8 mTScaled = mTotal * CST_MTSUN_SI;
    // Approximant seobApproximant = seobParams->seobApproximant;

    /* Vector of times without attachment */
    REAL8Vector timeVec;
    timeVec.length = retLen;
    timeVec.data = seobdynamics->tVec;
    REAL8 finalM, finalS = 0;

    REAL8Vector *chi1temp = NULL;
    REAL8Vector *chi2temp = NULL;
    REAL8Vector *hPlmRe = NULL;
    REAL8Vector *hPlmIm = NULL;
    SEOBLFrameVectors(&chi1temp, &chi2temp, seobvalues, m1, m2, seobParams->hParams->flagZframe);

    REAL8 chi1Lx = chi1temp->data[0];
    REAL8 chi1Ly = chi1temp->data[1];
    REAL8 chi1Lz = chi1temp->data[2];
    REAL8 chi2Lx = chi2temp->data[0];
    REAL8 chi2Ly = chi2temp->data[1];
    REAL8 chi2Lz = chi2temp->data[2];

    /* finalSpin interpolation is available only between -0.9996 and 0.9996 */
    /* Set finalSpin to +/- 0.9996 if it is out of this range */
    finalS = finalSpin;
    finalM = finalMass;

    if (finalS < -0.9996)
        finalS = -0.9996;
    if (finalS > 0.9996)
        finalS = 0.9996;

    // if (debug) {
    // XLAL_PRINT_INFO("In RD attachment: final mass = %e, final spin = %e, "
    //                 "total_mass = %e \n",
    //                 finalM, finalS, mTotal);
    // }

    /* Compute QNM frequencies */
    /* Generate 1 overtone (the 0th overtone) */
    // NOTE: this is redone internally in XLALSimIMREOBAttachFitRingdown

    *sigmaQNMlm0 = CreateCOMPLEX16Vector(nmodes);
    COMPLEX16Vector sigmaQNMlm0physicalVec;
    sigmaQNMlm0physicalVec.length = 1;
    COMPLEX16 sigmaQNMlm0physicalval = 0.;
    sigmaQNMlm0physicalVec.data = &sigmaQNMlm0physicalval;
    for ( nmode = 0; nmode < nmodes; nmode++) 
    {

        INT4 l = modes[nmode][0];
        INT4 m = modes[nmode][1];

        /* NOTE: XLALSimIMREOBGenerateQNMFreqV2Prec returns the complex frequency in
            * physical units... */
        if (XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec(&sigmaQNMlm0physicalVec, m1,
            m2, finalM, finalS, l, m,
            1) == CEV_FAILURE) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec for mode (l,m) = (%d,%d).", l, m);
            failed = 1;
            goto QUIT;
        }

        (*sigmaQNMlm0)->data[nmode] = mTScaled * sigmaQNMlm0physicalVec.data[0];

        // if (debug) {
        //     XLAL_PRINT_INFO("Complex QNM frequency: (%d,%d,0) = %.16e + I*%.16e\n", l,
        //                     m, creal((*sigmaQNMlm0)->data[nmode]),
        //                     cimag((*sigmaQNMlm0)->data[nmode]));
        // }
    }

    /* Find the time sample that is closest to tAttach */
    REAL8 timeNeartAttach = FindClosestValueInIncreasingVector(&timeVec, tAttach);
    // R.C: The attachment point is different for the 55 mode wrt the other modes
    // they are related by tAttach55 = tAttach -10M
    REAL8 timeNeartAttach55 =
        FindClosestValueInIncreasingVector(&timeVec, tAttach - 10.);
    /* This structure is inherited from the time when the attachment was done with
    * a comb */
    /* Only the value data[1] will be used -- the other two are ignored */
    REAL8Vector rdMatchPoint;
    REAL8 rdMatchPointdata[4] = {0.};
    rdMatchPoint.length = 3;
    rdMatchPoint.data = rdMatchPointdata;
    rdMatchPoint.data[0] = 0.; /* unused */
    rdMatchPoint.data[1] =
        timeNeartAttach;       /* this is the time closest to tAttach */
    rdMatchPoint.data[2] = 0.; /* unused */
    rdMatchPoint.data[3] = timeNeartAttach55; /* tAttach55 = tAttach -10M */

    /* Create real and imaginary part vectors for the mode with RD attached */
    if (!(hPlmRe = CreateREAL8Vector(retLen_RDattached))) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failed to allocate REAL8Vector hPlmRe.");
        failed = 1;
        goto QUIT;
    }
    if (!(hPlmIm = CreateREAL8Vector(retLen_RDattached))) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failed to allocate REAL8Vector hPlmIm.");
        failed = 1;
        goto QUIT;
    }
    memset(hPlmRe->data, 0, (hPlmRe->length) * sizeof(REAL8));
    memset(hPlmIm->data, 0, (hPlmIm->length) * sizeof(REAL8));

    /* This is used to keep track of the location of the max amplitude of the 22
    * mode */
    UINT indAmpMax = 0;

    /* Loop over modes */
    for ( nmode = 0; nmode < nmodes; nmode++) 
    {

        INT4 l = modes[nmode][0];
        INT4 m = modes[nmode][1];

        /* Get the relevant modes in lists */
        CAmpPhaseSequence *hPlm =
            SphHarmListCAmpPhaseSequence_GetMode(listhPlm, l, m)->campphase;
        CAmpPhaseSequence *hPlm_RDattached =
            SphHarmListCAmpPhaseSequence_GetMode(*listhPlm_RDattached, l, m)->campphase;

        COMPLEX16 hPlm_val = 0.;
        for ( i = 0; i < retLen; i++) 
        {
            hPlm_val = (hPlm->camp_real->data[i] + I * hPlm->camp_imag->data[i]) *
                        cexp(I * hPlm->phase->data[i]);
            hPlmRe->data[i] = creal(hPlm_val);
            hPlmIm->data[i] = cimag(hPlm_val);
        }

        /* NOTE: deltaT here is in physical units (s) */
        REAL8 deltaTseconds = deltaT * mTScaled;
        if (XLALSimIMREOBAttachFitRingdown(
            hPlmRe, hPlmIm, l, m, deltaTseconds, m1, m2, chi1Lx, chi1Ly, chi1Lz,
            chi2Lx, chi2Ly, chi2Lz, finalM, finalS, &timeVec, &rdMatchPoint,
            &indAmpMax) == CEV_FAILURE) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "XLALSimIMREOBAttachFitRingdown failed for mode (l,m) = (%d,%d).", l, m);
            failed = 1;
            goto QUIT;
        }

        /* Copy times in output, using deltaT */
        for ( i = 0; i < retLen; i++) 
        {
            hPlm_RDattached->xdata->data[i] = hPlm->xdata->data[i];
        }
        for ( i = 0; i < retLen_RDattached - retLen; i++) 
        {
            hPlm_RDattached->xdata->data[retLen + i] =
                hPlm_RDattached->xdata->data[retLen - 1] + i * deltaT;
        }

        /* Translate result in amp/phase, unwrap phase in place */
        for ( i = 0; i < retLen_RDattached; i++) 
        {
            hPlm_val = hPlmRe->data[i] + I * hPlmIm->data[i];
            hPlm_RDattached->camp_real->data[i] = cabs(hPlm_val);
            hPlm_RDattached->camp_imag->data[i] = 0.;
            hPlm_RDattached->phase->data[i] = carg(hPlm_val);
        }
        XLALREAL8VectorUnwrapAngle(hPlm_RDattached->phase, hPlm_RDattached->phase);
    }
QUIT:
    /* Cleanup */
    STRUCTFREE(hPlmRe, REAL8Vector);
    STRUCTFREE(hPlmIm, REAL8Vector);
    STRUCTFREE(chi1temp, REAL8Vector);
    STRUCTFREE(chi2temp, REAL8Vector);
    if (failed)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}

int SEOBPrecJoinTimeVector(
    REAL8Vector **tVecPmodes, /**<< Output: vector of times for P-modes
                                 (AdaS+HiS+RDpatch) */
    UINT *retLenPmodes,      /**<< Output: length of output vector of times for
                                 P-modes */
    REAL8 *tJoinHiS,          /**<< Output: first time >= tstartHiS */
    UINT *indexJoinHiS,      /**<< Output: first index >= tstartHiS */
    REAL8 *tJoinAttach,       /**<< Output: first time >= tAttach */
    UINT *indexJoinAttach,   /**<< Output: first index >= tAttach */
    UINT retLenHiSRDpatch,   /**<< Input: length of RD patch to be added at the
                                 end of HiS with the same constant sampling */
    REAL8 deltaTHiS,          /**<< Input: time step for the high sampling */
    REAL8 tstartHiS,          /**<< Input: time of start of HiS */
    REAL8 tAttach,            /**<< Input: time of attachment */
    SEOBPrecdynamics
        *seobdynamicsAdaS, /**<< Input: SEOB dynamics with adaptive-sampling */
    SEOBPrecdynamics
        *seobdynamicsHiS /**<< Input: SEOB dynamics with high-sampling */
) 
{
    /* Read from inputs */
    INT4 lenAdaS = seobdynamicsAdaS->length;
    REAL8 *tAdaS = seobdynamicsAdaS->tVec;
    INT4 lenHiS = seobdynamicsHiS->length;
    REAL8 *tHiS = seobdynamicsHiS->tVec;

    /* Determine how many time samples of AdaS to include */
    UINT iJoinHiS = 0, i;
    while ((iJoinHiS < lenAdaS - 1) && (tAdaS[iJoinHiS] < tstartHiS))
        iJoinHiS++;

    /* Total length of output and create output vector */
    INT4 lenPmodes = iJoinHiS + lenHiS + retLenHiSRDpatch;
    if (!((*tVecPmodes) = CreateREAL8Vector(lenPmodes))) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failed to allocate REAL8Vector tVecPmodes.");
        return CEV_FAILURE;
    }
    *retLenPmodes = lenPmodes;

    /* Copy time values for AdaS and HiS */
    memcpy(&((*tVecPmodes)->data[0]), tAdaS, iJoinHiS * sizeof(REAL8));
    memcpy(&((*tVecPmodes)->data[iJoinHiS]), tHiS, lenHiS * sizeof(REAL8));

    /* Set time values for RD patch */
    for ( i = 0; i < retLenHiSRDpatch; i++) 
    {
        (*tVecPmodes)->data[iJoinHiS + lenHiS + i] =
            tHiS[lenHiS - 1] + (i + 1) * deltaTHiS;
    }

    /* Determine how many time samples after tAttach */
    INT4 iJoinAttach = lenPmodes - 1;
    while ((iJoinAttach > 0) && ((*tVecPmodes)->data[iJoinAttach] > tAttach))
        iJoinAttach--;

    /* Output joining indices and times */
    *indexJoinHiS = iJoinHiS;
    *tJoinHiS = (*tVecPmodes)->data[iJoinHiS];
    *indexJoinAttach = iJoinAttach;
    *tJoinAttach = (*tVecPmodes)->data[iJoinAttach];

    return CEV_SUCCESS;
}


int SEOBPrecJoinDynamics(
    SEOBPrecdynamics **seobdynamicsJoined,     /**<< Output: pointer to joined dynamics */
    SEOBPrecdynamics *seobdynamics1, /**<< Input: first dynamics */
    SEOBPrecdynamics *seobdynamics2, /**<< Input: second dynamics */
    UINT indexJoin12, /**<< Input: index where to join the two dynamics */
    UINT indexEnd2    /**<< Input: index of the joined dynamics where to stop
                          dynamics 2 (excluded) */
) 
{
    UINT j;

    /* Lengths */
    INT lenJoined = indexEnd2;
    INT lendyn1Joined = indexJoin12;
    INT lendyn2Joined = indexEnd2 - indexJoin12;
    INT lendyn1 = seobdynamics1->length;
    INT lendyn2 = seobdynamics2->length;

    /* Initialize output dynamics */
    *seobdynamicsJoined = CreateSEOBPrecdynamics(lenJoined);

    /* Copy truncated dynamics 1 - v4PdynamicsVariables data fields */
    for (j = 0; j < v4PrecdynamicsVariables; j++) 
    {
        memcpy(&((*seobdynamicsJoined)->array->data[j * lenJoined]),
            &(seobdynamics1->array->data[j * lendyn1]),
            lendyn1Joined * sizeof(REAL8));
    }

    /* Copy truncated dynamics 2 - v4PdynamicsVariables data fields */
    for (j = 0; j < v4PrecdynamicsVariables; j++) 
    {
        memcpy(&((*seobdynamicsJoined)->array->data[j * lenJoined + lendyn1Joined]),
            &(seobdynamics2->array->data[j * lendyn2]),
            lendyn2Joined * sizeof(REAL8));
    }

    return CEV_SUCCESS;
}


int SEOBPrecEulerJ2PFromDynamics(
    REAL8Vector **alphaJ2P, /**<< Output: pointer to vector for alpha J2P */
    REAL8Vector **betaJ2P,  /**<< Output: pointer to vector for beta J2P */
    REAL8Vector **gammaJ2P, /**<< Output: pointer to vector for gamma J2P */
    REAL8Vector *e1J,       /**<< Input: unit Jframe vector e1J */
    REAL8Vector *e2J,       /**<< Input: unit Jframe vector e2J */
    REAL8Vector *e3J,       /**<< Input: unit Jframe vector e3J */
    UINT retLen, /**<< Input: total length of Euler angles data to be allocated
                     (length of P-modes) */
    UINT indexStop, /**<< Input: index where we stop the computation (excluded,
                        index of time of attachment) */
    SEOBPrecdynamics *seobdynamics, /**<<Input: SEOB dynamics (joined AdaS+HiS, up
                                   to tAttach) */
    SpinEOBParams *seobParams  /**<< SEOB params */
) 
{
    UINT i, j;
    UINT dynlength = seobdynamics->length;
    UINT SpinsAlmostAligned = seobParams->alignedSpins;

    /* Length of the subset where we want to compute the Euler angles from the
    * dynamics */
    UINT retLenDyn = indexStop;
    /* Check lengths -- note that indexStop is excluded */
    if (!((retLenDyn <= dynlength) && (dynlength <= retLen))) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "incompatible lengths.");
        return CEV_FAILURE;;
    }

    /* Create output vectors */
    if (!((*alphaJ2P) = CreateREAL8Vector(retLen))) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL,
            "failed to allocate REAL8Vector alphaJ2P.");
        return CEV_FAILURE;
    }
    if (!((*betaJ2P) = CreateREAL8Vector(retLen))) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failed to allocate REAL8Vector betaJ2P.");
        return CEV_FAILURE;
    }
    if (!((*gammaJ2P) = CreateREAL8Vector(retLen))) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failed to allocate REAL8Vector gammaJ2P.");
        return CEV_FAILURE;
    }
    memset((*alphaJ2P)->data, 0, retLen * sizeof(REAL8));
    memset((*betaJ2P)->data, 0, retLen * sizeof(REAL8));
    memset((*gammaJ2P)->data, 0, retLen * sizeof(REAL8));

    if (1) /* if spins are almost aligned, leave all Euler
                                angles to 0 */
    {
        /* Time vector of SEOB dynamics */
        REAL8 *tVec = seobdynamics->tVec;

        /* Local variables */
        REAL8 rvec[3] = {0, 0, 0};
        REAL8 pvec[3] = {0, 0, 0};
        REAL8 rdotvec[3] = {0, 0, 0};
        REAL8 Lhat[3] = {0, 0, 0};
        REAL8 LNhat[3] = {0, 0, 0};
        REAL8 Zframe[3] = {0, 0, 0};
        REAL8 e1PiniIbasis[3] = {0, 0, 0};
        REAL8 e2PiniIbasis[3] = {0, 0, 0};
        REAL8 e3PiniIbasis[3] = {0, 0, 0};
        REAL8 e1PiniJbasis[3] = {0, 0, 0};
        REAL8 e2PiniJbasis[3] = {0, 0, 0};

        /* Loop over time samples to compute alpha and beta, stopping at retLenDyn
        */
        for (i = 0; i < retLenDyn; i++) 
        {

        /* Read from the extended dynamics values */
            for (j = 0; j < 3; j++) 
            {
                rvec[j] = seobdynamics->array->data[(1 + j) * dynlength + i];
                pvec[j] = seobdynamics->array->data[(4 + j) * dynlength + i];
                rdotvec[j] = seobdynamics->array->data[(15 + j) * dynlength + i];
            }

            /* Compute Z-axis of the precessing frame, L or LN */
            if (seobParams->hParams->flagZframe == FLAG_SEOBNRv4P_ZFRAME_L) 
            {
                /* Compute Lhat */
                cross_product3d(rvec, pvec, Lhat);
                REAL8 Lhatnorm = sqrt(inner_product3d(Lhat, Lhat));
                for (j = 0; j < 3; j++) 
                {
                    Lhat[j] /= Lhatnorm;
                }
                memcpy(Zframe, Lhat, 3 * sizeof(REAL8));
            }
            else if (seobParams->hParams->flagZframe == FLAG_SEOBNRv4P_ZFRAME_LN) 
            {
                /* Compute LNhat */
                cross_product3d(rvec, rdotvec, LNhat);
                REAL8 LNhatnorm = sqrt(inner_product3d(LNhat, LNhat));
                for (j = 0; j < 3; j++) 
                {
                    LNhat[j] /= LNhatnorm;
                }
                memcpy(Zframe, LNhat, 3 * sizeof(REAL8));
            }
            else 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "flagZframe not recognized.");
                return CEV_FAILURE;
            }

            /* Compute Z projected in the J-frame (e1J,e2J,e3J) */
            REAL8 Ze1J = inner_product3d(Zframe, e1J->data);
            REAL8 Ze2J = inner_product3d(Zframe, e2J->data);
            REAL8 Ze3J = inner_product3d(Zframe, e3J->data);

            /* Get Euler angles alpha (to be unwrapped later) and beta */
            (*alphaJ2P)->data[i] = atan2(Ze2J, Ze1J);
            (*betaJ2P)->data[i] = acos(Ze3J);

            /* At initial time, compute the initial vectors (e1P, e2P, e3P) decomposed
            * in the frame (e1J, e2J, e3J) */
            /* This will be needed to set initial gamma angle */
            if (i == 0)
            {
                /* e3P is the Z axis of the precessing frame */
                memcpy(e3PiniIbasis, Zframe, 3 * sizeof(REAL8));
                /* e1P is the unit separation vector n */
                memcpy(e1PiniIbasis, rvec, 3 * sizeof(REAL8));
                REAL8 e1PiniIbasisnorm = sqrt(inner_product3d(e1PiniIbasis, e1PiniIbasis));
                for (j = 0; j < 3; j++) 
                {
                    e1PiniIbasis[j] /= e1PiniIbasisnorm;
                }
                /* e2P is obtained by completing the triad */
                cross_product3d(e3PiniIbasis, e1PiniIbasis, e2PiniIbasis);
                /* Components of vectors eP in the frame eJ */
                e1PiniJbasis[0] = inner_product3d(e1PiniIbasis, e1J->data);
                e1PiniJbasis[1] = inner_product3d(e1PiniIbasis, e2J->data);
                e1PiniJbasis[2] = inner_product3d(e1PiniIbasis, e3J->data);
                e2PiniJbasis[0] = inner_product3d(e2PiniIbasis, e1J->data);
                e2PiniJbasis[1] = inner_product3d(e2PiniIbasis, e2J->data);
                e2PiniJbasis[2] = inner_product3d(e2PiniIbasis, e3J->data);
            }
        }

        /* Unwrap alpha in-place */
        XLALREAL8VectorUnwrapAngle((*alphaJ2P), (*alphaJ2P));

        /* Compute gamma according to minimal rotation condition */
        /* gamma is set initially so that the initial P-frame reproduces the initial
        * (n, lambda, Lhat or LNhat) frame */
        REAL8 InitialGamma = atan2(e2PiniJbasis[2], -e1PiniJbasis[2]);

        // Integrate \dot{\alpha} \cos{\beta} to get the final Euler angle
        // Eq. 20 of PRD 89, 084006 (2014) [arXiv:1307.6232]

        // Here 1000 referes to the number of subintervals that can be used when
        // performing adaptive quadrature to compute the integral.
        // See
        // https://www.gnu.org/software/gsl/manual/html_node/QAG-adaptive-integration.html
        REAL8 precEulerresult = 0, precEulererror = 0;
        gsl_integration_workspace *precEulerw = gsl_integration_workspace_alloc(1000);
        gsl_function precEulerF;
        PrecEulerAnglesIntegration precEulerparams;
        gsl_spline *x_spline = gsl_spline_alloc(gsl_interp_cspline, retLenDyn);
        gsl_spline *y_spline = gsl_spline_alloc(gsl_interp_cspline, retLenDyn);
        gsl_interp_accel *x_acc = gsl_interp_accel_alloc();
        gsl_interp_accel *y_acc = gsl_interp_accel_alloc();
        gsl_spline_init(x_spline, seobdynamics->tVec, (*alphaJ2P)->data, retLenDyn);
        gsl_spline_init(y_spline, seobdynamics->tVec, (*betaJ2P)->data, retLenDyn);

        precEulerparams.alpha_spline = x_spline;
        precEulerparams.alpha_acc = x_acc;
        precEulerparams.beta_spline = y_spline;
        precEulerparams.beta_acc = y_acc;

        precEulerF.function = &f_alphadotcosi;
        precEulerF.params = &precEulerparams;

        for (i = 0; i < retLenDyn; i++) 
        {
            if (i == 0) 
            {
                (*gammaJ2P)->data[i] = InitialGamma;
            } 
            else 
            {
                gsl_integration_qags(&precEulerF, tVec[i - 1], tVec[i], 1e-9, 1e-9,
                                    1000, precEulerw, &precEulerresult,
                                    &precEulererror);
                (*gammaJ2P)->data[i] = (*gammaJ2P)->data[i - 1] + precEulerresult;
            }
        }
        gsl_integration_workspace_free(precEulerw);
        gsl_spline_free(x_spline);
        gsl_spline_free(y_spline);
        gsl_interp_accel_free(x_acc);
        gsl_interp_accel_free(y_acc);

        /* Unwrap gamma in-place -- note that with integration, this should not be
        * necessary */
        XLALREAL8VectorUnwrapAngle((*gammaJ2P), (*gammaJ2P));
    }

    return CEV_SUCCESS;
}
