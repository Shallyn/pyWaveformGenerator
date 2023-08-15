/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "pCore.h"
#include "myLog.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_multiroots.h>

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

/**
 * This function finds the value in a vector that is closest to an input value.
 * Assumes the input vector is increasing (typically, times or frequencies of
 * series). Purely for convenience.
 */
REAL8 FindClosestValueInIncreasingVector(
    REAL8Vector *vec, /**<< Input: monotonically increasing vector */
    REAL8 value       /**<< Input: value to look for */
) 
{
  UINT index = FindClosestIndex(vec, value);
  return vec->data[index];
}

REAL8 SEOBCalculatetplspin(REAL8 m1, REAL8 m2, REAL8 eta, REAL8 chi1dotZ, REAL8 chi2dotZ)
{
    REAL8 chiS, chiA, tplspin;
    chiS = SEOBCalculateChiS(chi1dotZ, chi2dotZ);
    chiA = SEOBCalculateChiA(chi1dotZ, chi2dotZ);
    tplspin = (1. - 2. * eta) * chiS + (m1 - m2) / (m1 + m2) * chiA;
    return tplspin;
}

static INT XLALEOBSpinPrecCalcSEOBHCoeffConstants(REAL8 eta, SEOBHCoeffConstants *ret)
{
    if (!ret)
        return CEV_FAILURE;
    REAL8 c20  = 1.712;
    REAL8 c21  = -1.803949138004582;
    REAL8 c22  = -39.77229225266885;
    REAL8 c23  = 103.16588921239249;
    REAL8 coeffsKK =  c20 + c21*eta + c22*eta*eta + c23*eta*eta*eta;
    REAL8 m1PlusEtaKK = -1. + eta*coeffsKK;
    REAL8 tmp4=(m1PlusEtaKK*m1PlusEtaKK*m1PlusEtaKK);
    REAL8 tmp10=m1PlusEtaKK*m1PlusEtaKK;
    REAL8 tmp7=(m1PlusEtaKK*m1PlusEtaKK*m1PlusEtaKK*m1PlusEtaKK);
    REAL8 tmp6=coeffsKK*coeffsKK;
    REAL8 tmp16=pow(m1PlusEtaKK,5.);
    REAL8 tmp19=pow(m1PlusEtaKK,6.);
    REAL8 tmp18=(coeffsKK*coeffsKK*coeffsKK);
    REAL8 tmp23=4.*coeffsKK*tmp7;
    REAL8 tmp34=pow(m1PlusEtaKK,7.);
    REAL8 tmp28=CST_PI*CST_PI;
    REAL8 tmp32=16.*coeffsKK*tmp16;
    REAL8 tmp37=pow(m1PlusEtaKK,8.);
    REAL8 tmp36=(coeffsKK*coeffsKK*coeffsKK*coeffsKK);
    REAL8 tmp64=pow(m1PlusEtaKK,9.);
    REAL8 kC0=0.+4.*coeffsKK*tmp4+2.*tmp6*tmp7;
    REAL8 kC1=1.*coeffsKK*tmp10-coeffsKK*tmp4;
    REAL8 kC2=0.+2.*tmp10-1.3333333333333335*tmp18*tmp19-8.*tmp16*tmp6-8.*coeffsKK*tmp7;
    REAL8 kC3=tmp23-2.*coeffsKK*tmp4+2.*tmp16*tmp6-2.*tmp6*tmp7;
    REAL8 kC4=0.+31.333333333333332*tmp10-1.28125*tmp10*tmp28+tmp32+8.*tmp18*tmp34+0.6666666666666666*tmp36*tmp37-4.*tmp4+24.*tmp19*tmp6-4.*coeffsKK*tmp7;
    REAL8 kC5=-12.*coeffsKK*tmp16+2.*tmp18*tmp19+tmp23-2.*tmp18*tmp34+8.*tmp16*tmp6-12.*tmp19*tmp6;
    REAL8 kC6=1.*coeffsKK*tmp16-tmp16*tmp6+0.5*tmp19*tmp6-coeffsKK*tmp7+0.5*tmp6*tmp7;
    REAL8 kC7=-35.12753102199746*tmp10+25.6*CST_GAMMA*tmp10-32.*coeffsKK*tmp19+4.443359375*tmp10*tmp28+tmp32-32.*tmp18*tmp37-62.666666666666664*tmp4+2.5625*tmp28*tmp4+4.*tmp19*tmp6-64.*tmp34*tmp6-5.333333333333334*tmp36*tmp64+8.*tmp7-62.666666666666664*coeffsKK*tmp7+2.5625*coeffsKK*tmp28*tmp7-0.2666666666666661*pow(coeffsKK,5.)*pow(m1PlusEtaKK,10.);
    REAL8 kC8=-10.*coeffsKK*tmp16+32.*coeffsKK*tmp19-12.*tmp18*tmp34+16.*tmp18*tmp37-1.3333333333333337*tmp36*tmp37-24.*tmp19*tmp6+48.*tmp34*tmp6+1.3333333333333337*tmp36*tmp64-2.*tmp7+2.*coeffsKK*tmp7;
    REAL8 kC9=4.*coeffsKK*tmp16-6.*coeffsKK*tmp19-tmp18*tmp19+2.*tmp18*tmp34-tmp18*tmp37-2.*tmp16*tmp6+8.*tmp19*tmp6-6.*tmp34*tmp6;
    ret->a0k2 = kC0;
    ret->a1k2 = kC1;
    ret->a0k3 = kC2;
    ret->a1k3 = kC3;
    ret->a0k4 = kC4;
    ret->a1k4 = kC5;
    ret->a2k4 = kC6;
    ret->a0k5 = kC7;
    ret->a1k5 = kC8;
    ret->a2k5 = kC9;
    return CEV_SUCCESS;
}

// Calc...
static REAL8 CalPN_calculateAmpDot(REAL8 re, REAL8 im, REAL8 reDot, REAL8 imDot)
{
    return (re * reDot + im * imDot) / (sqrt(re*re + im*im));
}

static REAL8 CalPN_calculateAmpDDot(REAL8 re, REAL8 im, REAL8 reDot, REAL8 imDot, REAL8 reDDot, REAL8 imDDot)
{
    return ((im*imDDot + pow(imDot,2))*pow(re,2) + pow(re,3)*reDDot + im*re*(im*reDDot - 2*imDot*reDot) + 
     pow(im,2)*(im*imDDot + pow(reDot,2)))/pow(pow(im,2) + pow(re,2),1.5);
}

static REAL8 CalPN_calculateOmega(REAL8 re, REAL8 im, REAL8 reDot, REAL8 imDot)
{
    return (re * imDot - im * reDot) / (re*re + im*im);
}

static REAL8 CalPN_calculateOmegaDot(REAL8 re, REAL8 im, REAL8 reDot, REAL8 imDot, REAL8 reDDot, REAL8 imDDot)
{
    return ((pow(im,2) + pow(re,2))*(imDDot*re - im*reDDot) - 2*(imDot*re - im*reDot)*(im*imDot + re*reDot))/pow(pow(im,2) + pow(re,2),2);
}
typedef
struct tagCalPNCoeffsRootParams
{
    REAL8 eAmp;
    REAL8 eAmpDot;
    REAL8 eAmpDDot;
    REAL8 eOmega;
    REAL8 eOmegaDot;
    COMPLEX16 rholm;
    COMPLEX16 rholmDot;
    COMPLEX16 rholmDDot;
    REAL8 nrAmp;
    REAL8 nrAmpDot;
    REAL8 nrAmpDDot;
    REAL8 nrOmega;
    REAL8 nrOmegaDot;
    REAL8 x0;
    REAL8 x0Dot;
    REAL8 x0DDot;
    REAL8 x1;
    REAL8 x1Dot;
    REAL8 x1DDot;
    REAL8 x2;
    REAL8 x2Dot;
    REAL8 x2DDot;
    int order;
}
CalPNCoeffsRootParams;

static int gslfuncFindCalPNCoeffs(const gsl_vector *x, void *params, gsl_vector *f)
{
    CalPNCoeffsRootParams *rootParams = (CalPNCoeffsRootParams *) params;
    REAL8 k1, k2, k3, k4, k5;
    REAL8 eAmp, eAmpDot, eAmpDDot, eOmega, eOmegaDot, nrAmp, nrOmega, nrAmpDot, nrAmpDDot, nrOmegaDot;
    REAL8 x0, x0Dot, x0DDot, x1, x1Dot, x1DDot, x2, x2Dot, x2DDot;
    REAL8 eobAmp, eobOmega, eobAmpDot, eobOmegaDot, eobAmpDDot;
    REAL8 fr, fi, frDot, fiDot, frDDot, fiDDot;
    REAL8 Rpart, Ipart, RpartDot, IpartDot, RpartDDot, IpartDDot;
    REAL8 Rtot, RtotDot, RtotDDot, Itot, ItotDot, ItotDDot;
    k1 = gsl_vector_get(x, 0);
    k2 = gsl_vector_get(x, 1);
    k3 = gsl_vector_get(x, 2);
    k4 = gsl_vector_get(x, 3);
    k5 = gsl_vector_get(x, 4);
    x0 = rootParams->x0;
    x0Dot = rootParams->x0Dot;
    x0DDot = rootParams->x0DDot;
    x1 = rootParams->x1;
    x1Dot = rootParams->x1Dot;
    x1DDot = rootParams->x1DDot;
    x2 = rootParams->x2;
    x2Dot = rootParams->x2Dot;
    x2DDot = rootParams->x2DDot;
    COMPLEX16 comb, combDot, combDDot;
    REAL8 pre, preDot, preDDot;
    int order = rootParams->order;
    comb = I*(k1*x0*x1 + k2*x1*x2) + k3*x0*x2 + k4*x0*x0 + k5*x2*x2;
    combDot = 2*k4*x0*x0Dot + I*k1*x0Dot*x1 + I*k1*x0*x1Dot + k3*x0Dot*x2 + I*k2*x1Dot*x2 + k3*x0*x2Dot + I*k2*x1*x2Dot + 2*k5*x2*x2Dot;
    combDDot = 2*k4*x0*x0DDot + 2*k4*pow(x0Dot,2) + I*k1*x0DDot*x1 + I*k1*x0*x1DDot + I*k1*x0Dot*x1Dot + k3*x0DDot*x2 + 
        I*k2*x1DDot*x2 + k3*x0*x2DDot + I*k2*x1*x2DDot + 2*k5*x2*x2DDot + 2*k3*x0Dot*x2Dot + I*k2*x1Dot*x2Dot + 
        2*k5*pow(x2Dot,2);
    // print_debug("comb = %f + i%f\n", comb);
    // print_debug("combDot = %f + i%f\n", combDot);
    // print_debug("combDDot = %f + i%f\n\n", combDDot);
    pre = pow(x0, order);
    preDot = order * pow(x0, order-1) * x0Dot;
    preDDot = order * (order - 1) * pow(x0, order-2) * x0Dot * x0Dot + order * pow(x0, order-1) * x0DDot;
    

    eAmp = rootParams->eAmp;
    eAmpDot = rootParams->eAmpDot;
    eAmpDDot = rootParams->eAmpDDot;
    eOmega = rootParams->eOmega;
    eOmegaDot = rootParams->eOmegaDot;

    fr = creal(rootParams->rholm);
    fi = cimag(rootParams->rholm);
    frDot = creal(rootParams->rholmDot);
    fiDot = cimag(rootParams->rholmDot);
    frDDot = creal(rootParams->rholmDDot);
    fiDDot = cimag(rootParams->rholmDDot);

    nrAmp = rootParams->nrAmp;
    nrAmpDot = rootParams->nrAmpDot;
    nrAmpDDot = rootParams->nrAmpDDot;

    Rpart = pre * creal(comb);
    Ipart = pre * cimag(comb);

    RpartDot = creal(pre * combDot + preDot * comb);
    IpartDot = cimag(pre * combDot + preDot * comb);

    RpartDDot = creal(pre * combDDot + 2*preDot*combDot + preDDot * comb);
    IpartDDot = cimag(pre * combDDot + 2*preDot*combDot + preDDot * comb);

    Rtot = fr + Rpart;
    Itot = fi + Ipart;

    RtotDot = frDot + RpartDot;
    ItotDot = fiDot + IpartDot;

    RtotDDot = frDDot + RpartDDot;
    ItotDDot = fiDDot + IpartDDot;

    nrOmega = rootParams->nrOmega;
    nrOmegaDot = rootParams->nrOmegaDot;
    eobAmp = eAmp * sqrt(Rtot*Rtot + Itot*Itot);
    eobOmega = eOmega + CalPN_calculateOmega(Rtot, Itot, RtotDot, ItotDot);
    eobAmpDot = eAmp * CalPN_calculateAmpDot(Rtot, Itot, RtotDot, ItotDot) + eAmpDot * sqrt(Rtot*Rtot + Itot*Itot);
    eobOmegaDot = eOmegaDot + CalPN_calculateOmegaDot(Rtot, Itot, RtotDot, ItotDot, RtotDDot, ItotDDot);
    eobAmpDDot = eAmp * CalPN_calculateAmpDDot(Rtot, Itot, RtotDot, ItotDot, RtotDDot, ItotDDot) + 
        2*eAmpDot * CalPN_calculateAmpDot(Rtot, Itot, RtotDot, ItotDot) + eAmpDDot * sqrt(Rtot*Rtot + Itot*Itot);
    // print_debug("IN: nra = %f, ea = %f, eoba = %f\n", nrAmp, eAmp, eobAmp);
    // print_debug("IN: nrOmega = %f, eO = %f, eobOmega = %f\n", nrOmega, eOmega, eobOmega);
    // print_debug("IN: nrAmpDot = %f, eAD = %f, eobAmpDot = %f\n", nrAmpDot, eAmpDot, eobAmpDot);
    // print_debug("IN: nrOmegaDot = %f, eOD = %f, eobOmegaDot = %f\n\n", nrOmegaDot, eOmegaDot, eobOmegaDot);
    // print_debug("IN: nrAmpDDot = %f, eADD = %f, eobAmpDDot = %f\n", nrAmpDDot, eAmpDDot, eobAmpDDot);
    // print_debug("RST: x = %f + i%f\n", yr, yi);
    gsl_vector_set( f, 0, nrAmp - eobAmp);
    gsl_vector_set( f, 1, nrOmega - eobOmega);
    gsl_vector_set( f, 2, nrAmpDot - eobAmpDot);
    gsl_vector_set( f, 3, nrOmegaDot - eobOmegaDot);
    gsl_vector_set( f, 4, nrAmpDDot - eobAmpDDot);
    return CEV_SUCCESS;
}

static INT FindCalPNCoeffs(COMPLEX16 hLM,
                           COMPLEX16 hLMDot,
                           COMPLEX16 hLMDDot,
                           COMPLEX16 rholm,
                           COMPLEX16 rholmDot,
                           COMPLEX16 rholmDDot,
                           REAL8 nrAmp,
                           REAL8 nrAmpDot,
                           REAL8 nrAmpDDot,
                           REAL8 nrOmega,
                           REAL8 nrOmegaDot,
                           REAL8 x0,
                           REAL8 dx0dt,
                           REAL8 dx0dt2,
                           REAL8 x1,
                           REAL8 dx1dt,
                           REAL8 dx1dt2,
                           REAL8 x2,
                           REAL8 dx2dt,
                           REAL8 dx2dt2,
                           INT order,
                           REAL8 initval,
                           REAL8 *out1,
                           REAL8 *out2,
                           REAL8 *out3,
                           REAL8 *out4,
                           REAL8 *out5)
{
    int i=0, gslStatus;
    const int maxIter = 300;
    const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver *rootSolver = NULL;
    gsl_vector *initValues  = NULL;
    gsl_vector *finalValues = NULL;
    gsl_multiroot_function rootFunction;
    CalPNCoeffsRootParams rootParams;
    REAL8 rst1, rst2, rst3, rst4, rst5;
    REAL8 eAmp, eAmpDot, eOmega;
    eAmp = cabs(hLM);
    eAmpDot = CalPN_calculateAmpDot(creal(hLM), cimag(hLM), creal(hLMDot), cimag(hLMDot));
    eOmega = CalPN_calculateOmega(creal(hLM), cimag(hLM), creal(hLMDot), cimag(hLMDot));
    rootParams.eAmp = eAmp;
    rootParams.eAmpDot = eAmpDot;
    rootParams.eAmpDDot = CalPN_calculateAmpDDot(creal(hLM), cimag(hLM), creal(hLMDot), cimag(hLMDot), creal(hLMDDot), cimag(hLMDDot));
    rootParams.eOmega = eOmega;
    rootParams.eOmegaDot = CalPN_calculateOmegaDot(creal(hLM), cimag(hLM), creal(hLMDot), cimag(hLMDot), creal(hLMDDot), cimag(hLMDDot));

    rootParams.rholm = rholm;
    rootParams.rholmDot = rholmDot;
    rootParams.rholmDDot = rholmDDot;
    REAL8 tgtAmp, tgtAmpDot, tgtOmega;
    tgtAmp = eAmp + ((nrAmp - eAmp > 0) ? 1. : -1) * (fabs(nrAmp - eAmp)*0.1 < eAmp*0.01 ? fabs(nrAmp-eAmp)*0.1 : eAmp*0.01);
    tgtAmpDot = eAmpDot + ((nrAmpDot - eAmpDot > 0) ? 1. : -1) * (fabs(nrAmpDot-eAmpDot)*0.5 < fabs(eAmpDot)*0.1 ? fabs(nrAmpDot-eAmpDot)*0.5 : fabs(eAmpDot)*0.1);
    tgtOmega = eOmega + ((nrOmega - eOmega > 0) ? 1. : -1.) * (fabs(nrOmega - eOmega)*0.2 < fabs(eOmega) * 0.1 ? fabs(nrOmega - eOmega)*0.2 : fabs(eOmega)*0.1);
    rootParams.nrAmp = tgtAmp;
    rootParams.nrAmpDot = tgtAmpDot;
    rootParams.nrAmpDDot = nrAmpDDot;
    rootParams.nrOmega = tgtOmega;
    rootParams.nrOmegaDot = nrOmegaDot;

    rootParams.x0 = x0;
    rootParams.x0Dot = dx0dt;
    rootParams.x0DDot = dx0dt2;
    rootParams.x1 = x1;
    rootParams.x1Dot = dx1dt;
    rootParams.x1DDot = dx1dt2;
    rootParams.x2 = x2;
    rootParams.x2Dot = dx2dt;
    rootParams.x2DDot = dx2dt2;

    rootParams.order = order;
    rootSolver = gsl_multiroot_fsolver_alloc( T, 5 );
    if ( !rootSolver )
    {
        return CEV_FAILURE;
    }

    initValues = gsl_vector_calloc( 5 );
    if ( !initValues )
    {
        gsl_multiroot_fsolver_free( rootSolver );
        return CEV_FAILURE;
    }
    gsl_vector_set( initValues, 0, 0.01);
    gsl_vector_set( initValues, 1, 0.01 );
    gsl_vector_set( initValues, 2, 0.01 );
    gsl_vector_set( initValues, 3, initval );
    gsl_vector_set( initValues, 4, 0.01 );

    rootFunction.f      = gslfuncFindCalPNCoeffs;
    rootFunction.n      = 5;
    rootFunction.params = &rootParams;
    gsl_multiroot_fsolver_set( rootSolver, &rootFunction, initValues );

    do
    {
        gslStatus = gsl_multiroot_fsolver_iterate( rootSolver );
        if ( gslStatus != GSL_SUCCESS )
        {
            print_warning( "Error in GSL iteration function!\n" );
            gsl_multiroot_fsolver_free( rootSolver );
            gsl_vector_free( initValues );
            return CEV_FAILURE;
        }
        
        gslStatus = gsl_multiroot_test_residual( rootSolver->f, 1.0e-10 );
        i++;
    }
    while ( gslStatus == GSL_CONTINUE && i <= maxIter );
    
    if ( i > maxIter && gslStatus != GSL_SUCCESS )
    {
        gsl_multiroot_fsolver_free( rootSolver );
        gsl_vector_free( initValues );
        return CEV_FAILURE;
    }
    
    finalValues = gsl_multiroot_fsolver_root( rootSolver );
    rst1 = gsl_vector_get( finalValues, 0 );
    rst2 = gsl_vector_get( finalValues, 1 );
    // rst1_i = 0;
    rst3 = gsl_vector_get( finalValues, 2 );
    rst4 = gsl_vector_get( finalValues, 3 );
    // rst2_r = 0;
    // rst2_i = 0;
    rst5 = gsl_vector_get(finalValues, 4);
    *out1 = rst1;
    *out2 = rst2;
    *out3 = rst3;
    *out4 = rst4;
    *out5 = rst5;
    // print_debug("FINAL: nra = %f, fiteobAmp = %f, nraDot = %f, fiteobAmpDot = %f\n", nrAmp, cabs(fithLM), nrAmpDot, ampDot);
    gsl_multiroot_fsolver_free( rootSolver );
    gsl_vector_free( initValues );
    // gsl_vector_free( finalValues );
    return CEV_SUCCESS;
}

static INT XLALSimIMREOBCalcCalibCoefficientHigherModesPrec (
               SpinEOBParams * params, /**Output **/
               const UINT modeL, /*<< Mode index L */
               const UINT modeM, /*<< Mode index M */
               SEOBdynamics *seobdynamics, /*<< Dynamics vector */
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
    HamVec.data = seobdynamics->hamVec;
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
        polarDynamics.data[2] = seobdynamics->polarprVec[i];
        polarDynamics.data[3] = seobdynamics->polarpphiVec[i];
        chi1dotZ = s1dotZ * mtot*mtot / (m1*m1);
        chi2dotZ = s2dotZ * mtot*mtot / (m2*m2);
        chiS = 0.5*(chi1dotZ+chi2dotZ);
        chiA = 0.5*(chi1dotZ-chi2dotZ);
        tplspin = (1.-2.*eta) * chiS + (m1 - m2)/(m1 + m2) * chiA;
        v = cbrt (orbOmegaVec.data[i]);
        if (CODE_VERSION == 3)
        {
            if (EccPrec_CalcSpinPrecFacWaveformCoefficients(
                    params->hCoeffs, m1, m2, eta, tplspin, chiS, chiA,
                    seobdynamics->chiSxVec[i], seobdynamics->chiSyVec[i], seobdynamics->chiSzVec[i],
                    seobdynamics->chiAxVec[i], seobdynamics->chiAyVec[i], seobdynamics->chiAzVec[i],
                    451) == CEV_FAILURE)
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in EccPrec_CalcSpinPrecFacWaveformCoefficients at step %d of the loop.", i);
                return CEV_FAILURE;
            }
        }
        else
        {
            if ( XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients( params->hCoeffs, m1, m2, eta, tplspin, chiS, chiA, SEOBWaveformVersion ) == CEV_FAILURE )
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients at step %d of the loop.", i);
                failed = 1; 
                goto QUIT;
            }
        }
        // if (CODE_VERSION == 0)
        // {
        if ( XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform( &(hLMVec->data[i]), &polarDynamics, &values, v, HamVec.data[i], modeL, modeM, params ) == CEV_FAILURE )
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
        // }
        // else
        // {
            dr = (seobdynamics->posVecx[i]*seobdynamics->velVecx[i] + 
                seobdynamics->posVecy[i]*seobdynamics->velVecy[i] + 
                seobdynamics->posVecz[i]*seobdynamics->velVecz[i]) / seobdynamics->polarrVec[i];
            ncrv = seobdynamics->polarrVec[i] * orbOmegaVec.data[i];
            drVec->data[i] = dr;
            ncrvVec->data[i] = ncrv;
            if ( XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2( &hLME, &polarDynamics, &values, v, dr, ncrv, seobdynamics->polarprDotVec[i], HamVec.data[i], modeL, modeM, params ) == CEV_FAILURE )
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2 at step %d of the loop.", i);
                failed = 1; 
                goto QUIT;
            }
        // if (isnan(params->hCoeffs->h21T2ff10))
        // {
        //     print_debug("ca = (%g, %g, %g), cs = (%g, %g, %g)\n", seobdynamics->chiSxVec[i], seobdynamics->chiSyVec[i], seobdynamics->chiSzVec[i],
        //             seobdynamics->chiAxVec[i], seobdynamics->chiAyVec[i], seobdynamics->chiAzVec[i]);
        // }

            if (XLALSimIMRSpinEOBGetAmplitudeResidualPrecV2 (&rholmpwrlE, &polarDynamics, &values, v, dr, ncrv, HamVec.data[i], modeL, modeM, params) == CEV_FAILURE)
            //RC: for the 21 and 55 mode rholmpwrlVec is always real. This is not true for the 33 mode. For this reason, in order to make this function general, we use a complex variable for it.
            {
                /* TODO: Clean-up */
                failed = 1; 
                goto QUIT;
            }
        // if (isnan(params->hCoeffs->h21T2ff10))
        // {
        //     print_debug("it become nan\n");
        // }
            // if (creal(rholmpwrlVec->data[i]) == 0.0)
            //     print_debug("[%d]get\n", i);
        // }
        rholmpwrlVecReal->data[i] = (REAL8)creal(rholmpwrlVec->data[i]);
        hLMdivrholmVec->data[i] = ((REAL8)cabs(hLMVec->data[i]))/fabs(rholmpwrlVecReal->data[i]);

        rholmpwrlEVecReal->data[i] = (REAL8) creal(rholmpwrlE);
        rholmpwrlEVecImag->data[i] = (REAL8) cimag(rholmpwrlE);
        hLMdivrholmEVec->data[i] = ((REAL8)cabs(hLME) / fabs(rholmpwrlEVecReal->data[i]));

        hLMrealVec->data[i] = (REAL8) creal(hLME/rholmpwrlE);
        hLMimagVec->data[i] = (REAL8) cimag(hLME/rholmpwrlE);
        // print_debug("[%d]")
        // if (rholmpwrlVecReal->data[i] == 0.0)
            // print_debug("[%d]rholmPwrlReal = %.3e + i%.3e, v = %g, dr = %g, ncrv = %g \n", 
            //     i, creal(rholmpwrlVec->data[i]), cimag(rholmpwrlVec->data[i]), v, dr, ncrv);
        // if (timeVec->data[i] > 0.9*timewavePeak)
        //     print_debug("t[%d] = %g, rholmpwrlEVec = %.3e + i%.3e\n",
        //             i, timeVec->data[i], rholmpwrlEVecReal->data[i], rholmpwrlEVecImag->data[i]);
            // print_debug("t[%d] = %g, rholmpwrlEVec = %.3e + i%.3e, hLM = %.3e + i%.3e,r = %g, v = %g, dr = %g, ncrv = %g, ham = %g \n", 
            //     i, timeVec->data[i], rholmpwrlEVecReal->data[i], rholmpwrlEVecImag->data[i], creal(hLME), cimag(hLME), seobdynamics->polarrVec[i], v, dr, ncrv, HamVec.data[i]);
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

    status = gsl_spline_init (spline, timeVec->data, hLMdivrholmEVec->data, hLMdivrholmEVec->length);
    if (status != GSL_SUCCESS)
    {
        PRINT_LOG_INFO (LOG_CRITICAL, "gsl_spline_init failed.");
        failed = 1; 
        goto QUIT;
    }
    gsl_interp_accel_reset (acc);
    hLMdivrholmEAttachmentPoint = gsl_spline_eval (spline, timewavePeak, acc);


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
    rholmENRAttachmentPoint = nra/hLMdivrholmEAttachmentPoint;
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

    REAL8 rholmErealDot, rholmEimagDot, rholmErealDDot, rholmEimagDDot;
    status = gsl_spline_init (spline, timeVec->data, rholmpwrlEVecReal->data, rholmpwrlEVecReal->length);
    if (status != GSL_SUCCESS)
    {
        PRINT_LOG_INFO (LOG_CRITICAL, "gsl_spline_init failed.");
        failed = 1; 
        goto QUIT;
    }
    gsl_interp_accel_reset (acc);
    /* rho_lm^l + f_lm before the calibration parameter is set */
    rholmERealBeforeCalAP = gsl_spline_eval (spline, timewavePeak, acc);
    rholmErealDot = gsl_spline_eval_deriv (spline, timewavePeak, acc);
    rholmErealDDot = gsl_spline_eval_deriv2 (spline, timewavePeak, acc);

    status = gsl_spline_init (spline, timeVec->data, rholmpwrlEVecImag->data, rholmpwrlEVecReal->length);
    if (status != GSL_SUCCESS)
    {
        PRINT_LOG_INFO (LOG_CRITICAL, "gsl_spline_init failed.");
        failed = 1; 
        goto QUIT;
    }
    gsl_interp_accel_reset (acc);
    /* rho_lm^l + f_lm before the calibration parameter is set */
    rholmEImagBeforeCalAP = gsl_spline_eval (spline, timewavePeak, acc);
    rholmEimagDot = gsl_spline_eval_deriv (spline, timewavePeak, acc);
    rholmEimagDDot = gsl_spline_eval_deriv2 (spline, timewavePeak, acc);
    COMPLEX16 rholmEBeforeCalAP = rholmERealBeforeCalAP + I*rholmEImagBeforeCalAP;
    COMPLEX16 rholmEDot = rholmErealDot + I*rholmEimagDot;
    COMPLEX16 rholmEDDot = rholmErealDDot + I*rholmEimagDDot;

    status = gsl_spline_init (spline, timeVec->data, hLMrealVec->data, hLMdivrholmVec->length);
    if (status != GSL_SUCCESS)
    {
        PRINT_LOG_INFO (LOG_CRITICAL, "gsl_spline_init failed.");
        failed = 1; 
        goto QUIT;
    }
    gsl_interp_accel_reset (acc);
    hLMrealAP = gsl_spline_eval(spline, timewavePeak, acc);
    hLMrealAPDot = gsl_spline_eval_deriv(spline, timewavePeak, acc);
    hLMrealAPDDot = gsl_spline_eval_deriv2(spline, timewavePeak, acc);

    status = gsl_spline_init (spline, timeVec->data, hLMimagVec->data, hLMdivrholmVec->length);
    if (status != GSL_SUCCESS)
    {
        PRINT_LOG_INFO (LOG_CRITICAL, "gsl_spline_init failed.");
        failed = 1; 
        goto QUIT;
    }
    gsl_interp_accel_reset (acc);
    hLMimagAP = gsl_spline_eval(spline, timewavePeak, acc);
    hLMimagAPDot = gsl_spline_eval_deriv(spline, timewavePeak, acc);
    hLMimagAPDDot = gsl_spline_eval_deriv2(spline, timewavePeak, acc);
    // if(debugRC==1){
    //     FILE *timeomegapeak = NULL;
    //     timeomegapeak = fopen ("timeomegapeak.dat", "w");
    //     fprintf(timeomegapeak, "%.16f\n", timewavePeak);
    //     fclose(timeomegapeak);
    // }
    REAL8 x1AP, x1APDot, x1APDDot, x2AP, x2APDot, x2APDDot;
    status = gsl_spline_init (spline, timeVec->data, drVec->data, drVec->length);
    if (status != GSL_SUCCESS)
    {
        PRINT_LOG_INFO (LOG_CRITICAL, "gsl_spline_init failed.");
        failed = 1; 
        goto QUIT;
    }
    gsl_interp_accel_reset (acc);
    x1AP = gsl_spline_eval(spline, timewavePeak, acc);
    x1APDot = gsl_spline_eval_deriv(spline, timewavePeak, acc);
    x1APDDot = gsl_spline_eval_deriv2(spline, timewavePeak, acc);

    gsl_spline_init (spline, timeVec->data, ncrvVec->data, drVec->length);
    gsl_interp_accel_reset (acc);
    x2AP = gsl_spline_eval(spline, timewavePeak, acc);
    x2APDot = gsl_spline_eval_deriv(spline, timewavePeak, acc);
    x2APDDot = gsl_spline_eval_deriv2(spline, timewavePeak, acc);
    // print_debug("x1 = %f, x1Dot = %f, x1DDot = %f\n", x1AP, x1APDot, x1APDDot);
    // print_debug("x2 = %f, x2Dot = %f, x2DDot = %f\n", x2AP, x2APDot, x2APDDot);
    REAL8 PNCalCoeff1, PNCalCoeff2, PNCalCoeff3, PNCalCoeff4, PNCalCoeff5;
    COMPLEX16 initval;
    COMPLEX16 hLMPeak = hLMrealAP + I*hLMimagAP;
    COMPLEX16 hLMPeakDot = hLMrealAPDot + I*hLMimagAPDot;
    COMPLEX16 hLMPeakDDot = hLMrealAPDDot + I*hLMimagAPDDot;
    int retval;
    if ((modeL == 2) && (modeM == 1))
    {
        /* Here we compute ((rho_lm^l + f_lm + CalPar*omega^7/3)_NR - (rho_lm^l + f_lm)_EOB)/omega^7/3 to get CalPar.
        The factor rholmpwrlVecReal->data[0])/cabs(rholmpwrlVecReal->data[0]) is used to know the sign of the function (rho_lm^l + f_lm + CalPar*omega^7/3)_NR which is computed as absolute value */
        params->cal21 = (rholmNRAttachmentPoint - rholmBeforeCalibAttachmentPoint)/(pow(omegaAttachmentPoint,7./3.));
        initval = (rholmENRAttachmentPoint - rholmEBeforeCalAP)*(pow(rAttachmentPoint,12./2.));
        // initval = (rholmENRAttachmentPoint - rholmEBeforeCalAP) /(pow(omegaAttachmentPoint,7./3.));
        // retval = FindCalPNCoeffs(hLMPeak, hLMPeakDot, hLMPeakDDot, 
        //             rholmEBeforeCalAP, rholmEDot, rholmEDDot,
        //             nra, nraDot, nraDDot, nrOmega, nrOmegaDot, 
        //             x0AP, x0DotAP, x0DDotAP, 
        //             x1AP, x1APDot, x1APDDot, 
        //             x2AP, x2APDot, x2APDDot, 
        //             5, creal(initval), &PNCalCoeff1, &PNCalCoeff2, &PNCalCoeff3, &PNCalCoeff4, &PNCalCoeff5);
        params->cal21E = creal(initval);
        // if (retval != CEV_SUCCESS)
        //     params->cal21E = creal(initval);
        // else
        // {
            // params->cal21E1 = PNCalCoeff1;// * rholmEBeforeCalAP / (hLMrealAP + I*hLMimagAP);
            // params->cal21E2 = PNCalCoeff2;// * rholmEBeforeCalAP / (hLMrealAP + I*hLMimagAP);
            // params->cal21E3 = PNCalCoeff3;
            // params->cal21E = PNCalCoeff4;
            // params->cal21E4 = PNCalCoeff5;
            // print_debug("fitpms = %f, %f, %f, %f, %f, %f\n", PNCalCoeff1, PNCalCoeff2, PNCalCoeff3, PNCalCoeff4, PNCalCoeff5, params->hCoeffs->f21v6);
        // }
        // Test
        // COMPLEX16 fithLMPeak, fithLMDotPeak;
        // REAL8 fitALMDotPeak;
        // REAL8 Re, Im, ReDot, ImDot;
        // fithLMPeak = hLMPeak * (rholmEBeforeCalAP + pow(x0AP, 7) * params->cal21E);
        // fithLMPeak = hLMPeak + pow(x0AP, 7)*PNCalCoeff1 + pow(x0AP, 8) * PNCalCoeff2;
        // fithLMDotPeak = hLMrealAPDot + I*hLMimagAPDot + 7*pow(x0AP, 6)*x0DotAP*PNCalCoeff1 + 8*pow(x0AP, 7)*x0DotAP * PNCalCoeff2;
        // fitALMDotPeak = (cimag(fithLMPeak) * cimag(fithLMDotPeak) + creal(fithLMPeak) * creal(fithLMDotPeak)) / cabs(fithLMPeak);

        // print_debug("nra = %f, eoba = %f, nraDot = %f, eobaDot = %f\n", nra, cabs(fithLMPeak), nraDot, fitALMDotPeak);
                            // print_debug("cal21 = %.16f, nra = %.3f\n", params->cal21, nra);
    }
    if ((modeL == 5) && (modeM == 5))
    {
        /* Here we compute ((rho_lm^l + f_lm + CalPar*omega^7/3)_NR - (rho_lm^l + f_lm)_EOB)/omega^7/3 to get CalPar.
        The factor rholmpwrlVecReal->data[0])/cabs(rholmpwrlVecReal->data[0]) is used to know the sign of the function (rho_lm^l + f_lm + CalPar*omega^7/3)_NR which is computed as absolute value */
        params->cal55 = (rholmNRAttachmentPoint - rholmBeforeCalibAttachmentPoint)/(pow(omegaAttachmentPoint,5./3.));
                            //printf("params->cal55 = %.16f\n",params->cal55);
    }
    //if (isnan(params->cal21) || isnan(params->cal21E1) || isnan(params->cal21E2) || isnan(params->cal21E) || isnan(params->cal55))
    if (isnan(creal(params->cal21)) || isnan(cimag(params->cal21)) || 
        isnan(params->cal21E1) || isnan(params->cal21E2) || isnan(params->cal21E) || 
        isnan(creal(params->cal55)) || isnan(cimag(params->cal55)))
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "(l,m) = (%d, %d) calibration get nan", modeL, modeM);
        PRINT_LOG_INFO(LOG_DEBUG, "rholmENRAttachmentPoint = %g, rholmEBeforeCalAP = %g + i %g\n", rholmENRAttachmentPoint, creal(rholmEBeforeCalAP), cimag(rholmEBeforeCalAP));
        PRINT_LOG_INFO(LOG_DEBUG, "cal21 = %g + i %g, cal21E = %g, rAttachmentPoint = %g\n", creal(params->cal21), cimag(params->cal21), params->cal21E, rAttachmentPoint);
        PRINT_LOG_INFO(LOG_DEBUG, "rholmNRAttachmentPoint = %g, rholmBeforeCalibAttachmentPoint = %g, omegaAttachmentPoint = %g",
            rholmNRAttachmentPoint, rholmBeforeCalibAttachmentPoint, omegaAttachmentPoint);
        PRINT_LOG_INFO(LOG_DEBUG, "nra = %g, hLMdivrholmAttachmentPoint = %g", nra, hLMdivrholmAttachmentPoint);
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


static INT XLALSimIMREOBSACalcCalibCoefficientHigherModes (
               SpinEOBParams * params, /**Output **/
               const UINT modeL, /*<< Mode index L */
               const UINT modeM, /*<< Mode index M */
               SEOBSAdynamics *seobdynamics, /*<< Dynamics vector */
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
    REAL8Vector polarDynamics;
    // REAL8 valuesdata[14] = {0.};
    REAL8 polarDynamicsdata[4] = {0.};
    polarDynamics.length = 4;
                    // values.length = 14;
    polarDynamics.data = polarDynamicsdata;
                    // values.data = valuesdata;
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
    REAL8 spin1z_omegaPeak = params->chi1;
    REAL8 spin2z_omegaPeak = params->chi2;
    REAL8 chiS_omegaPeak = 0.5*(spin1z_omegaPeak+ spin2z_omegaPeak);
    REAL8 chiA_omegaPeak = 0.5*(spin1z_omegaPeak-spin2z_omegaPeak);

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

    orbOmegaVec.data = seobdynamics->dphiVec;
    HamVec.data = seobdynamics->HVec;
    rVec.data = seobdynamics->rVec;

    GSL_START;
    /* Stuff for interpolating function */
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
    int status;
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
    chi1dotZ = params->chi1;
    chi2dotZ = params->chi2;
    chiS = 0.5*(chi1dotZ+chi2dotZ);
    chiA = 0.5*(chi1dotZ-chi2dotZ);
    tplspin = (1.-2.*eta) * chiS + (m1 - m2)/(m1 + m2) * chiA;
    if ( XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients( params->hCoeffs, m1, m2, eta, tplspin, chiS, chiA, SEOBWaveformVersion ) == CEV_FAILURE )
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients at step %d of the loop.", i);
        failed = 1; 
        goto QUIT;
    }

    for(i=0; i<orbOmegaVec.length; i++)
    {
        // print_out("%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n", 
        //     timeVec->data[i], seobdynamics->rVec[i], seobdynamics->phiVec[i], seobdynamics->prTVec[i], seobdynamics->pphiVec[i]);
        // for ( j=0; j<14; j++) 
        // {
        //     values.data[j] = seobdynamics->array->data[i + (j+1)*retLen];
        // }
        // s1dotZ = seobdynamics->s1dotZVec[i];
        // s2dotZ = seobdynamics->s2dotZVec[i];
        polarDynamics.data[0] = seobdynamics->rVec[i];
        polarDynamics.data[1] = seobdynamics->phiVec[i];
        polarDynamics.data[2] = seobdynamics->prTVec[i];
        polarDynamics.data[3] = seobdynamics->pphiVec[i];
        v = cbrt (orbOmegaVec.data[i]);
        // if (CODE_VERSION == 0)
        // {
            // DEBUG_START;
        if ( XLALSimIMRSpinEOBGetSASpinFactorizedWaveform( &(hLMVec->data[i]), &polarDynamics, v, HamVec.data[i], modeL, modeM, params ) == CEV_FAILURE )
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform at step %d of the loop.", i);
            failed = 1; 
            goto QUIT;
        }
        // DEBUG_END;

// print_out("%.16e\t%.16e\t%.16e\n", 
//     timeVec->data[i], hLMVec->data[i], hLMVec->data[i]);
// print_out("%.16e\t%.16e\t%.16e\n", 
//     timeVec->data[i], v, HamVec.data[i]);

        if (XLALSimIMRSpinEOBGetAmplitudeResidualPrec (&(rholmpwrlVec->data[i]), v, HamVec.data[i], modeL, modeM, params) == CEV_FAILURE)
        //RC: for the 21 and 55 mode rholmpwrlVec is always real. This is not true for the 33 mode. For this reason, in order to make this function general, we use a complex variable for it.
        {
            /* TODO: Clean-up */
            failed = 1; 
            goto QUIT;
        }
        // }
        // else
        // {
            // dr = (seobdynamics->posVecx[i]*seobdynamics->velVecx[i] + 
            //     seobdynamics->posVecy[i]*seobdynamics->velVecy[i] + 
            //     seobdynamics->posVecz[i]*seobdynamics->velVecz[i]) / seobdynamics->polarrVec[i];
            dr = seobdynamics->drVec[i];
            // ncrv = seobdynamics->polarrVec[i] * orbOmegaVec.data[i];
            ncrv = seobdynamics->rVec[i] * orbOmegaVec.data[i];
            drVec->data[i] = dr;
            ncrvVec->data[i] = ncrv;
            if ( XLALSimIMRSpinEOBGetSASpinFactorizedWaveformV2( &hLME, &polarDynamics, v, dr, ncrv, 0, HamVec.data[i], modeL, modeM, params ) == CEV_FAILURE )
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2 at step %d of the loop.", i);
                failed = 1; 
                goto QUIT;
            }
        // if (isnan(params->hCoeffs->h21T2ff10))
        // {
        //     print_debug("ca = (%g, %g, %g), cs = (%g, %g, %g)\n", seobdynamics->chiSxVec[i], seobdynamics->chiSyVec[i], seobdynamics->chiSzVec[i],
        //             seobdynamics->chiAxVec[i], seobdynamics->chiAyVec[i], seobdynamics->chiAzVec[i]);
        // }

            if (XLALSimIMRSpinEOBSAGetAmplitudeResidualPrecV2 (&rholmpwrlE, &polarDynamics, v, dr, ncrv, HamVec.data[i], modeL, modeM, params) == CEV_FAILURE)
            //RC: for the 21 and 55 mode rholmpwrlVec is always real. This is not true for the 33 mode. For this reason, in order to make this function general, we use a complex variable for it.
            {
                /* TODO: Clean-up */
                failed = 1; 
                goto QUIT;
            }
        // if (isnan(params->hCoeffs->h21T2ff10))
        // {
        //     print_debug("it become nan\n");
        // }
            // if (creal(rholmpwrlVec->data[i]) == 0.0)
            //     print_debug("[%d]get\n", i);
        // }
        rholmpwrlVecReal->data[i] = (REAL8)creal(rholmpwrlVec->data[i]);
        hLMdivrholmVec->data[i] = ((REAL8)cabs(hLMVec->data[i]))/fabs(rholmpwrlVecReal->data[i]);

        rholmpwrlEVecReal->data[i] = (REAL8) creal(rholmpwrlE);
        rholmpwrlEVecImag->data[i] = (REAL8) cimag(rholmpwrlE);
        hLMdivrholmEVec->data[i] = ((REAL8)cabs(hLME) / fabs(rholmpwrlEVecReal->data[i]));

        hLMrealVec->data[i] = (REAL8) creal(hLME/rholmpwrlE);
        hLMimagVec->data[i] = (REAL8) cimag(hLME/rholmpwrlE);
        // print_debug("[%d]");
        // if (rholmpwrlVecReal->data[i] == 0.0)
        //     print_debug("[%d]rholmPwrlReal = %.3e + i%.3e, v = %g, dr = %g, ncrv = %g \n", 
        //         i, creal(rholmpwrlVec->data[i]), cimag(rholmpwrlVec->data[i]), v, dr, ncrv);
        // if (timeVec->data[i] > 0.9*timewavePeak)
        //     print_debug("t[%d] = %g, rholmpwrlEVec = %.3e + i%.3e\n",
        //             i, timeVec->data[i], rholmpwrlEVecReal->data[i], rholmpwrlEVecImag->data[i]);
            // print_debug("t[%d] = %g, rholmpwrlEVec = %.3e + i%.3e, hLM = %.3e + i%.3e,r = %g, v = %g, dr = %g, ncrv = %g, ham = %g \n", 
            //     i, timeVec->data[i], rholmpwrlEVecReal->data[i], rholmpwrlEVecImag->data[i], creal(hLME), cimag(hLME), seobdynamics->polarrVec[i], v, dr, ncrv, HamVec.data[i]);
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

    status = gsl_spline_init (spline, timeVec->data, hLMdivrholmEVec->data, hLMdivrholmEVec->length);
    if (status != GSL_SUCCESS)
    {
        PRINT_LOG_INFO (LOG_CRITICAL, "gsl_spline_init failed.");
        failed = 1; 
        goto QUIT;
    }
    gsl_interp_accel_reset (acc);
    hLMdivrholmEAttachmentPoint = gsl_spline_eval (spline, timewavePeak, acc);


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
    rholmENRAttachmentPoint = nra/hLMdivrholmEAttachmentPoint;
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
    REAL8 rholmErealDot, rholmEimagDot, rholmErealDDot, rholmEimagDDot;
    status = gsl_spline_init (spline, timeVec->data, rholmpwrlEVecReal->data, rholmpwrlEVecReal->length);
    if (status != GSL_SUCCESS)
    {
        PRINT_LOG_INFO (LOG_CRITICAL, "gsl_spline_init failed.");
        failed = 1; 
        goto QUIT;
    }
    gsl_interp_accel_reset (acc);
    /* rho_lm^l + f_lm before the calibration parameter is set */
    rholmERealBeforeCalAP = gsl_spline_eval (spline, timewavePeak, acc);
    rholmErealDot = gsl_spline_eval_deriv (spline, timewavePeak, acc);
    rholmErealDDot = gsl_spline_eval_deriv2 (spline, timewavePeak, acc);

    status = gsl_spline_init (spline, timeVec->data, rholmpwrlEVecImag->data, rholmpwrlEVecReal->length);
    if (status != GSL_SUCCESS)
    {
        PRINT_LOG_INFO (LOG_CRITICAL, "gsl_spline_init failed.");
        failed = 1; 
        goto QUIT;
    }
    gsl_interp_accel_reset (acc);
    /* rho_lm^l + f_lm before the calibration parameter is set */
    rholmEImagBeforeCalAP = gsl_spline_eval (spline, timewavePeak, acc);
    rholmEimagDot = gsl_spline_eval_deriv (spline, timewavePeak, acc);
    rholmEimagDDot = gsl_spline_eval_deriv2 (spline, timewavePeak, acc);
    COMPLEX16 rholmEBeforeCalAP = rholmERealBeforeCalAP + I*rholmEImagBeforeCalAP;
    COMPLEX16 rholmEDot = rholmErealDot + I*rholmEimagDot;
    COMPLEX16 rholmEDDot = rholmErealDDot + I*rholmEimagDDot;

    status = gsl_spline_init (spline, timeVec->data, hLMrealVec->data, hLMdivrholmVec->length);
    if (status != GSL_SUCCESS)
    {
        PRINT_LOG_INFO (LOG_CRITICAL, "gsl_spline_init failed.");
        failed = 1; 
        goto QUIT;
    }
    gsl_interp_accel_reset (acc);
    hLMrealAP = gsl_spline_eval(spline, timewavePeak, acc);
    hLMrealAPDot = gsl_spline_eval_deriv(spline, timewavePeak, acc);
    hLMrealAPDDot = gsl_spline_eval_deriv2(spline, timewavePeak, acc);

    status = gsl_spline_init (spline, timeVec->data, hLMimagVec->data, hLMdivrholmVec->length);
    if (status != GSL_SUCCESS)
    {
        PRINT_LOG_INFO (LOG_CRITICAL, "gsl_spline_init failed.");
        failed = 1; 
        goto QUIT;
    }
    gsl_interp_accel_reset (acc);
    hLMimagAP = gsl_spline_eval(spline, timewavePeak, acc);
    hLMimagAPDot = gsl_spline_eval_deriv(spline, timewavePeak, acc);
    hLMimagAPDDot = gsl_spline_eval_deriv2(spline, timewavePeak, acc);
    // if(debugRC==1){
    //     FILE *timeomegapeak = NULL;
    //     timeomegapeak = fopen ("timeomegapeak.dat", "w");
    //     fprintf(timeomegapeak, "%.16f\n", timewavePeak);
    //     fclose(timeomegapeak);
    // }
    REAL8 x1AP, x1APDot, x1APDDot, x2AP, x2APDot, x2APDDot;
    status = gsl_spline_init (spline, timeVec->data, drVec->data, drVec->length);
    if (status != GSL_SUCCESS)
    {
        PRINT_LOG_INFO (LOG_CRITICAL, "gsl_spline_init failed.");
        failed = 1; 
        goto QUIT;
    }
    gsl_interp_accel_reset (acc);
    x1AP = gsl_spline_eval(spline, timewavePeak, acc);
    x1APDot = gsl_spline_eval_deriv(spline, timewavePeak, acc);
    x1APDDot = gsl_spline_eval_deriv2(spline, timewavePeak, acc);

    status = gsl_spline_init (spline, timeVec->data, ncrvVec->data, drVec->length);
    if (status != GSL_SUCCESS)
    {
        PRINT_LOG_INFO (LOG_CRITICAL, "gsl_spline_init failed.");
        failed = 1; 
        goto QUIT;
    }
    gsl_interp_accel_reset (acc);
    x2AP = gsl_spline_eval(spline, timewavePeak, acc);
    x2APDot = gsl_spline_eval_deriv(spline, timewavePeak, acc);
    x2APDDot = gsl_spline_eval_deriv2(spline, timewavePeak, acc);
    // print_debug("x1 = %f, x1Dot = %f, x1DDot = %f\n", x1AP, x1APDot, x1APDDot);
    // print_debug("x2 = %f, x2Dot = %f, x2DDot = %f\n", x2AP, x2APDot, x2APDDot);
    REAL8 PNCalCoeff1, PNCalCoeff2, PNCalCoeff3, PNCalCoeff4, PNCalCoeff5;
    COMPLEX16 initval;
    COMPLEX16 hLMPeak = hLMrealAP + I*hLMimagAP;
    COMPLEX16 hLMPeakDot = hLMrealAPDot + I*hLMimagAPDot;
    COMPLEX16 hLMPeakDDot = hLMrealAPDDot + I*hLMimagAPDDot;
    int retval;
    if ((modeL == 2) && (modeM == 1))
    {
        /* Here we compute ((rho_lm^l + f_lm + CalPar*omega^7/3)_NR - (rho_lm^l + f_lm)_EOB)/omega^7/3 to get CalPar.
        The factor rholmpwrlVecReal->data[0])/cabs(rholmpwrlVecReal->data[0]) is used to know the sign of the function (rho_lm^l + f_lm + CalPar*omega^7/3)_NR which is computed as absolute value */
        params->cal21 = (rholmNRAttachmentPoint - rholmBeforeCalibAttachmentPoint)/(pow(omegaAttachmentPoint,7./3.));
        initval = (rholmENRAttachmentPoint - rholmEBeforeCalAP)*(pow(rAttachmentPoint,12./2.));
        // initval = (rholmENRAttachmentPoint - rholmEBeforeCalAP) /(pow(omegaAttachmentPoint,7./3.));
        // retval = FindCalPNCoeffs(hLMPeak, hLMPeakDot, hLMPeakDDot, 
        //             rholmEBeforeCalAP, rholmEDot, rholmEDDot,
        //             nra, nraDot, nraDDot, nrOmega, nrOmegaDot, 
        //             x0AP, x0DotAP, x0DDotAP, 
        //             x1AP, x1APDot, x1APDDot, 
        //             x2AP, x2APDot, x2APDDot, 
        //             5, creal(initval), &PNCalCoeff1, &PNCalCoeff2, &PNCalCoeff3, &PNCalCoeff4, &PNCalCoeff5);
        params->cal21E = creal(initval);
        // if (retval != CEV_SUCCESS)
        //     params->cal21E = creal(initval);
        // else
        // {
            // params->cal21E1 = PNCalCoeff1;// * rholmEBeforeCalAP / (hLMrealAP + I*hLMimagAP);
            // params->cal21E2 = PNCalCoeff2;// * rholmEBeforeCalAP / (hLMrealAP + I*hLMimagAP);
            // params->cal21E3 = PNCalCoeff3;
            // params->cal21E = PNCalCoeff4;
            // params->cal21E4 = PNCalCoeff5;
            // print_debug("fitpms = %f, %f, %f, %f, %f, %f\n", PNCalCoeff1, PNCalCoeff2, PNCalCoeff3, PNCalCoeff4, PNCalCoeff5, params->hCoeffs->f21v6);
        // }
        // Test
        // COMPLEX16 fithLMPeak, fithLMDotPeak;
        // REAL8 fitALMDotPeak;
        // REAL8 Re, Im, ReDot, ImDot;
        // fithLMPeak = hLMPeak * (rholmEBeforeCalAP + pow(x0AP, 7) * params->cal21E);
        // fithLMPeak = hLMPeak + pow(x0AP, 7)*PNCalCoeff1 + pow(x0AP, 8) * PNCalCoeff2;
        // fithLMDotPeak = hLMrealAPDot + I*hLMimagAPDot + 7*pow(x0AP, 6)*x0DotAP*PNCalCoeff1 + 8*pow(x0AP, 7)*x0DotAP * PNCalCoeff2;
        // fitALMDotPeak = (cimag(fithLMPeak) * cimag(fithLMDotPeak) + creal(fithLMPeak) * creal(fithLMDotPeak)) / cabs(fithLMPeak);

        // print_debug("nra = %f, eoba = %f, nraDot = %f, eobaDot = %f\n", nra, cabs(fithLMPeak), nraDot, fitALMDotPeak);
                            // print_debug("cal21 = %.16f, nra = %.3f\n", params->cal21, nra);
    }
    if ((modeL == 5) && (modeM == 5))
    {
        /* Here we compute ((rho_lm^l + f_lm + CalPar*omega^7/3)_NR - (rho_lm^l + f_lm)_EOB)/omega^7/3 to get CalPar.
        The factor rholmpwrlVecReal->data[0])/cabs(rholmpwrlVecReal->data[0]) is used to know the sign of the function (rho_lm^l + f_lm + CalPar*omega^7/3)_NR which is computed as absolute value */
        params->cal55 = (rholmNRAttachmentPoint - rholmBeforeCalibAttachmentPoint)/(pow(omegaAttachmentPoint,5./3.));
                            //printf("params->cal55 = %.16f\n",params->cal55);
    }
    //if (isnan(params->cal21) || isnan(params->cal21E1) || isnan(params->cal21E2) || isnan(params->cal21E) || isnan(params->cal55))
    if (isnan(creal(params->cal21)) || isnan(cimag(params->cal21)) || 
        isnan(params->cal21E1) || isnan(params->cal21E2) || isnan(params->cal21E) || 
        isnan(creal(params->cal55)) || isnan(cimag(params->cal55)))
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "(l,m) = (%d, %d) calibration get nan", modeL, modeM);
        PRINT_LOG_INFO(LOG_DEBUG, "rholmENRAttachmentPoint = %g, rholmEBeforeCalAP = %g + i %g\n", rholmENRAttachmentPoint, creal(rholmEBeforeCalAP), cimag(rholmEBeforeCalAP));
        PRINT_LOG_INFO(LOG_DEBUG, "cal21 = %g + i %g, cal21E = %g, rAttachmentPoint = %g\n", creal(params->cal21), cimag(params->cal21), params->cal21E, rAttachmentPoint);
        PRINT_LOG_INFO(LOG_DEBUG, "rholmNRAttachmentPoint = %g, rholmBeforeCalibAttachmentPoint = %g, omegaAttachmentPoint = %g",
            rholmNRAttachmentPoint, rholmBeforeCalibAttachmentPoint, omegaAttachmentPoint);
        PRINT_LOG_INFO(LOG_DEBUG, "nra = %g, hLMdivrholmAttachmentPoint = %g", nra, hLMdivrholmAttachmentPoint);
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

/**
 * This function generates a waveform mode for a given SEOB dynamics.
 */
// NOTE: as is written here, the step
// XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients in the loop will be repeated
// across modes -- would be more efficient to loop on modes inside the loop on
// times
static int SEOBCalculatehlmAmpPhase(
    CAmpPhaseSequence *
        *hlm, /**<< Output: hlm in complex amplitude / phase form */
    INT4 l,   /**<< Input: mode index l */
    INT4 m,   /**<< Input: mode index m */
    SEOBdynamics *seobdynamics, /**<< Input: SEOB dynamics */
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
    if (includeNQC == 0) 
    {
        if (((l == 2) && (m == 1)) || ((l == 5) && (m == 5))) 
        {
            if (XLALSimIMREOBCalcCalibCoefficientHigherModesPrec(
                    seobParams, l, m, seobdynamics,
                    tPeakOmega - seobdynamics->tVec[0], m1, m2,
                    deltaT) == CEV_FAILURE) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBCalcCalibCoefficientHigherModesPrec.");
                return CEV_FAILURE;
            }
        }
    }

    /* Loop to compute compute amplitude and phase of the hlm mode */
    REAL8 s1dotZ, s2dotZ, chiS, chiA, tplspin, t, omega, ham, v;
    UINT i, j;
    // CHAR fdebug[256];
    // sprintf(fdebug, "fdebug_new_%d%d.dat", l, m);
    // FILE *out = NULL;
    // if (!includeNQC)
    // {
    //     print_debug("dump to %s\n", fdebug);
    //     // print_debug("code version = %d\n", CODE_VERSION);
    //     out = fopen(fdebug, "w");
    // }
    for (i = 0; i < retLen; i++) 
    {
        /* Compute waveform coefficients */
        t = seobdynamics->tVec[i];
        omega = seobdynamics->omegaVec[i];
        ham = seobdynamics->hamVec[i];
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

        if (CODE_VERSION == 3)
        {
            if (EccPrec_CalcSpinPrecFacWaveformCoefficients(
                    seobParams->hCoeffs, m1, m2, eta, tplspin, chiS, chiA,
                    seobdynamics->chiSxVec[i], seobdynamics->chiSyVec[i], seobdynamics->chiSzVec[i],
                    seobdynamics->chiAxVec[i], seobdynamics->chiAyVec[i], seobdynamics->chiAzVec[i],
                    451) == CEV_FAILURE) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in EccPrec_CalcSpinPrecFacWaveformCoefficients at step %d of the loop.", i);
                return CEV_FAILURE;
            }
        }
        else
        {
            if (XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(
                    seobParams->hCoeffs, m1, m2, eta, tplspin, chiS, chiA,
                    SpinAlignedEOBversionWaveform) == CEV_FAILURE) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients at step %d of the loop.", i);
                return CEV_FAILURE;
            }
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
        polarDynamics.data[2] = seobdynamics->polarprVec[i];
        polarDynamics.data[3] = seobdynamics->polarpphiVec[i];
        v = cbrt(omega);
        COMPLEX16 hlm_val = 0.;
        if (CODE_VERSION == 0)
        {
            if (XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform(
                    &hlm_val, &polarDynamics, &values, v, ham, l, m, seobParams) != CEV_SUCCESS) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform at step %d of the loop.", i);
                return CEV_FAILURE;
            }
        }
        else
        {
            dr = (seobdynamics->posVecx[i]*seobdynamics->velVecx[i] + 
                seobdynamics->posVecy[i]*seobdynamics->velVecy[i] + 
                seobdynamics->posVecz[i]*seobdynamics->velVecz[i]) / seobdynamics->polarrVec[i];
            ncrv = seobdynamics->polarrVec[i] * omega;

            if (XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2(
                    &hlm_val, &polarDynamics, &values, v, dr, ncrv, seobdynamics->polarprDotVec[i], ham, l, m, seobParams) != CEV_SUCCESS) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2 at step %d of the loop.", i);
                return CEV_FAILURE;
            }
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
        // if (!includeNQC)
        //     fprintf(out, "%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n",
        //         t, creal(hlm_val), cimag(hlm_val), v, omega, ham);
    }
//   if (!includeNQC)
//     fclose(out);
    /* Unwrap the phase vector, in place */
    XLALREAL8VectorUnwrapAngle((*hlm)->phase, (*hlm)->phase);

    return CEV_SUCCESS;
}

static int SEOBSACalculatehlmAmpPhase(
    CAmpPhaseSequence *
        *hlm, /**<< Output: hlm in complex amplitude / phase form */
    INT4 l,   /**<< Input: mode index l */
    INT4 m,   /**<< Input: mode index m */
    SEOBSAdynamics *seobdynamics, /**<< Input: SEOB dynamics */
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
    // UINT SpinAlignedEOBversion = seobParams->SpinAlignedEOBversion;
    UINT SpinAlignedEOBversion = 4;
    UINT SpinAlignedEOBversionWaveform; // RC: I use this different variable
                                        // because the PN terms in the waveform
                                        // are different from those in the flux

    /* Length of dynamics data and sampling step */
    UINT retLen = seobdynamics->length;
    // PRINT_LOG_INFO(LOG_DEBUG, "l,m = (%d, %d), deltaT = %.16e, retLen = %d\n", l, m, deltaT, retLen);

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
    if (includeNQC == 0) 
    {
        if (((l == 2) && (m == 1)) || ((l == 5) && (m == 5))) 
        {
            if (XLALSimIMREOBSACalcCalibCoefficientHigherModes(
                    seobParams, l, m, seobdynamics,
                    tPeakOmega - seobdynamics->tVec[0], m1, m2,
                    deltaT) == CEV_FAILURE) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBSACalcCalibCoefficientHigherModes.");
                return CEV_FAILURE;
            }
        }
    }

    /* Loop to compute compute amplitude and phase of the hlm mode */
    REAL8 s1dotZ, s2dotZ, chiS, chiA, tplspin, t, omega, ham, v;
    // s1dotZ = seobdynamics->s1dotZVec[i];
    // s2dotZ = seobdynamics->s2dotZVec[i];
    REAL8 chi1dotZ = seobParams->chi1;
    REAL8 chi2dotZ = seobParams->chi2;
    chiS = SEOBCalculateChiS(chi1dotZ, chi2dotZ);
    chiA = SEOBCalculateChiA(chi1dotZ, chi2dotZ);

    tplspin = SEOBCalculatetplspin(m1, m2, eta, s1dotZ, s2dotZ);
    SpinAlignedEOBversionWaveform = 451;

    if (XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(
            seobParams->hCoeffs, m1, m2, eta, tplspin, chiS, chiA,
            SpinAlignedEOBversionWaveform) == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients.");
        return CEV_FAILURE;
    }
    seobParams->hCoeffs->f21v7c = seobParams->cal21;
    seobParams->hCoeffs->f21v7cEff00 = seobParams->cal21E;
    seobParams->hCoeffs->f21v7cEff10 = seobParams->cal21E1;
    seobParams->hCoeffs->f21v7cEff11 = seobParams->cal21E2;
    seobParams->hCoeffs->f21v7cEff01 = seobParams->cal21E3;
    seobParams->hCoeffs->f21v7cEff02 = seobParams->cal21E4;
    seobParams->hCoeffs->f55v5c = seobParams->cal55;

    UINT i, j;
    for (i = 0; i < retLen; i++) 
    {
        /* Compute waveform coefficients */
        t = seobdynamics->tVec[i];
        omega = seobdynamics->dphiVec[i];
        ham = seobdynamics->HVec[i];
        // print_debug("f21v7c = %.16f\n", seobParams->hCoeffs->f21v7c);
        /* Dynamics, polar dynamics, omega */
        // for (j = 0; j < 14; j++) 
        //     values.data[j] = seobdynamics->array->data[i + (j + 1) * retLen];
        polarDynamics.data[0] = seobdynamics->rVec[i];
        polarDynamics.data[1] = seobdynamics->phiVec[i];
        polarDynamics.data[2] = seobdynamics->prTVec[i];
        polarDynamics.data[3] = seobdynamics->pphiVec[i];
        v = cbrt(omega);
        COMPLEX16 hlm_val = 0.;
        if (CODE_VERSION == 0)
        {
            if (XLALSimIMRSpinEOBGetSASpinFactorizedWaveform(
                    &hlm_val, &polarDynamics, v, ham, l, m, seobParams) != CEV_SUCCESS) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform at step %d of the loop.", i);
                return CEV_FAILURE;
            }
        }
        else
        {
            // dr = (seobdynamics->posVecx[i]*seobdynamics->velVecx[i] + 
            //     seobdynamics->posVecy[i]*seobdynamics->velVecy[i] + 
            //     seobdynamics->posVecz[i]*seobdynamics->velVecz[i]) / seobdynamics->polarrVec[i];
            dr = seobdynamics->drVec[i];
            ncrv = seobdynamics->rVec[i] * omega;

            if (XLALSimIMRSpinEOBGetSASpinFactorizedWaveformV2(
                    &hlm_val, &polarDynamics, v, dr, ncrv, 0, ham, l, m, seobParams) != CEV_SUCCESS) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2 at step %d of the loop.", i);
                return CEV_FAILURE;
            }
        }
        /* NQC correction */
        COMPLEX16 factor_nqc = 1.;
        if (includeNQC) 
        {
            if (XLALSimIMRSpinEOBSANonQCCorrection(&factor_nqc, &polarDynamics, omega, t, 
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
        // if (i==0 || i==1)
        // {
        //     print_debug("[%d,%d]%.16e\t%.16e\t%.16e\n", l, m, t, creal(hlmNQC), cimag(hlmNQC));
        //     // print_debug("t = %.16e, omega = %.16e, ham = %.16e, r = %.16e, phi = %.16e, prT = %.16e, pphi = %.16e\n",
        //     //     t, omega, ham, seobdynamics->rVec[i], seobdynamics->phiVec[i], seobdynamics->prTVec[i], seobdynamics->pphiVec[i]);
        // }
        // if (IS_DEBUG)
        //     print_out("%.16e\t%.16e\t%.16e\n", t, creal(hlmNQC), cimag(hlmNQC));
            // print_out("%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n", t, 
            //     seobdynamics->rVec[i], seobdynamics->phiVec[i],
            //     seobdynamics->prTVec[i], seobdynamics->pphiVec[i],
            //     seobdynamics->drVec[i], seobdynamics->dphiVec[i],
            //     seobdynamics->dprTVec[i], seobdynamics->dpphiVec[i],
            //     seobdynamics->HVec[i]);
    }

    /* Unwrap the phase vector, in place */
    XLALREAL8VectorUnwrapAngle((*hlm)->phase, (*hlm)->phase);

    return CEV_SUCCESS;
}

static int SEOBCalculatehlmAmpPhase_noNQC(
    CAmpPhaseSequence *
        *hlm, /**<< Output: hlm in complex amplitude / phase form */
    INT4 l,   /**<< Input: mode index l */
    INT4 m,   /**<< Input: mode index m */
    SEOBdynamics *seobdynamics, /**<< Input: SEOB dynamics */
    SpinEOBParams *seobParams  /**<< SEOB params */
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


    /* Loop to compute compute amplitude and phase of the hlm mode */
    REAL8 s1dotZ, s2dotZ, chiS, chiA, tplspin, t, omega, ham, v;
    UINT i, j;
    for (i = 0; i < retLen; i++) 
    {
        /* Compute waveform coefficients */
        t = seobdynamics->tVec[i];
        omega = seobdynamics->omegaVec[i];
        ham = seobdynamics->hamVec[i];
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

        if (CODE_VERSION == 3)
        {
            if (EccPrec_CalcSpinPrecFacWaveformCoefficients(
                    seobParams->hCoeffs, m1, m2, eta, tplspin, chiS, chiA,
                    seobdynamics->chiSxVec[i], seobdynamics->chiSyVec[i], seobdynamics->chiSzVec[i],
                    seobdynamics->chiAxVec[i], seobdynamics->chiAyVec[i], seobdynamics->chiAzVec[i],
                    451) == CEV_FAILURE) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in EccPrec_CalcSpinPrecFacWaveformCoefficients at step %d of the loop.", i);
                return CEV_FAILURE;
            }
        }
        else
        {
            if (XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(
                    seobParams->hCoeffs, m1, m2, eta, tplspin, chiS, chiA,
                    SpinAlignedEOBversionWaveform) == CEV_FAILURE) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients at step %d of the loop.", i);
                return CEV_FAILURE;
            }
        }
        seobParams->hCoeffs->f21v7c = 0.0;
        seobParams->hCoeffs->f21v7cEff00 = 0.0;
        seobParams->hCoeffs->f21v7cEff10 = 0.0;
        seobParams->hCoeffs->f21v7cEff11 = 0.0;
        seobParams->hCoeffs->f21v7cEff01 = 0.0;
        seobParams->hCoeffs->f21v7cEff02 = 0.0;
        seobParams->hCoeffs->f55v5c = 0.0;
        // print_debug("f21v7c = %.16f\n", seobParams->hCoeffs->f21v7c);
        /* Dynamics, polar dynamics, omega */
        for (j = 0; j < 14; j++) 
            values.data[j] = seobdynamics->array->data[i + (j + 1) * retLen];
        polarDynamics.data[0] = seobdynamics->polarrVec[i];
        polarDynamics.data[1] = seobdynamics->polarphiVec[i];
        polarDynamics.data[2] = seobdynamics->polarprVec[i];
        polarDynamics.data[3] = seobdynamics->polarpphiVec[i];
        v = cbrt(omega);
        COMPLEX16 hlm_val = 0.;
        if (CODE_VERSION == 0)
        {
            if (XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform(
                    &hlm_val, &polarDynamics, &values, v, ham, l, m, seobParams) != CEV_SUCCESS) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform at step %d of the loop.", i);
                return CEV_FAILURE;
            }
        }
        else
        {
            dr = (seobdynamics->posVecx[i]*seobdynamics->velVecx[i] + 
                seobdynamics->posVecy[i]*seobdynamics->velVecy[i] + 
                seobdynamics->posVecz[i]*seobdynamics->velVecz[i]) / seobdynamics->polarrVec[i];
            ncrv = seobdynamics->polarrVec[i] * omega;

            if (XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2(
                    &hlm_val, &polarDynamics, &values, v, dr, ncrv, seobdynamics->polarprDotVec[i], ham, l, m, seobParams) != CEV_SUCCESS) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2 at step %d of the loop.", i);
                return CEV_FAILURE;
            }
        }
        /* Result and output */
        COMPLEX16 hlmNQC = hlm_val;
        (*hlm)->xdata->data[i] = t; /* Copy times */
        (*hlm)->camp_real->data[i] = creal(hlmNQC);
        (*hlm)->camp_imag->data[i] = cimag(hlmNQC); /* We use only real amplitudes */
        (*hlm)->phase->data[i] = carg(hlmNQC);
    }

    /* Unwrap the phase vector, in place */
    XLALREAL8VectorUnwrapAngle((*hlm)->phase, (*hlm)->phase);

    return CEV_SUCCESS;
}

/**
 * This function converts a spin-aligned dynamics as output by the Runge-Kutta
 * integrator to a generic-spin dynamics. Spin-aligned dynamics format: t, r,
 * phi, pr, pphi Generic-spin dynamics format: t, x, y, z, px, py, pz, s1x, s1y,
 * s1z, s2x, s2y, s2z, phiMod, phiDMod
 */
static int SEOBConvertSpinAlignedDynamicsToGenericSpins(
    REAL8Array **dynamics, /**<< Output: pointer to array for the generic-spin
                              dynamics */
    REAL8Array *dynamics_spinaligned, /**<< Input: array for the aligned-spin
                                         dynamics */
    UINT retLen,                     /**<< Input: length of dynamics */
    REAL8 chi1, /**<< Input: spin 1 aligned component (dimensionless) */
    REAL8 chi2, /**<< Input: spin 2 aligned component (dimensionless) */
    SpinEOBParams *seobParams /**<< SEOB params */
) 
{
    UINT i;

    /* Masses */
    REAL8 m1 = seobParams->m1;
    REAL8 m2 = seobParams->m2;
    REAL8 mTotal = m1 + m2;

    /* Create output dynamics */
    *dynamics = CreateREAL8Array(2, 15, retLen);

    /* Convert the spin-aligned dynamics to a generic-spins dynamics */
    REAL8Vector tVec, rVec, phiVec, prVec, pPhiVec;
    tVec.length = rVec.length = phiVec.length = prVec.length = pPhiVec.length =
        retLen;
    tVec.data = dynamics_spinaligned->data;
    rVec.data = dynamics_spinaligned->data + retLen;
    phiVec.data = dynamics_spinaligned->data + 2 * retLen;
    prVec.data = dynamics_spinaligned->data + 3 * retLen;
    pPhiVec.data = dynamics_spinaligned->data + 4 * retLen;
    for (i = 0; i < retLen; i++) 
    {
        (*dynamics)->data[i] = tVec.data[i];
        (*dynamics)->data[retLen + i] = rVec.data[i] * cos(phiVec.data[i]);
        (*dynamics)->data[2 * retLen + i] = rVec.data[i] * sin(phiVec.data[i]);
        (*dynamics)->data[3 * retLen + i] = 0.;
        (*dynamics)->data[4 * retLen + i] =
            prVec.data[i] * cos(phiVec.data[i]) -
            pPhiVec.data[i] / rVec.data[i] * sin(phiVec.data[i]);
        (*dynamics)->data[5 * retLen + i] =
            prVec.data[i] * sin(phiVec.data[i]) +
            pPhiVec.data[i] / rVec.data[i] * cos(phiVec.data[i]);
        (*dynamics)->data[6 * retLen + i] = 0.;
        (*dynamics)->data[7 * retLen + i] = 0.;
        (*dynamics)->data[8 * retLen + i] = 0.;
        (*dynamics)->data[9 * retLen + i] = chi1 * (m1 * m1 / mTotal / mTotal);
        (*dynamics)->data[10 * retLen + i] = 0.;
        (*dynamics)->data[11 * retLen + i] = 0.;
        (*dynamics)->data[12 * retLen + i] = chi2 * (m2 * m2 / mTotal / mTotal);
        (*dynamics)->data[13 * retLen + i] = phiVec.data[i];
        (*dynamics)->data[14 * retLen + i] = 0.;
    }

    return CEV_SUCCESS;
}


/**
 * Stopping condition for the regular resolution SEOBNRv1/2 orbital evolution
 * -- stop when reaching max orbital frequency in strong field.
 * At each test,
 * if omega starts to decrease, return 1 to stop evolution;
 * if not, update omega with current value and return GSL_SUCCESS to continue
 * evolution.
 */
static int XLALEOBSpinPrecAlignedStopCondition(
    double t,       /**< UNUSED */
    const double values[], /**< dynamical variable values */
    double dvalues[],      /**< dynamical variable time derivative values */
    void *funcParams       /**< physical parameters */
) 
{
    REAL8 omega, r;
    SpinEOBParams *params = (SpinEOBParams *)funcParams;

    r = values[0];
    omega = dvalues[1];
    if (r < 6. && omega < params->omega) 
    {
        return 1;
    }
    params->omega = omega;
    return GSL_SUCCESS;
}

REAL8 g_h_rISCO = 6.;
REAL8 get_h_rISCO()
{
    return g_h_rISCO;
}

void set_h_rISCO(REAL8 rISCO)
{
    if (rISCO < 3.)
        PRINT_LOG_INFO(LOG_WARNING, "r_ISCO = %fM is too small!", rISCO);
    g_h_rISCO = rISCO;
    return;
}

static int EOBHighSRStopConditionEcc(double t,
                                  const double values[],
                                  double dvalues[],
                                  void *funcParams)
{
    REAL8 r = values[0];
    SpinEOBParams *params = (SpinEOBParams *) funcParams;
    if(r > g_h_rISCO)
    {
        params->omegaPeaked = 0;
        return GSL_SUCCESS;
    }
    REAL8 eta, ham, omega;
    eta = params->eta;
    REAL8Vector xVec, pVec;
    REAL8 xVal[3] = {0,0,0}, pVal[3] = {0,0,0};
    xVal[0] = values[0];
    pVal[0] = values[2];
    pVal[1] = values[3] / values[0];
    xVec.length = pVec.length = 3;
    xVec.data = xVal;
    pVec.data = pVal;
    omega = dvalues[1];
    ham = EOBHamiltonian(eta, &xVec, &pVec,
                    params->s1Vec, params->s2Vec,
                    params->sigmaKerr, params->sigmaStar,
                    params->tortoise, params->seobCoeffs);

    UINT counter = params->omegaPeaked;
    if (omega < params->omega)
    {
        params->omegaPeaked = counter + 1;
    }

    if( isnan(ham) || isnan (dvalues[3]) || isnan (dvalues[2])
        || isnan (dvalues[1]) || isnan (dvalues[0]) || 
        (dvalues[2] >= 0.&& dvalues[0] >=0.) || params->omegaPeaked >= 5 )
    {
        if (dvalues[2] >= 0.&& dvalues[0] >=0. && r>g_h_rISCO-1.)
        {
            params->omegaPeaked = 0;
            return GSL_SUCCESS;
        }
        // if (isnan(ham) || isnan (dvalues[3]) || isnan (dvalues[2])
        // || isnan (dvalues[1]) || isnan (dvalues[0]))
        //     print_debug("The integration stopped because of nan (%f/)\n", r);
        // else if (dvalues[2] >= 0.&& dvalues[0] >=0. )
        //     print_debug("The integration stopped because of dr&ddr > 0 (%f)\n", r);
        // else
        //     print_debug("The integration stopped because of omegaPeaked > 5 (%f)\n", r);
        // print_log("ham = %e, dvalues = (%e, %e, %e, %e)\n", ham, dvalues[0], dvalues[1], dvalues[2], dvalues[3]);
        return 1;
    }
    params->omega = omega;
    return GSL_SUCCESS;
}

/**
 * Stopping condition for the high resolution SEOBNRv4.
 */
static int XLALSpinPrecAlignedHiSRStopCondition(
    double t,              /**< UNUSED */
    const double values[], /**< dynamical variable values */
    double dvalues[],       /**< dynamical variable time derivative values */
    void *funcParams /**< physical parameters */
) 
{
    REAL8 omega, r;
    UINT counter;
    SpinEOBParams *params = (SpinEOBParams *)funcParams;
    r = values[0];
    omega = dvalues[1];
    counter = params->omegaPeaked;

    if (r < 6. && omega < params->omega) 
    {
        params->omegaPeaked = counter + 1;
    }
    if (dvalues[2] >= 0. || params->omegaPeaked == 5 ||
        isnan(dvalues[3]) || isnan(dvalues[2]) || isnan(dvalues[1]) ||
        isnan(dvalues[0])) 
    {
        return 1;
    }
    params->omega = omega;
    return GSL_SUCCESS;
}

/**
 * Stopping conditions for dynamics integration for SEOBNRv4P
 */
static int
XLALEOBSpinPrecStopConditionBasedOnPR(double t, 
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
    //     XLAL_PRINT_INFO("XLALEOBSpinPrecStopConditionBasedOnPR:: r = %e %e\n",
    //                     sqrt(r2), omega);
    // }
    // if (debugPK) {
    //     XLAL_PRINT_INFO(
    //         "XLALEOBSpinPrecStopConditionBasedOnPR:: values = %e %e %e %e %e %e\n",
    //         values[6], values[7], values[8], values[9], values[10], values[11]);
    // }
    // if (debugPK) {
    //     XLAL_PRINT_INFO(
    //         "XLALEOBSpinPrecStopConditionBasedOnPR:: dvalues = %e %e %e %e %e %e\n",
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

SpinEOBParams *CreateSpinEOBParams(REAL8 m1, REAL8 m2, 
                                   REAL8 s1x, REAL8 s1y, 
                                   REAL8 s1z, REAL8 s2x, 
                                   REAL8 s2y, REAL8 s2z,
                                   REAL8 e0,
                                   HyperParams *params)
{
    SpinEOBParams *seobParams = NULL;
    SpinEOBHCoeffs *seobCoeffs = NULL;
    SpinEOBHSACoeffs *saCoeffs = NULL;
    SEOBHCoeffConstants *seobCoeffConsts = NULL;
    HyperParams *hParams = NULL;
    NewtonMultipolePrefixes *prefixes = NULL;
    FacWaveformCoeffs *hCoeffs = NULL;
    EOBNonQCCoeffs *nqcCoeffs = NULL;
    EccCorrectionCoeffs *eccCoeffs = NULL;

    REAL8Vector *s1Vec = NULL;
    REAL8Vector *s2Vec = NULL;
    REAL8Vector *J0Vec = NULL;
    REAL8Vector *sigmaStar = NULL;
    REAL8Vector *sigmaKerr = NULL;
    REAL8Vector *chi1Vec = NULL;
    REAL8Vector *chi2Vec = NULL;

    seobParams = (SpinEOBParams *) MYCalloc(1, sizeof(SpinEOBParams));
    if (!seobParams) goto QUIT;
    saCoeffs = (SpinEOBHSACoeffs *) MYCalloc(1, sizeof(SpinEOBHSACoeffs));
    if (!saCoeffs) goto QUIT;
    seobCoeffs = (SpinEOBHCoeffs *) MYCalloc(1, sizeof(SpinEOBHCoeffs));
    if (!seobCoeffs) goto QUIT;
    seobCoeffConsts = (SEOBHCoeffConstants *) MYCalloc(1, sizeof(SEOBHCoeffConstants));
    if (!seobCoeffConsts) goto QUIT;
    hParams = (HyperParams *) MYCalloc(1, sizeof(HyperParams));
    if (!hParams) goto QUIT;
    prefixes = (NewtonMultipolePrefixes *) MYCalloc(1, sizeof(NewtonMultipolePrefixes));
    if (!prefixes) goto QUIT;
    hCoeffs = (FacWaveformCoeffs *) MYCalloc(1, sizeof(FacWaveformCoeffs));
    if (!hCoeffs) goto QUIT;
    nqcCoeffs = (EOBNonQCCoeffs *) MYCalloc(1, sizeof(EOBNonQCCoeffs));
    if (!nqcCoeffs) goto QUIT;
    eccCoeffs = (EccCorrectionCoeffs *)MYCalloc(1, sizeof(EccCorrectionCoeffs));
    if (!eccCoeffs) goto QUIT;

    if (params)
        memcpy(hParams, params, sizeof(HyperParams));
    s1Vec = CreateREAL8Vector(3);
    if (!s1Vec) goto QUIT;
    s2Vec = CreateREAL8Vector(3);
    if (!s2Vec) goto QUIT;
    J0Vec = CreateREAL8Vector(3);
    if (!J0Vec) goto QUIT;
    sigmaStar = CreateREAL8Vector(3);
    if (!sigmaStar) goto QUIT;
    sigmaKerr = CreateREAL8Vector(3);
    if (!sigmaKerr) goto QUIT;
    // Here m1 + m2 = 1
    REAL8 eta, q, mtotal;
    mtotal = m1 + m2;
    q = m1/m2;
    eta = q / (q+1.) / (q+1.);
    if (eta > 0.25)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Somehow the mass ratio eta = %.16e > 0.25 (q = %.16e), which implies that double precision floating-point numbers are not enough.", eta, q);
        eta = 0.25;
    }
    s1Vec->data[0] = s1x * m1 * m1 / mtotal / mtotal;
    s1Vec->data[1] = s1y * m1 * m1 / mtotal / mtotal;
    s1Vec->data[2] = s1z * m1 * m1 / mtotal / mtotal;
    s2Vec->data[0] = s2x * m2 * m2 / mtotal / mtotal;
    s2Vec->data[1] = s2y * m2 * m2 / mtotal / mtotal;
    s2Vec->data[2] = s2z * m2 * m2 / mtotal / mtotal;
    J0Vec->data[0] = s1Vec->data[0] + s2Vec->data[0];
    J0Vec->data[1] = s1Vec->data[1] + s2Vec->data[1];
    J0Vec->data[2] = s1Vec->data[2] + s2Vec->data[2];

    EOBCalculateSigmaStar(sigmaStar, m1, m2, s1Vec, s2Vec);
    EOBCalculateSigmaKerr(sigmaKerr, s1Vec, s2Vec);
    seobParams->eta = eta;
    seobParams->m1 = m1;
    seobParams->m2 = m2;
    seobParams->s1Vec = s1Vec;
    seobParams->s2Vec = s2Vec;
    seobParams->J0Vec = J0Vec;
    seobParams->sigmaKerr = sigmaKerr;
    seobParams->sigmaStar = sigmaStar;
    seobParams->chi1 = s1z;
    seobParams->chi2 = s2z;
    if (XLALEOBSpinPrecCalcSEOBHCoeffConstants(eta, seobCoeffConsts) != CEV_SUCCESS)
        goto QUIT;
    seobParams->seobCoeffConsts = seobCoeffConsts;
    REAL8 a;
    a = sqrt(inner_product3d(sigmaKerr->data, sigmaKerr->data));
    seobParams->a = a;
    REAL8 Lhat[3] = {0.0, 0.0, 1.0};
    REAL8 tempS1_p = inner_product3d(s1Vec->data, Lhat);
    REAL8 tempS2_p = inner_product3d(s2Vec->data, Lhat);
    REAL8 S1_perp[3] = {0, 0, 0};
    REAL8 S2_perp[3] = {0, 0, 0};
    INT jj;
    for (jj = 0; jj < 3; jj++) 
    {
        S1_perp[jj] = s1Vec->data[jj] - tempS1_p * Lhat[jj];
        S2_perp[jj] = s2Vec->data[jj] - tempS2_p * Lhat[jj];
    }
    REAL8 S_con = 0.0;
    if (a > 1e-6)
    {
        S_con = inner_product3d(sigmaKerr->data, Lhat);
        S_con /= (1 - 2 * eta);
        S_con += (inner_product3d(S1_perp, sigmaKerr->data) +
                inner_product3d(S2_perp, sigmaKerr->data)) /
                a / (1 - 2 * eta) / 2.;
    }
    if (XLALSimIMRCalculateSpinPrecEOBHCoeffs_v2(seobCoeffs, eta, a, S_con, 4, hParams) != CEV_SUCCESS)
        goto QUIT;
    seobParams->seobCoeffs = seobCoeffs;
    CalculateSpinEOBHSACoeffs(m1, m2, s1z, s2z, saCoeffs);
    seobParams->saCoeffs = saCoeffs; // only used in SA code... 2022.9.15
    if (XLALSimIMREOBComputeNewtonMultipolePrefixes(prefixes, m1, m2) != CEV_SUCCESS)
        goto QUIT;
    CalculateEccCorrectionCoeffs(eta, s1z, s2z, eccCoeffs); // only used in SA code... 2022.9.14
    seobParams->prefixes = prefixes;
    seobParams->hCoeffs = hCoeffs;
    seobParams->hParams = hParams;
    seobParams->nqcCoeffs = nqcCoeffs;
    seobParams->eccCoeffs = eccCoeffs;
    seobParams->tortoise = 1;
    seobParams->use_hm = FALSE;
    REAL8 EPS_ALIGN = 1e-4;
    if (sqrt(s1x*s1x + s1y*s1y) < EPS_ALIGN && sqrt(s2x*s2x + s2y*s2y) < EPS_ALIGN)
    {
        PRINT_LOG_INFO(LOG_INFO, "Almost spinAligned.");
        seobParams->alignedSpins = TRUE;
    }
    else
        seobParams->alignedSpins = FALSE;
    seobParams->ignoreflux = FALSE;
    seobParams->e0 = e0;
    seobParams->p0 = params->sl_p;
    seobParams->x0 = params->x0;
    seobParams->cal21E = 0.;
    seobParams->cal21E1 = 0.;
    seobParams->cal21E2 = 0.;
    seobParams->cal21E3 = 0.;
    seobParams->cal21E4 = 0.;
    seobParams->hCoeffs->f21v7cEff00 = 0.0;
    seobParams->hCoeffs->f21v7cEff10 = 0.0;
    seobParams->hCoeffs->f21v7cEff11 = 0.0;
    seobParams->hCoeffs->f21v7cEff01 = 0.0;
    seobParams->hCoeffs->f21v7cEff02 = 0.0;
    return seobParams;
QUIT:
    if (seobParams)
        STRUCTFREE(seobParams, SpinEOBParams);
    return NULL;
}

INT SEOBInitialConditions(REAL8Vector *ICvalues,
                          REAL8 MfMin,
                          REAL8 ecc,
                          SpinEOBParams *seobParams)
{
    PRINT_LOG_INFO(LOG_INFO, "Set initial conditions");
    if (!seobParams)
        return CEV_FAILURE;
    INT j;
    memset((ICvalues)->data, 0, ((ICvalues)->length) * sizeof(REAL8));
    REAL8 eta = seobParams->eta;
    REAL8 m1, m2, mTotal;
    m1 = seobParams->m1;
    m2 = seobParams->m2;
    mTotal = m1 + m2;
    REAL8 fMin = MfMin / (m1 + m2) / CST_MTSUN_SI;
    REAL8 mSpin1data[3] = {0., 0., 0.};
    REAL8 mSpin2data[3] = {0., 0., 0.};
    /* Check given params */
    // if (seobParams->hParams)
    // {
    //     REAL8 d_ini = seobParams->hParams->d_ini; //initial seperation
    //     REAL8 pr_ini = seobParams->hParams->pr_ini;
    //     REAL8 pphi_ini = seobParams->hParams->pphi_ini;
    //     REAL8 ptheta_ini = seobParams->hParams->ptheta_ini;
    //     REAL8 xSph[3] = {d_ini, 0., 0.};
    //     REAL8 pSph[3] = {pr_ini, ptheta_ini, pphi_ini};
    //     REAL8 xCart[3] = {0,0,0};
    //     REAL8 pCart[3] = {0,0,0};
    //     if (d_ini > 0. && pphi_ini > 0.)
    //     {
    //         if (ptheta_ini > 0)
    //             seobParams->alignedSpins = TRUE;
    //         else
    //             pSph[1] = 0.0;
    //         SphericalToCartesian(xCart, pCart, xSph, pSph);
    //         memcpy( ICvalues->data, xCart, sizeof(xCart) );
    //         memcpy( ICvalues->data+3, pCart, sizeof(pCart) );
    //         memcpy( ICvalues->data+6, seobParams->s1Vec->data, sizeof(mSpin1data) );
    //         memcpy( ICvalues->data+9, seobParams->s2Vec->data, sizeof(mSpin1data) );
    //         return CEV_SUCCESS;
    //     }
    // }

    if (seobParams->alignedSpins)
    {
        REAL8 chi1dotZ = seobParams->chi1;
        REAL8 chi2dotZ = seobParams->chi2;
        REAL8 chiS = SEOBCalculateChiS(chi1dotZ, chi2dotZ);
        REAL8 chiA = SEOBCalculateChiA(chi1dotZ, chi2dotZ);
        REAL8 tplspin = SEOBCalculatetplspin(m1, m2, eta, chi1dotZ, chi2dotZ);
        if (XLALSimIMREOBCalcSpinFacWaveformCoefficients(seobParams->hCoeffs, seobParams, tplspin, chiS, chiA) != CEV_SUCCESS)
            return CEV_FAILURE;
        mSpin1data[2] = chi1dotZ * m1 * m1;
        mSpin2data[2] = chi2dotZ * m2 * m2;
        // print_debug("m1 = %g, m2 = %g, fMin = %g\n", m1, m2, fMin);
        // print_debug("mSpin1data = (%g, %g, %g)\n", mSpin1data[0], mSpin1data[1], mSpin1data[2]);
        // print_debug("mSpin2data = (%g, %g, %g)\n", mSpin2data[0], mSpin2data[1], mSpin2data[2]);
        if (EOBInitialConditions(ICvalues, m1, m2, fMin, ecc, 0, mSpin1data, mSpin2data, seobParams) != CEV_SUCCESS)
            return CEV_FAILURE;
    }
    else
    {
        for (j=0; j<3; j++)
        {
            mSpin1data[j] = seobParams->s1Vec->data[j] * mTotal * mTotal;
            mSpin2data[j] = seobParams->s2Vec->data[j] * mTotal * mTotal;
        }
        if (EOBInitialConditionsPrec(ICvalues, m1, m2, fMin, ecc, 0, mSpin1data, mSpin2data, seobParams) != CEV_SUCCESS)
            return CEV_FAILURE;
    }
    return CEV_SUCCESS;
}

INT SEOBInitialConditions_Conserve(REAL8Vector *ICvalues,
                          REAL8 MfMin,
                          REAL8 ecc,
                          SpinEOBParams *seobParams)
{
    PRINT_LOG_INFO(LOG_INFO, "Set conserved initial conditions");
    if (!seobParams)
        return CEV_FAILURE;
    INT j;
    memset((ICvalues)->data, 0, ((ICvalues)->length) * sizeof(REAL8));
    REAL8 eta = seobParams->eta;
    REAL8 m1, m2, mTotal;
    m1 = seobParams->m1;
    m2 = seobParams->m2;
    mTotal = m1 + m2;
    REAL8 fMin = MfMin / (m1 + m2) / CST_MTSUN_SI;
    REAL8 mSpin1data[3] = {0., 0., 0.};
    REAL8 mSpin2data[3] = {0., 0., 0.};
    /* Check given params */
    // if (seobParams->hParams)
    // {
    //     REAL8 d_ini = seobParams->hParams->d_ini; //initial seperation
    //     REAL8 pr_ini = seobParams->hParams->pr_ini;
    //     REAL8 pphi_ini = seobParams->hParams->pphi_ini;
    //     REAL8 ptheta_ini = seobParams->hParams->ptheta_ini;
    //     REAL8 xSph[3] = {d_ini, 0., 0.};
    //     REAL8 pSph[3] = {pr_ini, ptheta_ini, pphi_ini};
    //     REAL8 xCart[3] = {0,0,0};
    //     REAL8 pCart[3] = {0,0,0};
    //     if (d_ini > 0. && pphi_ini > 0.)
    //     {
    //         if (ptheta_ini > 0)
    //             seobParams->alignedSpins = TRUE;
    //         else
    //             pSph[1] = 0.0;
    //         SphericalToCartesian(xCart, pCart, xSph, pSph);
    //         memcpy( ICvalues->data, xCart, sizeof(xCart) );
    //         memcpy( ICvalues->data+3, pCart, sizeof(pCart) );
    //         memcpy( ICvalues->data+6, seobParams->s1Vec->data, sizeof(mSpin1data) );
    //         memcpy( ICvalues->data+9, seobParams->s2Vec->data, sizeof(mSpin1data) );
    //         return CEV_SUCCESS;
    //     }
    // }

    REAL8 tcons;
    tcons = GET_CONSERV_TIME;
    if (tcons < 0)
    {
        if (seobParams->alignedSpins)
        {
            REAL8 chi1dotZ = seobParams->chi1;
            REAL8 chi2dotZ = seobParams->chi2;
            REAL8 chiS = SEOBCalculateChiS(chi1dotZ, chi2dotZ);
            REAL8 chiA = SEOBCalculateChiA(chi1dotZ, chi2dotZ);
            REAL8 tplspin = SEOBCalculatetplspin(m1, m2, eta, chi1dotZ, chi2dotZ);
            if (XLALSimIMREOBCalcSpinFacWaveformCoefficients(seobParams->hCoeffs, seobParams, tplspin, chiS, chiA) != CEV_SUCCESS)
                return CEV_FAILURE;
            mSpin1data[2] = chi1dotZ * m1 * m1;
            mSpin2data[2] = chi2dotZ * m2 * m2;
            // print_debug("m1 = %g, m2 = %g, fMin = %g\n", m1, m2, fMin);
            // print_debug("mSpin1data = (%g, %g, %g)\n", mSpin1data[0], mSpin1data[1], mSpin1data[2]);
            // print_debug("mSpin2data = (%g, %g, %g)\n", mSpin2data[0], mSpin2data[1], mSpin2data[2]);
            // if (EOBInitialConditions(ICvalues, m1, m2, fMin, ecc, 0, mSpin1data, mSpin2data, seobParams) != CEV_SUCCESS)
            //     return CEV_FAILURE;
            if (EOBInitialConditionsPrec(ICvalues, m1, m2, fMin, ecc, 0, mSpin1data, mSpin2data, seobParams) != CEV_SUCCESS)
                return CEV_FAILURE;
        }
        else
        {
            for (j=0; j<3; j++)
            {
                mSpin1data[j] = seobParams->s1Vec->data[j] * mTotal * mTotal;
                mSpin2data[j] = seobParams->s2Vec->data[j] * mTotal * mTotal;
            }
            if (EOBInitialConditionsPrec(ICvalues, m1, m2, fMin, ecc, 0, mSpin1data, mSpin2data, seobParams) != CEV_SUCCESS)
                return CEV_FAILURE;
        }
    } else {
        if (seobParams->alignedSpins)
        {
            REAL8 chi1dotZ = seobParams->chi1;
            REAL8 chi2dotZ = seobParams->chi2;
            REAL8 chiS = SEOBCalculateChiS(chi1dotZ, chi2dotZ);
            REAL8 chiA = SEOBCalculateChiA(chi1dotZ, chi2dotZ);
            REAL8 tplspin = SEOBCalculatetplspin(m1, m2, eta, chi1dotZ, chi2dotZ);
            if (XLALSimIMREOBCalcSpinFacWaveformCoefficients(seobParams->hCoeffs, seobParams, tplspin, chiS, chiA) != CEV_SUCCESS)
                return CEV_FAILURE;
            mSpin1data[2] = chi1dotZ * m1 * m1;
            mSpin2data[2] = chi2dotZ * m2 * m2;
            // print_debug("m1 = %g, m2 = %g, fMin = %g\n", m1, m2, fMin);
            // print_debug("mSpin1data = (%g, %g, %g)\n", mSpin1data[0], mSpin1data[1], mSpin1data[2]);
            // print_debug("mSpin2data = (%g, %g, %g)\n", mSpin2data[0], mSpin2data[1], mSpin2data[2]);
            // if (EOBInitialConditions(ICvalues, m1, m2, fMin, ecc, 0, mSpin1data, mSpin2data, seobParams) != CEV_SUCCESS)
            //     return CEV_FAILURE;
            if (seobParams->p0 > 0)
            {
                if (EOBInitialConditionsPrec_epi(ICvalues, m1, m2, seobParams->p0, ecc, seobParams->x0, mSpin1data, mSpin2data, seobParams) != CEV_SUCCESS)
                    return CEV_FAILURE;
            } else if (EOBInitialConditionsPrec_Conserve(ICvalues, m1, m2, fMin, ecc, 0, mSpin1data, mSpin2data, seobParams) != CEV_SUCCESS)
                return CEV_FAILURE;
        }
        else
        {
            for (j=0; j<3; j++)
            {
                mSpin1data[j] = seobParams->s1Vec->data[j] * mTotal * mTotal;
                mSpin2data[j] = seobParams->s2Vec->data[j] * mTotal * mTotal;
            }
            // print_debug("here, s2x = %.16e\n", mSpin2data[0]);
            if (seobParams->p0 > 0)
            {
                if (EOBInitialConditionsPrec_epi(ICvalues, m1, m2, seobParams->p0, ecc, seobParams->x0, mSpin1data, mSpin2data, seobParams) != CEV_SUCCESS)
                    return CEV_FAILURE;
            } else
                if (EOBInitialConditionsPrec_Conserve(ICvalues, m1, m2, fMin, ecc, 0, mSpin1data, mSpin2data, seobParams) != CEV_SUCCESS)
                    return CEV_FAILURE;
        }
    }

    return CEV_SUCCESS;
}


/* Check that the ringdown frequency of the highest ell mode is less than the
 * Nyquist frequency */
int XLALEOBCheckNyquistFrequency(REAL8 m1, REAL8 m2, REAL8 spin1[3],
                                 REAL8 spin2[3], REAL8 deltaT, UINT ell_max) 
{
    UINT mode_highest_freqL = ell_max;
    UINT mode_highest_freqM = ell_max;
    /* Ringdown freq used to check the sample rate */
    COMPLEX16Vector modefreqVec;
    COMPLEX16 modeFreq;
    modefreqVec.length = 1;
    modefreqVec.data = &modeFreq;

    if (XLALSimIMREOBGenerateQNMFreqV2Prec(&modefreqVec, m1, m2, spin1, spin2,
                                            mode_highest_freqL, mode_highest_freqM,
                                            1) == CEV_FAILURE) 
    {
        return CEV_FAILURE;
    }
    if (deltaT > CST_PI / creal(modeFreq)) {
        PRINT_LOG_INFO(LOG_CRITICAL, "Ringdown frequency > Nyquist");
        return CEV_FAILURE;
    }
// print_debug("deltaT = %.16e, modeFreq = %.16e\n", deltaT, CST_PI /creal(modeFreq));
    return CEV_SUCCESS;
}


INT SEOBIntegrateDynamics(REAL8Array **dynamics,
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
    UINT nb_Hamiltonian_variables_spinsaligned = 4;

    REAL8Vector *values = CreateREAL8Vector(nb_Hamiltonian_variables);
    if (!values) {failed = 1; goto QUIT;}
    memcpy(values->data, ICvalues->data, values->length * sizeof(REAL8));

    REAL8Vector *values_spinaligned = CreateREAL8Vector(nb_Hamiltonian_variables_spinsaligned);
    if (!values_spinaligned) {failed = 1; goto QUIT;}
    memset(values_spinaligned->data, 0, values_spinaligned->length * sizeof(REAL8));

    ARKIntegrator *integrator = NULL;
    if (seobParams->alignedSpins)
    {
        // Spin Aligned
        REAL8 temp_r = sqrt(ICvalues->data[0] * ICvalues->data[0] +
                            ICvalues->data[1] * ICvalues->data[1] +
                            ICvalues->data[2] * ICvalues->data[2]);
        REAL8 temp_phi = ICvalues->data[12];

        values_spinaligned->data[0] = temp_r;   // General form of r
        values_spinaligned->data[1] = temp_phi; // phi
        values_spinaligned->data[2] = 
            ICvalues->data[3] * cos(temp_phi) +
            ICvalues->data[4] * sin(temp_phi); // p_r^*
        values_spinaligned->data[3] =
            temp_r * (ICvalues->data[4] * cos(temp_phi) -
                    ICvalues->data[3] * sin(temp_phi)); // p_phi

        /* We have to use different stopping conditions depending
        we are in the low-sampling or high-sampling portion
        of the waveform. We can tell this apart because for the low-sampling (or
        ada sampling) we always start at t=0
        */
        if (tstart > 0) 
        {
            // High sampling
            if (seobParams->e0 != 0)
                integrator = XLALAdaptiveRungeKutta4Init(
                    nb_Hamiltonian_variables_spinsaligned, XLALSpinAlignedHcapDerivative,
                    EOBHighSRStopConditionEcc, EPS_ABS, EPS_REL);
            else
                integrator = XLALAdaptiveRungeKutta4Init(
                    nb_Hamiltonian_variables_spinsaligned, XLALSpinAlignedHcapDerivative,
                    XLALSpinPrecAlignedHiSRStopCondition, EPS_ABS, EPS_REL);
        } else {
            // Low sampling
            integrator = XLALAdaptiveRungeKutta4Init(
                nb_Hamiltonian_variables_spinsaligned, XLALSpinAlignedHcapDerivative,
                XLALEOBSpinPrecAlignedStopCondition, EPS_ABS, EPS_REL);
        }
    }
    else
    {
        // Prec
        if (tstart > 0) 
        {
            integrator = XLALAdaptiveRungeKutta4Init(
                nb_Hamiltonian_variables, XLALSpinPrecHcapNumericalDerivative,
                XLALEOBSpinPrecStopConditionBasedOnPR, EPS_ABS, EPS_REL);
        } else {
            integrator = XLALAdaptiveRungeKutta4Init(
                nb_Hamiltonian_variables, XLALSpinPrecHcapNumericalDerivative,
                XLALEOBSpinPrecStopConditionBasedOnPR, EPS_ABS, EPS_REL);
        }
    }

    if (!integrator) {failed = 1; goto QUIT;}

    /* Ensure that integration stops ONLY when the stopping condition is True */
    integrator->stopontestonly = 1;
    /* When this option is set to 0, the integration can be exceedingly slow for
    * spin-aligned systems */
    integrator->retries = 1;
    /* Computing the dynamical evolution of the system */
    if (seobParams->alignedSpins)
    {
        // Spin Aligned
        // flagConstantSampling = seobParams->alignedSpins
        // EOBversion = 2
        PRINT_LOG_INFO(LOG_INFO, "Adaptive RungeKutta4 No Interpolation");
        if (!flagConstantSampling)
            retLen = XLALAdaptiveRungeKutta4NoInterpolate(
                integrator, seobParams, values_spinaligned->data, 0., tend - tstart,
                deltaT, deltaT_min, &dynamics_spinaligned);
        else
            retLen = XLALAdaptiveRungeKutta4(integrator, seobParams, 
                values_spinaligned->data, 0., tend - tstart, 
                deltaT, &dynamics_spinaligned);

        if (retLen < 0)
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Integration Failed!!!");
            failed = 1;
            goto QUIT;
        }
        /* Convert the spin-aligned dynamics to a generic-spins dynamics */
        PRINT_LOG_INFO(LOG_INFO, "Convert Spin Aligned Dynamics to Generic Spins");
        status = SEOBConvertSpinAlignedDynamicsToGenericSpins(
            dynamics, dynamics_spinaligned, retLen, 
            seobParams->chi1, seobParams->chi2, seobParams);
        if (status != CEV_SUCCESS)
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Failure in dynamics transformation.");
            failed = 1;
            goto QUIT;
        }
    } 
    else 
    {
        // Prec
        if (!flagConstantSampling)
            retLen = XLALAdaptiveRungeKutta4NoInterpolate(integrator, seobParams, 
                values->data, 0., tend-tstart, 
                deltaT, deltaT_min, dynamics);
        else
            retLen = XLALAdaptiveRungeKutta4(integrator, seobParams, values->data, 0.,
                                            tend - tstart, deltaT, dynamics);
        if (retLen < 0)
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Integration Failed!!!");
            failed = 1;
            goto QUIT;
        }
    }
    PRINT_LOG_INFO(LOG_INFO, "Integration End");
    // NOTE: functions like XLALAdaptiveRungeKutta4 would give nans if the times
    // do not start at 0 -- we have to adjust the starting time after integration
    /* Adjust starting time */
    for ( i = 0; i < retLen; i++)
        (*dynamics)->data[i] += tstart;

QUIT:
    STRUCTFREE(values, REAL8Vector);
    STRUCTFREE(values_spinaligned, REAL8Vector);
    STRUCTFREE(integrator, ARKIntegrator);
    STRUCTFREE(dynamics_spinaligned, REAL8Array);
    *retLenOut = retLen;
    if (failed)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}


INT SEOBIntegrateDynamics_SA(REAL8Array **dynamics,
                          INT *retLenOut,
                          REAL8Vector *ICvalues,
                          REAL8 EPS_ABS,
                          REAL8 EPS_REL,
                          REAL8 deltaT,
                          REAL8 deltaT_min,
                          REAL8 tstart,
                          REAL8 tend ,
                          SpinEOBParams *seobParams)
{
    INT retLen;
    UINT i;
    REAL8Array *dynamics_spinaligned = NULL;
    INT status, failed = 0;
    /* Dimensions of vectors of dynamical variables to be integrated */
    // UINT nb_Hamiltonian_variables = 14;
    UINT nb_Hamiltonian_variables_spinsaligned = 4;

    REAL8Vector *values_spinaligned = CreateREAL8Vector(nb_Hamiltonian_variables_spinsaligned);
    if (!values_spinaligned) {failed = 1; goto QUIT;}
    memcpy(values_spinaligned->data, ICvalues->data, values_spinaligned->length * sizeof(REAL8));

    ARKIntegrator *integrator = NULL;

    // REAL8 temp_r = sqrt(ICvalues->data[0] * ICvalues->data[0] +
    //                     ICvalues->data[1] * ICvalues->data[1] +
    //                     ICvalues->data[2] * ICvalues->data[2]);
    // REAL8 temp_phi = ICvalues->data[12];

    // values_spinaligned->data[0] = temp_r;   // General form of r
    // values_spinaligned->data[1] = temp_phi; // phi
    // values_spinaligned->data[2] = 
    //     ICvalues->data[3] * cos(temp_phi) +
    //     ICvalues->data[4] * sin(temp_phi); // p_r^*
    // values_spinaligned->data[3] =
    //     temp_r * (ICvalues->data[4] * cos(temp_phi) -
    //             ICvalues->data[3] * sin(temp_phi)); // p_phi

    /* We have to use different stopping conditions depending
    we are in the low-sampling or high-sampling portion
    of the waveform. We can tell this apart because for the low-sampling (or
    ada sampling) we always start at t=0 */
    if (tstart > 0) 
    {
        // High sampling
        if (seobParams->e0 != 0)
            integrator = XLALAdaptiveRungeKutta4Init(
                nb_Hamiltonian_variables_spinsaligned, XLALSpinAlignedHcapDerivative_SA,
                EOBHighSRStopConditionEcc, EPS_ABS, EPS_REL);
        else
            integrator = XLALAdaptiveRungeKutta4Init(
                nb_Hamiltonian_variables_spinsaligned, XLALSpinAlignedHcapDerivative_SA,
                XLALSpinPrecAlignedHiSRStopCondition, EPS_ABS, EPS_REL);
    } else {
        // Low sampling
        integrator = XLALAdaptiveRungeKutta4Init(
            nb_Hamiltonian_variables_spinsaligned, XLALSpinAlignedHcapDerivative_SA,
            XLALEOBSpinPrecAlignedStopCondition, EPS_ABS, EPS_REL);
    }

    if (!integrator) {failed = 1; goto QUIT;}

    /* Ensure that integration stops ONLY when the stopping condition is True */
    integrator->stopontestonly = 1;
    /* When this option is set to 0, the integration can be exceedingly slow for
    * spin-aligned systems */
    integrator->retries = 1;
    /* Computing the dynamical evolution of the system */
    // Spin Aligned
    // flagConstantSampling = seobParams->alignedSpins
    // EOBversion = 2
    if (tstart > 0) 
    {
        PRINT_LOG_INFO(LOG_INFO, "Adaptive RungeKutta4");
            retLen = XLALAdaptiveRungeKutta4WithDeriv(
                integrator, seobParams, values_spinaligned->data, 0., tend - tstart,
                deltaT, &dynamics_spinaligned);
    } else {
        PRINT_LOG_INFO(LOG_INFO, "Adaptive RungeKutta4 No Interpolation");
            retLen = XLALAdaptiveRungeKutta4NoInterpolateWithDeriv(
                integrator, seobParams, values_spinaligned->data, 0., tend - tstart,
                deltaT, deltaT_min, &dynamics_spinaligned);
            // retLen = XLALAdaptiveRungeKutta4NoInterpolate(
            //     integrator, seobParams, values_spinaligned->data, 0., tend-tstart,
            //     deltaT, deltaT_min, &dynamics_spinaligned);
    }
    
    if (retLen < 0)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Integration Failed!!!");
        failed = 1;
        goto QUIT;
    }
    PRINT_LOG_INFO(LOG_INFO, "Integration End, retLen = %d", retLen);

    // NOTE: functions like XLALAdaptiveRungeKutta4 would give nans if the times
    // do not start at 0 -- we have to adjust the starting time after integration
    /* Adjust starting time */
    for ( i = 0; i < retLen; i++)
        dynamics_spinaligned->data[i] += tstart;
    *dynamics = dynamics_spinaligned;
QUIT:
    STRUCTFREE(values_spinaligned, REAL8Vector);
    STRUCTFREE(integrator, ARKIntegrator);
    *retLenOut = retLen;
    if (failed)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}

INT SEOBIntegrateDynamics_adaptive(REAL8Array **dynamics,
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
    UINT nb_Hamiltonian_variables_spinsaligned = 4;

    REAL8Vector *values = CreateREAL8Vector(nb_Hamiltonian_variables);
    if (!values) {failed = 1; goto QUIT;}
    memcpy(values->data, ICvalues->data, values->length * sizeof(REAL8));

    REAL8Vector *values_spinaligned = CreateREAL8Vector(nb_Hamiltonian_variables_spinsaligned);
    if (!values_spinaligned) {failed = 1; goto QUIT;}
    memset(values_spinaligned->data, 0, values_spinaligned->length * sizeof(REAL8));

    ARKIntegrator *integrator = NULL;
    if (seobParams->alignedSpins)
    {
        // Spin Aligned
        REAL8 temp_r = sqrt(ICvalues->data[0] * ICvalues->data[0] +
                            ICvalues->data[1] * ICvalues->data[1] +
                            ICvalues->data[2] * ICvalues->data[2]);
        REAL8 temp_phi = ICvalues->data[12];

        values_spinaligned->data[0] = temp_r;   // General form of r
        values_spinaligned->data[1] = temp_phi; // phi
        values_spinaligned->data[2] = 
            ICvalues->data[3] * cos(temp_phi) +
            ICvalues->data[4] * sin(temp_phi); // p_r^*
        values_spinaligned->data[3] =
            temp_r * (ICvalues->data[4] * cos(temp_phi) -
                    ICvalues->data[3] * sin(temp_phi)); // p_phi

        /* We have to use different stopping conditions depending
        we are in the low-sampling or high-sampling portion
        of the waveform. We can tell this apart because for the low-sampling (or
        ada sampling) we always start at t=0
        */
        if (tstart > 0) 
        {
            // High sampling
            if (seobParams->e0 != 0)
                integrator = XLALAdaptiveRungeKutta4Init(
                    nb_Hamiltonian_variables_spinsaligned, XLALSpinAlignedHcapDerivative,
                    EOBHighSRStopConditionEcc, EPS_ABS, EPS_REL);
            else
                integrator = XLALAdaptiveRungeKutta4Init(
                    nb_Hamiltonian_variables_spinsaligned, XLALSpinAlignedHcapDerivative,
                    XLALSpinPrecAlignedHiSRStopCondition, EPS_ABS, EPS_REL);
        } else {
            // Low sampling
            integrator = XLALAdaptiveRungeKutta4Init(
                nb_Hamiltonian_variables_spinsaligned, XLALSpinAlignedHcapDerivative,
                XLALEOBSpinPrecAlignedStopCondition, EPS_ABS, EPS_REL);
        }
    }
    else
    {
        // Prec
        if (tstart > 0) 
        {
            integrator = XLALAdaptiveRungeKutta4Init(
                nb_Hamiltonian_variables, XLALSpinPrecHcapNumericalDerivative,
                XLALEOBSpinPrecStopConditionBasedOnPR, EPS_ABS, EPS_REL);
        } else {
            integrator = XLALAdaptiveRungeKutta4Init(
                nb_Hamiltonian_variables, XLALSpinPrecHcapNumericalDerivative,
                XLALEOBSpinPrecStopConditionBasedOnPR, EPS_ABS, EPS_REL);
        }
    }

    if (!integrator) {failed = 1; goto QUIT;}

    /* Ensure that integration stops ONLY when the stopping condition is True */
    integrator->stopontestonly = 1;
    /* When this option is set to 0, the integration can be exceedingly slow for
    * spin-aligned systems */
    integrator->retries = 1;
    /* Computing the dynamical evolution of the system */
    if (seobParams->alignedSpins)
    {
        // Spin Aligned
        // flagConstantSampling = seobParams->alignedSpins
        // EOBversion = 2
        PRINT_LOG_INFO(LOG_INFO, "Adaptive RungeKutta4 No Interpolation");
        // if (!flagConstantSampling)
            retLen = XLALAdaptiveRungeKutta4NoInterpolate(
                integrator, seobParams, values_spinaligned->data, 0., tend - tstart,
                deltaT, deltaT_min, &dynamics_spinaligned);
        // else
            // retLen = XLALAdaptiveRungeKutta4(integrator, seobParams, 
            //     values_spinaligned->data, 0., tend - tstart, 
            //     deltaT, &dynamics_spinaligned);

        if (retLen < 0)
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Integration Failed!!!");
            failed = 1;
            goto QUIT;
        }
        /* Convert the spin-aligned dynamics to a generic-spins dynamics */
        PRINT_LOG_INFO(LOG_INFO, "Convert Spin Aligned Dynamics to Generic Spins");
        status = SEOBConvertSpinAlignedDynamicsToGenericSpins(
            dynamics, dynamics_spinaligned, retLen, 
            seobParams->chi1, seobParams->chi2, seobParams);
        if (status != CEV_SUCCESS)
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Failure in dynamics transformation.");
            failed = 1;
            goto QUIT;
        }
    } 
    else 
    {
        // Prec
        if (!flagConstantSampling)
            retLen = XLALAdaptiveRungeKutta4NoInterpolate(integrator, seobParams, 
                values->data, 0., tend-tstart, 
                deltaT, deltaT_min, dynamics);
        else
            retLen = XLALAdaptiveRungeKutta4(integrator, seobParams, values->data, 0.,
                                            tend - tstart, deltaT, dynamics);
        if (retLen < 0)
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Integration Failed!!!");
            failed = 1;
            goto QUIT;
        }
    }
    PRINT_LOG_INFO(LOG_INFO, "Integration End");
    // NOTE: functions like XLALAdaptiveRungeKutta4 would give nans if the times
    // do not start at 0 -- we have to adjust the starting time after integration
    /* Adjust starting time */
    for ( i = 0; i < retLen; i++)
        (*dynamics)->data[i] += tstart;

QUIT:
    STRUCTFREE(values, REAL8Vector);
    STRUCTFREE(values_spinaligned, REAL8Vector);
    STRUCTFREE(integrator, ARKIntegrator);
    STRUCTFREE(dynamics_spinaligned, REAL8Array);
    *retLenOut = retLen;
    if (failed)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}

INT SEOBIntegrateDynamics_Conserve(REAL8Array **dynamics,
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
    UINT nb_Hamiltonian_variables_spinsaligned = 4;
    REAL8 time_end;
    if (tend < 0)
        time_end = 1000.;
    else
        time_end = tend;
    // print_debug("tend = %.16e\n", tend);
    REAL8Vector *values = CreateREAL8Vector(nb_Hamiltonian_variables);
    if (!values) {failed = 1; goto QUIT;}
    memcpy(values->data, ICvalues->data, values->length * sizeof(REAL8));

    REAL8Vector *values_spinaligned = CreateREAL8Vector(nb_Hamiltonian_variables_spinsaligned);
    if (!values_spinaligned) {failed = 1; goto QUIT;}
    memset(values_spinaligned->data, 0, values_spinaligned->length * sizeof(REAL8));

    ARKIntegrator *integrator = NULL;
    if (seobParams->alignedSpins)
    {
        // Spin Aligned
        REAL8 temp_r = sqrt(ICvalues->data[0] * ICvalues->data[0] +
                            ICvalues->data[1] * ICvalues->data[1] +
                            ICvalues->data[2] * ICvalues->data[2]);
        REAL8 temp_phi = ICvalues->data[12];

        values_spinaligned->data[0] = temp_r;   // General form of r
        values_spinaligned->data[1] = temp_phi; // phi
        values_spinaligned->data[2] = 
            ICvalues->data[3] * cos(temp_phi) +
            ICvalues->data[4] * sin(temp_phi); // p_r^*
        values_spinaligned->data[3] =
            temp_r * (ICvalues->data[4] * cos(temp_phi) -
                    ICvalues->data[3] * sin(temp_phi)); // p_phi

        /* We have to use different stopping conditions depending
        we are in the low-sampling or high-sampling portion
        of the waveform. We can tell this apart because for the low-sampling (or
        ada sampling) we always start at t=0
        */
        if (CONSERVE_FLAG == 0) 
        {
            // High sampling
            // if (seobParams->e0 != 0)
            //     integrator = XLALAdaptiveRungeKutta4Init(
            //         nb_Hamiltonian_variables_spinsaligned, XLALSpinAlignedHcapDerivative,
            //         EOBHighSRStopConditionEcc, EPS_ABS, EPS_REL);
            // else
            PRINT_LOG_INFO(LOG_INFO, "Set unconserved dynamics : XLALSpinAlignedHcapDerivative");
            integrator = XLALAdaptiveRungeKutta4Init(
                nb_Hamiltonian_variables_spinsaligned, XLALSpinAlignedHcapDerivative,
                XLALEOBSpinPrecAlignedStopCondition, EPS_ABS, EPS_REL);
        } else {
            // Low sampling
            PRINT_LOG_INFO(LOG_INFO, "Set conserved dynamics : XLALSpinAlignedHcapDerivative_Conserve");
            integrator = XLALAdaptiveRungeKutta4Init(
                nb_Hamiltonian_variables_spinsaligned, XLALSpinAlignedHcapDerivative_Conserve,
                XLALEOBSpinPrecAlignedStopCondition, EPS_ABS, EPS_REL);
        }
    }
    else
    {
        // Prec
        // print_debug("prec...");
        if (CONSERVE_FLAG == 0) 
        {
            PRINT_LOG_INFO(LOG_INFO, "Set unconserved prec dynamics : XLALSpinPrecHcapNumericalDerivative");
            PRINT_LOG_INFO(LOG_INFO, "Set unconserved dynamics");
            integrator = XLALAdaptiveRungeKutta4Init(
                nb_Hamiltonian_variables, XLALSpinPrecHcapNumericalDerivative,
                XLALEOBSpinPrecStopConditionBasedOnPR, EPS_ABS, EPS_REL);
        } else {
            PRINT_LOG_INFO(LOG_INFO, "Set conserved prec dynamics : XLALSpinPrecHcapNumericalDerivative_Conserve");
            integrator = XLALAdaptiveRungeKutta4Init(
                nb_Hamiltonian_variables, XLALSpinPrecHcapNumericalDerivative_Conserve,
                XLALEOBSpinPrecStopConditionBasedOnPR, EPS_ABS, EPS_REL);
        }
    }

    if (!integrator) {failed = 1; goto QUIT;}

    /* The integration stops when the end time reached */
    integrator->stopontestonly = 0;

    /* When this option is set to 0, the integration can be exceedingly slow for
    * spin-aligned systems */
    integrator->retries = 1;
    /* Computing the dynamical evolution of the system */
    if (seobParams->alignedSpins)
    {
        // Spin Aligned
        // flagConstantSampling = seobParams->alignedSpins
        // EOBversion = 2
        PRINT_LOG_INFO(LOG_INFO, "Adaptive RungeKutta4 No Interpolation");
        if (!flagConstantSampling)
        {
            PRINT_LOG_INFO(LOG_INFO, "Run Integrator : XLALAdaptiveRungeKutta4NoInterpolate, tend = %f", time_end - tstart);
            retLen = XLALAdaptiveRungeKutta4NoInterpolate(
                integrator, seobParams, values_spinaligned->data, 0., time_end - tstart,
                deltaT, deltaT_min, &dynamics_spinaligned);
        }
        else
        {
            PRINT_LOG_INFO(LOG_INFO, "Run Integrator : XLALAdaptiveRungeKutta4, tend = %f", time_end - tstart);
            retLen = XLALAdaptiveRungeKutta4(integrator, seobParams, 
                values_spinaligned->data, 0., time_end - tstart, 
                deltaT, &dynamics_spinaligned);
        }

        if (retLen < 0)
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Integration Failed!!!");
            failed = 1;
            goto QUIT;
        }
        /* Convert the spin-aligned dynamics to a generic-spins dynamics */
        PRINT_LOG_INFO(LOG_INFO, "Convert Spin Aligned Dynamics to Generic Spins");
        status = SEOBConvertSpinAlignedDynamicsToGenericSpins(
            dynamics, dynamics_spinaligned, retLen, 
            seobParams->chi1, seobParams->chi2, seobParams);
        if (status != CEV_SUCCESS)
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Failure in dynamics transformation.");
            failed = 1;
            goto QUIT;
        }
    } 
    else 
    {
        // Prec
        if (!flagConstantSampling)
        {
            PRINT_LOG_INFO(LOG_INFO, "Run Integrator : XLALAdaptiveRungeKutta4NoInterpolate, tend = %f", time_end - tstart);
            retLen = XLALAdaptiveRungeKutta4NoInterpolate(integrator, seobParams, 
                values->data, 0., time_end-tstart, 
                deltaT, deltaT_min, dynamics);
        }
        else
        {
            PRINT_LOG_INFO(LOG_INFO, "Run Integrator : XLALAdaptiveRungeKutta4, tend = %f", time_end - tstart);
            retLen = XLALAdaptiveRungeKutta4(integrator, seobParams, values->data, 0.,
                                            time_end - tstart, deltaT, dynamics);
        }
        if (retLen < 0)
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Integration Failed!!!");
            failed = 1;
            goto QUIT;
        }
    }
    PRINT_LOG_INFO(LOG_INFO, "Integration End");
    // NOTE: functions like XLALAdaptiveRungeKutta4 would give nans if the times
    // do not start at 0 -- we have to adjust the starting time after integration
    /* Adjust starting time */
    for ( i = 0; i < retLen; i++)
        (*dynamics)->data[i] += tstart;

QUIT:
    STRUCTFREE(values, REAL8Vector);
    STRUCTFREE(values_spinaligned, REAL8Vector);
    STRUCTFREE(integrator, ARKIntegrator);
    STRUCTFREE(dynamics_spinaligned, REAL8Array);
    *retLenOut = retLen;
    if (failed)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}

#if 1
INT TortoiseConvertMatrix(SpinEOBParams *params,
                          REAL8 values[], REAL8 dvalues[],
                          REAL8 pVecRet[])
{
    int i, j;
    REAL8 eta = params->eta;
    REAL8 mass1 = params->m1;
    REAL8 mass2 = params->m2;
    REAL8 tmpP[3]= {0.}, rMag, rMag2, prT;
    REAL8 u, u2, u3, u4, u5, w2, a2;
    REAL8 D, m1PlusetaKK, bulk, logTerms, deltaU, deltaT, deltaR;
    REAL8 eobD_r, deltaU_u, deltaU_r, deltaT_r;
    REAL8 dcsi, csi;
    REAL8Vector rVec, pVec;
    REAL8 rData[3] = {0.}, pData[3] = {0.};
    REAL8Vector s1, s2, s1norm, s2norm, sKerr, sStar;
    REAL8       s1Data[3]= {0.}, s2Data[3]= {0.}, s1DataNorm[3]= {0.}, s2DataNorm[3]= {0.};
    REAL8       sKerrData[3]= {0.}, sStarData[3]= {0.};
    REAL8 Tmatrix[3][3]= {{0.}}, invTmatrix[3][3]= {{0.}}, dTijdXk[3][3][3]= {{{0.}}};
    SpinEOBHCoeffs *coeffs = (SpinEOBHCoeffs*) params->seobCoeffs;
    rVec.length = pVec.length = 3;
    rVec.data   = rData;
    pVec.data   = pData;
    memcpy( rData, values, sizeof(rData) );
    memcpy( pData, values+3, sizeof(pData) );
    s1.length = s2.length = s1norm.length = s2norm.length = 3;
    s1.data = s1Data;
    s2.data = s2Data;
    s1norm.data = s1DataNorm;
    s2norm.data = s2DataNorm;

    memcpy( s1Data, values+6, 3*sizeof(REAL8) );
    memcpy( s2Data, values+9, 3*sizeof(REAL8) );
    memcpy( s1DataNorm, values+6, 3*sizeof(REAL8) );
    memcpy( s2DataNorm, values+9, 3*sizeof(REAL8) );
    sKerr.length = 3;
    sKerr.data   = sKerrData;
    XLALSimIMRSpinEOBCalculateSigmaKerr( &sKerr, mass1, mass2, &s1, &s2 );

    sStar.length = 3;
    sStar.data   = sStarData;
    XLALSimIMRSpinEOBCalculateSigmaStar( &sStar, mass1, mass2, &s1, &s2 );
    REAL8 a = sqrt(sKerr.data[0]*sKerr.data[0] + sKerr.data[1]*sKerr.data[1]
        + sKerr.data[2]*sKerr.data[2]);

    rMag = sqrt(rData[0]*rData[0] + rData[1]*rData[1] + rData[2]*rData[2]);
    prT = pData[0]*(rData[0]/rMag) + pData[1]*(rData[1]/rMag)
                    + pData[2]*(rData[2]/rMag);

    rMag2 = rMag * rMag;
    u  = 1./rMag;
    u2 = u*u;
    u3 = u2*u;
    u4 = u2*u2;
    u5 = u4*u;
    a2 = a*a;
    w2 = rMag2 + a2;
    /* Eq. 5.83 of BB1, inverse */
    D = 1. + log(1. + 6.*eta*u2 + 2.*(26. - 3.*eta)*eta*u3);
    eobD_r =  (u2/(D*D))*(12.*eta*u + 6.*(26. - 3.*eta)*eta*u2)/(1.
        + 6.*eta*u2 + 2.*(26. - 3.*eta)*eta*u3);
    m1PlusetaKK = -1. + eta * coeffs->KK;
    /* Eq. 5.75 of BB1 */
    bulk = 1./(m1PlusetaKK*m1PlusetaKK) + (2.*u)/m1PlusetaKK + a2*u2;
    /* Eq. 5.73 of BB1 */
    logTerms = 1. + eta*coeffs->k0 + eta*log(fabs(1. + coeffs->k1*u
    + coeffs->k2*u2 + coeffs->k3*u3 + coeffs->k4*u4
    + coeffs->k5*u5 + coeffs->k5l*u5*log(u)));
    /* Eq. 5.73 of BB1 */
    deltaU = bulk*logTerms;
    deltaU = fabs(deltaU);

    /* Eq. 5.71 of BB1 */
    deltaT = rMag2*deltaU;
    /* ddeltaU/du */
    deltaU_u = 2.*(1./m1PlusetaKK + a2*u)*logTerms +
    bulk * (eta*(coeffs->k1 + u*(2.*coeffs->k2 + u*(3.*coeffs->k3
    + u*(4.*coeffs->k4 + 5.*(coeffs->k5+coeffs->k5l*log(u))*u)))))
    / (1. + coeffs->k1*u + coeffs->k2*u2 + coeffs->k3*u3
    + coeffs->k4*u4 + (coeffs->k5+coeffs->k5l*log(u))*u5);
    deltaU_r = -u2 * deltaU_u;
    /* Eq. 5.38 of BB1 */
    deltaR = deltaT*D;
    if ( params->tortoise )
        csi = sqrt( fabs(deltaT * deltaR) )/ w2; /* Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
    else
        csi = 1.0;

    for( i = 0; i < 3; i++ )
    {
        tmpP[i] = pData[i] - (rData[i]/rMag) * prT * (csi-1.)/csi;
    }

    for( i = 0; i < 3; i++ )
        for( j = 0; j <= i; j++ )
        {
            Tmatrix[i][j] = Tmatrix[j][i] = (rData[i]*rData[j]/rMag2)
                                        * (csi - 1.);

            invTmatrix[i][j] = invTmatrix[j][i] =
                    - (csi - 1)/csi * (rData[i]*rData[j]/rMag2);

            if( i==j ){
                Tmatrix[i][j]++;
                invTmatrix[i][j]++;  }
        }
    pVecRet[0] = tmpP[0];
    pVecRet[1] = tmpP[1];
    pVecRet[2] = tmpP[2];
    return CEV_SUCCESS;
}
#endif

INT SEOBComputeExtendedSEOBSAdynamics(SEOBSAdynamics **seobsadynamics,
                                      REAL8Array *dynamics,
                                      INT length,
                                      SpinEOBParams *seobParams)
{
    SEOBSAdynamics *ret = NULL;
    INT status, failed = 0;
    UINT jj;
    ret = (SEOBSAdynamics *) MYMalloc(sizeof(SEOBSAdynamics));
    if (!ret)
        return CEV_FAILURE;
    ret->length = length;
    ret->array = dynamics;
    if (!ret->array)
        return CEV_FAILURE;
    ret->tVec = ret->array->data;

    ret->rVec = ret->array->data + length;
    ret->phiVec = ret->array->data + 2*length;
    ret->prTVec = ret->array->data + 3*length;
    ret->pphiVec = ret->array->data + 4*length;

    ret->drVec = ret->array->data + 5*length;
    ret->dphiVec = ret->array->data + 6*length;
    ret->HVec = ret->array->data + 7*length;
    
    ret->dprTVec = ret->array->data + 8*length;
    ret->dpphiVec = ret->array->data + 9*length;
    *seobsadynamics = ret;
    return CEV_SUCCESS;
}

INT SEOBComputeExtendedSEOBdynamics(SEOBdynamics **seobdynamics,
                                    REAL8Array *dynamics,
                                    INT retLen,
                                    SpinEOBParams *seobParams)
{
    SEOBdynamics *seobdyn = NULL;
    INT status, failed = 0;
    UINT jj;
    seobdyn = CreateSEOBdynamics(retLen);
    // flagSEOBNRv4P_Zframe = 0
    // flagSEOBNRv4P_hamiltonian_derivative = 1
    /* Local variables */
    REAL8 rvec[3] = {0, 0, 0};
    REAL8 pvec[3] = {0, 0, 0};
    REAL8 spin1vec[3] = {0, 0, 0}; /* in units of mTotal^2 */
    REAL8 spin2vec[3] = {0, 0, 0}; /* in units of mTotal^2 */
    REAL8 rdotvec[3] = {0, 0, 0};
    REAL8 rcrossrdot[3] = {0, 0, 0};
    REAL8 rcrossp[3] = {0, 0, 0};
    REAL8 LNhat[3] = {0, 0, 0};
    REAL8 Lhat[3] = {0, 0, 0};
    REAL8 polarr, polarphi, polarpr, polarpphi, omega, s1dotZ, s2dotZ, ham;
    REAL8 nhat[3] = {0,0,0};
    /* Allocate temporary vectors values, dvalues */
    REAL8Vector *values = NULL;
    REAL8Vector *dvalues = NULL;
    values = CreateREAL8Vector(14);
    dvalues = CreateREAL8Vector(14);
    memset(values->data, 0, (values->length) * sizeof(REAL8));
    memset(dvalues->data, 0, (dvalues->length) * sizeof(REAL8));

    REAL8 m1 = seobParams->m1;
    REAL8 m2 = seobParams->m2;
    REAL8 mtot = m1 + m2;
    REAL8 m1sq = m1*m1/(m1+m2)/(m1+m2);
    REAL8 m2sq = m2*m2/(m1+m2)/(m1+m2);
    REAL8 eta = seobParams->eta;
    SpinEOBHCoeffs *seobCoeffs = seobParams->seobCoeffs;
    /* Copying directly the dynamics data in seobdynamics - 15 vectors of length
    * retLen */
    memcpy(seobdyn->array->data, dynamics->data, 15 * retLen * sizeof(REAL8));
    REAL8 c1x, c1y, c1z, c2x, c2y, c2z;
    REAL8 yhat[3] = {0,0,0};
    REAL8 magY;
    REAL8Vector *sigmaStar = NULL, *sigmaKerr = NULL;
    sigmaStar = CreateREAL8Vector(3);
    sigmaKerr = CreateREAL8Vector(3);
    memset(sigmaStar->data, 0, 3 * sizeof(REAL8));
    memset(sigmaKerr->data, 0, 3 * sizeof(REAL8));
    REAL8 flux, dr, ncrv;
	REAL8Vector	polarDynamics, cartDynamics;
	REAL8		polData[4];
    REAL8       cartData[12];
    polarDynamics.length = 4;
    polarDynamics.data = polData;
    cartDynamics.length = 12;
    cartDynamics.data = cartData;

    // DEBUG: RR force with Schott correction
    RRForceCoeffs coeffsFr;
    RRForceCoeffs coeffsFf;
    REAL8 Ff, Fr, psq;
    CalculateRRForceCoeffs(&coeffsFf, &coeffsFr, seobParams);
    /* Loop to compute the derived quantities from the dynamics */
    UINT i, j;
    for (i=0; i<retLen; i++)
    {
        for (j = 0; j < 14; j++) 
        {
            values->data[j] = dynamics->data[i + (j + 1) * retLen];
        }
        status = XLALSpinPrecHcapRvecDerivative(0, values->data, dvalues->data,
                                         (void *)seobParams);
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
        for (j=0; j<3; j++)
        {
            rvec[j] = values->data[j];
            pvec[j] = values->data[3 + j];
            spin1vec[j] = values->data[6 + j];
            spin2vec[j] = values->data[9 + j];
            rdotvec[j] = dvalues->data[j];
        }
        cross_product3d(rvec, pvec, rcrossp);
        cross_product3d(rvec, rdotvec, rcrossrdot);
        REAL8 rcrossrdotNorm = sqrt(inner_product3d(rcrossrdot, rcrossrdot));
        for (j = 0; j < 3; j++) 
            LNhat[j] = rcrossrdot[j] / rcrossrdotNorm;
        /* Polar dynamics */
        polarr = sqrt(inner_product3d(rvec, rvec));
        polarpr = inner_product3d(rvec, pvec) / polarr;
        polarphi = values->data[12] + values->data[13];
        // REAL8 chi1dotZ = s1dotZ * mtot * mtot / (m1 * m1);
        // REAL8 chi2dotZ = s2dotZ * mtot * mtot / (m2 * m2);
        // FIX ME
        cross_product3d(LNhat, rvec, yhat);
        magY = sqrt(yhat[0]*yhat[0] + yhat[1]*yhat[1] + yhat[2]*yhat[2]);
        // S1 = sqrt((spin1vec[0]*spin1vec[0] + spin1vec[1]*spin1vec[1] + spin1vec[2]*spin1vec[2]))/m1sq;
        // S2 = sqrt((spin2vec[0]*spin2vec[0] + spin2vec[1]*spin2vec[1] + spin2vec[2]*spin2vec[2]))/m2sq;

        c1x = inner_product3d(spin1vec, rvec) / polarr / m1sq;
        c1z = inner_product3d(spin1vec, LNhat) / m1sq;
        // c1y = sqrt( (spin1vec[0]*spin1vec[0] + spin1vec[1]*spin1vec[1] + spin1vec[2]*spin1vec[2])/(m1sq*m1sq) - c1x*c1x - c1z*c1z);
        c1y = inner_product3d(spin1vec, yhat) / m1sq / magY;

        c2x = inner_product3d(spin2vec, rvec) / polarr / m2sq;
        c2z = inner_product3d(spin2vec, LNhat) / m2sq;
        // c2y = sqrt( (spin2vec[0]*spin2vec[0] + spin2vec[1]*spin2vec[1] + spin2vec[2]*spin2vec[2])/(m2sq*m2sq) - c2x*c2x - c2z*c2z);
        c2y = inner_product3d(spin2vec, yhat) / m2sq / magY;

        seobdyn->chiAxVec[i] = 0.5*(c1x - c2x);
        seobdyn->chiAyVec[i] = 0.5*(c1y - c2y);
        seobdyn->chiAzVec[i] = 0.5*(c1z - c2z);

        seobdyn->chiSxVec[i] = 0.5*(c1x + c2x);
        seobdyn->chiSyVec[i] = 0.5*(c1y + c2y);
        seobdyn->chiSzVec[i] = 0.5*(c1z + c2z);
        // if(isnan(seobdyn->chiAyVec[i]) || isnan(seobdyn->chiSyVec[i]))
        //     print_debug("get r = %g, chi1 = (%g, %g, %g), chi2 = (%g, %g, %g)\n", polarr,
        //         c1x, c1y, c1z, c2x, c2y, c2z);

        REAL8 magL = sqrt(inner_product3d(rcrossp, rcrossp));
        for (j = 0; j < 3; j++)
            Lhat[j] = rcrossp[j] / magL;
        polarpphi = magL;

        /* Computing omega */
        omega = rcrossrdotNorm / (polarr * polarr);
        if (seobParams->hParams->flagZframe == FLAG_SEOBNRv4P_ZFRAME_L)
        {
            s1dotZ = inner_product3d(spin1vec, Lhat);
            s2dotZ = inner_product3d(spin2vec, Lhat);

        } else if (seobParams->hParams->flagZframe == FLAG_SEOBNRv4P_ZFRAME_LN)
        {
            s1dotZ = inner_product3d(spin1vec, LNhat);
            s2dotZ = inner_product3d(spin2vec, LNhat);
        }
        /* Compute Hamiltonian */
        // UINT SpinAlignedEOBversion = seobParams->seobCoeffs->SpinAlignedEOBversion;
        REAL8Vector cartPosVec, cartMomVec, s1Vec, s2Vec;
        cartPosVec.length = cartMomVec.length = s1Vec.length = s2Vec.length = 3;
        cartPosVec.data = rvec;
        cartMomVec.data = pvec;
        s1Vec.data = spin1vec; /* in units of mTotal^2 */
        s2Vec.data = spin2vec; /* in units of mTotal^2 */
        EOBCalculateSigmaStar(sigmaStar, m1, m2, &s1Vec, &s2Vec);
        EOBCalculateSigmaKerr(sigmaKerr, &s1Vec, &s2Vec);

        // Compute the augmented spin used in the Hamiltonian calibration
        // coefficients. See LIGO-T1900601-v1.

        REAL8 tempS1_p = inner_product3d(s1Vec.data, Lhat);
        REAL8 tempS2_p = inner_product3d(s2Vec.data, Lhat);
        REAL8 S1_perp[3] = {0, 0, 0};
        REAL8 S2_perp[3] = {0, 0, 0};
        for ( jj = 0; jj < 3; jj++) 
        {
            S1_perp[jj] = spin1vec[jj] - tempS1_p * Lhat[jj];
            S2_perp[jj] = spin2vec[jj] - tempS2_p * Lhat[jj];
        }

        REAL8 sKerr_norm =
            sqrt(inner_product3d(sigmaKerr->data, sigmaKerr->data));
        REAL8 S_con = 0.0;
        if (sKerr_norm > 1e-6) 
        {
            S_con = sigmaKerr->data[0] * Lhat[0] + sigmaKerr->data[1] * Lhat[1] +
                    sigmaKerr->data[2] * Lhat[2];
            S_con /= (1 - 2 * eta);
            S_con += (inner_product3d(S1_perp, sigmaKerr->data) +
                        inner_product3d(S2_perp, sigmaKerr->data)) /
                    sKerr_norm / (1 - 2 * eta) / 2.;
        }

        REAL8 a = sqrt(inner_product3d(sigmaKerr->data, sigmaKerr->data));
        status = XLALSimIMRCalculateSpinPrecEOBHCoeffs_v2(seobCoeffs, eta, a, S_con, 4, seobParams->hParams);
        if (status != CEV_SUCCESS) {failed = 1;goto QUIT;}
        ham = XLALSimIMRSpinPrecEOBHamiltonian(eta, &cartPosVec, &cartMomVec,
                                            &s1Vec, &s2Vec, sigmaKerr, sigmaStar,
                                            seobParams->tortoise, seobCoeffs, seobParams->hParams);
        /* Output values in seobdynamics */
        seobdyn->velVecx[i] = rdotvec[0];
        seobdyn->velVecy[i] = rdotvec[1];
        seobdyn->velVecz[i] = rdotvec[2];
        seobdyn->polarrVec[i] = polarr;
        seobdyn->polarphiVec[i] = polarphi;
        seobdyn->polarprVec[i] = polarpr;
        seobdyn->polarpphiVec[i] = polarpphi;
        seobdyn->omegaVec[i] = omega;
        seobdyn->s1dotZVec[i] = s1dotZ;
        seobdyn->s2dotZVec[i] = s2dotZ;
        seobdyn->hamVec[i] = ham;

        dr = (rvec[0]*rdotvec[0] + rvec[1]*rdotvec[1] + rvec[2]*rdotvec[2]) / polarr;
        ncrv = omega * polarr;

        memcpy(cartData, values->data, 12*sizeof(REAL8));
        polData[0] = polarr;
        polData[1] = 0;
        polData[2] = polarpr;
        polData[3] = polarpphi;
        REAL8 cFr, cFf;
        REAL8 c_prDot, c_pr;
#if 0
        flux = XLALInspiralPrecSpinFactorizedFlux(&polarDynamics, &cartDynamics,
                        NULL, omega, dr, ncrv, seobParams, ham, 8, 4);
        // print_debug("%.16e\n", flux);
        // // DEBUG: Fr and Ff
        REAL8 tmpP[3];
        TortoiseConvertMatrix(seobParams, values->data, dvalues->data, tmpP);
        c_pr = (tmpP[0] * rvec[0] + tmpP[1] * rvec[1] + tmpP[2] * rvec[2]) / polarr;
        c_prDot = -dr*c_pr/polarr + (rvec[0]*dvalues->data[3] + rvec[1]*dvalues->data[4] + rvec[2]*dvalues->data[5] + 
            rdotvec[0]*tmpP[0] + rdotvec[1]*tmpP[1] + rdotvec[2]*tmpP[2])/polarr - flux*c_pr/eta/polarpphi;
        // CalculateEccCorrectionToFlux(eta, c1z, c2z, polarr, dr, c_prDot, &cFr, &cFf);
        // CalculateEccCorrectionToFluxV2(eta, c1z, c2z, polarr, polarpr, c_prDot, &cFr, &cFf);
        CalculateEccCorrectionToFluxV3(eta, c1z, c2z, polarr, dr, c_prDot, &cFr, &cFf, seobParams->e0);
        // CalculateEccCorrectionToFluxV4(eta, c1z, c2z, polarr, polarpr, c_prDot, &cFr, &cFf);
        // CalculateRRForceSpinCoeffs(&coeffsFf, &coeffsFr, seobParams->m1, seobParams->m2, c1z, c2z);
        // psq = pvec[0]*pvec[0] + pvec[1]*pvec[1] + pvec[2]*pvec[2];
        // CalculateRRForce(&coeffsFf, &coeffsFr, &Ff, &Fr, polarr, polarpr, psq, polarpphi);
#else
        cFr = 0.0;
        cFf = 0.0;
        flux = 0.0;
        c_prDot = 0.0;
#endif
        seobdyn->FrVec[i] = dr;
        seobdyn->FfVec[i] = polarpr;
        seobdyn->polarprDotVec[i] = omega;
        seobdyn->fluxVec[i] = flux/eta;
        // print_debug("%e\t%e\t%e\t%e\n", cFr, cFf, c_prDot, polarr * c_prDot);
    }

QUIT:
    STRUCTFREE(values, REAL8Vector);
    STRUCTFREE(dvalues, REAL8Vector);
    STRUCTFREE(sigmaStar, REAL8Vector);
    STRUCTFREE(sigmaKerr, REAL8Vector);
    if (failed)
    {
        STRUCTFREE(seobdyn, SEOBdynamics);
        return CEV_FAILURE;
    }
    *seobdynamics = seobdyn;
    return CEV_SUCCESS;
}


INT SEOBComputeExtendedSEOBdynamics_Conserve(SEOBdynamics **seobdynamics,
                                    REAL8Array *dynamics,
                                    INT retLen,
                                    SpinEOBParams *seobParams)
{
    SEOBdynamics *seobdyn = NULL;
    INT status, failed = 0;
    UINT jj;
    seobdyn = CreateSEOBdynamics(retLen);
    // flagSEOBNRv4P_Zframe = 0
    // flagSEOBNRv4P_hamiltonian_derivative = 1
    /* Local variables */
    REAL8 rvec[3] = {0, 0, 0};
    REAL8 pvec[3] = {0, 0, 0};
    REAL8 spin1vec[3] = {0, 0, 0}; /* in units of mTotal^2 */
    REAL8 spin2vec[3] = {0, 0, 0}; /* in units of mTotal^2 */
    REAL8 rdotvec[3] = {0, 0, 0};
    REAL8 rcrossrdot[3] = {0, 0, 0};
    REAL8 rcrossp[3] = {0, 0, 0};
    REAL8 LNhat[3] = {0, 0, 0};
    REAL8 Lhat[3] = {0, 0, 0};
    REAL8 polarr, polarphi, polarpr, polarpphi, omega, s1dotZ, s2dotZ, ham;
    REAL8 nhat[3] = {0,0,0};
    /* Allocate temporary vectors values, dvalues */
    REAL8Vector *values = NULL;
    REAL8Vector *dvalues = NULL;
    values = CreateREAL8Vector(14);
    dvalues = CreateREAL8Vector(14);
    memset(values->data, 0, (values->length) * sizeof(REAL8));
    memset(dvalues->data, 0, (dvalues->length) * sizeof(REAL8));

    REAL8 m1 = seobParams->m1;
    REAL8 m2 = seobParams->m2;
    REAL8 mtot = m1 + m2;
    REAL8 m1sq = m1*m1/(m1+m2)/(m1+m2);
    REAL8 m2sq = m2*m2/(m1+m2)/(m1+m2);
    REAL8 eta = seobParams->eta;
    SpinEOBHCoeffs *seobCoeffs = seobParams->seobCoeffs;
    /* Copying directly the dynamics data in seobdynamics - 15 vectors of length
    * retLen */
    memcpy(seobdyn->array->data, dynamics->data, 15 * retLen * sizeof(REAL8));
    REAL8 c1x, c1y, c1z, c2x, c2y, c2z;
    REAL8 yhat[3] = {0,0,0};
    REAL8 magY;
    REAL8Vector *sigmaStar = NULL, *sigmaKerr = NULL;
    sigmaStar = CreateREAL8Vector(3);
    sigmaKerr = CreateREAL8Vector(3);
    memset(sigmaStar->data, 0, 3 * sizeof(REAL8));
    memset(sigmaKerr->data, 0, 3 * sizeof(REAL8));
    REAL8 flux, dr, ncrv;
	REAL8Vector	polarDynamics, cartDynamics;
	REAL8		polData[4];
    REAL8       cartData[12];
    polarDynamics.length = 4;
    polarDynamics.data = polData;
    cartDynamics.length = 12;
    cartDynamics.data = cartData;

    // DEBUG: RR force with Schott correction
    RRForceCoeffs coeffsFr;
    RRForceCoeffs coeffsFf;
    REAL8 Ff, Fr, psq;
    CalculateRRForceCoeffs(&coeffsFf, &coeffsFr, seobParams);
    /* Loop to compute the derived quantities from the dynamics */
    UINT i, j;
    for (i=0; i<retLen; i++)
    {
        for (j = 0; j < 14; j++) 
        {
            values->data[j] = dynamics->data[i + (j + 1) * retLen];
        }
        status = XLALSpinPrecHcapRvecDerivative(0, values->data, dvalues->data,
                                         (void *)seobParams);
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
        for (j=0; j<3; j++)
        {
            rvec[j] = values->data[j];
            pvec[j] = values->data[3 + j];
            spin1vec[j] = values->data[6 + j];
            spin2vec[j] = values->data[9 + j];
            rdotvec[j] = dvalues->data[j];
        }
        cross_product3d(rvec, pvec, rcrossp);
        cross_product3d(rvec, rdotvec, rcrossrdot);
        REAL8 rcrossrdotNorm = sqrt(inner_product3d(rcrossrdot, rcrossrdot));
        for (j = 0; j < 3; j++) 
            LNhat[j] = rcrossrdot[j] / rcrossrdotNorm;
        /* Polar dynamics */
        polarr = sqrt(inner_product3d(rvec, rvec));
        polarpr = inner_product3d(rvec, pvec) / polarr;
        polarphi = values->data[12] + values->data[13];
        // REAL8 chi1dotZ = s1dotZ * mtot * mtot / (m1 * m1);
        // REAL8 chi2dotZ = s2dotZ * mtot * mtot / (m2 * m2);
        // FIX ME
        cross_product3d(LNhat, rvec, yhat);
        magY = sqrt(yhat[0]*yhat[0] + yhat[1]*yhat[1] + yhat[2]*yhat[2]);
        // S1 = sqrt((spin1vec[0]*spin1vec[0] + spin1vec[1]*spin1vec[1] + spin1vec[2]*spin1vec[2]))/m1sq;
        // S2 = sqrt((spin2vec[0]*spin2vec[0] + spin2vec[1]*spin2vec[1] + spin2vec[2]*spin2vec[2]))/m2sq;

        c1x = inner_product3d(spin1vec, rvec) / polarr / m1sq;
        c1z = inner_product3d(spin1vec, LNhat) / m1sq;
        // c1y = sqrt( (spin1vec[0]*spin1vec[0] + spin1vec[1]*spin1vec[1] + spin1vec[2]*spin1vec[2])/(m1sq*m1sq) - c1x*c1x - c1z*c1z);
        c1y = inner_product3d(spin1vec, yhat) / m1sq / magY;

        c2x = inner_product3d(spin2vec, rvec) / polarr / m2sq;
        c2z = inner_product3d(spin2vec, LNhat) / m2sq;
        // c2y = sqrt( (spin2vec[0]*spin2vec[0] + spin2vec[1]*spin2vec[1] + spin2vec[2]*spin2vec[2])/(m2sq*m2sq) - c2x*c2x - c2z*c2z);
        c2y = inner_product3d(spin2vec, yhat) / m2sq / magY;

        seobdyn->chiAxVec[i] = 0.5*(c1x - c2x);
        seobdyn->chiAyVec[i] = 0.5*(c1y - c2y);
        seobdyn->chiAzVec[i] = 0.5*(c1z - c2z);

        seobdyn->chiSxVec[i] = 0.5*(c1x + c2x);
        seobdyn->chiSyVec[i] = 0.5*(c1y + c2y);
        seobdyn->chiSzVec[i] = 0.5*(c1z + c2z);
        // if(isnan(seobdyn->chiAyVec[i]) || isnan(seobdyn->chiSyVec[i]))
        //     print_debug("get r = %g, chi1 = (%g, %g, %g), chi2 = (%g, %g, %g)\n", polarr,
        //         c1x, c1y, c1z, c2x, c2y, c2z);

        REAL8 magL = sqrt(inner_product3d(rcrossp, rcrossp));
        for (j = 0; j < 3; j++)
            Lhat[j] = rcrossp[j] / magL;
        polarpphi = magL;

        /* Computing omega */
        omega = rcrossrdotNorm / (polarr * polarr);
        if (seobParams->hParams->flagZframe == FLAG_SEOBNRv4P_ZFRAME_L)
        {
            s1dotZ = inner_product3d(spin1vec, Lhat);
            s2dotZ = inner_product3d(spin2vec, Lhat);

        } else if (seobParams->hParams->flagZframe == FLAG_SEOBNRv4P_ZFRAME_LN)
        {
            s1dotZ = inner_product3d(spin1vec, LNhat);
            s2dotZ = inner_product3d(spin2vec, LNhat);
        }
        /* Compute Hamiltonian */
        // UINT SpinAlignedEOBversion = seobParams->seobCoeffs->SpinAlignedEOBversion;
        REAL8Vector cartPosVec, cartMomVec, s1Vec, s2Vec;
        cartPosVec.length = cartMomVec.length = s1Vec.length = s2Vec.length = 3;
        cartPosVec.data = rvec;
        cartMomVec.data = pvec;
        s1Vec.data = spin1vec; /* in units of mTotal^2 */
        s2Vec.data = spin2vec; /* in units of mTotal^2 */
        EOBCalculateSigmaStar(sigmaStar, m1, m2, &s1Vec, &s2Vec);
        EOBCalculateSigmaKerr(sigmaKerr, &s1Vec, &s2Vec);

        // Compute the augmented spin used in the Hamiltonian calibration
        // coefficients. See LIGO-T1900601-v1.

        REAL8 tempS1_p = inner_product3d(s1Vec.data, Lhat);
        REAL8 tempS2_p = inner_product3d(s2Vec.data, Lhat);
        REAL8 S1_perp[3] = {0, 0, 0};
        REAL8 S2_perp[3] = {0, 0, 0};
        for ( jj = 0; jj < 3; jj++) 
        {
            S1_perp[jj] = spin1vec[jj] - tempS1_p * Lhat[jj];
            S2_perp[jj] = spin2vec[jj] - tempS2_p * Lhat[jj];
        }

        REAL8 sKerr_norm =
            sqrt(inner_product3d(sigmaKerr->data, sigmaKerr->data));
        REAL8 S_con = 0.0;
        if (sKerr_norm > 1e-6) 
        {
            S_con = sigmaKerr->data[0] * Lhat[0] + sigmaKerr->data[1] * Lhat[1] +
                    sigmaKerr->data[2] * Lhat[2];
            S_con /= (1 - 2 * eta);
            S_con += (inner_product3d(S1_perp, sigmaKerr->data) +
                        inner_product3d(S2_perp, sigmaKerr->data)) /
                    sKerr_norm / (1 - 2 * eta) / 2.;
        }

        REAL8 a = sqrt(inner_product3d(sigmaKerr->data, sigmaKerr->data));
        status = XLALSimIMRCalculateSpinPrecEOBHCoeffs_v2(seobCoeffs, eta, a, S_con, 4, seobParams->hParams);
        if (status != CEV_SUCCESS) {failed = 1;goto QUIT;}
        ham = XLALSimIMRSpinPrecEOBHamiltonian(eta, &cartPosVec, &cartMomVec,
                                            &s1Vec, &s2Vec, sigmaKerr, sigmaStar,
                                            seobParams->tortoise, seobCoeffs, seobParams->hParams);
        /* Output values in seobdynamics */
        seobdyn->velVecx[i] = rdotvec[0];
        seobdyn->velVecy[i] = rdotvec[1];
        seobdyn->velVecz[i] = rdotvec[2];
        seobdyn->polarrVec[i] = polarr;
        seobdyn->polarphiVec[i] = polarphi;
        seobdyn->polarprVec[i] = polarpr;
        seobdyn->polarpphiVec[i] = polarpphi;
        seobdyn->omegaVec[i] = omega;
        seobdyn->s1dotZVec[i] = s1dotZ;
        seobdyn->s2dotZVec[i] = s2dotZ;
        seobdyn->hamVec[i] = ham;

        dr = (rvec[0]*rdotvec[0] + rvec[1]*rdotvec[1] + rvec[2]*rdotvec[2]) / polarr;
        ncrv = omega * polarr;

        memcpy(cartData, values->data, 12*sizeof(REAL8));
        polData[0] = polarr;
        polData[1] = 0;
        polData[2] = polarpr;
        polData[3] = polarpphi;
        REAL8 cFr, cFf;
        REAL8 c_prDot, c_pr;
#if 1
        flux = XLALInspiralPrecSpinFactorizedFlux(&polarDynamics, &cartDynamics,
                        NULL, omega, dr, ncrv, seobParams, ham, 8, 4);
        // print_debug("%.16e\n", flux);
        // // DEBUG: Fr and Ff
        REAL8 tmpP[3];
        TortoiseConvertMatrix(seobParams, values->data, dvalues->data, tmpP);
        c_pr = (tmpP[0] * rvec[0] + tmpP[1] * rvec[1] + tmpP[2] * rvec[2]) / polarr;
        c_prDot = -dr*c_pr/polarr + (rvec[0]*dvalues->data[3] + rvec[1]*dvalues->data[4] + rvec[2]*dvalues->data[5] + 
            rdotvec[0]*tmpP[0] + rdotvec[1]*tmpP[1] + rdotvec[2]*tmpP[2])/polarr - flux*c_pr/eta/polarpphi;
        // CalculateEccCorrectionToFlux(eta, c1z, c2z, polarr, dr, c_prDot, &cFr, &cFf);
        // CalculateEccCorrectionToFluxV2(eta, c1z, c2z, polarr, polarpr, c_prDot, &cFr, &cFf);
        CalculateEccCorrectionToFluxV3(eta, c1z, c2z, polarr, dr, c_prDot, &cFr, &cFf, seobParams->e0);
        // prDot = dvalues[2] - ( values[2] / values[3] / csi ) * flux / omega;
        // CalculateEccCorrectionToFluxV3(eta, params.params->chi1, params.params->chi1, r, dvalues[0], prDot, &cFr, &cFf, params.params->e0);
        // CalculateEccCorrectionToFluxV4(eta, c1z, c2z, polarr, polarpr, c_prDot, &cFr, &cFf);
        // CalculateRRForceSpinCoeffs(&coeffsFf, &coeffsFr, seobParams->m1, seobParams->m2, c1z, c2z);
        // psq = pvec[0]*pvec[0] + pvec[1]*pvec[1] + pvec[2]*pvec[2];
        // CalculateRRForce(&coeffsFf, &coeffsFr, &Ff, &Fr, polarr, polarpr, psq, polarpphi);
#else
        cFr = 0.0;
        cFf = 0.0;
        flux = 0.0;
        c_prDot = 0.0;
#endif
        seobdyn->FrVec[i] = cFr;
        seobdyn->FfVec[i] = cFf;
        seobdyn->polarprDotVec[i] = omega;
        seobdyn->fluxVec[i] = flux/eta;
        // print_debug("%e\t%e\t%e\t%e\n", cFr, cFf, c_prDot, polarr * c_prDot);
    }
QUIT:
    STRUCTFREE(values, REAL8Vector);
    STRUCTFREE(dvalues, REAL8Vector);
    STRUCTFREE(sigmaStar, REAL8Vector);
    STRUCTFREE(sigmaKerr, REAL8Vector);
    if (failed)
    {
        STRUCTFREE(seobdyn, SEOBdynamics);
        return CEV_FAILURE;
    }
    *seobdynamics = seobdyn;
    return CEV_SUCCESS;
}


/**
 * This function computes all extended dynamics values at a given time by
 * interpolating the dynamics array. We build a cubic spline limited to +- 20
 * samples on each side of the time of interest.
 */
int SEOBInterpolateDynamicsAtTime(
    REAL8Vector **seobdynamics_values, /**<< Output: pointer to vector for
                                          seobdynamics interpolated values */
    REAL8 t,                           /**<< Input: time at which to evaluate */
    SEOBdynamics *seobdynamics         /**<< Input: SEOB dynamics */
) 
{
    UINT j;
    /* Create output vector */
    if (!((*seobdynamics_values) = CreateREAL8Vector(v4PdynamicsVariables))) 
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

    /* Spline allocation */
    GSL_START;
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, interp_length);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    int status, is_failed = 0;
    /* Interpolate all quantities */
    (*seobdynamics_values)->data[0] = t;
    for ( j = 1; j < v4PdynamicsVariables; j++) {
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

int SEOBInterpolateSADynamicsAtTime(
    REAL8Vector **seobdynamics_values, /**<< Output: pointer to vector for
                                          seobdynamics interpolated values */
    REAL8 t,                           /**<< Input: time at which to evaluate */
    SEOBSAdynamics *seobdynamics         /**<< Input: SEOB dynamics */
)
{
    // turn off gsl error handler
    int is_failed = 0;
    UINT j;
    /* Create output vector */
    if (!((*seobdynamics_values) = CreateREAL8Vector(5))) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failed to allocate REAL8Vector ");
        return CEV_FAILURE;
    }
    memset((*seobdynamics_values)->data, 0, (*seobdynamics_values)->length * sizeof(REAL8));

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
    /* Spline allocation */
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, interp_length);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    GSL_START;
    /* Interpolate all quantities */
    int status;
    (*seobdynamics_values)->data[0] = t;
    for ( j = 1; j < 5; j++) {
        status = gsl_spline_init(spline, &(tVec[indexstart]),
                        &(seobdynamics->array->data[j * retLen + indexstart]),
                        interp_length);
        if (status != GSL_SUCCESS)
        {
            is_failed = 1;
            goto QUIT;
        }
        (*seobdynamics_values)->data[j] = gsl_spline_eval(spline, t, acc);
        // status = gsl_spline_eval_e(spline, t, acc, &((*seobdynamics_values)->data[j]));
        // print_debug("status2 = %d\n", status);
        // if (status != GSL_SUCCESS)
        // {
        //     gsl_spline_free(spline);
        //     gsl_interp_accel_free(acc);
        //     return CEV_FAILURE;
        // }
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

int SEOBLocateTimePeakOmega(
    REAL8 *tPeakOmega, /**<< Output: time of peak of Omega if found (see inside
        XLALSimLocateOmegaTime for what is returned otherwise) */
    INT *foundPeakOmega, /**<< Output: flag indicating wether tPeakOmega has been found */
    REAL8Array *dynamics,      /**<< Input: array for dynamics */
    SEOBdynamics *seobdynamics,       /**<< Input: SEOB dynamics object */
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

int SEOBSALocateTimePeakOmega(
    REAL8 *tPeakOmega, /**<< Output: time of peak of Omega if found (see inside
        XLALSimLocateOmegaTime for what is returned otherwise) */
    INT *foundPeakOmega, /**<< Output: flag indicating wether tPeakOmega has been found */
    SEOBSAdynamics *seobdynamics,       /**<< Input: SEOB dynamics object */
    UINT retLen,                     /**<< Input: length of dynamics */
    SpinEOBParams *seobParams /**<< SEOB params */
)
{ 
    REAL8Vector tVec;
    tVec.length = retLen;
    tVec.data = seobdynamics->tVec;
    REAL8Vector omegaVec;
    omegaVec.length = retLen;
    omegaVec.data = seobdynamics->dphiVec;
    XLALEOBFindRobustPeak(tPeakOmega, &tVec, &omegaVec, 3);
    *foundPeakOmega = 1;
    // {
    //     int j;
    //     REAL8 *vector;
    //     for (int i=0; i<seobdynamics->length; i++)
    //     {
    //         print_out("%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n", 
    //             seobdynamics->tVec[i], 
    //             seobdynamics->rVec[i], seobdynamics->phiVec[i],
    //             seobdynamics->prTVec[i], seobdynamics->pphiVec[i],
    //             seobdynamics->drVec[i], seobdynamics->dphiVec[i],
    //             seobdynamics->dprTVec[i], seobdynamics->dpphiVec[i],
    //             seobdynamics->HVec[i]);
    //     }
    // }

    return CEV_SUCCESS;
}


int SEOBLFrameVectors(
    REAL8Vector **S1,        /**<<Output: S1 in L-n frame */
    REAL8Vector **S2,        /**<<Output: S2 in L-n frame */
    REAL8Vector *seobvalues, /**<<Input: vector of extended dynamics */
    REAL8 m1, /**<<Input: mass of the first object in solar masses */
    REAL8 m2, /**<<Input: mass of the second object in solar masses */
    INT flagZframe
) 
{

    if ((!S1) || (!S2) || (!seobvalues)) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Passed null pointers to SEOBLFrameVectors!");
        return CEV_FAILURE;
    }
    REAL8 mTotal = m1 + m2;
    // Scaled masses, useful to compute spins in sane units
    REAL8 m1_sc = m1 / mTotal;
    REAL8 m2_sc = m2 / mTotal;
    /* Create output vectors */
    if (!((*S1) = CreateREAL8Vector(3))) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failed to allocate REAL8Vector S1.");
        return CEV_FAILURE;
    }
    if (!((*S2) = CreateREAL8Vector(3))) {
        STRUCTFREE(*S1, REAL8Vector); // Free the memory above
        PRINT_LOG_INFO(LOG_CRITICAL, "failed to allocate REAL8Vector S2.");
        return CEV_FAILURE;
    }
    memset((*S1)->data, 0, 3 * sizeof(REAL8));
    memset((*S2)->data, 0, 3 * sizeof(REAL8));

    /* Local variables */
    REAL8 rvec[3] = {0, 0, 0};
    REAL8 drvec[3] = {0, 0, 0};
    REAL8 pvec[3] = {0, 0, 0};
    REAL8 spin1vec[3] = {0, 0, 0}; /* in units of mTotal^2 */
    REAL8 spin2vec[3] = {0, 0, 0}; /* in units of mTotal^2 */
    REAL8 crossp[3] = {0, 0, 0};
    REAL8 L_hat[3] = {0, 0, 0};
    REAL8 n_hat[3] = {0, 0, 0};
    REAL8 lambda_hat[3] = {0, 0, 0};
    int j, jj;
    /* Read from the extended dynamics values */
    for ( j = 0; j < 3; j++) 
    {
        rvec[j] = seobvalues->data[1 + j];
        drvec[j] = seobvalues->data[15 + j];
        pvec[j] = seobvalues->data[4 + j];
        spin1vec[j] = seobvalues->data[7 + j] / m1_sc / m1_sc;
        spin2vec[j] = seobvalues->data[10 + j] / m2_sc / m2_sc;
    }
    if (flagZframe == FLAG_SEOBNRv4P_ZFRAME_L) 
    {
        // Note: pvec is missing a factor of nu, but don't care about it's magnitide
        // anyway
        cross_product3d(rvec, pvec, crossp);
    } else if (flagZframe == FLAG_SEOBNRv4P_ZFRAME_LN) 
    {
        cross_product3d(rvec, drvec, crossp);
    }
    REAL8 Lmag = sqrt(inner_product3d(crossp, crossp));
    REAL8 sep = sqrt(inner_product3d(rvec, rvec));

    for ( jj = 0; jj < 3; jj++) 
    {
        L_hat[jj] = crossp[jj] / Lmag;
        n_hat[jj] = rvec[jj] / sep;
    }

    cross_product3d(L_hat, n_hat, lambda_hat);
    // Project onto the new frame
    (*S1)->data[0] = inner_product3d(spin1vec, n_hat);
    (*S2)->data[0] = inner_product3d(spin2vec, n_hat);

    (*S1)->data[1] = inner_product3d(spin1vec, lambda_hat);
    (*S2)->data[1] = inner_product3d(spin2vec, lambda_hat);

    (*S1)->data[2] = inner_product3d(spin1vec, L_hat);
    (*S2)->data[2] = inner_product3d(spin2vec, L_hat);
    return CEV_SUCCESS;
}


int SEOBJfromDynamics(
    REAL8Vector **J,         /**<< Output: pointer to vector J */
    REAL8Vector *seobvalues, /**<< Input: vector for extended dynamics values */
    SpinEOBParams *seobParams /**<< SEOB params */
) 
{
    UINT j;
    if ((!J) || (!seobvalues) || (!seobParams)) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Some pointers passed to SEOBJfromDynamics were null");
        return CEV_FAILURE;
    }
    /* Create output vector */
    if (!((*J) = CreateREAL8Vector(3))) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "XLAL Error failed to allocate REAL8Vector J");
        return CEV_FAILURE;
    }
    memset((*J)->data, 0, 3 * sizeof(REAL8));

    /* Masses */
    REAL8 eta = seobParams->eta;

    /* Local variables */
    REAL8 rvec[3] = {0, 0, 0};
    REAL8 pvec[3] = {0, 0, 0};
    REAL8 spin1vec[3] = {0, 0, 0}; /* in units of mTotal^2 */
    REAL8 spin2vec[3] = {0, 0, 0}; /* in units of mTotal^2 */
    REAL8 rcrossp[3] = {0, 0, 0};

    /* Read from the extended dynamics values */
    for (j = 0; j < 3; j++) 
    {
        rvec[j] = seobvalues->data[1 + j];
        pvec[j] = seobvalues->data[4 + j];
        spin1vec[j] = seobvalues->data[7 + j];
        spin2vec[j] = seobvalues->data[10 + j];
    }
    cross_product3d(rvec, pvec, rcrossp);

    /* Compute J - restoring the factor eta in L (p stored in the dynamics is
    * p/mu) */
    for (j = 0; j < 3; j++) 
    {
        (*J)->data[j] = eta * rcrossp[j] + spin1vec[j] + spin2vec[j];
    }

    return CEV_SUCCESS;
}


int SEOBLhatfromDynamics(
    REAL8Vector **L,         /**<< Output: pointer to vector L */
    REAL8Vector *seobvalues, /**<< Input: vector for extended dynamics values */
    SpinEOBParams *seobParams /**<< SEOB params */
)
{
    if ((!L) || (!seobvalues) || (!seobParams))
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Some pointers passed to SEOBLfromDynamics were null");
        return CEV_FAILURE;
    }
    /* Create output vector */
    if (!((*L) = CreateREAL8Vector(3)))
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "XLAL Error failed to allocate REAL8Vector L.");
        return CEV_FAILURE;
    }
    memset((*L)->data, 0, 3 * sizeof(REAL8));

    /* Local variables */
    REAL8 rvec[3] = {0, 0, 0};
    REAL8 drvec[3] = {0, 0, 0};
    REAL8 pvec[3] = {0, 0, 0};
    REAL8 crossp[3] = {0, 0, 0};
    int j, jj;
    /* Read from the extended dynamics values */
    for ( j = 0; j < 3; j++)
    {
        rvec[j] = seobvalues->data[1 + j];
        drvec[j] = seobvalues->data[15 + j];
        pvec[j] = seobvalues->data[4 + j];
    }
    if (seobParams->hParams->flagZframe == FLAG_SEOBNRv4P_ZFRAME_L)
    {
        // Note: pvec is missing a factor of nu, but don't care about it's magnitide
        // anyway
        cross_product3d(rvec, pvec, crossp);
    }
    else if (seobParams->hParams->flagZframe == FLAG_SEOBNRv4P_ZFRAME_LN)
    {
        cross_product3d(rvec, drvec, crossp);
    }
    REAL8 Lmag = sqrt(inner_product3d(crossp, crossp));

    for ( jj = 0; jj < 3; jj++)
    {
        (*L)->data[jj] = crossp[jj] / Lmag;
    }

    return CEV_SUCCESS;
}

int SEOBBuildJframeVectors(
    REAL8Vector *e1J, /**<< Output: vector for e1J, already allocated */
    REAL8Vector *e2J, /**<< Output: vector for e2J, already allocated */
    REAL8Vector *e3J, /**<< Output: vector for e3J, already allocated */
    REAL8Vector *JVec /**<< Input: vector J */
) 
{
    UINT j;

    /* Checking size and of input vectors */
    if ((!e1J) || (!e2J) || (!e3J) || (!JVec)) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "at least one input pointer is NULL.");
        return CEV_FAILURE;
    }
    if ((!(e1J->length == 3)) || (!(e2J->length == 3)) || (!(e2J->length == 3))) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "at least one input vector is not of length 3.");
        return CEV_FAILURE;
    }

    /* Set e3J to Jhat */
    REAL8 Jnorm = sqrt(inner_product3d(JVec->data, JVec->data));
    for (j = 0; j < 3; j++) 
    {
        e3J->data[j] = JVec->data[j] / Jnorm;
    }

    /* Set e1J to enforce the condition that x is in the plane (e1J, e3J) and
    * x.e1J>0 */
    /* Added a protection against the degenerate case where e3J, x are aligned:
    * let lambda = 1 - |ex.e3J|, measuring the alignment of e3J, x
    * for lambda < 1e-5 use y instead of x
    * for lambda > 1e-4 use x normally
    * for lambda in [1e-4, 1e-5] use a lambda-weighted combination of both
    * thresholds are arbitrary */
    REAL8 normfacx = 0.;
    REAL8 normfacy = 0.;
    REAL8 weightx = 0.;
    REAL8 weighty = 0.;
    REAL8 e1Jblendednorm = 0.;
    REAL8 exvec[3] = {1, 0, 0};
    REAL8 eyvec[3] = {0, 1, 0};
    REAL8 exdote3J = inner_product3d(exvec, e3J->data);
    REAL8 eydote3J = inner_product3d(eyvec, e3J->data);
    REAL8 lambda = 1. - fabs(exdote3J);
    if ((lambda < 0.) || (lambda > 1.)) {
        PRINT_LOG_INFO(LOG_CRITICAL, "Problem: lambda=1-|e3J.ex|=%g, should be in [0,1]", lambda);
        return CEV_FAILURE;
    }
    if (lambda > 1e-4) 
    {
        normfacx = 1. / sqrt(1. - exdote3J * exdote3J);
        for (j = 0; j < 3; j++)
            e1J->data[j] = (exvec[j] - exdote3J * e3J->data[j]) / normfacx;
    } 
    else if (lambda < 1e-5) 
    {
        normfacy = 1. / sqrt(1. - eydote3J * eydote3J);
        for (j = 0; j < 3; j++)
            e1J->data[j] = (eyvec[j] - eydote3J * e3J->data[j]) / normfacy;
    } 
    else 
    {
        weightx = (lambda - 1e-5) / (1e-4 - 1e-5);
        weighty = 1. - weightx;
        normfacx = 1. / sqrt(1. - exdote3J * exdote3J);
        normfacy = 1. / sqrt(1. - eydote3J * eydote3J);
        for (j = 0; j < 3; j++)
            e1J->data[j] = weightx * (exvec[j] - exdote3J * e3J->data[j]) / normfacx +
                            weighty * (eyvec[j] - eydote3J * e3J->data[j]) / normfacy;
        e1Jblendednorm = sqrt(inner_product3d(e1J->data, e1J->data));
        for (j = 0; j < 3; j++)
            e1J->data[j] /= e1Jblendednorm;
    }

    /* Get e2J = e3J * e1J */
    cross_product3d(e3J->data, e1J->data, e2J->data);

    /* Normally, vectors already of unit norm - we normalize again to eliminate
    * possible round-off error */
    REAL8 e1Jnorm = sqrt(inner_product3d(e1J->data, e1J->data));
    REAL8 e2Jnorm = sqrt(inner_product3d(e2J->data, e2J->data));
    REAL8 e3Jnorm = sqrt(inner_product3d(e3J->data, e3J->data));
    for (j = 0; j < 3; j++) 
    {
        e1J->data[j] /= e1Jnorm;
        e2J->data[j] /= e2Jnorm;
        e3J->data[j] /= e3Jnorm;
    }

    return CEV_SUCCESS;
}

int SEOBEulerI2JFromJframeVectors(
    REAL8 *alphaI2J,  /**<< Output: Euler angle alpha I2J */
    REAL8 *betaI2J,   /**<< Output: Euler angle beta I2J */
    REAL8 *gammaI2J,  /**<< Output: Euler angle gamma I2J */
    REAL8Vector *e1J, /**<< Input: unit Jframe vector e1J */
    REAL8Vector *e2J, /**<< Input: unit Jframe vector e2J */
    REAL8Vector *e3J  /**<< Input: unit Jframe vector e3J */
) 
{
    /* Active rotation matrix from frame (x,y,z) to frame (e1J,e2J,e3J) */
    /* The input vectors (eJ) are decomposed on the basis (x,y,z) */
    REAL8Array *R = CreateREAL8Array(2, 3, 3);
    if (!R) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Allocating the rotation matrix failed!");
        return CEV_FAILURE;
    }
    RotationMatrixActiveFromBasisVectors(R, e1J->data, e2J->data, e3J->data);

    /* Compute Euler angles in the Z-Y-Z convention */
    EulerAnglesZYZFromRotationMatrixActive(alphaI2J, betaI2J, gammaI2J, R);

    /* Cleanup */
    STRUCTFREE(R, REAL8Array);
    return CEV_SUCCESS;
}

int SEOBCalculateSphHarmListNQCCoefficientsV4(
    SphHarmListEOBNonQCCoeffs *
        *nqcCoeffsList, /**<< Output: non-quasi-circular coefficients as a list
                           for each mode */
    INT modes[][2],    /**<< Input: array of modes (l,m) */
    UINT nmodes,       /**<< Input: number of modes (l,m) */
    REAL8 tPeakOmega,   /**<< Input: time of peak of Omega */
    SEOBdynamics *seobdynamics,  /**<< Input: SEOB dynamics */
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
    pr.data = seobdynamics->polarprVec;
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
        SEOBCalculatehlmAmpPhase(&hlm, l, m, seobdynamics, nqcCoeffs, seobParams,
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
#if 0
        CHAR fname[256];
        sprintf(fname, "debug_hlmNQC_new_%d%d.dat", l, m);
        FILE *out = fopen( fname ,"w");
        for (INT ii=0; ii<hlm->xdata->length; ii++)
        {
            fprintf(out, "%.16e\t%.16e\t%.16e\n",
            deltaT*ii- tPeakOmegaFromStartDyn,
            hlm_amp->data[ii], hlm_phase->data[ii]);
        }
        fclose(out);
#endif
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

int SEOBSACalculateSphHarmListNQCCoefficientsV4(
    SphHarmListEOBNonQCCoeffs *
        *nqcCoeffsList, /**<< Output: non-quasi-circular coefficients as a list
                           for each mode */
    INT modes[][2],    /**<< Input: array of modes (l,m) */
    UINT nmodes,       /**<< Input: number of modes (l,m) */
    REAL8 tPeakOmega,   /**<< Input: time of peak of Omega */
    SEOBSAdynamics *seobdynamics,  /**<< Input: SEOB dynamics */
    SpinEOBParams *seobParams
)
{
    // PRINT_LOG_INFO(LOG_INFO, "Calculate NQC Coeffs");
    INT failed = 0;
    /* Masses */
    REAL8 m1 = seobParams->m1;
    REAL8 m2 = seobParams->m2;
    seobParams->tPeakOmega = tPeakOmega;
    seobParams->spin1z_omegaPeak = seobParams->chi1;
    seobParams->spin2z_omegaPeak = seobParams->chi2;
    // REAL8 mtot = m1 + m2;
    // UINT4 SpinAlignedEOBversion =
    // seobParams->seobCoeffs->SpinAlignedEOBversion;

    /* Length of dynamics data and sampling step */
    UINT retLen = seobdynamics->length;
    REAL8 deltaT = seobdynamics->tVec[1] - seobdynamics->tVec[0];

    /* Create vectors from dynamics */
    REAL8Vector r, pr, omega;
    r.length = pr.length = omega.length = retLen;
    r.data = seobdynamics->rVec;
    pr.data = seobdynamics->prTVec;
    omega.data = seobdynamics->dphiVec;

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

    REAL8 chi1dotZfinal = seobParams->chi1;
    REAL8 chi2dotZfinal = seobParams->chi2;
    REAL8 chiSfinal = SEOBCalculateChiS(chi1dotZfinal, chi2dotZfinal);
    REAL8 chiAfinal = SEOBCalculateChiA(chi1dotZfinal, chi2dotZfinal);
    REAL8 q = m1/m2;
    //printf("chiA = %.16f\n",chiAfinal);
    /* Time elapsed from the start of the dynamics to tPeakOmega */
    REAL8 tPeakOmegaFromStartDyn = tPeakOmega - seobdynamics->tVec[0];
PRINT_LOG_INFO(LOG_DEBUG,"tPeak = %.16e, t0 = %.16e, tend = %.16e, tPeakOmegaFromStartDyn = %.16e", 
        tPeakOmega, seobdynamics->tVec[0], seobdynamics->tVec[retLen-1]-seobdynamics->tVec[0], tPeakOmegaFromStartDyn);
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
        if (q<1.0005 && (m % 2 != 0) && (fabs(chiAfinal) < 0.15)) 
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
            // SEOBCalculatehlmAmpPhase(&hlm, l, m, seobdynamics, nqcCoeffs, seobParams,
            //                         includeNQC);
            // print_debug("(%d,%d)calc SA ampphase\n", l, m);
            // if (l==2 && m==1)
            //     DEBUG_START;
            SEOBSACalculatehlmAmpPhase(&hlm, l, m, seobdynamics, nqcCoeffs, seobParams,
                                    includeNQC);
            // if (l==2 && m==1)
            //     DEBUG_END;
            // print_debug("(%d,%d)calc SA ampphase...done\n", l, m);
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
            // if (l==2 && m==2)
            //     DEBUG_START;
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
            // if (l==2 && m==2)
            //     DEBUG_END;
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

int SEOBCalculateSphHarmListhlmAmpPhase(
    SphHarmListCAmpPhaseSequence **listhlm,               /**<< Output: list of modes for hlm */
    INT modes[][2],            /**<< Input: array of modes (l,m) */
    UINT nmodes,               /**<< Input: number of modes (l,m) */
    SEOBdynamics *seobdynamics, /**<< Input: SEOB dynamics */
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

        EOBNonQCCoeffs *nqcCoeffs =
            SphHarmListEOBNonQCCoeffs_GetMode(listnqcCoeffs, l, m)->nqcCoeffs;

        CAmpPhaseSequence *hlm = NULL;
        SEOBCalculatehlmAmpPhase(&hlm, l, m, seobdynamics, nqcCoeffs, seobParams,
                                    includeNQC);

        SphHarmListCAmpPhaseSequence_AddMode(listhlm, hlm, l, m);
    }

    return CEV_SUCCESS;
}

int SEOBSACalculateSphHarmListhlmAmpPhase(
    SphHarmListCAmpPhaseSequence **listhlm,               /**<< Output: list of modes for hlm */
    INT modes[][2],            /**<< Input: array of modes (l,m) */
    UINT nmodes,               /**<< Input: number of modes (l,m) */
    SEOBSAdynamics *seobdynamics, /**<< Input: SEOB dynamics */
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

        EOBNonQCCoeffs *nqcCoeffs =
            SphHarmListEOBNonQCCoeffs_GetMode(listnqcCoeffs, l, m)->nqcCoeffs;

        CAmpPhaseSequence *hlm = NULL;
        SEOBSACalculatehlmAmpPhase(&hlm, l, m, seobdynamics, nqcCoeffs, seobParams,
                                    includeNQC);

        SphHarmListCAmpPhaseSequence_AddMode(listhlm, hlm, l, m);
    }

    return CEV_SUCCESS;
}

int SEOBCalculateSphHarmListhlmAmpPhase_noNQC(
    SphHarmListCAmpPhaseSequence **listhlm,               /**<< Output: list of modes for hlm */
    INT modes[][2],            /**<< Input: array of modes (l,m) */
    UINT nmodes,               /**<< Input: number of modes (l,m) */
    SEOBdynamics *seobdynamics, /**<< Input: SEOB dynamics */
    SpinEOBParams *seobParams                /**<< SEOB params */
) 
{
    /* Read version of SEOB to be used */
    // UINT SpinAlignedEOBversion = seobParams->seobCoeffs->SpinAlignedEOBversion;
    UINT SpinAlignedEOBversion = 4;

    /* Flag for inclusion of NQC, useful for higher modes that have no NQC
    * implemented for some SpinAlignedEOBversion */
    UINT includeNQC = 0, nmode;
    EOBNonQCCoeffs *nqcCoeffs = NULL;

    /* Loop over modes */
    for ( nmode = 0; nmode < nmodes; nmode++) 
    {

        INT l = modes[nmode][0];
        INT m = modes[nmode][1];

        CAmpPhaseSequence *hlm = NULL;
        SEOBCalculatehlmAmpPhase_noNQC(&hlm, l, m, seobdynamics, seobParams);

        SphHarmListCAmpPhaseSequence_AddMode(listhlm, hlm, l, m);
    }

    return CEV_SUCCESS;
}

int SEOBGetFinalSpinMass(
    REAL8
        *finalMass, /**<< Output: final mass computed from fit (scaled by M) */
    REAL8 *
        finalSpin, /**<< Output: final spin computed from fit (dimensionless) */
    REAL8Vector *seobvalues, /**<< Input: vector for dynamics values at time of
                                peak of omega */
    SpinEOBParams *seobParams /**<< Input: SEOB params */
) 
{
    int failed = 0;
    /* Masses */
    REAL8 m1 = seobParams->m1;
    REAL8 m2 = seobParams->m2;

    // Approximant seobApproximant = seobParams->seobApproximant;

    /* Compute components of vectors chi1, chi2 from seobvalues at tPeakOmega in
    * the L-frame */

    REAL8Vector *chi1temp = NULL;
    REAL8Vector *chi2temp = NULL;
    if (SEOBLFrameVectors(&chi1temp, &chi2temp, seobvalues, m1, m2, seobParams->hParams->flagZframe) ==
        CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in SEOBLFrameVectors.");
        failed = 1;
        goto QUIT;
    }

    REAL8 chi1L[3] = {chi1temp->data[0], chi1temp->data[1], chi1temp->data[2]};
    REAL8 chi2L[3] = {chi2temp->data[0], chi2temp->data[1], chi2temp->data[2]};
    if (XLALSimIMREOBFinalMassSpinPrec(finalMass, finalSpin, m1, m2, chi1L, chi2L) == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBFinalMassSpinPrec.");
        failed = 1;
        goto QUIT;
    }
// print_debug("chi1L = %.16e, chi2L = %.16e, m1 = %.16e, m2 = %.16e, finalmass = %.16e, finalspin = %.16e\n", 
//     chi1L[2], chi2L[2], m1, m2, *finalMass, *finalSpin);

QUIT:
    STRUCTFREE(chi1temp, REAL8Vector);
    STRUCTFREE(chi2temp, REAL8Vector);
    if (failed)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}

int SEOBSAGetFinalSpinMass(
    REAL8
        *finalMass, /**<< Output: final mass computed from fit (scaled by M) */
    REAL8 *
        finalSpin, /**<< Output: final spin computed from fit (dimensionless) */
    SpinEOBParams *seobParams /**<< Input: SEOB params */
) 
{
    int failed = 0;
    /* Masses */
    REAL8 m1 = seobParams->m1;
    REAL8 m2 = seobParams->m2;

    // Approximant seobApproximant = seobParams->seobApproximant;

    /* Compute components of vectors chi1, chi2 from seobvalues at tPeakOmega in
    * the L-frame */

    // REAL8Vector *chi1temp = NULL;
    // REAL8Vector *chi2temp = NULL;
    // if (SEOBLFrameVectors(&chi1temp, &chi2temp, seobvalues, m1, m2, seobParams->hParams->flagZframe) ==
    //     CEV_FAILURE) 
    // {
    //     PRINT_LOG_INFO(LOG_CRITICAL, "failure in SEOBLFrameVectors.");
    //     failed = 1;
    //     goto QUIT;
    // }

    REAL8 chi1L[3] = {0, 0, seobParams->chi1};
    REAL8 chi2L[3] = {0, 0, seobParams->chi2};
    if (XLALSimIMREOBFinalMassSpinPrec(finalMass, finalSpin, m1, m2, chi1L, chi2L) == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBFinalMassSpinPrec.");
        failed = 1;
        goto QUIT;
    }
// print_debug("chi1L = %.16e, chi2L = %.16e, m1 = %.16e, m2 = %.16e, finalmass = %.16e, finalspin = %.16e\n", 
//     chi1L[2], chi2L[2], m1, m2, *finalMass, *finalSpin);

QUIT:
    if (failed)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}



INT SEOBAttachRDToSphHarmListhPlm(
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
    SEOBdynamics *seobdynamics,            /**<< Input: SEOB dynamics */
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

INT SEOBSAAttachRDToSphHarmListhPlm(
    SphHarmListCAmpPhaseSequence **listhPlm_RDattached, /**<< Output: list of extended modes hlm with RD
                                 attached */
    COMPLEX16Vector **sigmaQNMlm0, /**<< Output: list of QNM complex frequency for modes lm,
                         0th overtone (dimensionless) */
    INT4 modes[][2],  /**<< Input: array of modes (l,m) */
    UINT nmodes,     /**<< Input: number of modes (l,m) */
    REAL8 finalMass,  /**<< Input: final mass computed from fit (scaled by M) */
    REAL8 finalSpin, /**<< Input: final spin computed from fit (dimensionless) */
    SphHarmListCAmpPhaseSequence *listhPlm, /**<< Input: list of modes hlm */
    REAL8 deltaT,                           /**<< Input: time step */
    UINT retLen,        /**<< Input: length of the input modes and dynamics */
    UINT retLenRDPatch, /**<< Input: length of the ringdown patch */
    REAL8 tAttach,       /**<< Input: time of attachment */
    // REAL8 tStart, /**<< Input: starting time (of the HiS part) */
    REAL8Vector *seobvalues, /**<< Input: vector for dynamics values at time of
                                peak of omega */
    SEOBSAdynamics *seobdynamics,            /**<< Input: SEOB dynamics */
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

    // REAL8Vector *chi1temp = NULL;
    // REAL8Vector *chi2temp = NULL;
    REAL8Vector *hPlmRe = NULL;
    REAL8Vector *hPlmIm = NULL;
    // SEOBLFrameVectors(&chi1temp, &chi2temp, seobvalues, m1, m2, seobParams->hParams->flagZframe);

    REAL8 chi1Lx = 0;
    REAL8 chi1Ly = 0;
    REAL8 chi1Lz = seobParams->chi1;
    REAL8 chi2Lx = 0;
    REAL8 chi2Ly = 0;
    REAL8 chi2Lz = seobParams->chi2;

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
    // STRUCTFREE(chi1temp, REAL8Vector);
    // STRUCTFREE(chi2temp, REAL8Vector);
    if (failed)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}


int SEOBJoinTimeVector(
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
    SEOBdynamics
        *seobdynamicsAdaS, /**<< Input: SEOB dynamics with adaptive-sampling */
    SEOBdynamics
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

int SEOBSAJoinTimeVector(
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
    SEOBSAdynamics
        *seobdynamicsAdaS, /**<< Input: SEOB dynamics with adaptive-sampling */
    SEOBSAdynamics
        *seobdynamicsHiS /**<< Input: SEOB dynamics with high-sampling */
) 
{
    /* Read from inputs */
    INT4 lenAdaS = seobdynamicsAdaS->length;
    REAL8 *tAdaS = seobdynamicsAdaS->tVec;
    INT4 lenHiS = seobdynamicsHiS->length;
    REAL8 *tHiS = seobdynamicsHiS->tVec;
// print_debug("lenAdaS = %d, lenHiS = %d, retLenHiSRDpatch = %d\n", lenAdaS, lenHiS, retLenHiSRDpatch);
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

    // print_debug("------test seobdynamicsAdaS---------\n");

    // INT j = 0;
    // print_debug("%d, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e\n", j, 
    //     seobdynamicsAdaS->tVec[j], seobdynamicsAdaS->rVec[j], seobdynamicsAdaS->phiVec[j],
    //     seobdynamicsAdaS->prTVec[j], seobdynamicsAdaS->pphiVec[j], seobdynamicsAdaS->drVec[j],
    //     seobdynamicsAdaS->dphiVec[j], seobdynamicsAdaS->dprTVec[j], seobdynamicsAdaS->dpphiVec[j],
    //     seobdynamicsAdaS->HVec[j]);

    return CEV_SUCCESS;
}


int SEOBJoinDynamics(
    SEOBdynamics **seobdynamicsJoined,     /**<< Output: pointer to joined dynamics */
    SEOBdynamics *seobdynamics1, /**<< Input: first dynamics */
    SEOBdynamics *seobdynamics2, /**<< Input: second dynamics */
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
    *seobdynamicsJoined = CreateSEOBdynamics(lenJoined);

    /* Copy truncated dynamics 1 - v4PdynamicsVariables data fields */
    for (j = 0; j < v4PdynamicsVariables; j++) 
    {
        memcpy(&((*seobdynamicsJoined)->array->data[j * lenJoined]),
            &(seobdynamics1->array->data[j * lendyn1]),
            lendyn1Joined * sizeof(REAL8));
    }

    /* Copy truncated dynamics 2 - v4PdynamicsVariables data fields */
    for (j = 0; j < v4PdynamicsVariables; j++) 
    {
        memcpy(&((*seobdynamicsJoined)->array->data[j * lenJoined + lendyn1Joined]),
            &(seobdynamics2->array->data[j * lendyn2]),
            lendyn2Joined * sizeof(REAL8));
    }

    return CEV_SUCCESS;
}

int SEOBSAJoinDynamics(
    SEOBSAdynamics **seobdynamicsJoined,     /**<< Output: pointer to joined dynamics */
    SEOBSAdynamics *seobdynamics1, /**<< Input: first dynamics */
    SEOBSAdynamics *seobdynamics2, /**<< Input: second dynamics */
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
    // print_debug("lenJoined = %d\n", lenJoined);
    // print_debug("lendyn1Joined = %d\n", lendyn1Joined);
    // print_debug("lendyn2Joined = %d\n", lendyn2Joined);
    // print_debug("lendyn1 = %d\n", lendyn1);
    // print_debug("lendyn2 = %d\n", lendyn2);

    /* Initialize output dynamics */
    // print_debug("Initialize output dynamics\n");
    *seobdynamicsJoined = CreateSEOBSAdynamics(lenJoined);
    // for (j=0; j<lendyn1; j++)
    j = 0;
    // print_debug("%d, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e\n", j, 
    //     seobdynamics1->tVec[j], seobdynamics1->rVec[j], seobdynamics1->phiVec[j],
    //     seobdynamics1->prTVec[j], seobdynamics1->pphiVec[j], seobdynamics1->drVec[j],
    //     seobdynamics1->dphiVec[j], seobdynamics1->dprTVec[j], seobdynamics1->dpphiVec[j],
    //     seobdynamics1->HVec[j]);
    // j = lendyn1-1;
    // print_debug("%d, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e\n", j, 
    //     seobdynamics1->tVec[j], seobdynamics1->rVec[j], seobdynamics1->phiVec[j],
    //     seobdynamics1->prTVec[j], seobdynamics1->pphiVec[j], seobdynamics1->drVec[j],
    //     seobdynamics1->dphiVec[j], seobdynamics1->dprTVec[j], seobdynamics1->dpphiVec[j],
    //     seobdynamics1->HVec[j]);

    // print_debug("--------------------------\n");
    // j = 0;
    // print_debug("%d, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e\n", j, 
    //     seobdynamics2->tVec[j], seobdynamics2->rVec[j], seobdynamics2->phiVec[j],
    //     seobdynamics2->prTVec[j], seobdynamics2->pphiVec[j], seobdynamics2->drVec[j],
    //     seobdynamics2->dphiVec[j], seobdynamics2->dprTVec[j], seobdynamics2->dpphiVec[j],
    //     seobdynamics2->HVec[j]);
    // j = lendyn2-1;
    // print_debug("%d, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e\n", j, 
    //     seobdynamics2->tVec[j], seobdynamics2->rVec[j], seobdynamics2->phiVec[j],
    //     seobdynamics2->prTVec[j], seobdynamics2->pphiVec[j], seobdynamics2->drVec[j],
    //     seobdynamics2->dphiVec[j], seobdynamics2->dprTVec[j], seobdynamics2->dpphiVec[j],
    //     seobdynamics2->HVec[j]);

    /* Copy truncated dynamics 1 - v4PdynamicsVariables data fields */
    // print_debug("Copy truncated dynamics 1 - v4PdynamicsVariables data fields\n");
    for (j = 0; j < v4SAdynamicsVariables; j++) 
    {
        // print_debug("Copy truncated dynamics 1 - v4PdynamicsVariables data fields: %d\n", j);
        memcpy(&((*seobdynamicsJoined)->array->data[j * lenJoined]),
            &(seobdynamics1->array->data[j * lendyn1]),
            lendyn1Joined * sizeof(REAL8));
    }

    /* Copy truncated dynamics 2 - v4PdynamicsVariables data fields */
    for (j = 0; j < v4SAdynamicsVariables; j++) 
    {
        // print_debug("Copy truncated dynamics 2 - v4PdynamicsVariables data fields: %d\n", j);
        memcpy(&((*seobdynamicsJoined)->array->data[j * lenJoined + lendyn1Joined]),
            &(seobdynamics2->array->data[j * lendyn2]),
            lendyn2Joined * sizeof(REAL8));
    }
    // print_debug("done\n");
    return CEV_SUCCESS;
}


int SEOBSAAttachAdaSandHiSRTimeVector(REAL8Vector **tVecPmodes, 
    SEOBSAdynamics *seobdynamicsAdaS, 
    SEOBSAdynamics *seobdynamicsHiS, 
    REAL8 tstartHiS)
{
    INT4 lenAdaS = seobdynamicsAdaS->length;
    REAL8 *tAdaS = seobdynamicsAdaS->tVec;
    INT4 lenHiS = seobdynamicsHiS->length;
    REAL8 *tHiS = seobdynamicsHiS->tVec;
    UINT iJoinHiS = 0, i;
    while ((iJoinHiS < lenAdaS - 1) && (tAdaS[iJoinHiS] < tstartHiS))
        iJoinHiS++;
    INT lenPmodes = iJoinHiS + lenHiS;
    if (!((*tVecPmodes) = CreateREAL8Vector(lenPmodes))) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failed to allocate REAL8Vector tVecPmodes.");
        return CEV_FAILURE;
    }
    /* Copy time values for AdaS and HiS */
    memcpy(&((*tVecPmodes)->data[0]), tAdaS, iJoinHiS * sizeof(REAL8));
    memcpy(&((*tVecPmodes)->data[iJoinHiS]), tHiS, lenHiS * sizeof(REAL8));
    return CEV_SUCCESS;
}

int SEOBAmplitudePeakFromAmp22Amp21(
    REAL8 *tPeak,                           /**<< Output: time of peak */
    UINT *indexPeak,                       /**<< Output: index of peak */
    SphHarmListCAmpPhaseSequence *listhPlm, /**<< Input: list of modes hlm */
    INT4 modes[][2],                        /**<< Input: array of modes (l,m) */
    UINT nmodes,     /**<< Input: number of modes (l,m) */
    REAL8Vector *tVec /**<< Input: vector of times */
) 
{
    UINT length = tVec->length;

    /* Looking form modes 22 and 21 */
    UINT found22 = 0, found21 = 0, nmode;
    for ( nmode = 0; nmode < nmodes; nmode++) 
    {
        if (modes[nmode][0] == 2 && modes[nmode][1] == 2)
        found22 = 1;
        if (modes[nmode][0] == 2 && modes[nmode][1] == 1)
        found21 = 1;
    }
    if ((!found22)) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "mode 22 not found.");
        return CEV_FAILURE;
    }

    CAmpPhaseSequence *hP22 = NULL;
    CAmpPhaseSequence *hP21 = NULL;

    hP22 = SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 2, 2)->campphase;
    if (found21)
        hP21 = SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 2, 1)->campphase;

    /* Check lengths */
    if ((!(hP22->xdata->length == length)) ||
        (found21 && !(hP21->xdata->length == length))) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "lengths of input amplitude and time REAL8Vector do not match.");
        return CEV_FAILURE;
    }

    /* Look for max combined amplitude - done discretely, no interpolation */
    UINT indexMax = 0, i;
    REAL8 AsquareMax = 0.;
    REAL8 A22_real = 0., A22_imag = 0., A21_real = 0., A21_imag = 0.,
            Asquare = 0.;
    for ( i = 0; i < length; i++) 
    {
        A22_real = hP22->camp_real->data[i];
        A22_imag = hP22->camp_imag->data[i];
        Asquare = A22_real * A22_real + A22_imag * A22_imag;
        if (found21) 
        {
            A21_real = hP21->camp_real->data[i];
            A21_imag = hP21->camp_imag->data[i];
            Asquare += A21_real * A21_real + A21_imag * A21_imag;
        }
        if (Asquare > AsquareMax) 
        {
            AsquareMax = Asquare;
            indexMax = i;
        }
    }

    /* Output */
    *indexPeak = indexMax;
    *tPeak = tVec->data[indexMax];

    return CEV_SUCCESS;
}

int SEOBJoinSphHarmListhlm(
    SphHarmListCAmpPhaseSequence **listhlm_joined, /**<< Output: list of joined modes */
    SphHarmListCAmpPhaseSequence *listhlm_1, /**<< Input: list of modes 1 */
    SphHarmListCAmpPhaseSequence *listhlm_2, /**<< Input: list of modes 2 */
    INT4 modes[][2],  /**<< Input: array of modes (l,m) */
    UINT nmodes,     /**<< Input: number of modes (l,m) */
    UINT indexJoin12 /**<< Input: index where to join the two dynamics */
) 
{
    /* Loop over modes */
    UINT nmode, i;
    for ( nmode = 0; nmode < nmodes; nmode++) 
    {

        INT4 l = modes[nmode][0];
        INT4 m = modes[nmode][1];

        /* Get the relevant modes in lists */
        CAmpPhaseSequence *hlm_1 =
            SphHarmListCAmpPhaseSequence_GetMode(listhlm_1, l, m)->campphase;
        CAmpPhaseSequence *hlm_2 =
            SphHarmListCAmpPhaseSequence_GetMode(listhlm_2, l, m)->campphase;

        /* Lengths */
        UINT len2 = hlm_2->xdata->length;
        INT lenJoined1 = indexJoin12;
        INT lenJoined = lenJoined1 + len2;

        /* Real and imaginary part vectors for the mode with RD attached */
        CAmpPhaseSequence *hlm_joined = NULL;
        if (CAmpPhaseSequence_Init(&hlm_joined, lenJoined) == CEV_FAILURE) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "failed to allocate CAmpPhaseSequence hlm_joined for mode (l,m) = (%d,%d).", l, m);
            return CEV_FAILURE;
        }

        /* Copy data, stopping part 1 at the specified index */
        memcpy(&(hlm_joined->xdata->data[0]), hlm_1->xdata->data,
                lenJoined1 * sizeof(REAL8));
        memcpy(&(hlm_joined->camp_real->data[0]), hlm_1->camp_real->data,
                lenJoined1 * sizeof(REAL8));
        memcpy(&(hlm_joined->camp_imag->data[0]), hlm_1->camp_imag->data,
                lenJoined1 * sizeof(REAL8));
        memcpy(&(hlm_joined->phase->data[0]), hlm_1->phase->data,
                lenJoined1 * sizeof(REAL8));
        /* Copy data, for part 2 starting from the specified index */
        memcpy(&(hlm_joined->xdata->data[lenJoined1]), hlm_2->xdata->data,
                len2 * sizeof(REAL8));
        memcpy(&(hlm_joined->camp_real->data[lenJoined1]), hlm_2->camp_real->data,
                len2 * sizeof(REAL8));
        memcpy(&(hlm_joined->camp_imag->data[lenJoined1]), hlm_2->camp_imag->data,
                len2 * sizeof(REAL8));
        memcpy(&(hlm_joined->phase->data[lenJoined1]), hlm_2->phase->data,
                len2 * sizeof(REAL8));
        /* Adjust for a 2kpi-shift in the phase - use the difference between the
            * last included sample of 1 and the first sample of 2 */
        REAL8 phase_diff =
            hlm_2->phase->data[0] - hlm_1->phase->data[lenJoined1 - 1];
        REAL8 shift_2kpi =
            -floor((phase_diff + CST_PI) / (2 * CST_PI)) * (2 * CST_PI);
        for ( i = 0; i < len2; i++) {
            hlm_joined->phase->data[lenJoined1 + i] += shift_2kpi;
        }

        SphHarmListCAmpPhaseSequence_AddMode(listhlm_joined, hlm_joined, l, m);
    }

    return CEV_SUCCESS;
}

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

int SEOBEulerJ2PFromDynamics(
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
    SEOBdynamics *seobdynamics, /**<<Input: SEOB dynamics (joined AdaS+HiS, up
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
        return CEV_FAILURE;
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

    if (!SpinsAlmostAligned) /* if spins are almost aligned, leave all Euler
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
        GSL_START;
        gsl_spline *x_spline = gsl_spline_alloc(gsl_interp_cspline, retLenDyn);
        gsl_spline *y_spline = gsl_spline_alloc(gsl_interp_cspline, retLenDyn);
        int status;
        gsl_interp_accel *x_acc = gsl_interp_accel_alloc();
        gsl_interp_accel *y_acc = gsl_interp_accel_alloc();
        status = gsl_spline_init(x_spline, seobdynamics->tVec, (*alphaJ2P)->data, retLenDyn);
        if (status != GSL_SUCCESS)
        {
            gsl_integration_workspace_free(precEulerw);
            gsl_spline_free(x_spline);
            gsl_spline_free(y_spline);
            gsl_interp_accel_free(x_acc);
            gsl_interp_accel_free(y_acc);
            GSL_END;
            return CEV_FAILURE;
        }
        status = gsl_spline_init(y_spline, seobdynamics->tVec, (*betaJ2P)->data, retLenDyn);
        if (status != GSL_SUCCESS)
        {
            gsl_integration_workspace_free(precEulerw);
            gsl_spline_free(x_spline);
            gsl_spline_free(y_spline);
            gsl_interp_accel_free(x_acc);
            gsl_interp_accel_free(y_acc);
            return CEV_FAILURE;
        }
        GSL_END;
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


int SEOBEulerJ2PPostMergerExtension(
    REAL8Vector
        *alphaJ2P, /**<< Output: vector for alpha J2P, already allocated */
    REAL8Vector
        *betaJ2P, /**<< Output: vector for beta J2P, already allocated */
    REAL8Vector
        *gammaJ2P, /**<< Output: vector for gamma J2P, already allocated */
    COMPLEX16
        sigmaQNM220, /**<< Input: complex frequency for QNM 22, 0th overtone */
    COMPLEX16
        sigmaQNM210, /**<< Input: complex frequency for QNM 21, 0th overtone */
    REAL8Vector *tVec, /**<< Input: time vector for Euler angles data (length of
                          P-modes) */
    UINT retLen,      /**<< Input: total length of Euler angles data (length of
                          P-modes) */
    UINT indexStart, /**<< Input: index where we start the extension (included,
                         index of time of attachment) */
    SpinEOBParams *seobParams, /**<< SEOB params */
    INT flip /** << a flag of whether to flip the sign of the precession frequency
              */
) 
{
    UINT i;
    UINT SpinsAlmostAligned = seobParams->alignedSpins;

    if (!SpinsAlmostAligned) 
    {

        /* Initial values, that were computed from the dynamics */
        REAL8 timeAttach = tVec->data[indexStart - 1];
        REAL8 alphaAttach = alphaJ2P->data[indexStart - 1];
        REAL8 betaAttach = betaJ2P->data[indexStart - 1];
        REAL8 gammaAttach = gammaJ2P->data[indexStart - 1];

        if (1/*flagEulerextension == FLAG_SEOBNRv4P_EULEREXT_QNM_SIMPLE_PRECESSION*/) 
        {
            /* Precession rate around final J */
            REAL8 omegaQNM220 = creal(sigmaQNM220);
            REAL8 omegaQNM210 = creal(sigmaQNM210);
            REAL8 precRate = omegaQNM220 - omegaQNM210;
            // flip is either 1 or -1. This is needed because we want the sign of the precession frequency
            // to be correct even when the projected final spin is negative, at which point the QNMs
            // flip sign
            precRate *= flip;
            REAL8 cosbetaAttach = cos(betaAttach);
            for (i = indexStart; i < retLen; i++) 
            {
                alphaJ2P->data[i] =
                    alphaAttach + (tVec->data[i] - timeAttach) * precRate;
                betaJ2P->data[i] = betaAttach;
                gammaJ2P->data[i] = gammaAttach - cosbetaAttach *
                                                        (tVec->data[i] - timeAttach) *
                                                        precRate;
            }
        }
        else if (0/*flagEulerextension == FLAG_SEOBNRv4P_EULEREXT_CONSTANT*/) 
        {
            for (i = indexStart; i < retLen; i++) 
            {
                alphaJ2P->data[i] = alphaAttach;
                betaJ2P->data[i] = betaAttach;
                gammaJ2P->data[i] = gammaAttach;
            }
        }
        else 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "flagEulerextension not recognized.");
            return CEV_FAILURE;
        }

    }
    else /* if spins are almost aligned, set all Euler angles to 0 */
    {
        memset(&(alphaJ2P->data[indexStart]), 0,
                (retLen - indexStart) * sizeof(REAL8));
        memset(&(betaJ2P->data[indexStart]), 0,
                (retLen - indexStart) * sizeof(REAL8));
        memset(&(gammaJ2P->data[indexStart]), 0,
                (retLen - indexStart) * sizeof(REAL8));
    }

    return CEV_SUCCESS;
}

int SEOBRotateInterpolatehJlmReImFromSphHarmListhPlmAmpPhase(

    SphHarmTimeSeries **hJlm, /**<< Output: hJlm time series, will contain
                                 complex values on fixed sampling */
    INT modes[][2],          /**<< Input: array of modes (l,m) */
    UINT nmodes,             /**<< Input: number of modes (l,m) */
    INT modes_lmax,          /**<< Input: maximum value of l in modes (l,m) */
    REAL8 deltaT,             /**<< Input: time step for the hJlm timeseries */
    UINT retLenTS, /**<< Input: number of samples for the hJlm timeseries */
    REAL8Vector *tVecPmodes, /**<< Input: irregular time vector on which the
                                hPlm and Euler angles are given */
    SphHarmListCAmpPhaseSequence
        *listhPlm,         /**<< Input: list of P-frame modes hPlm */
    REAL8Vector *alphaJ2P, /**<< Input: vector for Euler angle alpha J2P */
    REAL8Vector *betaJ2P,  /**<< Input: vector for Euler angle beta J2P */
    REAL8Vector *gammaJ2P /**<< Input: vector for Euler angle gamma J2P */
) 
{
    INT flagSymmetrizehPlminusm = 1;
    UINT i;
    INT l, m, mp;
    REAL8 t = 0., camp_real = 0., camp_imag = 0., phase = 0., alpha = 0., beta = 0., gamma = 0.;
    UINT retLenP = tVecPmodes->length;

    /* Create output time vector */
    REAL8Vector *timeTS = CreateREAL8Vector(retLenTS);
    REAL8 *timeTSdata = timeTS->data;
    for (i = 0; i < retLenTS; i++) {
        timeTSdata[i] = i * deltaT;
    }

    /* Create output list of timeseries, with all (l,m) up to modes_lmax */
    *hJlm = NULL;
    // LIGOTimeGPS tGPS = LIGOTIMEGPSZERO;
    REAL8 tGPS = 0;
    char mode_string[32];
    for (l = 2; l <= modes_lmax; l++) 
    {
        for (m = -l; m <= l; m++) 
        {
            // sprintf(mode_string, "H_%d%d", l, m);
            // COMPLEX16TimeSeries *hJlm_TS = XLALCreateCOMPLEX16TimeSeries(
            //     mode_string, &tGPS, 0., deltaT, &lalStrainUnit, retLenTS);
            COMPLEX16TimeSeries *hJlm_TS = CreateCOMPLEX16TimeSeries(tGPS, deltaT, retLenTS);
            memset(hJlm_TS->data->data, 0, retLenTS * sizeof(COMPLEX16));

            /* Note: with the AddMode function, data is copied over */
            *hJlm = XLALSphHarmTimeSeriesAddMode(*hJlm, hJlm_TS, l, m);

            /* Data has been copied over, we need to destroy */
            STRUCTFREE(hJlm_TS, COMPLEX16TimeSeries);
        }
    }

    /* Set time data */
    XLALSphHarmTimeSeriesSetTData(*hJlm, timeTS);

    /* Create working space for amp/phase of Wigner coefficient */
    REAL8Vector *Dlmmp_amp = CreateREAL8Vector(retLenP);
    REAL8Vector *Dlmmp_phase = CreateREAL8Vector(retLenP);
    REAL8Vector *Dlmminusmp_amp = CreateREAL8Vector(retLenP);
    REAL8Vector *Dlmminusmp_phase = CreateREAL8Vector(retLenP);
    memset(Dlmmp_amp->data, 0, (Dlmmp_amp->length) * sizeof(REAL8));
    memset(Dlmmp_phase->data, 0, (Dlmmp_phase->length) * sizeof(REAL8));
    memset(Dlmminusmp_amp->data, 0, (Dlmminusmp_amp->length) * sizeof(REAL8));
    memset(Dlmminusmp_phase->data, 0, (Dlmminusmp_phase->length) * sizeof(REAL8));

    /* Interpolating splines for Wigner coefficients */
    gsl_spline *spline_Dlmmp_amp = gsl_spline_alloc(gsl_interp_cspline, retLenP);
    gsl_spline *spline_Dlmmp_phase =
        gsl_spline_alloc(gsl_interp_cspline, retLenP);
    gsl_spline *spline_Dlmminusmp_amp =
        gsl_spline_alloc(gsl_interp_cspline, retLenP);
    gsl_spline *spline_Dlmminusmp_phase =
        gsl_spline_alloc(gsl_interp_cspline, retLenP);
    gsl_interp_accel *accel_Dlmmp_amp = gsl_interp_accel_alloc();
    gsl_interp_accel *accel_Dlmmp_phase = gsl_interp_accel_alloc();
    gsl_interp_accel *accel_Dlmminusmp_amp = gsl_interp_accel_alloc();
    gsl_interp_accel *accel_Dlmminusmp_phase = gsl_interp_accel_alloc();

    /* Interpolating splines for hPlm modes */
    GSL_START;
    int status, is_failed = 0;
    gsl_spline *spline_camp_real = gsl_spline_alloc(gsl_interp_cspline, retLenP);
    gsl_spline *spline_camp_imag = gsl_spline_alloc(gsl_interp_cspline, retLenP);
    gsl_spline *spline_phase = gsl_spline_alloc(gsl_interp_cspline, retLenP);
    gsl_interp_accel *accel_camp_real = gsl_interp_accel_alloc();
    gsl_interp_accel *accel_camp_imag = gsl_interp_accel_alloc();
    gsl_interp_accel *accel_phase = gsl_interp_accel_alloc();

    /* Interpolate P-frame modes hPlm on the time samples needed for the output as
    * a time series */
    SphHarmListCAmpPhaseSequence *listhPlm_TS = NULL;
    /* Loop over modes */
    UINT nmode;
    for ( nmode = 0; nmode < nmodes; nmode++) 
    {

        l = modes[nmode][0];
        m = modes[nmode][1];

        CAmpPhaseSequence *hPlm =
            SphHarmListCAmpPhaseSequence_GetMode(listhPlm, l, m)->campphase;

        CAmpPhaseSequence *hPlm_TS = NULL;
        if (CAmpPhaseSequence_Init(&hPlm_TS, retLenTS) == CEV_FAILURE) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "failure in CAmpPhaseSequence_Init for mode (l,m) = (%d,%d).", l, m);
            return CEV_FAILURE;
        }
        status = gsl_spline_init(spline_camp_real, tVecPmodes->data, hPlm->camp_real->data,
                        retLenP);
        if (status != GSL_SUCCESS)
        {
            is_failed = 1;
            goto QUIT;
        }
        status = gsl_spline_init(spline_camp_imag, tVecPmodes->data, hPlm->camp_imag->data,
                        retLenP);
        if (status != GSL_SUCCESS)
        {
            is_failed = 1;
            goto QUIT;
        }
        status = gsl_spline_init(spline_phase, tVecPmodes->data, hPlm->phase->data, retLenP);
        if (status != GSL_SUCCESS)
        {
            is_failed = 1;
            goto QUIT;
        }

        COMPLEX16 hPlm_val = 0.;
        REAL8 *hPlmTS_tdata = hPlm_TS->xdata->data;
        REAL8 *hPlmTS_camprealdata = hPlm_TS->camp_real->data;
        REAL8 *hPlmTS_campimagdata = hPlm_TS->camp_imag->data;
        REAL8 *hPlmTS_phasedata = hPlm_TS->phase->data;
        for (i = 0; i < retLenTS; i++) 
        {
            t = timeTSdata[i];
            /* Here we include a possible imaginary part of the complex envelope, but
            * at the moment it is simply 0 (only real part is used) */
            camp_real = gsl_spline_eval(spline_camp_real, t, accel_camp_real);
            camp_imag = gsl_spline_eval(spline_camp_imag, t, accel_camp_imag);
            phase = gsl_spline_eval(spline_phase, t, accel_phase);
            hPlm_val = (camp_real + I * camp_imag) * cexp(I * phase);
            /* We output the interpolated value for the mode hPlm in Re/Im form,
            * setting the phase to 0 */
            hPlmTS_tdata[i] = t;
            hPlmTS_camprealdata[i] = creal(hPlm_val);
            hPlmTS_campimagdata[i] = cimag(hPlm_val);
            hPlmTS_phasedata[i] = 0.;
        }

        SphHarmListCAmpPhaseSequence_AddMode(&listhPlm_TS, hPlm_TS, l, m);
    }

    /* Main computation */
    /* hJlm = \sum_mp Dlmmpstar hPlmp */

    REAL8 *Dlmmp_amp_data = Dlmmp_amp->data;
    REAL8 *Dlmmp_phase_data = Dlmmp_phase->data;
    REAL8 *Dlmminusmp_amp_data = Dlmminusmp_amp->data;
    REAL8 *Dlmminusmp_phase_data = Dlmminusmp_phase->data;
    REAL8 *alphadata = alphaJ2P->data;
    REAL8 *betadata = betaJ2P->data;
    REAL8 *gammadata = gammaJ2P->data;
    COMPLEX16TimeSeries *hJlmmode = NULL;
    COMPLEX16 *hJlmmode_data = NULL;
    REAL8 *hPlmp_campreal_data;
    REAL8 *hPlmp_campimag_data;
    COMPLEX16 hPlmp_val;
    COMPLEX16 Dlmmp_amp_val, Dlmmp_phase_val, Dlmmp_val;
    COMPLEX16 Dlmminusmp_amp_val, Dlmminusmp_phase_val, Dlmminusmp_val;

    /* Loop on l */
    for (l = 2; l <= modes_lmax; l++) 
    {

        /* Loop on m */
        for (m = -l; m <= l; m++) 
        {

            /* Get hJlm mode */
            hJlmmode = XLALSphHarmTimeSeriesGetMode(*hJlm, l, m);
            hJlmmode_data = hJlmmode->data->data;

            /* Loop on the modes hPlmp */
            for ( nmode = 0; nmode < nmodes; nmode++) 
            {

                /* Select modes with the same l */
                if (modes[nmode][0] != l)
                    continue;

                mp = modes[nmode][1];

                /* We do not allow mp<=0 in the P-frame modes hPlmp */
                /* We will use hPl-mp = (-1)^l hPlmp* */
                if (mp <= 0) 
                {
                    PRINT_LOG_INFO(LOG_CRITICAL, "mode (l,mp) = (%d,%d) is not allowed as mp<=0.", l, mp);
                    return CEV_FAILURE;
                }

                /* Get interpolated mode hPlmp */
                CAmpPhaseSequence *hPlmp_TS =
                    SphHarmListCAmpPhaseSequence_GetMode(listhPlm_TS, l, mp)->campphase;
                hPlmp_campreal_data = hPlmp_TS->camp_real->data;
                hPlmp_campimag_data = hPlmp_TS->camp_imag->data;

                /* Compute Wigner coefficients amp/phase on the same input sampling as
                    * the P-frame modes */
                for (i = 0; i < retLenP; i++) 
                {
                    alpha = alphadata[i];
                    beta = betadata[i];
                    gamma = gammadata[i];
                    Dlmmp_amp_data[i] = SEOBWignerDAmp(l, m, mp, beta);
                    Dlmmp_phase_data[i] = SEOBWignerDPhase(m, mp, alpha, gamma);
                    if (flagSymmetrizehPlminusm) 
                    {
                        Dlmminusmp_amp_data[i] = SEOBWignerDAmp(l, m, -mp, beta);
                        Dlmminusmp_phase_data[i] = SEOBWignerDPhase(m, -mp, alpha, gamma);
                    }
                }

                /* Interpolate amplitude/phase of the Wigner coefficient, add
                    * contribution to hJlm mode */
                status = gsl_spline_init(spline_Dlmmp_amp, tVecPmodes->data, Dlmmp_amp->data,
                                retLenP);
                if (status != GSL_SUCCESS)
                {
                    is_failed = 1;
                    goto QUIT;
                }
                status = gsl_spline_init(spline_Dlmmp_phase, tVecPmodes->data, Dlmmp_phase->data,
                                retLenP);
                if (status != GSL_SUCCESS)
                {
                    is_failed = 1;
                    goto QUIT;
                }
                if (flagSymmetrizehPlminusm) 
                {
                    status = gsl_spline_init(spline_Dlmminusmp_amp, tVecPmodes->data,
                                    Dlmminusmp_amp->data, retLenP);
                    if (status != GSL_SUCCESS)
                    {
                        is_failed = 1;
                        goto QUIT;
                    }
                    status = gsl_spline_init(spline_Dlmminusmp_phase, tVecPmodes->data,
                                    Dlmminusmp_phase->data, retLenP);
                    if (status != GSL_SUCCESS)
                    {
                        is_failed = 1;
                        goto QUIT;
                    }
                }
                for (i = 0; i < retLenTS; i++) 
                {
                    t = timeTSdata[i];
                    Dlmmp_amp_val = gsl_spline_eval(spline_Dlmmp_amp, t, accel_Dlmmp_amp);
                    Dlmmp_phase_val =
                        gsl_spline_eval(spline_Dlmmp_phase, t, accel_Dlmmp_phase);
                    Dlmmp_val =
                        Dlmmp_amp_val *
                        cexp(-I * Dlmmp_phase_val); /* mind the conjugation Dlmmpstar */
                    hPlmp_val = hPlmp_campreal_data[i] + I * hPlmp_campimag_data[i];
                    hJlmmode_data[i] += Dlmmp_val * hPlmp_val;
                    if (flagSymmetrizehPlminusm) 
                    {
                        Dlmminusmp_amp_val =
                            gsl_spline_eval(spline_Dlmminusmp_amp, t, accel_Dlmminusmp_amp);
                        Dlmminusmp_phase_val = gsl_spline_eval(spline_Dlmminusmp_phase, t,
                                                                accel_Dlmminusmp_phase);
                        Dlmminusmp_val =
                            Dlmminusmp_amp_val *
                            cexp(-I * Dlmminusmp_phase_val); /* mind the conjugation
                                                                Dlmminusmpstar */
                        hJlmmode_data[i] += Dlmminusmp_val * pow(-1, l) * conj(hPlmp_val);
                    }
                }
            }
        }
    }
#if 0
    {
        CAmpPhaseSequence *h22 = SphHarmListCAmpPhaseSequence_GetMode(listhPlm_TS, 2, 2)->campphase;
        CAmpPhaseSequence *h21 = SphHarmListCAmpPhaseSequence_GetMode(listhPlm_TS, 2, 1)->campphase;
        CAmpPhaseSequence *h33 = SphHarmListCAmpPhaseSequence_GetMode(listhPlm_TS, 3, 3)->campphase;
        CAmpPhaseSequence *h44 = SphHarmListCAmpPhaseSequence_GetMode(listhPlm_TS, 4, 4)->campphase;
        CAmpPhaseSequence *h55 = SphHarmListCAmpPhaseSequence_GetMode(listhPlm_TS, 5, 5)->campphase;
        char fout[STR_COMM_SIZE];
        strncpy(fout, "tmp_debug_ecc.dat", STR_COMM_SIZE);
        FILE *out = fopen(fout, "w");
        for (i = 0; i < retLenTS; i++)
        {
            t = timeTSdata[i];
            fprintf(out, "%.16e\t%.16e\t%.16e\n", 
                t, h22->camp_real->data[i],
                   h22->camp_imag->data[i]);
        }
        fclose(out);
    }
#endif
QUIT:
    /* Cleanup */
    STRUCTFREE(listhPlm_TS, SphHarmListCAmpPhaseSequence);
    gsl_spline_free(spline_camp_real);
    gsl_spline_free(spline_camp_imag);
    gsl_spline_free(spline_phase);
    gsl_interp_accel_free(accel_camp_real);
    gsl_interp_accel_free(accel_camp_imag);
    gsl_interp_accel_free(accel_phase);
    gsl_spline_free(spline_Dlmmp_amp);
    gsl_spline_free(spline_Dlmmp_phase);
    gsl_spline_free(spline_Dlmminusmp_amp);
    gsl_spline_free(spline_Dlmminusmp_phase);
    gsl_interp_accel_free(accel_Dlmmp_amp);
    gsl_interp_accel_free(accel_Dlmmp_phase);
    gsl_interp_accel_free(accel_Dlmminusmp_amp);
    gsl_interp_accel_free(accel_Dlmminusmp_phase);
    GSL_END;
    STRUCTFREE(Dlmmp_amp, REAL8Vector);
    STRUCTFREE(Dlmmp_phase, REAL8Vector);
    STRUCTFREE(Dlmminusmp_amp, REAL8Vector);
    STRUCTFREE(Dlmminusmp_phase, REAL8Vector);
    if (is_failed)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}

int SEOBPlmListToSphHarmTimeSeries(

    SphHarmTimeSeries **hJlm, /**<< Output: hJlm time series, will contain
                                 complex values on fixed sampling */
    INT modes[][2],          /**<< Input: array of modes (l,m) */
    UINT nmodes,             /**<< Input: number of modes (l,m) */
    INT modes_lmax,          /**<< Input: maximum value of l in modes (l,m) */
    REAL8 deltaT,             /**<< Input: time step for the hJlm timeseries */
    UINT retLenTS, /**<< Input: number of samples for the hJlm timeseries */
    SphHarmListCAmpPhaseSequence
        *listhPlm,
    REAL8Vector *alphaJ2P, /**<< Input: vector for Euler angle alpha J2P */
    REAL8Vector *betaJ2P,  /**<< Input: vector for Euler angle beta J2P */
    REAL8Vector *gammaJ2P /**<< Input: vector for Euler angle gamma J2P */
)
{
    INT flagSymmetrizehPlminusm = 1;
    UINT i;
    INT l, m, mp;
    REAL8 t = 0., camp_real = 0., camp_imag = 0., phase = 0., alpha = 0., beta = 0., gamma = 0.;
    UINT retLenP = retLenTS;

    /* Create output time vector */
    REAL8Vector *timeTS = CreateREAL8Vector(retLenTS);
    REAL8 *timeTSdata = timeTS->data;
    for (i = 0; i < retLenTS; i++) {
        timeTSdata[i] = i * deltaT;
    }

    /* Create output list of timeseries, with all (l,m) up to modes_lmax */
    *hJlm = NULL;
    // LIGOTimeGPS tGPS = LIGOTIMEGPSZERO;
    REAL8 tGPS = 0;
    char mode_string[32];
    for (l = 2; l <= modes_lmax; l++) 
    {
        for (m = -l; m <= l; m++) 
        {
            // sprintf(mode_string, "H_%d%d", l, m);
            // COMPLEX16TimeSeries *hJlm_TS = XLALCreateCOMPLEX16TimeSeries(
            //     mode_string, &tGPS, 0., deltaT, &lalStrainUnit, retLenTS);
            COMPLEX16TimeSeries *hJlm_TS = CreateCOMPLEX16TimeSeries(tGPS, deltaT, retLenTS);
            memset(hJlm_TS->data->data, 0, retLenTS * sizeof(COMPLEX16));

            /* Note: with the AddMode function, data is copied over */
            *hJlm = XLALSphHarmTimeSeriesAddMode(*hJlm, hJlm_TS, l, m);

            /* Data has been copied over, we need to destroy */
            STRUCTFREE(hJlm_TS, COMPLEX16TimeSeries);
        }
    }

    /* Set time data */
    XLALSphHarmTimeSeriesSetTData(*hJlm, timeTS);

    /* Create working space for amp/phase of Wigner coefficient */
    REAL8Vector *Dlmmp_amp = CreateREAL8Vector(retLenP);
    REAL8Vector *Dlmmp_phase = CreateREAL8Vector(retLenP);
    REAL8Vector *Dlmminusmp_amp = CreateREAL8Vector(retLenP);
    REAL8Vector *Dlmminusmp_phase = CreateREAL8Vector(retLenP);
    memset(Dlmmp_amp->data, 0, (Dlmmp_amp->length) * sizeof(REAL8));
    memset(Dlmmp_phase->data, 0, (Dlmmp_phase->length) * sizeof(REAL8));
    memset(Dlmminusmp_amp->data, 0, (Dlmminusmp_amp->length) * sizeof(REAL8));
    memset(Dlmminusmp_phase->data, 0, (Dlmminusmp_phase->length) * sizeof(REAL8));

    /* Interpolating splines for Wigner coefficients */
    // gsl_spline *spline_Dlmmp_amp = gsl_spline_alloc(gsl_interp_cspline, retLenP);
    // gsl_spline *spline_Dlmmp_phase =
    //     gsl_spline_alloc(gsl_interp_cspline, retLenP);
    // gsl_spline *spline_Dlmminusmp_amp =
    //     gsl_spline_alloc(gsl_interp_cspline, retLenP);
    // gsl_spline *spline_Dlmminusmp_phase =
    //     gsl_spline_alloc(gsl_interp_cspline, retLenP);
    // gsl_interp_accel *accel_Dlmmp_amp = gsl_interp_accel_alloc();
    // gsl_interp_accel *accel_Dlmmp_phase = gsl_interp_accel_alloc();
    // gsl_interp_accel *accel_Dlmminusmp_amp = gsl_interp_accel_alloc();
    // gsl_interp_accel *accel_Dlmminusmp_phase = gsl_interp_accel_alloc();

    /* Interpolating splines for hPlm modes */
    // gsl_spline *spline_camp_real = gsl_spline_alloc(gsl_interp_cspline, retLenP);
    // gsl_spline *spline_camp_imag = gsl_spline_alloc(gsl_interp_cspline, retLenP);
    // gsl_spline *spline_phase = gsl_spline_alloc(gsl_interp_cspline, retLenP);
    // gsl_interp_accel *accel_camp_real = gsl_interp_accel_alloc();
    // gsl_interp_accel *accel_camp_imag = gsl_interp_accel_alloc();
    // gsl_interp_accel *accel_phase = gsl_interp_accel_alloc();

    UINT nmode;

    /* Main computation */
    /* hJlm = \sum_mp Dlmmpstar hPlmp */

    REAL8 *Dlmmp_amp_data = Dlmmp_amp->data;
    REAL8 *Dlmmp_phase_data = Dlmmp_phase->data;
    REAL8 *Dlmminusmp_amp_data = Dlmminusmp_amp->data;
    REAL8 *Dlmminusmp_phase_data = Dlmminusmp_phase->data;
    REAL8 *alphadata = alphaJ2P->data;
    REAL8 *betadata = betaJ2P->data;
    REAL8 *gammadata = gammaJ2P->data;
    COMPLEX16TimeSeries *hJlmmode = NULL;
    COMPLEX16 *hJlmmode_data = NULL;
    REAL8 *hPlmp_campreal_data;
    REAL8 *hPlmp_campimag_data;
    COMPLEX16 hPlmp_val;
    COMPLEX16 Dlmmp_amp_val, Dlmmp_phase_val, Dlmmp_val;
    COMPLEX16 Dlmminusmp_amp_val, Dlmminusmp_phase_val, Dlmminusmp_val;

    // for ( nmode = 0; nmode < nmodes; nmode++) 
    // {
    //     /* Select modes with the same l */
    //     // if (modes[nmode][0] != l)
    //     //     continue;
    //     l = modes[nmode][0];
    //     m = modes[nmode][1];
    //     hJlmmode = XLALSphHarmTimeSeriesGetMode(*hJlm, l, m);
    //     hJlmmode_data = hJlmmode->data->data;
    //     CAmpPhaseSequence *hPlmp_TS =
    //         SphHarmListCAmpPhaseSequence_GetMode(listhPlm, l, m)->campphase;
    //     hPlmp_campreal_data = hPlmp_TS->camp_real->data;
    //     hPlmp_campimag_data = hPlmp_TS->camp_imag->data;
    //     for (i = 0; i < retLenTS; i++) 
    //     {
    //         hPlmp_val = hPlmp_campreal_data[i] + I * hPlmp_campimag_data[i];
    //         hJlmmode_data[i] = hPlmp_val;
    //     }
    // }

#if 1
    /* Loop on l */
    for (l = 2; l <= modes_lmax; l++) 
    {

        /* Loop on m */
        for (m = -l; m <= l; m++) 
        {
            /* Get hJlm mode */
            hJlmmode = XLALSphHarmTimeSeriesGetMode(*hJlm, l, m);
            hJlmmode_data = hJlmmode->data->data;

            /* Get interpolated mode hPlmp */
            // CAmpPhaseSequence *hPlmp_TS =
            //     SphHarmListCAmpPhaseSequence_GetMode(listhPlm, l, m)->campphase;
            // hPlmp_campreal_data = hPlmp_TS->camp_real->data;
            // hPlmp_campimag_data = hPlmp_TS->camp_imag->data;
            // for (i = 0; i < retLenTS; i++) 
            // {
            //     hPlmp_val = hPlmp_campreal_data[i] + I * hPlmp_campimag_data[i];
            //     hJlmmode_data[i] = hPlmp_val;
            // }
                // print_debug("l,m = (%d, %d)\n", l, m);

            /* Loop on the modes hPlmp */
            for ( nmode = 0; nmode < nmodes; nmode++) 
            {
                /* Select modes with the same l */
                if (modes[nmode][0] != l)
                    continue;

                mp = modes[nmode][1];

                /* We do not allow mp<=0 in the P-frame modes hPlmp */
                /* We will use hPl-mp = (-1)^l hPlmp* */
                if (mp <= 0) 
                {
                    PRINT_LOG_INFO(LOG_CRITICAL, "mode (l,mp) = (%d,%d) is not allowed as mp<=0.", l, mp);
                    return CEV_FAILURE;
                }
                /* Get interpolated mode hPlmp */
                CAmpPhaseSequence *hPlmp_TS =
                    SphHarmListCAmpPhaseSequence_GetMode(listhPlm, l, mp)->campphase;
                hPlmp_campreal_data = hPlmp_TS->camp_real->data;
                hPlmp_campimag_data = hPlmp_TS->camp_imag->data;

                /* Compute Wigner coefficients amp/phase on the same input sampling as
                    * the P-frame modes */
                for (i = 0; i < retLenP; i++) 
                {
                    alpha = alphadata[i];
                    beta = betadata[i];
                    gamma = gammadata[i];
                    Dlmmp_amp_data[i] = SEOBWignerDAmp(l, m, mp, beta);
                    Dlmmp_phase_data[i] = SEOBWignerDPhase(m, mp, alpha, gamma);
                    if (flagSymmetrizehPlminusm) 
                    {
                        Dlmminusmp_amp_data[i] = SEOBWignerDAmp(l, m, -mp, beta);
                        Dlmminusmp_phase_data[i] = SEOBWignerDPhase(m, -mp, alpha, gamma);
                    }
                }

                /* Interpolate amplitude/phase of the Wigner coefficient, add
                    * contribution to hJlm mode */
                for (i = 0; i < retLenTS; i++) 
                {
                    t = timeTSdata[i];
                    // Dlmmp_amp_val = gsl_spline_eval(spline_Dlmmp_amp, t, accel_Dlmmp_amp);
                    // Dlmmp_phase_val =
                    //     gsl_spline_eval(spline_Dlmmp_phase, t, accel_Dlmmp_phase);
                    Dlmmp_amp_val = Dlmmp_amp_data[i];
                    Dlmmp_phase_val = Dlmmp_phase_data[i];
                    Dlmmp_val = Dlmmp_amp_val * cexp(-I * Dlmmp_phase_val); /* mind the conjugation Dlmmpstar */
                    hPlmp_val = hPlmp_campreal_data[i] + I * hPlmp_campimag_data[i];
                    // print_debug("[%d]hPlmp_val = (%.16e, %.16e)\n", creal(hPlmp_val), cimag(hPlmp_val));
                    // hPlmp_val = 0.0;
                    hJlmmode_data[i] += Dlmmp_val * hPlmp_val;
                    if (flagSymmetrizehPlminusm) 
                    {
                        // Dlmminusmp_amp_val =
                        //     gsl_spline_eval(spline_Dlmminusmp_amp, t, accel_Dlmminusmp_amp);
                        // Dlmminusmp_phase_val = gsl_spline_eval(spline_Dlmminusmp_phase, t,
                        //                                         accel_Dlmminusmp_phase);
                        Dlmminusmp_amp_val = Dlmminusmp_amp_data[i];
                        Dlmminusmp_phase_val = Dlmminusmp_phase_data[i];
                        Dlmminusmp_val = Dlmminusmp_amp_val * cexp(-I * Dlmminusmp_phase_val); 
                        /* mind the conjugation Dlmminusmpstar */
                        hJlmmode_data[i] += Dlmminusmp_val * pow(-1, l) * conj(hPlmp_val);
                    }
                }
            }
        }
    }
#endif

    /* Cleanup */
    // gsl_spline_free(spline_camp_real);
    // gsl_spline_free(spline_camp_imag);
    // gsl_spline_free(spline_phase);
    // gsl_interp_accel_free(accel_camp_real);
    // gsl_interp_accel_free(accel_camp_imag);
    // gsl_interp_accel_free(accel_phase);
    // gsl_spline_free(spline_Dlmmp_amp);
    // gsl_spline_free(spline_Dlmmp_phase);
    // gsl_spline_free(spline_Dlmminusmp_amp);
    // gsl_spline_free(spline_Dlmminusmp_phase);
    // gsl_interp_accel_free(accel_Dlmmp_amp);
    // gsl_interp_accel_free(accel_Dlmmp_phase);
    // gsl_interp_accel_free(accel_Dlmminusmp_amp);
    // gsl_interp_accel_free(accel_Dlmminusmp_phase);
    STRUCTFREE(Dlmmp_amp, REAL8Vector);
    STRUCTFREE(Dlmmp_phase, REAL8Vector);
    STRUCTFREE(Dlmminusmp_amp, REAL8Vector);
    STRUCTFREE(Dlmminusmp_phase, REAL8Vector);

    return CEV_SUCCESS;
}

int SEOBRotatehIlmFromhJlm(
    SphHarmTimeSeries **hIlm, /**<< Output: hIlm time series, complex values on
                                 fixed sampling */
    SphHarmTimeSeries *hJlm,  /**<< Output: hJlm time series, complex values on
                                 fixed sampling */
    INT modes_lmax,          /**<< Input: maximum value of l in modes (l,m) */
    REAL8 alphaI2J,           /**<< Input: Euler angle alpha I->J */
    REAL8 betaI2J,            /**<< Input: Euler angle beta I->J */
    REAL8 gammaI2J,           /**<< Input: Euler angle gamma I->J */
    REAL8 deltaT              /**<< Input: time step, necessary to initialize new timeseries */
) 
{
    UINT i;
    INT l, m, mp;
    REAL8 amp_wigner = 0., phase_wigner = 0.;
    COMPLEX16 D_wigner = 0.;
    UINT retLen = hJlm->tdata->length;
    REAL8 *tJdata = hJlm->tdata->data;

    /* Copy time vector */
    REAL8Vector *tI = CreateREAL8Vector(retLen);
    memcpy(tI->data, tJdata, retLen * sizeof(REAL8));

    /* Create output list of timeseries, with all (l,m) up to modes_lmax */
    *hIlm = NULL;
    // LIGOTimeGPS tGPS = LIGOTIMEGPSZERO;
    REAL8 tGPS = 0;
    char mode_string[32];
    for (l = 2; l <= modes_lmax; l++) 
    {
        for (m = -l; m <= l; m++) 
        {
            // sprintf(mode_string, "H_%d%d", l, m);
            // COMPLEX16TimeSeries *hIlm_TS = CreateCOMPLEX16TimeSeries(
            //     mode_string, &tGPS, 0., deltaT, &lalStrainUnit, retLen);
            COMPLEX16TimeSeries *hIlm_TS = CreateCOMPLEX16TimeSeries(tGPS, deltaT, retLen);
            memset(hIlm_TS->data->data, 0, retLen * sizeof(COMPLEX16));

            /* Note: with the AddMode function, data is copied over */
            *hIlm = XLALSphHarmTimeSeriesAddMode(*hIlm, hIlm_TS, l, m);

            /* Data has been copied over, we need to destroy */
            STRUCTFREE(hIlm_TS, COMPLEX16TimeSeries);
        }
    }

    /* Set time data */
    XLALSphHarmTimeSeriesSetTData(*hIlm, tI);

    /* Main computation */
    /* hIlm = \sum_mp Dlmpm hJlmp */

    COMPLEX16TimeSeries *hIlmmode = NULL;
    COMPLEX16TimeSeries *hJlmpmode = NULL;
    COMPLEX16 *hIlmmode_data = NULL;
    COMPLEX16 *hJlmpmode_data = NULL;

    /* Loop on l */
    for (l = 2; l <= modes_lmax; l++) 
    {
        /* Loop on m */
        for (m = -l; m <= l; m++) 
        {
            /* Get hJlm mode */
            hIlmmode = XLALSphHarmTimeSeriesGetMode(*hIlm, l, m);
            hIlmmode_data = hIlmmode->data->data;

                /* Loop on mp - exclude value 0, since hPl0=0 in our approximation */
            for (mp = -l; mp <= l; mp++) 
            {
                /* Get hJlm mode */
                hJlmpmode = XLALSphHarmTimeSeriesGetMode(hJlm, l, mp);
                hJlmpmode_data = hJlmpmode->data->data;

                /* Compute constant Wigner coefficient */
                amp_wigner = SEOBWignerDAmp(l, m, mp, betaI2J);
                phase_wigner = SEOBWignerDPhase(m, mp, alphaI2J, gammaI2J);
                D_wigner = amp_wigner *
                            cexp(-I * phase_wigner); /* mind the conjugation Dlmmpstar */

                /* Evaluate mode contribution */
                for (i = 0; i < retLen; i++) 
                {
                    hIlmmode_data[i] += D_wigner * hJlmpmode_data[i];
                }
            }
        }
    }

    return CEV_SUCCESS;
}

int SEOBComputehplushcrossFromhIlm(
    REAL8TimeSeries *hplusTS, /**<< Output: time series for hplus, already created */
    REAL8TimeSeries *hcrossTS,   /**<< Output: time series for hplus, already created */
    INT modes_lmax, /**<< Input: maximum value of l */
    SphHarmTimeSeries *hIlm,  /**<< Input: list with time series for each mode hIlm */
    REAL8 amp0, /**<< Input: amplitude prefactor */
    REAL8 inc,  /**<< Input: inclination */
    REAL8 phi,   /**<< Input: phase */
    INT is_only22
) 
{
    INT l, m;
    UINT i;
    /* hplus, hcross */
    REAL8 *hplusdata = hplusTS->data->data;
    memset(hplusdata, 0, hplusTS->data->length * sizeof(REAL8));
    REAL8 *hcrossdata = hcrossTS->data->data;
    memset(hcrossdata, 0, hplusTS->data->length * sizeof(REAL8));
    COMPLEX16 hpc_contrib = 0.;

    /* Loop over modes */
    for (l = 2; l <= modes_lmax; l++) 
    {
        for (m = -l; m <= l; m++) 
        {
            if (is_only22 && l != 2 && (m!=2 || m!=-2))
                continue;
            /* Compute sYlm */
            COMPLEX16 sYlm;
            // COMPLEX16 sYlm =
            //     XLALSpinWeightedSphericalHarmonic(inc, CST_PI / 2. - phi, -2, l, m);
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - phi, -2, l, m, &sYlm);
            /* Get mode hIlm */
            COMPLEX16TimeSeries *hIlmTS = XLALSphHarmTimeSeriesGetMode(hIlm, l, m);
            COMPLEX16 *hIlm_data = hIlmTS->data->data;

            for ( i = 0; i < hplusTS->data->length; i++) 
            {
                hpc_contrib = sYlm * hIlm_data[i];
                hplusdata[i] += amp0 * creal(hpc_contrib);
                hcrossdata[i] += -amp0 * cimag(hpc_contrib);
            }
        }
    }

    return CEV_SUCCESS;
}

INT SEOBCalculateNQCWindowFactorsFromDyn(SEOBdynamics *dyn,
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
        pr_n = dyn->polarprVec[i];
        pr_o = dyn->polarprVec[i+1];
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

INT SEOBSACalculateNQCWindowFactorsFromDyn(SEOBSAdynamics *dyn,
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
        if (dyn->rVec[i] < tThresh) continue;
        if (is_first)
        {
            // dr_n = (dyn->posVecx[i]*dyn->velVecx[i] + 
            //     dyn->posVecy[i]*dyn->velVecy[i] + 
            //     dyn->posVecz[i]*dyn->velVecz[i]) / dyn->polarrVec[i];
            dr_n = dyn->drVec[i];
            // dr_o = (dyn->posVecx[i+1]*dyn->velVecx[i+1] + 
            //     dyn->posVecy[i+1]*dyn->velVecy[i+1] + 
            //     dyn->posVecz[i+1]*dyn->velVecz[i+1]) / dyn->polarrVec[i+1];
            dr_o = dyn->drVec[i+1];
            if (dr_n * dr_o < 0)
            {
                ip1 = i;
                ip2 = i;
                break;
            }
            continue;
        }
        // pr_n = dyn->polarprVec[i];
        pr_n = dyn->prTVec[i];
        // pr_o = dyn->polarprVec[i+1];
        pr_o = dyn->prTVec[i+1];
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

INT dbg_CalculateNQCTimeSeries(SEOBdynamics *seobdynamics, 
                                SpinEOBParams *seobParams, 
                                SphHarmListEOBNonQCCoeffs *listnqcCoeffs,
                                INT ModeL, 
                                INT ModeM, 
                                COMPLEX16TimeSeries **ret)
{
    EOBNonQCCoeffs *nqcCoeffs =
        SphHarmListEOBNonQCCoeffs_GetMode(listnqcCoeffs, ModeL, ModeM)->nqcCoeffs;
    // print_debug("nqc: a1 = %e, a2 = %e, a3 = %e\n", nqcCoeffs->a1, nqcCoeffs->a2, nqcCoeffs->a3);
    // print_debug("nqc: b1 = %e, b2 = %e\n", nqcCoeffs->b1, nqcCoeffs->b2);
    COMPLEX16 factor_nqc = 1.;
    REAL8Vector values;
    REAL8 valuesdata[14] = {0.};
    values.length = 14;
    values.data = valuesdata;
    UINT retLen = seobdynamics->length;
    REAL8 deltaT = seobdynamics->tVec[1] - seobdynamics->tVec[0];
    REAL8 s1dotZ, s2dotZ, chiS, chiA, tplspin, t, omega, ham, v;
    COMPLEX16TimeSeries *nqcSeries = CreateCOMPLEX16TimeSeries(0.0, deltaT, retLen);
    UINT i, j;
    for (i = 0; i < retLen; i++) 
    {
        t = seobdynamics->tVec[i];
        omega = seobdynamics->omegaVec[i];
        for (j = 0; j < 14; j++) 
            values.data[j] = seobdynamics->array->data[i + (j + 1) * retLen];
        if (XLALSimIMRSpinEOBNonQCCorrection(&factor_nqc, &values, omega, t, 
                                seobParams->tWind, seobParams->wWind,
                                nqcCoeffs) == CEV_FAILURE) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL,
                "failure in XLALSimIMRSpinEOBNonQCCorrection at step %d of the loop.", i);
            return CEV_FAILURE;
        }
        nqcSeries->data->data[i] = factor_nqc;
    }
    *ret = nqcSeries;
    return CEV_SUCCESS;
}

INT dbg_CalculateWaveformFromDynamicsAdaS(SEOBdynamics *seobdynamics, 
                                          SpinEOBParams *seobParams, 
                                          INT ModeL, 
                                          INT ModeM, 
                                          COMPLEX16TimeSeries **out_new, 
                                          COMPLEX16TimeSeries **out_old)
{
    /* Check that the input double pointer are not NULL */
    /* Masses */
    REAL8 m1 = seobParams->m1;
    REAL8 m2 = seobParams->m2;
    REAL8 deltaT = seobdynamics->tVec[1] - seobdynamics->tVec[0];
    REAL8 eta = seobParams->eta;
    REAL8 mtot = m1 + m2;
    REAL8 dr, ncrv;
    // UINT SpinAlignedEOBversion = seobParams->SpinAlignedEOBversion;
    UINT SpinAlignedEOBversion = 4;
    UINT SpinAlignedEOBversionWaveform; // RC: I use this different variable
                                        // because the PN terms in the waveform
                                        // are different from those in the flux

    /* Length of dynamics data and sampling step */
    UINT retLen = seobdynamics->length;
    COMPLEX16 hLM1, hLM2;
    COMPLEX16TimeSeries *hLMVec1 = NULL;
    COMPLEX16TimeSeries *hLMVec2 = NULL;
    hLMVec1 = CreateCOMPLEX16TimeSeries(0.0, deltaT, retLen);
    hLMVec2 = CreateCOMPLEX16TimeSeries(0.0, deltaT, retLen);
    /* Workspace vectors */
    REAL8Vector values, polarDynamics;
    REAL8 valuesdata[14] = {0.};
    REAL8 polarDynamicsdata[4] = {0.};
    values.length = 14;
    polarDynamics.length = 4;
    values.data = valuesdata;
    polarDynamics.data = polarDynamicsdata;
    REAL8 tPeakOmega = seobParams->tPeakOmega;

    /* Loop to compute compute amplitude and phase of the hlm mode */
    REAL8 s1dotZ, s2dotZ, chiS, chiA, tplspin, t, omega, ham, v;
    UINT i, j;
    for (i = 0; i < retLen; i++) 
    {
        /* Compute waveform coefficients */
        t = seobdynamics->tVec[i];
        omega = seobdynamics->omegaVec[i];
        ham = seobdynamics->hamVec[i];
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

        if (CODE_VERSION == 3)
        {
            if (EccPrec_CalcSpinPrecFacWaveformCoefficients(
                    seobParams->hCoeffs, m1, m2, eta, tplspin, chiS, chiA,
                    seobdynamics->chiSxVec[i], seobdynamics->chiSyVec[i], seobdynamics->chiSzVec[i],
                    seobdynamics->chiAxVec[i], seobdynamics->chiAyVec[i], seobdynamics->chiAzVec[i],
                    451) == CEV_FAILURE) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in EccPrec_CalcSpinPrecFacWaveformCoefficients at step %d of the loop.", i);
                return CEV_FAILURE;
            }
        }
        else
        {
            if (XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(
                    seobParams->hCoeffs, m1, m2, eta, tplspin, chiS, chiA,
                    SpinAlignedEOBversionWaveform) == CEV_FAILURE) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients at step %d of the loop.", i);
                return CEV_FAILURE;
            }
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
        polarDynamics.data[2] = seobdynamics->polarprVec[i];
        polarDynamics.data[3] = seobdynamics->polarpphiVec[i];
        v = cbrt(omega);
        if (XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform(
                &hLM1, &polarDynamics, &values, v, ham, ModeL, ModeM, seobParams) != CEV_SUCCESS) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform at step %d of the loop.", i);
            return CEV_FAILURE;
        }
        dr = (seobdynamics->posVecx[i]*seobdynamics->velVecx[i] + 
            seobdynamics->posVecy[i]*seobdynamics->velVecy[i] + 
            seobdynamics->posVecz[i]*seobdynamics->velVecz[i]) / seobdynamics->polarrVec[i];
        ncrv = seobdynamics->polarrVec[i] * omega;

        
        if (XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2(
                &hLM2, &polarDynamics, &values, v, dr, ncrv, seobdynamics->polarprDotVec[i], ham, ModeL, ModeM, seobParams) != CEV_SUCCESS) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2 at step %d of the loop.", i);
            return CEV_FAILURE;
        }
        hLMVec1->data->data[i] = hLM2;
        hLMVec2->data->data[i] = hLM1;
    }
    *out_new = hLMVec1;
    *out_old = hLMVec2;
    return CEV_SUCCESS;
}

INT CheckStopCondition(SEOBdynamics *seobdynamics, 
                       SpinEOBParams *seobParams, 
                       SphHarmListEOBNonQCCoeffs *listnqcCoeffs,
                       REAL8 tPeak)
{
    INT ModeL = 2;
    INT ModeM = 2;
    EOBNonQCCoeffs *nqcCoeffs =
        SphHarmListEOBNonQCCoeffs_GetMode(listnqcCoeffs, ModeL, ModeM)->nqcCoeffs;
    /* Check that the input double pointer are not NULL */
    /* Masses */
    REAL8 m1 = seobParams->m1;
    REAL8 m2 = seobParams->m2;
    REAL8 deltaT = seobdynamics->tVec[1] - seobdynamics->tVec[0];
    REAL8 eta = seobParams->eta;
    REAL8 mtot = m1 + m2;
    REAL8 dr, ncrv;
    // UINT SpinAlignedEOBversion = seobParams->SpinAlignedEOBversion;
    UINT SpinAlignedEOBversion = 4;
    UINT SpinAlignedEOBversionWaveform; // RC: I use this different variable
                                        // because the PN terms in the waveform
                                        // are different from those in the flux

    /* Length of dynamics data and sampling step */
    UINT retLen = seobdynamics->length;
    COMPLEX16 hLM;
    // REAL8Vector *amphLMVec = NULL;
    // amphLMVec = CreateREAL8Vector(retLen);
    // REAL8Vector *tVec = CreateREAL8Vector(retLen);
    /* Workspace vectors */
    REAL8Vector values, polarDynamics;
    REAL8 valuesdata[14] = {0.};
    REAL8 polarDynamicsdata[4] = {0.};
    values.length = 14;
    polarDynamics.length = 4;
    values.data = valuesdata;
    polarDynamics.data = polarDynamicsdata;
    REAL8 tPeakOmega = seobParams->tPeakOmega;
    COMPLEX16 factor_nqc;
    /* Loop to compute compute amplitude and phase of the hlm mode */
    REAL8 s1dotZ, s2dotZ, chiS, chiA, tplspin, t, omega, ham, v;
    REAL8 hPeak = 0.0, NPeak = 0.0;;
    INT idx_hPeak = 0, idx_NPeak = 0;
    UINT i, j;
    for (i = 0; i < retLen; i++) 
    {
        /* Compute waveform coefficients */
        t = seobdynamics->tVec[i];
        // tVec->data[i] = i*deltaT;
        omega = seobdynamics->omegaVec[i];
        ham = seobdynamics->hamVec[i];
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

        if (CODE_VERSION == 3)
        {
            if (EccPrec_CalcSpinPrecFacWaveformCoefficients(
                    seobParams->hCoeffs, m1, m2, eta, tplspin, chiS, chiA,
                    seobdynamics->chiSxVec[i], seobdynamics->chiSyVec[i], seobdynamics->chiSzVec[i],
                    seobdynamics->chiAxVec[i], seobdynamics->chiAyVec[i], seobdynamics->chiAzVec[i],
                    451) == CEV_FAILURE) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in EccPrec_CalcSpinPrecFacWaveformCoefficients at step %d of the loop.", i);
                return CEV_FAILURE;
            }
        }
        else
        {
            if (XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(
                    seobParams->hCoeffs, m1, m2, eta, tplspin, chiS, chiA,
                    SpinAlignedEOBversionWaveform) == CEV_FAILURE) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients at step %d of the loop.", i);
                return CEV_FAILURE;
            }
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
        polarDynamics.data[2] = seobdynamics->polarprVec[i];
        polarDynamics.data[3] = seobdynamics->polarpphiVec[i];
        v = cbrt(omega);
        if (CODE_VERSION == 0)
        {
            if (XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform(
                    &hLM, &polarDynamics, &values, v, ham, ModeL, ModeM, seobParams) != CEV_SUCCESS) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform at step %d of the loop.", i);
                return CEV_FAILURE;
            }
        } else {
            dr = (seobdynamics->posVecx[i]*seobdynamics->velVecx[i] + 
                seobdynamics->posVecy[i]*seobdynamics->velVecy[i] + 
                seobdynamics->posVecz[i]*seobdynamics->velVecz[i]) / seobdynamics->polarrVec[i];
            ncrv = seobdynamics->polarrVec[i] * omega;

            
            if (XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2(
                    &hLM, &polarDynamics, &values, v, dr, ncrv, seobdynamics->polarprDotVec[i], ham, ModeL, ModeM, seobParams) != CEV_SUCCESS) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2 at step %d of the loop.", i);
                return CEV_FAILURE;
            }
        }

        if (XLALSimIMRSpinEOBNonQCCorrection(&factor_nqc, &values, omega, t, 
                                seobParams->tWind, seobParams->wWind,
                                nqcCoeffs) == CEV_FAILURE) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL,
                "failure in XLALSimIMRSpinEOBNonQCCorrection at step %d of the loop.", i);
            return CEV_FAILURE;
        }

        hLM = hLM * factor_nqc;
        if (cabs(hLM) > hPeak)
        {
            hPeak = cabs(hLM);
            idx_hPeak = i;
        }
        if (cabs(factor_nqc) > NPeak)
        {
            NPeak = cabs(factor_nqc);
            idx_NPeak = i;
        }
    }
    INT idx_tPeak = (INT) (tPeak / deltaT);
    INT idx_max = retLen-1;
    // print_debug("di = %d, dNi = %d\n", idx_hPeak-idx_tPeak, idx_NPeak - retLen+1);
    if (idx_hPeak < idx_tPeak + 10 && idx_NPeak < idx_max + 10)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}

INT CheckStopConditionSA(SEOBSAdynamics *seobdynamics, 
                       SpinEOBParams *seobParams, 
                       SphHarmListEOBNonQCCoeffs *listnqcCoeffs,
                       REAL8 tPeak)
{
    INT ModeL = 2;
    INT ModeM = 2;
    EOBNonQCCoeffs *nqcCoeffs =
        SphHarmListEOBNonQCCoeffs_GetMode(listnqcCoeffs, ModeL, ModeM)->nqcCoeffs;
    /* Check that the input double pointer are not NULL */
    /* Masses */
    REAL8 m1 = seobParams->m1;
    REAL8 m2 = seobParams->m2;
    REAL8 deltaT = seobdynamics->tVec[1] - seobdynamics->tVec[0];
    REAL8 eta = seobParams->eta;
    REAL8 mtot = m1 + m2;
    REAL8 dr, ncrv;
    // UINT SpinAlignedEOBversion = seobParams->SpinAlignedEOBversion;
    UINT SpinAlignedEOBversion = 4;
    // UINT SpinAlignedEOBversionWaveform; // RC: I use this different variable
    //                                     // because the PN terms in the waveform
    //                                     // are different from those in the flux

    /* Length of dynamics data and sampling step */
    UINT retLen = seobdynamics->length;
    COMPLEX16 hLM;
    // REAL8Vector *amphLMVec = NULL;
    // amphLMVec = CreateREAL8Vector(retLen);
    // REAL8Vector *tVec = CreateREAL8Vector(retLen);
    /* Workspace vectors */
    REAL8Vector values, polarDynamics;
    REAL8 valuesdata[14] = {0.};
    REAL8 polarDynamicsdata[4] = {0.};
    values.length = 14;
    polarDynamics.length = 4;
    values.data = valuesdata;
    polarDynamics.data = polarDynamicsdata;
    REAL8 tPeakOmega = seobParams->tPeakOmega;
    COMPLEX16 factor_nqc;
    /* Loop to compute compute amplitude and phase of the hlm mode */
    REAL8 s1dotZ, s2dotZ, chiS, chiA, tplspin, t, omega, ham, v;
    REAL8 hPeak = 0.0, NPeak = 0.0;;
    INT idx_hPeak = 0, idx_NPeak = 0;

    // s1dotZ = seobdynamics->s1dotZVec[i];
    // s2dotZ = seobdynamics->s2dotZVec[i];
    // REAL8 chi1dotZ = s1dotZ * mtot * mtot / (m1 * m1);
    // REAL8 chi2dotZ = s2dotZ * mtot * mtot / (m2 * m2);
    REAL8 chi1dotZ = seobParams->chi1;
    REAL8 chi2dotZ = seobParams->chi2;
    chiS = SEOBCalculateChiS(chi1dotZ, chi2dotZ);
    chiA = SEOBCalculateChiA(chi1dotZ, chi2dotZ);

    tplspin = SEOBCalculatetplspin(m1, m2, eta, s1dotZ, s2dotZ);

    if (XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(
            seobParams->hCoeffs, m1, m2, eta, tplspin, chiS, chiA,
            451) == CEV_FAILURE)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients.");
        return CEV_FAILURE;
    }

    seobParams->hCoeffs->f21v7c = seobParams->cal21;
    seobParams->hCoeffs->f21v7cEff00 = seobParams->cal21E;
    seobParams->hCoeffs->f21v7cEff10 = seobParams->cal21E1;
    seobParams->hCoeffs->f21v7cEff11 = seobParams->cal21E2;
    seobParams->hCoeffs->f21v7cEff01 = seobParams->cal21E3;
    seobParams->hCoeffs->f21v7cEff02 = seobParams->cal21E4;
    seobParams->hCoeffs->f55v5c = seobParams->cal55;

    UINT i, j;
    for (i = 0; i < retLen; i++) 
    {
        /* Compute waveform coefficients */
        t = seobdynamics->tVec[i];
        // tVec->data[i] = i*deltaT;
        omega = seobdynamics->dphiVec[i];
        ham = seobdynamics->HVec[i];

        // print_debug("f21v7c = %.16f\n", seobParams->hCoeffs->f21v7c);
        /* Dynamics, polar dynamics, omega */
        // for (j = 0; j < 14; j++) 
        //     values.data[j] = seobdynamics->array->data[i + (j + 1) * retLen];
        polarDynamics.data[0] = seobdynamics->rVec[i];
        polarDynamics.data[1] = seobdynamics->phiVec[i];
        polarDynamics.data[2] = seobdynamics->prTVec[i];
        polarDynamics.data[3] = seobdynamics->pphiVec[i];
        v = cbrt(omega);
        if (CODE_VERSION == 0)
        {
            if (XLALSimIMRSpinEOBGetSASpinFactorizedWaveform(
                    &hLM, &polarDynamics, v, ham, ModeL, ModeM, seobParams) != CEV_SUCCESS) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform at step %d of the loop.", i);
                return CEV_FAILURE;
            }
        } else {
            // dr = (seobdynamics->posVecx[i]*seobdynamics->velVecx[i] + 
            //     seobdynamics->posVecy[i]*seobdynamics->velVecy[i] + 
            //     seobdynamics->posVecz[i]*seobdynamics->velVecz[i]) / seobdynamics->polarrVec[i];
            dr = seobdynamics->drVec[i];
            ncrv = seobdynamics->rVec[i] * omega;

            
            if (XLALSimIMRSpinEOBGetSASpinFactorizedWaveformV2(
                    &hLM, &polarDynamics, v, dr, ncrv, 0, ham, ModeL, ModeM, seobParams) != CEV_SUCCESS) 
            {
                PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2 at step %d of the loop.", i);
                return CEV_FAILURE;
            }
        }

        if (XLALSimIMRSpinEOBSANonQCCorrection(&factor_nqc, &polarDynamics, omega, t, 
                                seobParams->tWind, seobParams->wWind,
                                nqcCoeffs) == CEV_FAILURE) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL,
                "failure in XLALSimIMRSpinEOBNonQCCorrection at step %d of the loop.", i);
            return CEV_FAILURE;
        }

        hLM = hLM * factor_nqc;
        if (cabs(hLM) > hPeak)
        {
            hPeak = cabs(hLM);
            idx_hPeak = i;
        }
        if (cabs(factor_nqc) > NPeak)
        {
            NPeak = cabs(factor_nqc);
            idx_NPeak = i;
        }
    }
    INT idx_tPeak = (INT) (tPeak / deltaT);
    INT idx_max = retLen-1;
    // print_debug("di = %d, dNi = %d\n", idx_hPeak-idx_tPeak, idx_NPeak - retLen+1);
    if (idx_hPeak < idx_tPeak + 10 && idx_NPeak < idx_max + 10)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}


/*--------------------------------------------------------------*/
/*                                                              */
/*                                                              */
/*                                                              */
/*                           egw init                           */
/*                                                              */
/*                                                              */
/*                                                              */
/*--------------------------------------------------------------*/

static int XLALEOBSpinPrecAlignedStopCondition_egw(
    double t,       /**< UNUSED */
    const double values[], /**< dynamical variable values */
    double dvalues[],      /**< dynamical variable time derivative values */
    void *funcParams       /**< physical parameters */
) 
{
    REAL8 omega, r;
    SpinEOBParams *params = (SpinEOBParams *)funcParams;

    if (values[1] > 3.*CST_PI) 
    {
        return 1;
    }
    return GSL_SUCCESS;
}

INT SEOBIntegrateDynamics_egw(REAL8Array **dynamics,
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
    UINT nb_Hamiltonian_variables_spinsaligned = 4;
    // print_debug("tend = %.16e\n", tend);
    REAL8Vector *values = CreateREAL8Vector(nb_Hamiltonian_variables);
    if (!values) {failed = 1; goto QUIT;}
    memcpy(values->data, ICvalues->data, values->length * sizeof(REAL8));

    REAL8Vector *values_spinaligned = CreateREAL8Vector(nb_Hamiltonian_variables_spinsaligned);
    if (!values_spinaligned) {failed = 1; goto QUIT;}
    memset(values_spinaligned->data, 0, values_spinaligned->length * sizeof(REAL8));

    ARKIntegrator *integrator = NULL;
    if (seobParams->alignedSpins)
    {
        // Spin Aligned
        REAL8 temp_r = sqrt(ICvalues->data[0] * ICvalues->data[0] +
                            ICvalues->data[1] * ICvalues->data[1] +
                            ICvalues->data[2] * ICvalues->data[2]);
        REAL8 temp_phi = ICvalues->data[12];

        values_spinaligned->data[0] = temp_r;   // General form of r
        values_spinaligned->data[1] = temp_phi; // phi
        values_spinaligned->data[2] = 
            ICvalues->data[3] * cos(temp_phi) +
            ICvalues->data[4] * sin(temp_phi); // p_r^*
        values_spinaligned->data[3] =
            temp_r * (ICvalues->data[4] * cos(temp_phi) -
                    ICvalues->data[3] * sin(temp_phi)); // p_phi

        /* We have to use different stopping conditions depending
        we are in the low-sampling or high-sampling portion
        of the waveform. We can tell this apart because for the low-sampling (or
        ada sampling) we always start at t=0
        */
        integrator = XLALAdaptiveRungeKutta4Init(
            nb_Hamiltonian_variables_spinsaligned, XLALSpinAlignedHcapDerivative_SA,
            XLALEOBSpinPrecAlignedStopCondition_egw, EPS_ABS, EPS_REL);
    }
    else
    {
        // Prec
        // print_debug("prec...");
        integrator = XLALAdaptiveRungeKutta4Init(
            nb_Hamiltonian_variables, XLALSpinPrecHcapNumericalDerivative,
            XLALEOBSpinPrecStopConditionBasedOnPR, EPS_ABS, EPS_REL);
    }

    if (!integrator) {failed = 1; goto QUIT;}

    /* The integration stops when stop condition reached */
    integrator->stopontestonly = 1;

    /* When this option is set to 0, the integration can be exceedingly slow for
    * spin-aligned systems */
    integrator->retries = 1;
    /* Computing the dynamical evolution of the system */
    if (seobParams->alignedSpins)
    {
        // Spin Aligned
        // flagConstantSampling = seobParams->alignedSpins
        // EOBversion = 2
        PRINT_LOG_INFO(LOG_INFO, "Adaptive RungeKutta4 No Interpolation");
        retLen = XLALAdaptiveRungeKutta4NoInterpolate(
            integrator, seobParams, values_spinaligned->data, 0., 10.,
            deltaT, deltaT_min, &dynamics_spinaligned);

        if (retLen < 0)
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Integration Failed!!!");
            failed = 1;
            goto QUIT;
        }
        /* Convert the spin-aligned dynamics to a generic-spins dynamics */
        PRINT_LOG_INFO(LOG_INFO, "Convert Spin Aligned Dynamics to Generic Spins");
        status = SEOBConvertSpinAlignedDynamicsToGenericSpins(
            dynamics, dynamics_spinaligned, retLen, 
            seobParams->chi1, seobParams->chi2, seobParams);
        if (status != CEV_SUCCESS)
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Failure in dynamics transformation.");
            failed = 1;
            goto QUIT;
        }
    } 
    else 
    {
        // Prec
        retLen = XLALAdaptiveRungeKutta4NoInterpolate(integrator, seobParams, 
            values->data, 0., 10., 
            deltaT, deltaT_min, dynamics);
        if (retLen < 0)
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "Integration Failed!!!");
            failed = 1;
            goto QUIT;
        }
    }
    PRINT_LOG_INFO(LOG_INFO, "Integration End");
    // NOTE: functions like XLALAdaptiveRungeKutta4 would give nans if the times
    // do not start at 0 -- we have to adjust the starting time after integration
    /* Adjust starting time */
    for ( i = 0; i < retLen; i++)
        (*dynamics)->data[i] += tstart;

QUIT:
    STRUCTFREE(values, REAL8Vector);
    STRUCTFREE(values_spinaligned, REAL8Vector);
    STRUCTFREE(integrator, ARKIntegrator);
    STRUCTFREE(dynamics_spinaligned, REAL8Array);
    *retLenOut = retLen;
    if (failed)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}

static REAL8 find_dOmega0_sol_egw(gsl_spline *spline, gsl_interp_accel *acc, REAL8 t1, REAL8 t2)
{
    REAL8 ret, EPS = 1e-9;
    REAL8 tMid, dO, dO1, dO2;
    dO1 = gsl_spline_eval_deriv2(spline, t1, acc);
    if (fabs(dO1) < EPS)
        return gsl_spline_eval_deriv(spline, t1, acc);
    dO2 = gsl_spline_eval_deriv2(spline, t2, acc);
    if (fabs(dO2) < EPS)
        return gsl_spline_eval_deriv(spline, t2, acc);
    int i, imax;
    imax = 100;
    while(i < imax)
    {
        tMid = (t1 + t2)/2.;
        dO = gsl_spline_eval_deriv2(spline, tMid, acc);
        if (fabs(dO) < EPS)
            return gsl_spline_eval_deriv(spline, tMid, acc);
        if (dO * dO1 < 0)
            t2 = tMid;
        else
            t1 = tMid;
    }
    return gsl_spline_eval_deriv(spline, tMid, acc);
}

INT CalculateSAh22SeriesFromrpphi(REAL8 r, REAL8 pphi, SpinEOBParams *core,
    REAL8 *omega22, REAL8 *egw)
{
    // print_debug( "Values r = %.16e, pphi = %.16e\n", r, pphi );
    REAL8 EPS_REL = 1.0e-9;
    REAL8 EPS_ABS = 1.0e-10;
    INT retLenAdaS, failed = 0, status;
    SEOBdynamics *seobdynamicsAdaS = NULL;
    REAL8Array *dynamicsAdaS = NULL;
    REAL8Vector *ICvalues = NULL;
    COMPLEX16 h22;
    REAL8Vector *h22Phase = NULL;
    REAL8Vector *h22dOmega = NULL;
    gsl_spline *spline = NULL;
    gsl_interp_accel *acc = NULL;
    ICvalues = CreateREAL8Vector(14);

    // Set initial Conditions
    REAL8 xSph[3] = {r, 0., 0.};
    REAL8 pSph[3] = {0, 0, pphi};
    REAL8 xCart[3] = {0,0,0};
    REAL8 pCart[3] = {0,0,0};
    SphericalToCartesian(xCart, pCart, xSph, pSph);
    memcpy( ICvalues->data, xCart, sizeof(xCart) );
    memcpy( ICvalues->data+3, pCart, sizeof(pCart) );
    memcpy( ICvalues->data+6, core->s1Vec->data, sizeof(pCart) );
    memcpy( ICvalues->data+9, core->s2Vec->data, sizeof(pCart) );

    REAL8 deltaT, deltaT_min;
    REAL8 tthresh;
    deltaT = 0.5;
    deltaT_min = 8.0e-5;
    status = SEOBIntegrateDynamics_egw(&dynamicsAdaS, &retLenAdaS, 
        ICvalues, EPS_ABS, EPS_REL, deltaT, deltaT_min, 
        0, 10., core, core->alignedSpins);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}

    status = SEOBComputeExtendedSEOBdynamics_Conserve(&seobdynamicsAdaS, dynamicsAdaS, retLenAdaS, core);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}

    h22Phase = CreateREAL8Vector(retLenAdaS);
    h22dOmega = CreateREAL8Vector(retLenAdaS);
    REAL8 m1 = core->m1;
    REAL8 m2 = core->m2;
    REAL8 mtot = m1 + m2;
    REAL8 dr, ncrv, omega, ham, s1dotZ, s2dotZ;
    REAL8Vector values, polarDynamics;
    REAL8 valuesdata[14] = {0.};
    REAL8 polarDynamicsdata[4] = {0.};
    values.length = 14;
    polarDynamics.length = 4;
    values.data = valuesdata;
    polarDynamics.data = polarDynamicsdata;
    for (UINT i=0; i<retLenAdaS; i++)
    {
        omega = seobdynamicsAdaS->omegaVec[i];
        ham = seobdynamicsAdaS->hamVec[i];
        s1dotZ = seobdynamicsAdaS->s1dotZVec[i];
        s2dotZ = seobdynamicsAdaS->s2dotZVec[i];
        REAL8 chi1dotZ = s1dotZ * mtot * mtot / (m1 * m1);
        REAL8 chi2dotZ = s2dotZ * mtot * mtot / (m2 * m2);
        REAL8 chiS = SEOBCalculateChiS(chi1dotZ, chi2dotZ);
        REAL8 chiA = SEOBCalculateChiA(chi1dotZ, chi2dotZ);
        dr = (seobdynamicsAdaS->posVecx[i]*seobdynamicsAdaS->velVecx[i] + 
            seobdynamicsAdaS->posVecy[i]*seobdynamicsAdaS->velVecy[i] + 
            seobdynamicsAdaS->posVecz[i]*seobdynamicsAdaS->velVecz[i]) / seobdynamicsAdaS->polarrVec[i];
        ncrv = seobdynamicsAdaS->polarrVec[i] * omega;
        for (UINT j = 0; j < 14; j++) 
            values.data[j] = seobdynamicsAdaS->array->data[i + (j + 1) * retLenAdaS];
        polarDynamics.data[0] = seobdynamicsAdaS->polarrVec[i];
        polarDynamics.data[1] = seobdynamicsAdaS->polarphiVec[i];
        polarDynamics.data[2] = seobdynamicsAdaS->polarprVec[i];
        polarDynamics.data[3] = seobdynamicsAdaS->polarpphiVec[i];
        REAL8 v = cbrt(omega);
        if (XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2(
                &h22, &polarDynamics, &values, v, dr, ncrv, seobdynamicsAdaS->polarprDotVec[i], ham, 2, 2, core) != CEV_SUCCESS) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2 at step %d of the loop.", i);
            return CEV_FAILURE;
        }
        h22Phase->data[i] = carg(h22);
    }
    XLALREAL8VectorUnwrapAngle(h22Phase, h22Phase);
    REAL8 phi0 = h22Phase->data[0];
    for (int i = 0; i<retLenAdaS; i++)
    {
        h22Phase->data[i] = fabs(h22Phase->data[i] - phi0);
    }
    spline = gsl_spline_alloc(gsl_interp_cspline, retLenAdaS);
    acc = gsl_interp_accel_alloc();
    gsl_spline_init(spline, seobdynamicsAdaS->tVec, h22Phase->data, retLenAdaS);
    for (int i = 0; i<retLenAdaS; i++)
    {
        // print_out("%.16e\t%.16e\t%.16e\n", 
        //     seobdynamicsAdaS->tVec[i],
        //     seobdynamicsAdaS->polarrVec[i],
        //     h22Phase->data[i]);
        // h22Omega->data[i] = gsl_spline_eval_deriv(spline, seobdynamicsAdaS->tVec[i], acc);
        h22dOmega->data[i] = gsl_spline_eval_deriv2(spline, seobdynamicsAdaS->tVec[i], acc);
    }
    REAL8 omg22a, omg22p;
    REAL8 sqomg22a, sqomg22p, e22, Psi;
    int found0 = 0, ia, ip;
    for (int i = 0; i<retLenAdaS-1; i++)
    {
        if (h22dOmega->data[i] * h22dOmega->data[i+1] < 0)
        {
            if (!found0)
            {
                ia = i;
                found0++;
            } else {
                ip = i;
                break;
            }
        }
    }
    omg22a = find_dOmega0_sol_egw(spline, acc, seobdynamicsAdaS->tVec[ia], seobdynamicsAdaS->tVec[ia+1]);
    omg22p = find_dOmega0_sol_egw(spline, acc, seobdynamicsAdaS->tVec[ip], seobdynamicsAdaS->tVec[ip+1]);
    sqomg22a = sqrt(omg22a);
    sqomg22p = sqrt(omg22p);
    e22 = fabs(sqomg22a - sqomg22p) / (sqomg22a + sqomg22p);
    Psi = atan2(1. - e22*e22, 2.*e22);
    *omega22 = omg22a < omg22p ? omg22a : omg22p;
    *egw = cos(Psi/3.) - CST_SQ3*sin(Psi/3.);
    // print_debug("e22 = %.16e, egw = %.16e\nomg22a = %.16e, omg22p = %.16e, omg22In = %.16e\n",
    //     e22, *egw, omg22a, omg22p, CST_2PI*0.002);
    // print_debug("egw = %.16e\n", *egw);
QUIT:
    STRUCTFREE(dynamicsAdaS, REAL8Array);
    STRUCTFREE(seobdynamicsAdaS, SEOBdynamics);
    STRUCTFREE(ICvalues, REAL8Vector);
    STRUCTFREE(h22Phase, REAL8Vector);
    STRUCTFREE(h22dOmega, REAL8Vector);
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    if (failed)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}


INT SEOBInitialConditions_egw(REAL8Vector *ICvalues,
                            REAL8 MfMin,
                            REAL8 ecc,
                            SpinEOBParams *seobParams)
{
    PRINT_LOG_INFO(LOG_INFO, "Set initial conditions");
    if (!seobParams)
        return CEV_FAILURE;
    INT j;
    memset((ICvalues)->data, 0, ((ICvalues)->length) * sizeof(REAL8));
    REAL8 eta = seobParams->eta;
    REAL8 m1, m2, mTotal;
    m1 = seobParams->m1;
    m2 = seobParams->m2;
    mTotal = m1 + m2;
    REAL8 fMin = MfMin / (m1 + m2) / CST_MTSUN_SI;
    REAL8 mSpin1data[3] = {0., 0., 0.};
    REAL8 mSpin2data[3] = {0., 0., 0.};

    if (seobParams->alignedSpins)
    {
        REAL8 chi1dotZ = seobParams->chi1;
        REAL8 chi2dotZ = seobParams->chi2;
        REAL8 chiS = SEOBCalculateChiS(chi1dotZ, chi2dotZ);
        REAL8 chiA = SEOBCalculateChiA(chi1dotZ, chi2dotZ);
        REAL8 tplspin = SEOBCalculatetplspin(m1, m2, eta, chi1dotZ, chi2dotZ);
        if (XLALSimIMREOBCalcSpinFacWaveformCoefficients(seobParams->hCoeffs, seobParams, tplspin, chiS, chiA) != CEV_SUCCESS)
            return CEV_FAILURE;
        mSpin1data[2] = chi1dotZ * m1 * m1;
        mSpin2data[2] = chi2dotZ * m2 * m2;
        // print_debug("m1 = %g, m2 = %g, fMin = %g\n", m1, m2, fMin);
        // print_debug("mSpin1data = (%g, %g, %g)\n", mSpin1data[0], mSpin1data[1], mSpin1data[2]);
        // print_debug("mSpin2data = (%g, %g, %g)\n", mSpin2data[0], mSpin2data[1], mSpin2data[2]);
        if (EOBInitialConditionsSA_egw(ICvalues, m1, m2, fMin, ecc, 0, mSpin1data, mSpin2data, seobParams) != CEV_SUCCESS)
            return CEV_FAILURE;
    }
    else
    {
        for (j=0; j<3; j++)
        {
            mSpin1data[j] = seobParams->s1Vec->data[j] * mTotal * mTotal;
            mSpin2data[j] = seobParams->s2Vec->data[j] * mTotal * mTotal;
        }
        if (EOBInitialConditionsPrec(ICvalues, m1, m2, fMin, ecc, 0, mSpin1data, mSpin2data, seobParams) != CEV_SUCCESS)
            return CEV_FAILURE;
    }
    return CEV_SUCCESS;
}

#include <gsl/gsl_roots.h>
REAL8 calc_egw_from_e22(REAL8 e22)
{
    REAL8 Psi = atan2(1. - e22*e22, 2.*e22);
    return cos(Psi) - CST_SQ3*sin(Psi);
}

static double GSL_calc_egw_from_e22(const double x, void *params)
{
    REAL8 egw = *((REAL8*)params);
    REAL8 Psi = atan2(1. - x*x, 2.*x);
    return cos(Psi) - CST_SQ3*sin(Psi) - egw;
}

REAL8 calc_e22_from_egw(REAL8 egw)
{
    REAL8 e22;
    REAL8 egw_target = egw;
    gsl_function F;
    F.function = &GSL_calc_egw_from_e22;
    F.params = &egw_target;
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc (T);
    int iter = 0, max_iter = 100;
    REAL8 x_lo, x_hi, root;
    x_lo = 0.0;
    x_hi = 0.99;
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
    int status;
    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        root = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi, 1e-9, 1e-9);
    } while (status == GSL_CONTINUE && iter < max_iter);
    if ( iter > max_iter && status != GSL_SUCCESS )
    {
        gsl_root_fsolver_free( s );
        return CEV_FAILURE;
    }
    e22 = root;
    gsl_root_fsolver_free (s);
    return e22;
}

typedef struct {
    REAL8 rp;
    REAL8 rm;
    // SpinEOBParams *pms;
    HcapDerivParams *pms;
}DiffHamParams;

static double diffHamiltonianByrpm(const double x, /**<< Parameters requested by gsl root finder */
                       void *params        /**<< Spin EOB parameters */)
{
    DiffHamParams *pms = (DiffHamParams*)params;
    HcapDerivParams *hpms = pms->pms;
    REAL8         cartValues[6];
    hpms->values = cartValues;
    memset( cartValues, 0, sizeof( cartValues ) );

    REAL8 rp, rm;
    rp = pms->rp;
    rm = pms->rm;
    REAL8 Hp, Hm;
    cartValues[4] = x/rp; // pphi / r
    Hp = GSLSpinAlignedHamiltonianWrapper_SA(rp, hpms);
    cartValues[4] = x/rm; // pphi / r
    Hm = GSLSpinAlignedHamiltonianWrapper_SA(rm, hpms);
    return Hp-Hm;
}

INT find_SApphi_from_rpm(REAL8 rp, REAL8 rm, SpinEOBParams *core, REAL8 *pphi)
{
    print_debug("rp = %.16e, rm = %.16e\n", rp, rm);
    DiffHamParams rootpms;
    HcapDerivParams hcdpms;
    hcdpms.params  = core;
    hcdpms.varyParam = 0;
    rootpms.rp = rp;
    rootpms.rm = rm;
    rootpms.pms = (void*)&hcdpms;
    gsl_function F;
    F.function = &diffHamiltonianByrpm;
    F.params = &rootpms;
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc (T);
    int iter = 0, max_iter = 100;
    REAL8 x_lo, x_hi, root;
    x_lo = 0.0;
    x_hi = 3.*sqrt(2.*rp*rm / (rp + rm));
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
    int status;
    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        root = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi, 1e-8, 1e-7);
    } while (status == GSL_CONTINUE && iter < max_iter);
    if ( iter > max_iter && status != GSL_SUCCESS )
    {
        gsl_root_fsolver_free( s );
        return CEV_FAILURE;
    }
    *pphi = root;
    gsl_root_fsolver_free (s);
    return CEV_SUCCESS;
}

typedef struct {
    REAL8 *values; // len = 12
    SpinEOBParams *params;
}Rderiv2RootParams;

static double Rderiv2HamiltonianByr(const double x, /**<< Parameters requested by gsl root finder */
                       void *params        /**<< Spin EOB parameters */)
{
    Rderiv2RootParams *hpms = (Rderiv2RootParams*)params;
    REAL8 dHdr;
    hpms->values[4] = x / hpms->values[0]; // pphi
    dHdr = XLALSpinHcapNumDerivWRTParam( 0, hpms->values, hpms->params );
    return dHdr;
}

INT find_SACircpphi_from_r(REAL8 r, SpinEOBParams *core, REAL8 *pphi)
{
    Rderiv2RootParams rootpms;
    REAL8         sphValues[12];
    rootpms.values = sphValues;
    memset( sphValues, 0, sizeof( sphValues ) );
    sphValues[0] = r;
    for (int i=0; i<3; i++)
    {
        sphValues[i+6] = core->s1Vec->data[i] * core->m1 * core->m1;
        sphValues[i+9] = core->s2Vec->data[i] * core->m2 * core->m2;
    }
    rootpms.params  = core;

    gsl_function F;
    F.function = &Rderiv2HamiltonianByr;
    F.params = &rootpms;
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc (T);
    int iter = 0, max_iter = 100;
    REAL8 x_lo, x_hi, root;
    x_lo = 0.5*sqrt(r);
    x_hi = 3.*sqrt(r);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
    int status;
    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        root = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi, 1e-8, 1e-7);
    } while (status == GSL_CONTINUE && iter < max_iter);
    if ( iter > max_iter && status != GSL_SUCCESS )
    {
        gsl_root_fsolver_free( s );
        return CEV_FAILURE;
    }
    *pphi = root;
    gsl_root_fsolver_free (s);
    return CEV_SUCCESS;
}


static int XLALEOBSpinPrecAlignedStopCondition_egwPlus(
    double t,       /**< UNUSED */
    const double values[], /**< dynamical variable values */
    double dvalues[],      /**< dynamical variable time derivative values */
    void *funcParams       /**< physical parameters */
) 
{
    REAL8 omega, r;
    SpinEOBParams *params = (SpinEOBParams *)funcParams;

    if (values[1] > CST_PI/5.) 
    {
        return 1;
    }
    return GSL_SUCCESS;
}

static int XLALEOBSpinPrecAlignedStopCondition_egwMinus(
    double t,       /**< UNUSED */
    const double values[], /**< dynamical variable values */
    double dvalues[],      /**< dynamical variable time derivative values */
    void *funcParams       /**< physical parameters */
) 
{
    REAL8 omega, r;
    SpinEOBParams *params = (SpinEOBParams *)funcParams;

    if (values[1] < -CST_PI/5.) 
    {
        return 1;
    }
    return GSL_SUCCESS;
}

INT EvaluateOmega22SA_form_rpphi(REAL8 r, REAL8 pphi, SpinEOBParams *seobParams, REAL8 *omega22)
{
    INT retLenPlus, retLenMinus;
    UINT i;
    REAL8 EPS_ABS, EPS_REL;
    EPS_ABS = 1e-10;
    EPS_REL = 1e-10;
    REAL8 deltaT = 0.5, deltaT_min = 0.0001;
    REAL8Array *dynamics_Plus = NULL;
    REAL8Array *dynamics_Minus = NULL;
    REAL8Array *dynamics = NULL;
    SEOBdynamics *seobdynamicsAdaS = NULL;
    REAL8Vector *h22Phase = NULL;
    INT status, failed = 0;
    gsl_spline *spline = NULL;
    gsl_interp_accel *acc = NULL;
    /* Dimensions of vectors of dynamical variables to be integrated */
    UINT nb_Hamiltonian_variables_spinsaligned = 4;

    REAL8Vector *values_spinaligned = CreateREAL8Vector(nb_Hamiltonian_variables_spinsaligned);
    if (!values_spinaligned) {failed = 1; goto QUIT;}
    REAL8 initvals[4];
    memset(values_spinaligned->data, 0, values_spinaligned->length * sizeof(REAL8));

    ARKIntegrator *integrator = NULL;
    // Spin Aligned
    values_spinaligned->data[0] = r;   // General form of r
    values_spinaligned->data[1] = 0.0; // phi
    values_spinaligned->data[2] = 0.0; // p_r^*
    values_spinaligned->data[3] = pphi; // p_phi

    // Minus
    memcpy(initvals, values_spinaligned->data, 4*sizeof(REAL8));
    integrator = XLALAdaptiveRungeKutta4Init(
        nb_Hamiltonian_variables_spinsaligned, XLALSpinAlignedHcapDerivative_invConserve,
        XLALEOBSpinPrecAlignedStopCondition_egwMinus, EPS_ABS, EPS_REL);
    if (!integrator) {failed = 1; goto QUIT;}
    /* The integration stops when stop condition reached */
    integrator->stopontestonly = 1;
    integrator->retries = 1;
    retLenMinus = XLALAdaptiveRungeKutta4NoInterpolate(
        integrator, seobParams, initvals, 0., 10.,
        deltaT, deltaT_min, &dynamics_Minus);
    if (retLenMinus < 0)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Integration Failed!!!");
        failed = 1;
        goto QUIT;
    }
    REAL8 tAtPeak;
    tAtPeak = dynamics_Minus->data[retLenMinus-1];
    // plus
    integrator = XLALAdaptiveRungeKutta4Init(
        nb_Hamiltonian_variables_spinsaligned, XLALSpinAlignedHcapDerivative_Conserve,
        XLALEOBSpinPrecAlignedStopCondition_egwPlus, EPS_ABS, EPS_REL);
    if (!integrator) {failed = 1; goto QUIT;}
    /* The integration stops when stop condition reached */
    integrator->stopontestonly = 1;
    integrator->retries = 1;
    retLenPlus = XLALAdaptiveRungeKutta4NoInterpolate(
        integrator, seobParams, initvals, 0., 10.,
        deltaT, deltaT_min, &dynamics_Plus);
    if (retLenPlus < 0)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Integration Failed!!!");
        failed = 1;
        goto QUIT;
    }

    status = SEOBConvertSpinAlignedDynamicsToGenericSpins(
        &dynamics, dynamics_Plus, retLenPlus, 
        seobParams->chi1, seobParams->chi2, seobParams);
    if (status != CEV_SUCCESS)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Failure in dynamics transformation.");
        failed = 1;
        goto QUIT;
    }

    status = SEOBComputeExtendedSEOBdynamics_Conserve(&seobdynamicsAdaS, dynamics, retLenPlus, seobParams);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}

    h22Phase = CreateREAL8Vector(retLenPlus);
    COMPLEX16 h22;
    REAL8Vector values, polarDynamics;
    REAL8 valuesdata[14] = {0.};
    REAL8 polarDynamicsdata[4] = {0.};
    values.length = 14;
    polarDynamics.length = 4;
    values.data = valuesdata;
    polarDynamics.data = polarDynamicsdata;
    REAL8 m1 = seobParams->m1;
    REAL8 m2 = seobParams->m2;
    REAL8 mtot = m1 + m2;
    for (i=0; i<retLenPlus; i++)
    {
        REAL8 omega = seobdynamicsAdaS->omegaVec[i];
        REAL8 ham = seobdynamicsAdaS->hamVec[i];
        REAL8 s1dotZ = seobdynamicsAdaS->s1dotZVec[i];
        REAL8 s2dotZ = seobdynamicsAdaS->s2dotZVec[i];
        REAL8 chi1dotZ = s1dotZ * mtot * mtot / (m1 * m1);
        REAL8 chi2dotZ = s2dotZ * mtot * mtot / (m2 * m2);
        REAL8 chiS = SEOBCalculateChiS(chi1dotZ, chi2dotZ);
        REAL8 chiA = SEOBCalculateChiA(chi1dotZ, chi2dotZ);
        REAL8 dr = (seobdynamicsAdaS->posVecx[i]*seobdynamicsAdaS->velVecx[i] + 
            seobdynamicsAdaS->posVecy[i]*seobdynamicsAdaS->velVecy[i] + 
            seobdynamicsAdaS->posVecz[i]*seobdynamicsAdaS->velVecz[i]) / seobdynamicsAdaS->polarrVec[i];
        REAL8 ncrv = seobdynamicsAdaS->polarrVec[i] * omega;
        for (UINT j = 0; j < 14; j++) 
            values.data[j] = seobdynamicsAdaS->array->data[i + (j + 1) * retLenPlus];
        polarDynamics.data[0] = seobdynamicsAdaS->polarrVec[i];
        polarDynamics.data[1] = seobdynamicsAdaS->polarphiVec[i];
        polarDynamics.data[2] = seobdynamicsAdaS->polarprVec[i];
        polarDynamics.data[3] = seobdynamicsAdaS->polarpphiVec[i];
        REAL8 v = cbrt(omega);
        if (XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2(
                &h22, &polarDynamics, &values, v, dr, ncrv, seobdynamicsAdaS->polarprDotVec[i], ham, 2, 2, seobParams) != CEV_SUCCESS) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveformV2 at step %d of the loop.", i);
            return CEV_FAILURE;
        }
        h22Phase->data[i] = carg(h22);
    }

    XLALREAL8VectorUnwrapAngle(h22Phase, h22Phase);
    REAL8 phi0 = h22Phase->data[0];
    for (i = 0; i<retLenPlus; i++)
    {
        h22Phase->data[i] = fabs(h22Phase->data[i] - phi0);
    }
    spline = gsl_spline_alloc(gsl_interp_cspline, retLenPlus);
    acc = gsl_interp_accel_alloc();
    gsl_spline_init(spline, seobdynamicsAdaS->tVec, h22Phase->data, retLenPlus);
    for (i = 0; i<retLenPlus; i++)
    {
        print_out("%.16e\t%.16e\t%.16e\n", 
            seobdynamicsAdaS->tVec[i],
            seobdynamicsAdaS->omegaVec[i],
            h22Phase->data[i]);
        // h22Omega->data[i] = gsl_spline_eval_deriv(spline, seobdynamicsAdaS->tVec[i], acc);
    }
    *omega22 = gsl_spline_eval_deriv(spline, tAtPeak, acc);

QUIT:
    STRUCTFREE(values_spinaligned, REAL8Vector);
    STRUCTFREE(h22Phase, REAL8Vector);
    STRUCTFREE(integrator, ARKIntegrator);
    STRUCTFREE(dynamics_Plus, REAL8Array);
    STRUCTFREE(dynamics_Minus, REAL8Array);
    STRUCTFREE(dynamics, REAL8Array);
    STRUCTFREE(seobdynamicsAdaS, SEOBdynamics);
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    if (failed)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}

#if 0

INT CalculateEnergyFluxFromDynamics(SEOBdynamics *dyn,
                                    SpinEOBParams *seobpms,
                                    REAL8Vector *fluxVecOut)
{
    REAL8 flux;
    REAL8Vector *fluxVec;
    fluxVec = CreateREAL8Vector(dyn->length);
    memset(fluxVec->data, 0, fluxVec->length * sizeof(REAL8));
	REAL8Vector	polarDynamics, cartDynamics;
	REAL8		polData[4];
    REAL8       values[12], dvalues[12];
    polarDynamics.length = 4;
    REAL8 mass1 = seobpms->m1;
    REAL8 mass2 = seobpms->m2;
    REAL8 r, dr, ham;
    REAL8 rCrossV_x, rCrossV_y, rCrossV_z, omega;
    INT i;
    for (i=0; i<dyn->length; i++)
    {

        polarDynamics.data = polData;
        cartDynamics.data = values;
        values[0] = dyn->posVecx[i];
        values[1] = dyn->posVecy[i];
        values[2] = dyn->posVecz[i];
        values[3] = dyn->momVecx[i];
        values[4] = dyn->momVecy[i];
        values[5] = dyn->momVecz[i];
        values[6] = dyn->s1Vecx[i];
        values[7] = dyn->s1Vecy[i];
        values[8] = dyn->s1Vecz[i];
        values[9] = dyn->s2Vecx[i];
        values[10] = dyn->s2Vecy[i];
        values[11] = dyn->s2Vecz[i];
        dvalues[0] = dyn->velVecx[i];
        dvalues[1] = dyn->velVecy[i];
        dvalues[2] = dyn->velVecz[i];
        r = poldata[0] = dyn->polarrVec[i];
        poldata[1] = dyn->polarphiVec[i];
        poldata[2] = dyn->polarprVec[i];
        poldata[3] = dyn->polarpphiVec[i];     
        rCrossV_x = values[1] * dvalues[2] - values[2] * dvalues[1];
        rCrossV_y = values[2] * dvalues[0] - values[0] * dvalues[2];
        rCrossV_z = values[0] * dvalues[1] - values[1] * dvalues[0];

        omega = sqrt(rCrossV_x * rCrossV_x + rCrossV_y * rCrossV_y + rCrossV_z * rCrossV_z) / (r * r);
        dr = (values[0]*dvalues[0] + values[1]*dvalues[1] + values[2]*dvalues[2]) / r;
        ncrv = omega * r;

        flux = XLALInspiralPrecSpinFactorizedFlux(&polarDynamics, &cartDynamics,
                        nqcCoeffs, omega, dr, ncrv, params.params, ham / (mass1 + mass2), 8, 4);
    }
    return CEV_SUCCESS;
}
#endif