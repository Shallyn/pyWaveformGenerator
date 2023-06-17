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
#include "pDebug.h"

int calc_rpphi_from_e(REAL8 e, REAL8 omega, SpinEOBParams *params, REAL8 *r, REAL8 *pphi);

static INT debug_InitialCondition_egw(REAL8 m1,  REAL8 m2, 
        REAL8 s1x, REAL8 s1y, REAL8 s1z, 
        REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 phi0, REAL8 distance,
        REAL8 ecc, REAL8 f_min, REAL8 INdeltaT, REAL8 inc, HyperParams *hparams)
{
    PRINT_LOG_INFO(LOG_INFO, "DEBUG mode 1: test egw initial conditions");
    register INT i;
    INT failed = 0, this_step = 0, status = CEV_FAILURE;
    SpinEOBParams *core = NULL;
    REAL8Vector *ICvalues = NULL;
    REAL8Array *dynamicsAdaS = NULL;
    SEOBdynamics *seobdynamicsAdaS = NULL;

    REAL8 mTotal = m1 + m2;
    REAL8 mTScaled = mTotal * CST_MTSUN_SI;
    REAL8 EPS_ALIGN = 1.0e-4;
    REAL8 deltaT = INdeltaT / mTScaled;
    REAL8 tStepBack;

    if (m1 < 0 || m2 < 0) {failed = 1; goto QUIT;}
    if (m2 > m1)
        SWAP(m1, m2);

    core = CreateSpinEOBParams(m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, ecc, hparams);
    if (!core) {failed = 1; goto QUIT;}
    tStepBack = GET_MAX(tStepBack, core->hParams->tStepBack);

    ICvalues = CreateREAL8Vector(14);
    REAL8 MfMin = mTScaled * f_min;

    status = SEOBInitialConditions_egw(ICvalues, MfMin, ecc, core);
    PRINT_LOG_INFO(LOG_DEBUG, "initial conditions:");
    PRINT_LOG_INFO(LOG_DEBUG, "(x,y,z) = (%.16e %.16e %.16e)", 
            ICvalues->data[0], ICvalues->data[1], ICvalues->data[2]);
    PRINT_LOG_INFO(LOG_DEBUG, "(px,py,pz) = (%.16e %.16e %.16e)",
            ICvalues->data[3], ICvalues->data[4], ICvalues->data[5]);
    PRINT_LOG_INFO(LOG_DEBUG, "(S1x,S1y,S1z) = (%.16e %.16e %.16e)",
            ICvalues->data[6], ICvalues->data[7], ICvalues->data[8]);
    PRINT_LOG_INFO(LOG_DEBUG, "(S2x,S2y,S2z) = (%.16e %.16e %.16e)",
            ICvalues->data[9], ICvalues->data[10], ICvalues->data[11]);
    PRINT_LOG_INFO(LOG_DEBUG, "MfMin = %f, deltaT = %f\n", MfMin, deltaT);
    PRINT_LOG_INFO(LOG_DEBUG, "e0 = %f\n", ecc);

    INT retLenAdaS = 0;
    REAL8 tendAdaS = 20. / mTScaled;    
    REAL8 tstartAdaS = 0.;
    REAL8 deltaT_min = 8.0e-5;
    REAL8 EPS_ABS = 1.0e-8;
    REAL8 EPS_REL = 1.0e-8;
    if (core->alignedSpins)
    {
        EPS_REL = 1.0e-9;
        EPS_ABS = 1.0e-10;
    }
    REAL8 omega22, egw;
    print_debug("r =  %.16e, pphi = %.16e\n", ICvalues->data[0], ICvalues->data[0]*ICvalues->data[4]);
    CalculateSAh22SeriesFromrpphi(ICvalues->data[0], ICvalues->data[0]*ICvalues->data[4], core, &omega22, &egw);
    print_debug("egw =  %.16e\n", egw);
    // REAL8 pphi, pphiC;
    // REAL8 rm = ICvalues->data[0], rp;
    // rp = (rm + ecc*rm) / (1. - ecc);
    // print_debug("rp = %.16e, rm = %.16e\n", rp, rm);
    // find_SApphi_from_rpm(rp, rm, core, &pphi);
    // find_SACircpphi_from_r(rm, core, &pphiC);
    // print_debug("pphi = %.16e, pphiC = %.16e\n", pphi, pphiC);
    // // status = CalculateSAh22SeriesFromrpphi(rm, pphi, core, &omega22, &egw);
    // // if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    // // PRINT_LOG_INFO(LOG_DEBUG, "AdaS data length = %d", retLenAdaS);
    // EvaluateOmega22SA_form_rpphi(rm, pphi, core, &omega22);
    // print_debug("omega22 = %.16e\n", omega22);
    // REAL8 e22;
    // REAL8 egw_rec;
    // e22 = calc_e22_from_egw(ecc);
    // egw_rec = calc_egw_from_e22(e22);
    // print_debug("ecc = %.16e\n", ecc);
    // print_debug("e22 = %.16e\n", e22);
    // print_debug("egw_rec =  %.16e\n", egw_rec);
QUIT:
    STRUCTFREE(core, SpinEOBParams);
    STRUCTFREE(ICvalues, REAL8Vector);
    STRUCTFREE(dynamicsAdaS, REAL8Array);
    STRUCTFREE(seobdynamicsAdaS, SEOBdynamics);
    return CEV_SUCCESS;
}

INT choose_debug(INT debug_id,
    REAL8 m1,  REAL8 m2, 
    REAL8 s1x, REAL8 s1y, REAL8 s1z, 
    REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 phi0, REAL8 distance,
    REAL8 ecc, REAL8 f_min, REAL8 INdeltaT, REAL8 inc, HyperParams *hparams)
{
    switch(debug_id)
    {
        case 1:
        default:
            debug_InitialCondition_egw(m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, phi0, distance, ecc, f_min, INdeltaT, inc, hparams);
            break;
    }
    return CEV_SUCCESS;
}
