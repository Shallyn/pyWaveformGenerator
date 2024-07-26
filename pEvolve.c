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
#include "pPrec.h"
#include "pInitialConditionExact.h"
#include "myLog.h"

#define EFOLDS 150
INT SEOBConvertInitConditionFromGeneralToV1(REAL8Vector *ICvalues,
    SpinEOBParams *params, 
    REAL8 *ret_omega0,
    REAL8 *ret_e0);

REAL8 NQCWindow(REAL8 t, REAL8 t0, REAL8 W);
static int XLALEOBHighestInitialFreq(
    REAL8 *freqMinRad /**<< OUTPUT, lowest initial 22 mode frequency*/,
    REAL8 mTotal /**<< Total mass in units of solar masses */) 
{
    REAL8 mTScaled = mTotal * CST_MTSUN_SI;
    *freqMinRad = pow(10.5, -1.5) / (CST_PI * mTScaled);
    return CEV_SUCCESS;
}

// PrecEOBversion = 402
// SEOBNRv4P_SpinAlignedEOBversion = 4
// SEOBNRv4P_SymmetrizehPlminusm = 1
// SEOBNRv4P_euler_extension = FLAG_SEOBNRv4P_EULEREXT_QNM_SIMPLE_PRECESSION
// SEOBNRv4P_Zframe = FLAG_SEOBNRv4P_ZFRAME_L
INT evolve(REAL8 m1,  REAL8 m2, 
           REAL8 s1x, REAL8 s1y, REAL8 s1z, 
           REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 phi0, REAL8 distance,
           REAL8 ecc, REAL8 zeta, REAL8 xi, REAL8 f_min, REAL8 Mf_ref, REAL8 INdeltaT, REAL8 inc,
           INT is_only22,
           HyperParams *hparams, 
           REAL8TimeSeries **hPlusOut,
           REAL8TimeSeries **hCrossOut,
           SEOBCoreOutputs *all)
{


    register INT i;
    INT failed = 0, this_step = 0, status = CEV_FAILURE;
    SpinEOBParams *core = NULL;
    REAL8Vector *ICvalues = NULL;
    REAL8Vector *ICvaluesHiS = NULL;
    REAL8Vector *seobvalues_tstartHiS = NULL;
    REAL8Vector *seobvalues_tPeakOmega = NULL;
    REAL8Vector *seobvalues_test = NULL;
    REAL8Vector* m1rVec = NULL;
    REAL8Array *dynamicsAdaS = NULL;
    REAL8Array *dynamicsInverse = NULL;
    REAL8Array *dynamicsHiS = NULL;
    REAL8Vector *Jfinal = NULL;
    REAL8Vector *Lhatfinal = NULL;
    COMPLEX16Vector *sigmaQNMlm0 = NULL;
    SphHarmListEOBNonQCCoeffs *nqcCoeffsList = NULL;
    SphHarmListCAmpPhaseSequence *listhPlm_HiS = NULL;
    SphHarmListCAmpPhaseSequence *listhPlm_HiSRDpatch = NULL;
    SphHarmListCAmpPhaseSequence *listhPlm_AdaS = NULL;
    SphHarmListCAmpPhaseSequence *listhClm = NULL;
    SEOBdynamics *seobdynamicsAdaS = NULL;
    SEOBdynamics *seobdynamicsHiS = NULL;
    
    REAL8Vector *chi1L_tPeakOmega = NULL;
    REAL8Vector *chi2L_tPeakOmega = NULL;

    // Output
    REAL8Vector *tVecPmodes = NULL;
    SphHarmListCAmpPhaseSequence *listhPlm = NULL;
    SEOBdynamics *seobdynamicsAdaSHiS = NULL;
    REAL8Vector *alphaJ2P = NULL;
    REAL8Vector *betaJ2P = NULL;
    REAL8Vector *gammaJ2P = NULL;
    SphHarmTimeSeries *hIlm = NULL;
    SphHarmTimeSeries *hJlm = NULL;
    REAL8TimeSeries *hplusTS = NULL;
    REAL8TimeSeries *hcrossTS = NULL;

    UINT nmodes = 5;
    INT modes_lmax = 5;
    INT modes[5][2] = {{2,2}, {2,1}, {3,3}, {4,4}, {5,5}};
    // memset(modes, 0, 2 * nmodes * sizeof(INT));
    REAL8 mTotal = m1 + m2;
    REAL8 mTScaled = mTotal * CST_MTSUN_SI;
    REAL8 EPS_ALIGN = 1.0e-4;
    REAL8 deltaT = INdeltaT / mTScaled;
    REAL8 amp0 = mTotal * CST_MRSUN_SI / distance / 1e6 / CST_PC_SI;
    REAL8 tStepBack = 200.;
    INT flagZframe = FLAG_SEOBNRv4P_ZFRAME_L;
    hparams->flagZframe = flagZframe;


#if 0
DEBUG_START;
REAL8Vector dbg_xVec, dbg_pVec, dbg_s1Vec, dbg_s2Vec, dbg_sigKerr, dbg_sigStar;
REAL8 dbg_xdata[3], dbg_pdata[3], dbg_s1data[3], dbg_s2data[3], dbg_sigKdata[3], dbg_sigSdata[3];
dbg_xVec.length = dbg_pVec.length = dbg_s1Vec.length = dbg_s2Vec.length = dbg_sigKerr.length = dbg_sigStar.length = 3;
dbg_xVec.data = dbg_xdata;
dbg_pVec.data = dbg_pdata;
dbg_s1Vec.data = dbg_s1data;
dbg_s2Vec.data = dbg_s2data;
dbg_sigKerr.data = dbg_sigKdata;
dbg_sigStar.data = dbg_sigSdata;
// memcpy(dbg_xdata, ICvalues->data, 3*sizeof(REAL8));
// memcpy(dbg_pdata, ICvalues->data+3, 3*sizeof(REAL8));
// dbg_xdata[0] = 4.1659938868533658e1;
dbg_xdata[0] = 2.2786576810580382e+01;
dbg_xdata[1] = dbg_xdata[2] = 0.0;

// dbg_pdata[0] = -1e-2;
// dbg_pdata[1] = 1.5936732938778625e-1;
// dbg_pdata[2] = -2.6799608562176783e-12;
dbg_pdata[0] = -2.0223298742475641e-04;
dbg_pdata[1] = 2.2970285283699207e-01;
dbg_pdata[2] = 0.0;
// memcpy(dbg_s1data, ICvalues->data+6, 3*sizeof(REAL8));
// memcpy(dbg_s2data, ICvalues->data+9, 3*sizeof(REAL8));
// dbg_s1data[0] = 0.0;
// dbg_s1data[1] = -1.3888888888888892e-01;
// dbg_s1data[2] = 6.2500000000000011e-01;
dbg_s1data[0] = dbg_s1data[1] = 0.0;
dbg_s1data[2] = -5.2679500520291378e-01;

// dbg_s2data[0] = 2.7777777777777779e-03;
// dbg_s2data[1] = 0.0;
// dbg_s2data[2] = 2.2222222222222227e-02;
dbg_s2data[0] = dbg_s2data[1] = 0.0;
dbg_s2data[2] = -3.0343392299687830e-02;

REAL8 dbg_m1, dbg_m2, dbg_eta, dbg_a, dbg_S_con;
dbg_m1 = 48.3871;
dbg_m2 = 11.6129;
dbg_eta = dbg_m1 * dbg_m2 / (dbg_m1 + dbg_m2) / (dbg_m1 + dbg_m2);
EOBCalculateSigmaStar(&dbg_sigStar, dbg_m1, dbg_m2, &dbg_s1Vec, &dbg_s2Vec);
EOBCalculateSigmaKerr(&dbg_sigKerr, &dbg_s1Vec, &dbg_s2Vec);

SpinEOBHCoeffs dbg_seobCoeffs;
memset(&dbg_seobCoeffs, 0, sizeof(dbg_seobCoeffs));

  REAL8 dbg_Lhat[3] = {0.0, 0.0, 1.0}; // Not quite true but should be very close
  REAL8 dbg_tempS1_p = inner_product3d(dbg_s1data, dbg_Lhat);
  REAL8 dbg_tempS2_p = inner_product3d(dbg_s2data, dbg_Lhat);
  REAL8 dbg_S1_perp[3] = {0, 0, 0};
  REAL8 dbg_S2_perp[3] = {0, 0, 0};
  for (UINT jj = 0; jj < 3; jj++) {
    dbg_S1_perp[jj] = dbg_s1data[jj] - dbg_tempS1_p * dbg_Lhat[jj];
    dbg_S2_perp[jj] = dbg_s2data[jj] - dbg_tempS2_p * dbg_Lhat[jj];
  }
  REAL8 dbg_sKerr_norm = sqrt(inner_product3d(dbg_sigKerr.data, dbg_sigKerr.data));
  dbg_S_con = 0.0;
  if (dbg_sKerr_norm > 1e-6) {
    dbg_S_con = dbg_sigKerr.data[0] * dbg_Lhat[0] + dbg_sigKerr.data[1] * dbg_Lhat[1] +
            dbg_sigKerr.data[2] * dbg_Lhat[2];
    dbg_S_con /= (1 - 2 * dbg_eta);
    dbg_S_con += (inner_product3d(dbg_S1_perp, dbg_sigKerr.data) +
              inner_product3d(dbg_S2_perp, dbg_sigKerr.data)) /
              dbg_sKerr_norm / (1 - 2 * dbg_eta) / 2.;
  }

dbg_a = sqrt(inner_product3d(dbg_sigKerr.data, dbg_sigKerr.data));
XLALSimIMRCalculateSpinPrecEOBHCoeffs_v2(&dbg_seobCoeffs, dbg_eta, dbg_a, dbg_S_con, 4, NULL);

REAL8 dbg_Ham = XLALSimIMRSpinPrecEOBHamiltonian(dbg_eta, &dbg_xVec, &dbg_pVec, &dbg_s1Vec, &dbg_s2Vec, &dbg_sigKerr, &dbg_sigStar, 1, &dbg_seobCoeffs, NULL);
DEBUG_END;
if(1) {failed = 1; goto QUIT;}
#endif

    /*
    *
    *       Check params
    * 
    */
    PRINT_LOG_INFO(LOG_INFO, "Input Params: phi = %.16e\n\tINdeltaT = %.16e\n\tm1 = %.16e\n\tm2 = %.16e\n\tfMin = %.16e\n\tinc = %.16e\n\tchi1x = %.16e\n\tchi1y = %.16e\n\tchi1z = %.16e\n\tchi2x = %.16e\n\tchi2y = %.16e\n\tchi2z = %.16e", 
        phi0, INdeltaT, m1, m2, f_min, inc, s1x, s1y, s1z, s2x, s2y, s2z);

    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Check params.", this_step);
    if (m1 < 0 || m2 < 0) {failed = 1; goto QUIT;}
    if (m2 > m1)
        SWAP(m1, m2);
    if (m1/m2 > 100) {failed = 1; goto QUIT;}
    if (sqrt(s1x*s1x+s1y*s1y+s1z*s1z) > 1. ||
        sqrt(s2x*s2x+s2y*s2y+s2z*s2z) > 1)
    {failed = 1; goto QUIT;}

    REAL8 freqMinRad = 0;
    /*Compute the highest initial frequency of the 22 mode */
    XLALEOBHighestInitialFreq(&freqMinRad, mTotal);
    if (f_min > freqMinRad)
    {
        PRINT_LOG_INFO(LOG_WARNING, "Initial frequency is too high, the limit is %4.10f", freqMinRad);
        // failed = 1;
        // goto QUIT;
    }
    /* Check NyquistFrequency */
    UINT ell_max = 5;
    REAL8 INchi1[3] = {s1x, s1y, s1z};
    REAL8 INchi2[3] = {s2x, s2y, s2z};
    if (hparams->Mf_min > hparams->Mf_max && hparams->t_max < 0) {
        status = XLALEOBCheckNyquistFrequency(m1, m2, INchi1, INchi2, INdeltaT, ell_max);
        if (status != CEV_SUCCESS)
        {failed = 1; goto QUIT;}
    }
    /*
    *
    * 
    *       Initialization
    *
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Initialization.", this_step);
    core = CreateSpinEOBParams(m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, ecc, hparams);
    if (!core) {failed = 1; goto QUIT;}
    // if (1) {failed = 1; goto QUIT;}
    tStepBack = GET_MAX(tStepBack, core->hParams->tStepBack);
    /*
    *
    * 
    *       Solve Initial Conditions
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Solve Initial Conditions.", this_step);
    ICvalues = CreateREAL8Vector(14);
    REAL8 MfMin = mTScaled * f_min;
    // print_debug("f_min = %.16e, MfMin = %.16e, mTScaled = %.16e, mTotal = %.16e, CST_MTSUN_SI = %.16e\n",
    //      f_min, MfMin, mTScaled, mTotal, CST_MTSUN_SI);
    if (core->hParams && core->hParams->d_ini > 0.0)
    {
        REAL8 d_ini = core->hParams->d_ini; //initial seperation
        REAL8 pr_ini = core->hParams->pr_ini;
        REAL8 pphi_ini = core->hParams->pphi_ini;
        REAL8 ptheta_ini = core->hParams->ptheta_ini;
        REAL8 xSph[3] = {d_ini, 0., 0.};
        REAL8 pSph[3] = {pr_ini, ptheta_ini, pphi_ini};
        REAL8 xCart[3] = {0,0,0};
        REAL8 pCart[3] = {0,0,0};
        if (ptheta_ini > 0)
            core->alignedSpins = TRUE;
        SphericalToCartesian(xCart, pCart, xSph, pSph);
        memset(ICvalues->data, 0, ICvalues->length*sizeof(REAL8));
        memcpy(ICvalues->data, xCart, 3*sizeof(REAL8));
        memcpy(ICvalues->data+3, pCart, 3*sizeof(REAL8));
        memcpy(ICvalues->data+6, core->s1Vec->data, 3*sizeof(REAL8));
        memcpy(ICvalues->data+9, core->s2Vec->data, 3*sizeof(REAL8));
    } else {
        if (ecc > 0.0 && get_egw_flag()) {
            //status = SEOBInitialConditions_egw(ICvalues, MfMin, ecc, core);
            status = SEOBInitialConditions_e_anomaly(ICvalues, Mf_ref, ecc, zeta, xi, core);
            // if (status != CEV_SUCCESS && ecc < 0.05)
            // {
            //     status = SEOBInitialConditions(ICvalues, Mf_ref, ecc, core);
            // }
        } else 
            status = SEOBInitialConditions(ICvalues, Mf_ref, ecc, core);
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    }
    PRINT_LOG_INFO(LOG_DEBUG, "initial conditions:");
    PRINT_LOG_INFO(LOG_DEBUG, "(x,y,z) = (%.16e %.16e %.16e)", 
            ICvalues->data[0], ICvalues->data[1], ICvalues->data[2]);
    PRINT_LOG_INFO(LOG_DEBUG, "(px,py,pz) = (%.16e %.16e %.16e)",
            ICvalues->data[3], ICvalues->data[4], ICvalues->data[5]);
    PRINT_LOG_INFO(LOG_DEBUG, "(S1x,S1y,S1z) = (%.16e %.16e %.16e)",
            ICvalues->data[6], ICvalues->data[7], ICvalues->data[8]);
    PRINT_LOG_INFO(LOG_DEBUG, "(S2x,S2y,S2z) = (%.16e %.16e %.16e)",
            ICvalues->data[9], ICvalues->data[10], ICvalues->data[11]);
    PRINT_LOG_INFO(LOG_DEBUG, "(phi, phiD) = (%.16e %.16e)",
            ICvalues->data[12], ICvalues->data[13]);
    PRINT_LOG_INFO(LOG_DEBUG, "MfMin = %f, deltaT = %f\n", MfMin, deltaT);
    PRINT_LOG_INFO(LOG_DEBUG, "e0 = %f\n", ecc);
    /*
    *
    * 
    *       Evolve EOB trajectory with adaptive sampling 
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Evolve EOB trajectory with adaptive sampling.", this_step);
    INT retLenAdaS = 0, retLenInverse = 0;
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
    if (hparams->inEPS_REL > 0.)
        EPS_REL = hparams->inEPS_REL;
    if (hparams->inEPS_ABS > 0.)
        EPS_ABS = hparams->inEPS_ABS;
    // inv integrate dynamics
    if (MfMin < Mf_ref)
    {
        status = SEOBIntegrateDynamics_inverse(&dynamicsInverse, &retLenInverse, ICvalues, EPS_ABS, EPS_REL,
           deltaT, deltaT_min, tstartAdaS, tendAdaS, core, core->alignedSpins);
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    }
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    if (MfMin < hparams->Mf_max || hparams->t_max>0.0)
    {
        /**
         * @brief 20240304 by X.L.:
         *          if Mf_max > MfMin,
         *              let integrate stops at omega = pi*Mf_max
         */
        status = SEOBIntegrateDynamics_withfMax(&dynamicsAdaS, &retLenAdaS, 
            ICvalues, EPS_ABS, EPS_REL, 
            deltaT, deltaT_min, tstartAdaS, tendAdaS, core, core->alignedSpins);
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
        if (dynamicsInverse)
            SEOBConcactInverseDynToAdaSDyn(&dynamicsAdaS, dynamicsInverse, &retLenAdaS, retLenInverse);
        status = SEOBComputeExtendedSEOBdynamics(&seobdynamicsAdaSHiS, dynamicsAdaS, retLenAdaS, core);    
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
        if (MfMin > Mf_ref)
        {
            status = CutSEOBdynamics(&seobdynamicsAdaSHiS, MfMin);
            retLenAdaS = seobdynamicsAdaSHiS->length;
            if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
        }
        PRINT_LOG_INFO(LOG_DEBUG, "Get retLenAdaS = %d", retLenAdaS);
        tVecPmodes = CreateREAL8Vector(retLenAdaS);
        memcpy(tVecPmodes->data, seobdynamicsAdaSHiS->tVec, retLenAdaS*sizeof(REAL8));
        REAL8 tEndAtFMax = seobdynamicsAdaSHiS->tVec[retLenAdaS-1];
        PRINT_LOG_INFO(LOG_INFO, "Step %d_ Get final J/L/spins from HiS dynamics at peak of Omega, compute constant angles EulerI2J.", this_step);
        status = SEOBInterpolateDynamicsAtTime(&seobvalues_tPeakOmega, tEndAtFMax, seobdynamicsAdaSHiS);
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
        SEOBLFrameVectors(&chi1L_tPeakOmega, &chi2L_tPeakOmega, 
            seobvalues_tPeakOmega, m1, m2, core->hParams->flagZframe);
        PRINT_LOG_INFO(LOG_DEBUG, "chi1L_tPeakOmega = (%.16e, %.16e, %.16e)\n", 
            chi1L_tPeakOmega->data[0], chi1L_tPeakOmega->data[1], chi1L_tPeakOmega->data[2]);
        PRINT_LOG_INFO(LOG_DEBUG, "chi2L_tPeakOmega = (%.16e, %.16e, %.16e)\n",
            chi2L_tPeakOmega->data[0], chi2L_tPeakOmega->data[1], chi2L_tPeakOmega->data[2]);
        /* Compute final J from dynamics quantities */
        SEOBJfromDynamics(&Jfinal, seobvalues_tPeakOmega, core);
        /*Compute the L-hat vector. Note that it has unit norm */
        // SEOBLhatfromDynamics(&Lhatfinal, seobvalues_tPeakOmega, core);
        PRINT_LOG_INFO(LOG_DEBUG, "Jfinal = (%.16e, %.16e, %.16e)", Jfinal->data[0], Jfinal->data[1], Jfinal->data[2]);
        // PRINT_LOG_INFO(LOG_DEBUG, "Lhatfinal = (%.16e, %.16e, %.16e)", Lhatfinal->data[0], Lhatfinal->data[1], Lhatfinal->data[2]);

        REAL8Vector e1J_fmax, e2J_fmax, e3J_fmax;
        e1J_fmax.length = e2J_fmax.length = e3J_fmax.length = 3;
        REAL8 e1Jdata_fmax[3] = {0.};
        REAL8 e2Jdata_fmax[3] = {0.};
        REAL8 e3Jdata_fmax[3] = {0.};
        e1J_fmax.data = e1Jdata_fmax;
        e2J_fmax.data = e2Jdata_fmax;
        e3J_fmax.data = e3Jdata_fmax;
        SEOBBuildJframeVectors(&e1J_fmax, &e2J_fmax, &e3J_fmax, Jfinal);
        PRINT_LOG_INFO(LOG_DEBUG, "e1J = (%.16e, %.16e, %.16e)", e1J_fmax.data[0], e1J_fmax.data[1], e1J_fmax.data[2]);
        PRINT_LOG_INFO(LOG_DEBUG, "e2J = (%.16e, %.16e, %.16e)", e2J_fmax.data[0], e2J_fmax.data[1], e2J_fmax.data[2]);
        PRINT_LOG_INFO(LOG_DEBUG, "e2J = (%.16e, %.16e, %.16e)", e3J_fmax.data[0], e3J_fmax.data[1], e3J_fmax.data[2]);

        /* Compute Euler angles from initial I-frame to final-J-frame */
        /* Note: if spins are aligned, the function SEOBEulerI2JFromJframeVectors */
        /* becomes ill-defined - just keep these Euler angles to zero then */
        REAL8 alphaI2J_fmax = 0., betaI2J_fmax = 0., gammaI2J_fmax = 0.;
        if (!core->alignedSpins) 
            SEOBEulerI2JFromJframeVectors(&alphaI2J_fmax, &betaI2J_fmax, &gammaI2J_fmax, &e1J_fmax, &e2J_fmax, &e3J_fmax);

        SEOBCalculateSphHarmListhlmAmpPhase(&listhPlm_AdaS, modes, nmodes,
                                            seobdynamicsAdaSHiS, NULL,
                                            core, 0);
        PRINT_LOG_INFO(LOG_DEBUG, "Calculate hPlms done");
        status = SEOBEulerJ2PFromDynamics(&alphaJ2P, &betaJ2P, &gammaJ2P, 
                    &e1J_fmax, &e2J_fmax, &e3J_fmax,
                    retLenAdaS, retLenAdaS-1,
                    seobdynamicsAdaSHiS, core);
        PRINT_LOG_INFO(LOG_DEBUG, "Calculate Euler J2P done");
        UINT retLenTS_fmax = floor(((tVecPmodes)->data[retLenAdaS - 1] - (tVecPmodes)->data[0]) / deltaT);
        status = SEOBRotateInterpolatehJlmReImFromSphHarmListhPlmAmpPhase(
            &hJlm, &listhClm, modes, nmodes, modes_lmax, deltaT, retLenTS_fmax, tVecPmodes,
            listhPlm_AdaS, alphaJ2P, betaJ2P, gammaJ2P);
        if ( status == CEV_FAILURE) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "failure in SEOBRotateInterpolatehJlmReImFromSphHarmListhPlmAmpPhase for mode (l,m) = (2,2).");
            failed = 1;
            goto QUIT;
        }
        PRINT_LOG_INFO(LOG_DEBUG, "Calculate hJlm done");
        status = SEOBRotatehIlmFromhJlm(&hIlm, hJlm, modes_lmax, alphaI2J_fmax, betaI2J_fmax, gammaI2J_fmax, deltaT);
        if (status != CEV_SUCCESS)
        {failed = 1; goto QUIT;}
        PRINT_LOG_INFO(LOG_DEBUG, "Calculate hIlm done");
        hplusTS = CreateREAL8TimeSeries(-mTScaled * tEndAtFMax, INdeltaT, retLenTS_fmax);
        hcrossTS = CreateREAL8TimeSeries(-mTScaled * tEndAtFMax, INdeltaT, retLenTS_fmax);
        status = SEOBComputehplushcrossFromhIlm(hplusTS, hcrossTS, modes_lmax, hIlm, amp0,
            inc, phi0, is_only22);
        if (status != CEV_SUCCESS)
        {failed = 1; goto QUIT;}

        (*hPlusOut) = hplusTS;
        (*hCrossOut) = hcrossTS;
        seobdynamicsAdaSHiS->th22Peak = tEndAtFMax;
        all->dyn = seobdynamicsAdaSHiS;
        all->hLM = hIlm;
        all->Plm = listhClm;
        PRINT_LOG_INFO(LOG_DEBUG, "Finished");
        goto QUIT;
    }
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    status = SEOBIntegrateDynamics(&dynamicsAdaS, &retLenAdaS, 
        ICvalues, EPS_ABS, EPS_REL, 
        deltaT, deltaT_min, tstartAdaS, tendAdaS, core, core->alignedSpins);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    if (dynamicsInverse)
        SEOBConcactInverseDynToAdaSDyn(&dynamicsAdaS, dynamicsInverse, &retLenAdaS, retLenInverse);
    PRINT_LOG_INFO(LOG_DEBUG, "AdaS data length = %d", retLenAdaS);
    status = SEOBComputeExtendedSEOBdynamics(&seobdynamicsAdaS, dynamicsAdaS, retLenAdaS, core);    
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    if (MfMin > Mf_ref)
    {
        status = CutSEOBdynamics(&seobdynamicsAdaS, MfMin);
        retLenAdaS = seobdynamicsAdaS->length;
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    }
    m1rVec = CreateREAL8Vector(retLenAdaS);
    for (i = 0; i < retLenAdaS; i++)
    {
        m1rVec->data[i] = -1* seobdynamicsAdaS->polarrVec[i];
    }
    UINT index_6M = FindClosestIndex(m1rVec, -6.0);
    REAL8 time_6M = seobdynamicsAdaS->tVec[index_6M];
    tStepBack = GET_MAX(tStepBack, seobdynamicsAdaS->tVec[retLenAdaS-1] - time_6M + 10*deltaT);
#if 0
    FILE *out = fopen( "debug_dynamics_AdaS.dat","w");
    for (int ii=0; ii<retLenAdaS; ii++)
    {
        // t, x, y, z, px, py, pz
    fprintf(out, "%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n",
      seobdynamicsAdaS->tVec[ii],
      seobdynamicsAdaS->posVecx[ii],
      seobdynamicsAdaS->posVecy[ii],
      seobdynamicsAdaS->posVecz[ii],
      seobdynamicsAdaS->momVecx[ii],
      seobdynamicsAdaS->momVecy[ii],
      seobdynamicsAdaS->momVecz[ii],
      seobdynamicsAdaS->s1Vecx[ii],
      seobdynamicsAdaS->s1Vecy[ii],
      seobdynamicsAdaS->s1Vecz[ii],
      seobdynamicsAdaS->s2Vecx[ii],
      seobdynamicsAdaS->s2Vecy[ii],
      seobdynamicsAdaS->s2Vecz[ii]);
    }
    fclose(out);
#endif
    /*
    *
    * 
    *       (DEBUG MODE) Compare new waveform and old waveform
    * 
    * 
    */
    if (INPUT_DEBUG_FLAG == 1)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "DEBUG_FLAG[%d]:Compare new waveform and old waveform", INPUT_DEBUG_FLAG);
        INT dbg_LL;
        INT dbg_MM;
        COMPLEX16TimeSeries *dbg_hLM_new = NULL;
        COMPLEX16TimeSeries *dbg_hLM_old = NULL;
        CHAR dbg_file_dump_waveform[MAX_STR_LEN];
        FILE *dbg_fout;
        for (dbg_LL = 2; dbg_LL <= 4; dbg_LL++)
        {
            for (dbg_MM = 1; dbg_MM <= dbg_LL; dbg_MM++)
            {
                print_debug("Calculate Mode %d %d\n", dbg_LL, dbg_MM);
                // status = eobcore_object_CalculateWaveform_debug(core, dy, &hLM_new, &hLM_old, LL, MM);
                status = dbg_CalculateWaveformFromDynamicsAdaS(seobdynamicsAdaS, core, dbg_LL, dbg_MM, &dbg_hLM_new, &dbg_hLM_old);
                REAL8 dbg_dt = dbg_hLM_new->deltaT;
                // print_debug("deltaT = %f\n", dbg_dt);
                sprintf(dbg_file_dump_waveform, "debug_h%d%d.dat", dbg_LL, dbg_MM);
                dbg_fout = fopen(dbg_file_dump_waveform, "w");
                for(i=0; i < dbg_hLM_new->data->length; i++)
                {
                    fprintf(dbg_fout, "%.16e %.16e %.16e %.16e %.16e\n", i*dbg_dt, 
                        creal(dbg_hLM_new->data->data[i]), cimag(dbg_hLM_new->data->data[i]),
                        creal(dbg_hLM_old->data->data[i]), cimag(dbg_hLM_old->data->data[i]));
                }
                fclose(dbg_fout);
                DestroyCOMPLEX16TimeSeries(dbg_hLM_new);
                dbg_hLM_new = NULL;
                DestroyCOMPLEX16TimeSeries(dbg_hLM_old);
                dbg_hLM_old = NULL;
            }
        }
        failed = 1;
        goto QUIT;
    }
    /*
    *
    * 
    *       Step back and evolve EOB trajectory at high sampling rate (HiS)
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Step back and evolve EOB trajectory at high sampling rate (HiS).", this_step);
    if (!core->alignedSpins)
    {
        EPS_ABS = 1.0e-8;
        EPS_REL = 1.0e-8;
    } else {
        EPS_REL = 1.0e-9;
        EPS_ABS = 1.0e-10;
    }
    REAL8 deltaTHiS = 1. / 50; /* Fixed at 1/50M */
    REAL8 tstartHiSTarget = seobdynamicsAdaS->tVec[retLenAdaS - 1] - tStepBack;
    INT4 indexstartHiS = retLenAdaS - 1; /* index for the AdaS dynamics */
    while ((indexstartHiS > 0) && (seobdynamicsAdaS->tVec[indexstartHiS] > tstartHiSTarget))
        indexstartHiS--;
    REAL8 tstartHiS = seobdynamicsAdaS->tVec[indexstartHiS];
    PRINT_LOG_INFO(LOG_DEBUG, "tstartHiSTarget = %.16e\n\tindexstartHiS = %d\n\ttstartHiS = %.16e\n",
        tstartHiSTarget, indexstartHiS, tstartHiS);
    // PRINT_LOG_INFO(LOG_DEBUG, "Interpolate Dynamics At Time %f (%f)/%f", tstartHiS, tstartHiSTarget, seobdynamicsAdaS->tVec[retLenAdaS - 1]);
    status = SEOBInterpolateDynamicsAtTime(&seobvalues_tstartHiS, tstartHiS, seobdynamicsAdaS);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    ICvaluesHiS = CreateREAL8Vector(14);
    if (!ICvaluesHiS) {failed = 1; goto QUIT;}
    memcpy(ICvaluesHiS->data, &(seobvalues_tstartHiS->data[1]), 14 * sizeof(REAL8));
    /* Integrate again the dynamics with a constant high sampling rate */
    core->prev_dr = 0.; // This is used to check whether dr/dt is increasing
                            // in the stopping condition
    core->termination_reason =
        -999; // This is going to store the termination condition for the
                // high-sampling integration
    INT is_re_evolved = FALSE; 
    INT retLenHiS = 0;
    REAL8 tendHiS = tstartHiS + tStepBack; /* Seems it should be ignored anyway because of
                                    integrator->stopontestonly */
    REAL8 rISCO;
HISR:
    rISCO = get_h_rISCO();
    PRINT_LOG_INFO(LOG_DEBUG, "tStepBack = %g", tStepBack);
    PRINT_LOG_INFO(LOG_DEBUG, "Integrate Dynamics");
    status = SEOBIntegrateDynamics(&dynamicsHiS, &retLenHiS, 
        ICvaluesHiS, EPS_ABS, EPS_REL, 
        deltaTHiS, 0., tstartHiS, tendHiS, core, 1);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    PRINT_LOG_INFO(LOG_DEBUG, "HiSR data length = %d", retLenHiS);
    /* Compute derived quantities for the high-sampling dynamics */
    status = SEOBComputeExtendedSEOBdynamics(&seobdynamicsHiS, dynamicsHiS, retLenHiS, core);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}

    /* Find time of peak of omega for the High-sampling dynamics */
    INT foundPeakOmega = 0;
    REAL8 tPeakOmega = 0.;
    SEOBLocateTimePeakOmega(&tPeakOmega, &foundPeakOmega, dynamicsHiS,
                            seobdynamicsHiS, retLenHiS, core);
    PRINT_LOG_INFO(LOG_DEBUG, "Locate timePeak of omega = %f", tPeakOmega);

#if 0
    FILE *out2 = fopen( "debug_dynamics_newHiS.dat","w");
    for (INT ii=0; ii<retLenHiS; ii++)
    {
        // t, x, y, z, px, py, pz
        fprintf(out2, "%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n",
        seobdynamicsHiS->tVec[ii],
        seobdynamicsHiS->posVecx[ii],
        seobdynamicsHiS->posVecy[ii],
        seobdynamicsHiS->posVecz[ii],
        seobdynamicsHiS->momVecx[ii],
        seobdynamicsHiS->momVecy[ii],
        seobdynamicsHiS->momVecz[ii],
        seobdynamicsHiS->s1Vecx[ii],
        seobdynamicsHiS->s1Vecy[ii],
        seobdynamicsHiS->s1Vecz[ii],
        seobdynamicsHiS->s2Vecx[ii],
        seobdynamicsHiS->s2Vecy[ii],
        seobdynamicsHiS->s2Vecz[ii]);
    }
    fclose(out2);
#endif
    /*
    *
    * 
    *       Get final J/L/spins from HiS dynamics at peak of Omega, compute constant angles EulerI2J
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Get final J/L/spins from HiS dynamics at peak of Omega, compute constant angles EulerI2J.", this_step);
    status = SEOBInterpolateDynamicsAtTime(&seobvalues_tPeakOmega, tPeakOmega, seobdynamicsHiS);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}

    /* Interpolate the dynamics to r=10M. This is used as a
    feducial point to evaluate the final mass and spin fits */
    REAL8Vector timeVec;
    timeVec.length = retLenAdaS;
    timeVec.data = seobdynamicsAdaS->tVec;

    // FindClosestIndex requires an *increasing* function of time, so we use -r
    // instead of r
    UINT index_10M = FindClosestIndex(m1rVec, -10.0);
    REAL8 time_10M = timeVec.data[index_10M];
    status = SEOBInterpolateDynamicsAtTime(&seobvalues_test, time_10M, seobdynamicsAdaS);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}

    /* Compute the timeshift to get the attachment point */
    SEOBLFrameVectors(&chi1L_tPeakOmega, &chi2L_tPeakOmega, 
        seobvalues_tPeakOmega, m1, m2, core->hParams->flagZframe);
    PRINT_LOG_INFO(LOG_DEBUG, "chi1L_tPeakOmega = (%.16e, %.16e, %.16e)\n", 
        chi1L_tPeakOmega->data[0], chi1L_tPeakOmega->data[1], chi1L_tPeakOmega->data[2]);
    PRINT_LOG_INFO(LOG_DEBUG, "chi2L_tPeakOmega = (%.16e, %.16e, %.16e)\n",
        chi2L_tPeakOmega->data[0], chi2L_tPeakOmega->data[1], chi2L_tPeakOmega->data[2]);
    REAL8 Deltat22 = XLALSimIMREOBGetNRSpinPeakDeltaTv4(
        2, 2, m1, m2, chi1L_tPeakOmega->data[2], chi2L_tPeakOmega->data[2], core->hParams);

    /* Determine the time of attachment */
    REAL8 tAttach = tPeakOmega - Deltat22;
    PRINT_LOG_INFO(LOG_INFO, "tAttach = %.16e, Deltat22 = %.16e", tAttach, Deltat22);
    if (core->hParams->zero_dyncoaphase)
    {
        REAL8 de_phiD, de_phiM, de_phi;
        status = SetZeroPhaseAtTime(seobdynamicsHiS, tAttach, &de_phiD, &de_phiM, &de_phi);
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
        OrbitalPhaseReduce(seobdynamicsAdaS, de_phiD, de_phiM, de_phi);
    }
    //Compute NQC window factors
    if (ecc != 0.0)
    {
        status = SEOBCalculateNQCWindowFactorsFromDyn(seobdynamicsAdaS, tAttach, time_6M, tstartHiS, 6., 0, core);
        // core->wWind = 1.;
    }
    
    if (is_re_evolved)
    {
        PRINT_LOG_INFO(LOG_INFO, "re-calculating window parameters; rISCO = %f, rmin = %f", rISCO, seobdynamicsHiS->polarrVec[retLenHiS-1]);
        status = SEOBCalculateNQCWindowFactorsFromDyn(seobdynamicsHiS, tAttach, time_6M, tstartHiS, rISCO, 1, core);
    }

    /* Compute final J from dynamics quantities */
    SEOBJfromDynamics(&Jfinal, seobvalues_tPeakOmega, core);
    /*Compute the L-hat vector. Note that it has unit norm */
    SEOBLhatfromDynamics(&Lhatfinal, seobvalues_tPeakOmega, core);
    PRINT_LOG_INFO(LOG_DEBUG, "Jfinal = (%.16e, %.16e, %.16e)", Jfinal->data[0], Jfinal->data[1], Jfinal->data[2]);
    PRINT_LOG_INFO(LOG_DEBUG, "Lhatfinal = (%.16e, %.16e, %.16e)", Lhatfinal->data[0], Lhatfinal->data[1], Lhatfinal->data[2]);
    REAL8 Jmag = sqrt(inner_product3d(Jfinal->data, Jfinal->data));
    /* Cosine of the angle between L-hat and J. Needed to determine
    * the correct sign of the final spin
    */
    REAL8 cos_angle = inner_product3d(Jfinal->data, Lhatfinal->data) / Jmag;
    /* Compute final-J-frame unit vectors e1J, e2J, e3J=Jfinalhat */
    /* Convention: if (ex, ey, ez) is the initial I-frame, e1J chosen such that ex
    * is in the plane (e1J, e3J) and ex.e1J>0 */
    REAL8Vector e1J, e2J, e3J;
    e1J.length = e2J.length = e3J.length = 3;
    REAL8 e1Jdata[3] = {0.};
    REAL8 e2Jdata[3] = {0.};
    REAL8 e3Jdata[3] = {0.};
    e1J.data = e1Jdata;
    e2J.data = e2Jdata;
    e3J.data = e3Jdata;
    SEOBBuildJframeVectors(&e1J, &e2J, &e3J, Jfinal);
    PRINT_LOG_INFO(LOG_DEBUG, "e1J = (%.16e, %.16e, %.16e)", e1J.data[0], e1J.data[1], e1J.data[2]);
    PRINT_LOG_INFO(LOG_DEBUG, "e2J = (%.16e, %.16e, %.16e)", e2J.data[0], e2J.data[1], e2J.data[2]);
    PRINT_LOG_INFO(LOG_DEBUG, "e2J = (%.16e, %.16e, %.16e)", e3J.data[0], e3J.data[1], e3J.data[2]);

    /* Compute Euler angles from initial I-frame to final-J-frame */
    /* Note: if spins are aligned, the function SEOBEulerI2JFromJframeVectors */
    /* becomes ill-defined - just keep these Euler angles to zero then */
    REAL8 alphaI2J = 0., betaI2J = 0., gammaI2J = 0.;
    if (!core->alignedSpins) 
        SEOBEulerI2JFromJframeVectors(&alphaI2J, &betaI2J, &gammaI2J, &e1J, &e2J, &e3J);

    /*
    *
    * 
    *       Compute P-frame amp/phase for all modes on HiS and compute NQC Coeffs
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Compute P-frame amp/phase for all modes on HiS and compute NQC Coeffs.", this_step);
    // REAL8 wWind_old = core->wWind;
    status = SEOBCalculateSphHarmListNQCCoefficientsV4(
            &nqcCoeffsList, modes, nmodes, tPeakOmega, seobdynamicsHiS,
            core, chi1L_tPeakOmega, chi2L_tPeakOmega);

    if ( status != CEV_SUCCESS) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "NQC computation failed.");
        failed = 1;
        goto QUIT;
    }
    if (!is_re_evolved)
    {
        core->wWind = 1.;
    }
    // Check should we need to re-evolve HiSR EOB
// #if 1
    if (ecc != 0.0 && rISCO > 3 && CheckStopCondition(seobdynamicsHiS, core, nqcCoeffsList, tAttach - seobdynamicsHiS->tVec[0]) != CEV_SUCCESS)
    {
#if 1
        // re-evolve EOB
        SET_RISCO(rISCO-1.);
        // PRINT_LOG_INFO(LOG_INFO, "Re-evolve HiSR EOB...");
        // print_debug("Re-evolve EOB...");
        is_re_evolved = TRUE;
        STRUCTFREE(seobdynamicsHiS, SEOBdynamics);
        STRUCTFREE(Jfinal, REAL8Vector);
        STRUCTFREE(Lhatfinal, REAL8Vector);
        STRUCTFREE(seobvalues_tPeakOmega, REAL8Vector);
        STRUCTFREE(seobvalues_test, REAL8Vector);
        STRUCTFREE(chi1L_tPeakOmega, REAL8Vector);
        STRUCTFREE(chi2L_tPeakOmega, REAL8Vector);
        STRUCTFREE(dynamicsHiS, REAL8Array);
        // STRUCTFREE(seobvalues_tstartHiS, REAL8Vector);
        STRUCTFREE(nqcCoeffsList, SphHarmListEOBNonQCCoeffs);
        // STRUCTFREE(ICvaluesHiS, REAL8Vector);
        goto HISR;
#endif
    }
    /*
    *
    *       Compute P-frame amp/phase for all modes on HiS, now including NQC
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Compute P-frame amp/phase for all modes on HiS, now including NQC.", this_step);
    // /* We now include the NQC when generating modes */
    // flagNQC = 1;

    /* Compute amplitude and phase of the P-frame modes hPlm on high sampling,
    * with NQC */
    SEOBCalculateSphHarmListhlmAmpPhase(&listhPlm_HiS, modes, nmodes,
                                        seobdynamicsHiS, nqcCoeffsList,
                                        core, 1);

    /*
    *
    * 
    *       (DEBUG MODE) Compare new waveform and old waveform
    * 
    * 
    */
    if (INPUT_DEBUG_FLAG == 2)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "DEBUG_FLAG[%d]:Compare new waveform and old waveform with correction", INPUT_DEBUG_FLAG);
        print_debug("tAttach = %f\n", tAttach);
        INT dbg_LL;
        INT dbg_MM;
        COMPLEX16TimeSeries *dbg_hLM_new = NULL;
        COMPLEX16TimeSeries *dbg_hLM_old = NULL;
        CHAR dbg_file_dump_waveform[MAX_STR_LEN];
        FILE *dbg_fout;
        for (dbg_LL = 2; dbg_LL <= 4; dbg_LL++)
        {
            for (dbg_MM = 1; dbg_MM <= dbg_LL; dbg_MM++)
            {
                print_debug("Calculate Mode %d %d\n", dbg_LL, dbg_MM);
                // status = eobcore_object_CalculateWaveform_debug(core, dy, &hLM_new, &hLM_old, LL, MM);
                status = dbg_CalculateWaveformFromDynamicsAdaS(seobdynamicsAdaS, core, dbg_LL, dbg_MM, &dbg_hLM_new, &dbg_hLM_old);
                REAL8 dbg_dt = dbg_hLM_new->deltaT;
                sprintf(dbg_file_dump_waveform, "debug_h%d%d.dat", dbg_LL, dbg_MM);
                dbg_fout = fopen(dbg_file_dump_waveform, "w");
                for(i=0; i < dbg_hLM_new->data->length; i++)
                {
                    fprintf(dbg_fout, "%.16e %.16e %.16e %.16e %.16e\n", i*deltaT - tAttach, 
                        creal(dbg_hLM_new->data->data[i]), cimag(dbg_hLM_new->data->data[i]),
                        creal(dbg_hLM_old->data->data[i]), cimag(dbg_hLM_old->data->data[i]));
                }
                fclose(dbg_fout);
                DestroyCOMPLEX16TimeSeries(dbg_hLM_new);
                dbg_hLM_new = NULL;
                DestroyCOMPLEX16TimeSeries(dbg_hLM_old);
                dbg_hLM_old = NULL;
            }
        }
        failed = 1;
        goto QUIT;
    }

    if (INPUT_DEBUG_FLAG == 3)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "DEBUG_FLAG[%d]:Calculate NQC correction series", INPUT_DEBUG_FLAG);
        INT dbg_LL = 2;
        INT dbg_MM = 2;
        COMPLEX16TimeSeries *dbg_nqc = NULL;
        CHAR dbg_file_dump_nqc[MAX_STR_LEN];
        FILE *dbg_fout;
        sprintf(dbg_file_dump_nqc, "debug_nqc_h%d%d.dat", dbg_LL, dbg_MM);
        dbg_fout = fopen(dbg_file_dump_nqc, "w");
        print_debug("Calculate Mode %d %d\n", dbg_LL, dbg_MM);
        status = dbg_CalculateNQCTimeSeries(seobdynamicsHiS, core, nqcCoeffsList, dbg_LL, dbg_MM, &dbg_nqc);
        REAL8 dbg_dt = dbg_nqc->deltaT;
        REAL8 nWind;
        print_debug("tWind = %f, wWind = %f", core->tWind, core->wWind);
        for(i=0; i < dbg_nqc->data->length; i++)
        {
            nWind = NQCWindow(seobdynamicsHiS->tVec[i], core->tWind, core->wWind);
            // print_debug("nWind = %f\n", nWind);
            fprintf(dbg_fout, "%.16e %.16e %.16e %.16e\n", i*dbg_dt + seobdynamicsHiS->tVec[0] - tAttach, 
                creal(dbg_nqc->data->data[i]), cimag(dbg_nqc->data->data[i]), nWind);
        }
        fclose(dbg_fout);
        STRUCTFREE(dbg_nqc, COMPLEX16TimeSeries);
        // failed = 1;
        // goto QUIT;
    }
    /*
    *
    * 
    *       Attach RD to the P-frame waveform
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Attach RD to the P-frame waveform.", this_step);
    REAL8 finalMass = 0., finalSpin = 0.;
    status = SEOBGetFinalSpinMass(&finalMass, &finalSpin, seobvalues_test, core);
    if (status != CEV_SUCCESS) 
    {failed = 1; goto QUIT;}

    /* The function above returns only the magnitude of the spin.
    *  We pick the direction based on whether Lhat \cdot J is positive
    *  or negative */
    if (cos_angle < 0)
    {
        finalSpin *= -1;
    }
    /* finalSpin interpolation is available only between -0.9996 and 0.9996 */
    /* Set finalSpin to +/- 0.9996 if it is out of this range */
    if (finalSpin < -0.9996)
        finalSpin = -0.9996;
    if (finalSpin > 0.9996)
        finalSpin = 0.9996;

    PRINT_LOG_INFO(LOG_INFO, "final mass = %e, final spin = %e", finalMass, finalSpin);

    /* Estimate leading QNM , to determine how long the ringdown
    * patch will be */
    /* NOTE: XLALSimIMREOBGenerateQNMFreqV2Prec returns the complex frequency in
    * physical units... */
    COMPLEX16Vector sigmaQNM220estimatephysicalVec;
    COMPLEX16 sigmaQNM220estimatephysical = 0.;
    sigmaQNM220estimatephysicalVec.length = 1;
    sigmaQNM220estimatephysicalVec.data = &sigmaQNM220estimatephysical;
    status = XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec(&sigmaQNM220estimatephysicalVec, m1,
                                                m2, finalMass, finalSpin, 2,
                                                2, 1);
    if (status == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec for mode (l,m) = (2,2).");
        failed = 1;
        goto QUIT;
    }
    COMPLEX16 sigmaQNM220estimate = mTScaled * sigmaQNM220estimatephysical;

    /* Length of RD patch, 40 e-folds of decay of the estimated QNM220 */
    UINT retLenRDPatch =
    (UINT)ceil(EFOLDS / (cimag(sigmaQNM220estimate) * deltaTHiS));


    /* Attach RD to the P-frame modes */
    /* Vector holding the values of the 0-th overtone QNM complex frequencies for
    * the modes (l,m) */
    // NOTE: the QNM complex frequencies are computed inside
    // SEOBAttachRDToSphHarmListhPlm
    status = SEOBAttachRDToSphHarmListhPlm(
        &listhPlm_HiSRDpatch, &sigmaQNMlm0, modes, nmodes, finalMass, finalSpin,
        listhPlm_HiS, deltaTHiS, retLenHiS, retLenRDPatch, tAttach,
        seobvalues_tPeakOmega, seobdynamicsHiS, core);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}

    /*
    *
    * 
    *       Build the joined dynamics AdaS+HiS up to attachment, joined P-modes AdaS+HiS+RDpatch 
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Build the joined dynamics AdaS+HiS up to attachment, joined P-modes AdaS+HiS+RDpatch", this_step);

    /* Compute amplitude and phase of the P-frame modes hPlm on adaptive sampling,
    * with NQC */
    // flagNQC = 1;
    SEOBCalculateSphHarmListhlmAmpPhase(&listhPlm_AdaS, modes, nmodes,
                                    seobdynamicsAdaS, nqcCoeffsList,
                                    core, 1);

    /* Vector of times for the P-modes, that will be used for interpolation:
    * joining AdaS and HiS+RDpatch */
    UINT retLenPmodes = 0;
    /* First junction at indexAdaSHiS, tAdaSHiS */
    UINT indexJoinHiS = 0;
    REAL8 tJoinHiS = 0.;
    /* Second junction at seobdynamicsAdaS, tJoinAttach */
    UINT indexJoinAttach = 0;
    REAL8 tJoinAttach = 0.;
    /* Construct the joined vector of times (AdaS+HiS+RDpatch) and keep the
    * jonction indices and times */
    status = SEOBJoinTimeVector(&tVecPmodes, &retLenPmodes, &tJoinHiS, &indexJoinHiS,
                    &tJoinAttach, &indexJoinAttach, retLenRDPatch, deltaTHiS,
                    tstartHiS, tAttach, seobdynamicsAdaS, seobdynamicsHiS);
 

    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}
    PRINT_LOG_INFO(LOG_DEBUG, "tJoinAttach = %.16e, indexJoinAttach = %u", tJoinAttach, indexJoinAttach);
    /* Copy dynamics from AdaS<HiS and HiS<tAttach to form joined dynamics, ending
    * at the last time sample <tAttach */
    // NOTE: we cut the dynamics at tAttach, as we will extend the Euler
    // angles for t>=tAttach -- but we could also choose to finish at tPeakOmega
    // which is used for final-J and for the final mass/spin fit

    status = SEOBJoinDynamics(&seobdynamicsAdaSHiS, seobdynamicsAdaS, seobdynamicsHiS,
                indexJoinHiS, indexJoinAttach);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}

    /* Copy waveform modes from AdaS and HiS+RDpatch - adjusting 2pi-phase shift
    * at the junction point AdaS/HiS */
    status = SEOBJoinSphHarmListhlm(&listhPlm, listhPlm_AdaS, listhPlm_HiSRDpatch, modes,
                        nmodes, indexstartHiS);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}

    /* Get the time of the frame-invariant amplitude peak */
    REAL8 tPeak = 0;
    UINT indexPeak = 0;
    // NOTE: peak amplitude using l=2 only: h22 required, and h21 used if present
    status = SEOBAmplitudePeakFromAmp22Amp21(&tPeak, &indexPeak, listhPlm, modes, nmodes,
                                tVecPmodes);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}

    /*
    *
    * 
    *       Compute Euler angles J2P from AdaS and HiS dynamics up to attachment
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Compute Euler angles J2P from AdaS and HiS dynamics up to attachment", this_step);
    /* Compute Euler angles J2P from the dynamics before attachment point */
    /* If SpinsAlmostAligned, all Euler angles are set to 0 */
    status = SEOBEulerJ2PFromDynamics(&alphaJ2P, &betaJ2P, &gammaJ2P, 
                &e1J, &e2J, &e3J,
                retLenPmodes, indexJoinAttach,
                seobdynamicsAdaSHiS, core);
    if (status == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "SEOBEulerJ2PFromDynamics failed.");
        failed = 1;
        goto QUIT;
    }

    /*
    *
    * 
    *       Compute Euler angles J2P extension after attachment
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Compute Euler angles J2P extension after attachment", this_step);

    /* Compute Euler angles J2P according to the prescription flagEulerextension
    * after attachment point */
    /* If SpinsAlmostAligned, all Euler angles are set to 0 */
    /* NOTE: Regardless of the mode content of hPlm, the frame extension at the
    * moment is based on sigmaQNM22, sigmaQNM21 */
    COMPLEX16 sigmaQNM220 = 0., sigmaQNM210 = 0.;
    COMPLEX16Vector sigmaQNM220physicalVec, sigmaQNM210physicalVec;
    sigmaQNM220physicalVec.length = 1;
    sigmaQNM210physicalVec.length = 1;
    COMPLEX16 sigmaQNM220physicalval = 0.;
    COMPLEX16 sigmaQNM210physicalval = 0.;
    sigmaQNM220physicalVec.data = &sigmaQNM220physicalval;
    sigmaQNM210physicalVec.data = &sigmaQNM210physicalval;
    /* NOTE: XLALSimIMREOBGenerateQNMFreqV2Prec returns the complex frequency in
    * physical units... */
    status = XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec(&sigmaQNM220physicalVec, m1,
                                                m2, finalMass, finalSpin, 2,
                                                2, 1);
    if (status == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec for mode (l,m) = (2,2).");
        failed = 1;
        goto QUIT;
    }
    status = XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec(&sigmaQNM210physicalVec, m1,
                                                m2, finalMass, finalSpin, 2,
                                                1, 1);
    if (status == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec for mode (l,m) = (2,1).");
        failed = 1;
        goto QUIT;
    }
    sigmaQNM220 = mTScaled * sigmaQNM220physicalVec.data[0];
    sigmaQNM210 = mTScaled * sigmaQNM210physicalVec.data[0];
    INT flip = 1;
    if (cos_angle < 0)
        flip = -1;
    // flagEulerextension = 0
    status = SEOBEulerJ2PPostMergerExtension(
        alphaJ2P, betaJ2P, gammaJ2P, sigmaQNM220, sigmaQNM210, tVecPmodes,
        retLenPmodes, indexJoinAttach, core, flip);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}

    /*
    * 
    * 
    *       Compute modes hJlm on the output time series by rotating and interpolating the modes hPlm
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Compute modes hJlm on the output time series by rotating and interpolating the modes hPlm", this_step);

    /* Determine the length of the fixed-sampling output time series */
    UINT retLenTS = floor(((tVecPmodes)->data[retLenPmodes - 1] - (tVecPmodes)->data[0]) / deltaT);

#if 0
    CAmpPhaseSequence *dbg_h22 = SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 2, 2)->campphase;
    CAmpPhaseSequence *dbg_h21 = SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 2, 1)->campphase;
    CAmpPhaseSequence *dbg_h33 = SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 3, 3)->campphase;
    CAmpPhaseSequence *dbg_h44 = SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 4, 4)->campphase;
    CAmpPhaseSequence *dbg_h55 = SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 5, 5)->campphase;
    INT dbg_length;
    dbg_length = dbg_h22->xdata->length;
    FILE *outP = fopen( "debug_hlmP_new.dat","w");
    for (INT ii=0; ii<dbg_length; ii++)
    {
        // t, h22_r, h22_i, h21_r, h21_i, h33_r, h33_i, h44_r, h44_i, h55_r, h55_i
        fprintf(outP, "%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\n", dbg_h22->xdata->data[ii],
            dbg_h22->camp_real->data[ii]*cos(dbg_h22->phase->data[ii]), dbg_h22->camp_real->data[ii]*sin(dbg_h22->phase->data[ii]),
            dbg_h21->camp_real->data[ii]*cos(dbg_h21->phase->data[ii]), dbg_h21->camp_real->data[ii]*sin(dbg_h21->phase->data[ii]),
            dbg_h33->camp_real->data[ii]*cos(dbg_h33->phase->data[ii]), dbg_h33->camp_real->data[ii]*sin(dbg_h33->phase->data[ii]),
            dbg_h44->camp_real->data[ii]*cos(dbg_h44->phase->data[ii]), dbg_h44->camp_real->data[ii]*sin(dbg_h44->phase->data[ii]),
            dbg_h55->camp_real->data[ii]*cos(dbg_h55->phase->data[ii]), dbg_h55->camp_real->data[ii]*sin(dbg_h55->phase->data[ii]));
    }
    fclose(outP);
#endif

    /* Rotate waveform from P-frame to J-frame */
    // flagSymmetrizehPlminusm = 1
    status = SEOBRotateInterpolatehJlmReImFromSphHarmListhPlmAmpPhase(
        &hJlm, &listhClm, modes, nmodes, modes_lmax, deltaT, retLenTS, tVecPmodes,
        listhPlm, alphaJ2P, betaJ2P, gammaJ2P);
    if ( status == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in SEOBRotateInterpolatehJlmReImFromSphHarmListhPlmAmpPhase for mode (l,m) = (2,2).");
        failed = 1;
        goto QUIT;
    }
    /*
    *
    * 
    *       Rotate waveform from J-frame to the output I-frame on timeseries-sampling (constant Wigner coeffs)
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Rotate waveform from J-frame to the output I-frame on timeseries-sampling (constant Wigner coeffs)", this_step);

    /* Rotate waveform from J-frame to I-frame */
    status = SEOBRotatehIlmFromhJlm(&hIlm, hJlm, modes_lmax, alphaI2J, betaI2J, gammaI2J, deltaT);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}
    hIlm->tAttach = tAttach;

    /*
    *
    * 
    *       Compute h20 from I-frame waveform on timeseries sampling
    * 
    * 
    */
    // this_step++;
    // PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Compute h20 from I-frame waveform on timeseries sampling", this_step);

    // status = SEOBCalculateh20(hIlm);
    // {failed = 1; goto QUIT;}
#if 0
  COMPLEX16TimeSeries *dbg_h22 = XLALSphHarmTimeSeriesGetMode(hIlm, 2, 2);
  COMPLEX16TimeSeries *dbg_h21 = XLALSphHarmTimeSeriesGetMode(hIlm, 2, 1);
  COMPLEX16TimeSeries *dbg_h33 = XLALSphHarmTimeSeriesGetMode(hIlm, 3, 3);
  COMPLEX16TimeSeries *dbg_h44 = XLALSphHarmTimeSeriesGetMode(hIlm, 4, 4);
  COMPLEX16TimeSeries *dbg_h55 = XLALSphHarmTimeSeriesGetMode(hIlm, 5, 5);
  REAL8 dbg_deltaT = dbg_h22->deltaT;
  INT dbg_length = dbg_h22->data->length;
  FILE *outI = fopen( "debug_hlmI_new.dat","w");
  for(INT ii=0;ii<dbg_length;ii++)
  {
      fprintf(outI, "%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\n", ii*dbg_deltaT - tPeakOmega,
          creal(dbg_h22->data->data[ii]), cimag(dbg_h22->data->data[ii]),
          creal(dbg_h21->data->data[ii]), cimag(dbg_h21->data->data[ii]),
          creal(dbg_h33->data->data[ii]), cimag(dbg_h33->data->data[ii]),
          creal(dbg_h44->data->data[ii]), cimag(dbg_h44->data->data[ii]),
          creal(dbg_h55->data->data[ii]), cimag(dbg_h55->data->data[ii]));

  }
  fclose(outI);

#endif

    /*
    *
    * 
    *       Compute hplus, hcross from I-frame waveform on timeseries sampling
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Compute hplus, hcross from I-frame waveform on timeseries sampling", this_step);

    /* GPS time for output time series and modes */
    // LIGOTimeGPS tGPS = LIGOTIMEGPSZERO;
    // XLALGPSAdd(&tGPS, -mTScaled *
    //                     tPeak); /* tPeak converted back to dimensionfull time */
    REAL8 tGPS = -mTScaled * tPeak;
    /* Create output timeseries for hplus, hcross */
    /* Use the dimensionfull INdeltaT (s) as time step */
    // hplusTS = CreateREAL8TimeSeries("H_PLUS", &tGPS, 0.0, INdeltaT,
    //                                 &lalStrainUnit, retLenTS);
    // hcrossTS = CreateREAL8TimeSeries("H_CROSS", &tGPS, 0.0, INdeltaT,
    //                                 &lalStrainUnit, retLenTS);
    hplusTS = CreateREAL8TimeSeries(tGPS, INdeltaT, retLenTS);
    hcrossTS = CreateREAL8TimeSeries(tGPS, INdeltaT, retLenTS);
    /* Compute hplus, hcross from hIlm */
    // NOTE: azimuthal angle of the observer entering the -2Ylm is pi/2-phi
    // according to LAL conventions
    status = SEOBComputehplushcrossFromhIlm(hplusTS, hcrossTS, modes_lmax, hIlm, amp0,
        inc, phi0, is_only22);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}
#if 0
    FILE *outPC = fopen( "debug_hlmI_new.dat","w");
    INT dbg_length = hplusTS->data->length;
    for(INT ii=0;ii<dbg_length;ii++)
    {
        fprintf(outPC, "%.16e\t%.16e\t%.16e\n",
            ii*deltaT-tPeakOmega, hplusTS->data->data[ii], hcrossTS->data->data[ii]);
    }
    fclose(outPC);
#endif
    /*
    *
    * 
    *       Output and cleanup
    * 
    * 
    */

    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Output and cleanup", this_step);
    /* Point the output pointers to the relevant time series and return */
    (*hPlusOut) = hplusTS;
    (*hCrossOut) = hcrossTS;
    seobdynamicsAdaSHiS->th22Peak = tAttach;
    all->dyn = seobdynamicsAdaSHiS;
    all->hLM = hIlm;
    all->Plm = listhClm;
#if 0
    /* Output vector gathering quantities related to merger (similar to previous
    * AttachParams) */
    /* Format: tPeakOmega tAttach tPeak Jfinalx Jfinaly Jfinalz finalMassfit
    * finalSpinfit termination_reason [sigmaQNMlm0Re sigmaQNMlm0Im for lm in
    * modes] */
    /* NOTE: the size of this output vector depends on the number of modes, due to
    * the sigmaQNM */
    if (!((*mergerParams) = CreateREAL8Vector(9 + 2 * nmodes))) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failed to allocate REAL8Vector mergerParams.");
        failed = 1;
        goto QUIT;
    }
    (*mergerParams)->data[0] = tPeakOmega;
    (*mergerParams)->data[1] = tAttach;
    (*mergerParams)->data[2] = tPeak;
    (*mergerParams)->data[3] = Jfinal->data[0];
    (*mergerParams)->data[4] = Jfinal->data[1];
    (*mergerParams)->data[5] = Jfinal->data[2];
    (*mergerParams)->data[6] = finalMass;
    (*mergerParams)->data[7] = finalSpin;
    (*mergerParams)->data[8] = seobParams.termination_reason;
    for (UINT4 nmode = 0; nmode < nmodes; nmode++) {
    (*mergerParams)->data[9 + 2 * nmode] = creal(sigmaQNMlm0->data[nmode]);
    (*mergerParams)->data[9 + 2 * nmode + 1] = cimag(sigmaQNMlm0->data[nmode]);
    }


    /* Additional outputs */

    /* Dynamics */
    // NOTE: casting to REAL8Vector due to the SWIG wrapping
    *seobdynamicsAdaSVector =
    XLALCreateREAL8Vector(v4PdynamicsVariables * seobdynamicsAdaS->length);
    *seobdynamicsHiSVector =
    XLALCreateREAL8Vector(v4PdynamicsVariables * seobdynamicsHiS->length);
    *seobdynamicsAdaSHiSVector =
    XLALCreateREAL8Vector(v4PdynamicsVariables * seobdynamicsAdaSHiS->length);
    memcpy((*seobdynamicsAdaSVector)->data, seobdynamicsAdaS->array->data,
    (v4PdynamicsVariables * seobdynamicsAdaS->length) * sizeof(REAL8));
    memcpy((*seobdynamicsHiSVector)->data, seobdynamicsHiS->array->data,
    (v4PdynamicsVariables * seobdynamicsHiS->length) * sizeof(REAL8));
    memcpy((*seobdynamicsAdaSHiSVector)->data, seobdynamicsAdaSHiS->array->data,
    (v4PdynamicsVariables * seobdynamicsAdaSHiS->length) * sizeof(REAL8));

    /* Modes in the P-frame */
    // NOTE: casting to REAL8Vector due to the SWIG wrapping
    // NOTE: in the output, real amplitude instead of complex envelope
    /* (2,2) */
    *hP22_amp = XLALCreateREAL8Vector(retLenPmodes);
    *hP22_phase = XLALCreateREAL8Vector(retLenPmodes);
    SphHarmListCAmpPhaseSequence *hP22 =
    SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 2, 2);
    if (hP22 == NULL) {
    memset((*hP22_amp)->data, 0, retLenPmodes * sizeof(REAL8));
    memset((*hP22_phase)->data, 0, retLenPmodes * sizeof(REAL8));
    } else {
    memcpy((*hP22_amp)->data, hP22->campphase->camp_real->data,
        retLenPmodes * sizeof(REAL8));
    memcpy((*hP22_phase)->data, hP22->campphase->phase->data,
        retLenPmodes * sizeof(REAL8));
    }
    /* (2,1) */
    *hP21_amp = XLALCreateREAL8Vector(retLenPmodes);
    *hP21_phase = XLALCreateREAL8Vector(retLenPmodes);
    SphHarmListCAmpPhaseSequence *hP21 =
    SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 2, 1);
    if (hP21 == NULL) {
    memset((*hP21_amp)->data, 0, retLenPmodes * sizeof(REAL8));
    memset((*hP21_phase)->data, 0, retLenPmodes * sizeof(REAL8));
    } else {
    memcpy((*hP21_amp)->data, hP21->campphase->camp_real->data,
        retLenPmodes * sizeof(REAL8));
    memcpy((*hP21_phase)->data, hP21->campphase->phase->data,
        retLenPmodes * sizeof(REAL8));
    }
    /* (3,3) */
    *hP33_amp = XLALCreateREAL8Vector(retLenPmodes);
    *hP33_phase = XLALCreateREAL8Vector(retLenPmodes);
    SphHarmListCAmpPhaseSequence *hP33 =
    SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 3, 3);
    if (hP33 == NULL) {
    memset((*hP33_amp)->data, 0, retLenPmodes * sizeof(REAL8));
    memset((*hP33_phase)->data, 0, retLenPmodes * sizeof(REAL8));
    } else {
    memcpy((*hP33_amp)->data, hP33->campphase->camp_real->data,
        retLenPmodes * sizeof(REAL8));
    memcpy((*hP33_phase)->data, hP33->campphase->phase->data,
        retLenPmodes * sizeof(REAL8));
    }
    /* (4,4) */
    *hP44_amp = XLALCreateREAL8Vector(retLenPmodes);
    *hP44_phase = XLALCreateREAL8Vector(retLenPmodes);
    SphHarmListCAmpPhaseSequence *hP44 =
    SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 4, 4);
    if (hP44 == NULL) {
    memset((*hP44_amp)->data, 0, retLenPmodes * sizeof(REAL8));
    memset((*hP44_phase)->data, 0, retLenPmodes * sizeof(REAL8));
    } else {
    memcpy((*hP44_amp)->data, hP44->campphase->camp_real->data,
        retLenPmodes * sizeof(REAL8));
    memcpy((*hP44_phase)->data, hP44->campphase->phase->data,
        retLenPmodes * sizeof(REAL8));
    }
    /* (5,5) */
    *hP55_amp = XLALCreateREAL8Vector(retLenPmodes);
    *hP55_phase = XLALCreateREAL8Vector(retLenPmodes);
    SphHarmListCAmpPhaseSequence *hP55 =
    SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 5, 5);
    if (hP55 == NULL) {
    memset((*hP55_amp)->data, 0, retLenPmodes * sizeof(REAL8));
    memset((*hP55_phase)->data, 0, retLenPmodes * sizeof(REAL8));
    } else {
    memcpy((*hP55_amp)->data, hP55->campphase->camp_real->data,
        retLenPmodes * sizeof(REAL8));
    memcpy((*hP55_phase)->data, hP55->campphase->phase->data,
        retLenPmodes * sizeof(REAL8));
    }
#endif

QUIT:
    PRINT_LOG_INFO(LOG_INFO, "Program END\n\n");
    STRUCTFREE(core, SpinEOBParams);

    STRUCTFREE(dynamicsAdaS, REAL8Array);
    STRUCTFREE(dynamicsInverse, REAL8Array);
    STRUCTFREE(ICvalues, REAL8Vector);
    STRUCTFREE(seobdynamicsAdaS, SEOBdynamics);

    STRUCTFREE(dynamicsHiS, REAL8Array);
    STRUCTFREE(ICvaluesHiS, REAL8Vector);
    STRUCTFREE(seobvalues_tstartHiS, REAL8Vector);
    STRUCTFREE(seobdynamicsHiS, SEOBdynamics);

    STRUCTFREE(seobvalues_tPeakOmega, REAL8Vector);
    STRUCTFREE(seobvalues_test, REAL8Vector);
    STRUCTFREE(m1rVec,REAL8Vector);
    STRUCTFREE(Jfinal, REAL8Vector);
    STRUCTFREE(Lhatfinal, REAL8Vector);
    STRUCTFREE(sigmaQNMlm0, COMPLEX16Vector);

    STRUCTFREE(chi2L_tPeakOmega, REAL8Vector);
    STRUCTFREE(chi1L_tPeakOmega, REAL8Vector);
    STRUCTFREE(nqcCoeffsList, SphHarmListEOBNonQCCoeffs);

    STRUCTFREE(listhPlm_HiS, SphHarmListCAmpPhaseSequence);
    STRUCTFREE(listhPlm_HiSRDpatch, SphHarmListCAmpPhaseSequence);
    STRUCTFREE(listhPlm_AdaS, SphHarmListCAmpPhaseSequence);

    STRUCTFREE(tVecPmodes, REAL8Vector);
    STRUCTFREE(listhPlm, SphHarmListCAmpPhaseSequence);
    // STRUCTFREE(seobdynamicsAdaSHiS, SEOBdynamics);

    STRUCTFREE(alphaJ2P, REAL8Vector);
    STRUCTFREE(betaJ2P, REAL8Vector);
    STRUCTFREE(gammaJ2P, REAL8Vector);

    // STRUCTFREE(hIlm, SphHarmTimeSeries);
    STRUCTFREE(hJlm, SphHarmTimeSeries);
    if (failed)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Program abort at step %d\n", this_step);
        return CEV_FAILURE;
    }
    return CEV_SUCCESS;
}

INT evolve_conserv(REAL8 m1,  REAL8 m2, 
           REAL8 s1x, REAL8 s1y, REAL8 s1z, 
           REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 phi0, REAL8 distance,
           REAL8 ecc, REAL8 zeta, REAL8 xi, REAL8 f_min, REAL8 Mf_ref, REAL8 INdeltaT, REAL8 inc,
           HyperParams *hparams, 
           SEOBCoreOutputs *all)
{
    PRINT_LOG_INFO(LOG_INFO, "Conserve mode");
    register INT i;
    INT failed = 0, this_step = 0, status = CEV_FAILURE;
    SpinEOBParams *core = NULL;
    REAL8Vector *ICvalues = NULL;
    REAL8Array *dynamicsAdaS = NULL;
    REAL8Array *dynamicsInverse = NULL;
    REAL8Vector *Jfinal = NULL;
    REAL8Vector *Lhatfinal = NULL;
    COMPLEX16Vector *sigmaQNMlm0 = NULL;
    SphHarmListCAmpPhaseSequence *listhPlm_AdaS = NULL;
    SEOBdynamics *seobdynamicsAdaS = NULL;
    REAL8Vector *seobvalues_final = NULL;

    // Output
    REAL8Vector *tVecPmodes = NULL;
    // SphHarmListCAmpPhaseSequence *listhPlm = NULL;
    // SEOBdynamics *seobdynamicsAdaSHiS = NULL;
    REAL8Vector *alphaJ2P = NULL;
    REAL8Vector *betaJ2P = NULL;
    REAL8Vector *gammaJ2P = NULL;
    SphHarmTimeSeries *hIlm = NULL;
    SphHarmTimeSeries *hJlm = NULL;
    REAL8TimeSeries *hplusTS = NULL;
    REAL8TimeSeries *hcrossTS = NULL;

    UINT nmodes = 5;
    INT modes_lmax = 5;
    INT modes[5][2] = {{2,2}, {2,1}, {3,3}, {4,4}, {5,5}};
    // memset(modes, 0, 2 * nmodes * sizeof(INT));
    REAL8 mTotal = m1 + m2;
    REAL8 mTScaled = mTotal * CST_MTSUN_SI;
    REAL8 EPS_ALIGN = 1.0e-4;
    REAL8 deltaT = INdeltaT / mTScaled;
    PRINT_LOG_INFO(LOG_DEBUG, "deltaT/M = %.16e\n", deltaT);
    REAL8 amp0 = mTotal * CST_MRSUN_SI / distance / 1e6 / CST_PC_SI;
    REAL8 tStepBack = 200.;
    INT flagZframe = FLAG_SEOBNRv4P_ZFRAME_L;
    hparams->flagZframe = flagZframe;
    /*
    *
    *       Check params
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Check params.", this_step);
    if (m1 < 0 || m2 < 0) {failed = 1; goto QUIT;}
    if (m2 > m1)
        SWAP(m1, m2);
    // if (m1/m2 > 100) {failed = 1; goto QUIT;}
    if (sqrt(s1x*s1x+s1y*s1y+s1z*s1z) > 1. ||
        sqrt(s2x*s2x+s2y*s2y+s2z*s2z) > 1)
    {failed = 1; goto QUIT;}

    REAL8 freqMinRad = 0;
    /*Compute the highest initial frequency of the 22 mode */
    XLALEOBHighestInitialFreq(&freqMinRad, mTotal);
    if (f_min > freqMinRad)
    {
        PRINT_LOG_INFO(LOG_WARNING, "Initial frequency is too high, the limit is %4.10f", freqMinRad);
        // failed = 1;
        // goto QUIT;
    }
    /* Check NyquistFrequency */
    UINT ell_max = 4;
    REAL8 INchi1[3] = {s1x, s1y, s1z};
    REAL8 INchi2[3] = {s2x, s2y, s2z};
    if (hparams->Mf_min > hparams->Mf_max && hparams->t_max < 0) {
        status = XLALEOBCheckNyquistFrequency(m1, m2, INchi1, INchi2, INdeltaT, ell_max);
        if (status != CEV_SUCCESS)
        {failed = 1; goto QUIT;}
    }
    /*
    *
    * 
    *       Initialization
    *
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Initialization.", this_step);
    core = CreateSpinEOBParams(m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, ecc, hparams);
    if (!core) {failed = 1; goto QUIT;}
    tStepBack = GET_MAX(tStepBack, core->hParams->tStepBack);

    /*
    *
    * 
    *       Solve Initial Conditions
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Solve Initial Conditions (conserved).", this_step);
    ICvalues = CreateREAL8Vector(14);
    REAL8 MfMin = mTScaled * f_min;
    if (core->hParams && core->hParams->d_ini > 0.0)
    {
        REAL8 d_ini = core->hParams->d_ini; //initial seperation
        REAL8 pr_ini = core->hParams->pr_ini;
        REAL8 pphi_ini = core->hParams->pphi_ini;
        REAL8 ptheta_ini = core->hParams->ptheta_ini;
        REAL8 xSph[3] = {d_ini, 0., 0.};
        REAL8 pSph[3] = {pr_ini, ptheta_ini, pphi_ini};
        REAL8 xCart[3] = {0,0,0};
        REAL8 pCart[3] = {0,0,0};
        if (ptheta_ini > 0)
            core->alignedSpins = TRUE;
        SphericalToCartesian(xCart, pCart, xSph, pSph);
        memset(ICvalues->data, 0, ICvalues->length*sizeof(REAL8));
        memcpy(ICvalues->data, xCart, 3*sizeof(REAL8));
        memcpy(ICvalues->data+3, pCart, 3*sizeof(REAL8));
        memcpy(ICvalues->data+6, core->s1Vec->data, 3*sizeof(REAL8));
        memcpy(ICvalues->data+9, core->s2Vec->data, 3*sizeof(REAL8));
    } else {
        if (ecc > 0.0 && get_egw_flag()) {
            //status = SEOBInitialConditions_egw(ICvalues, MfMin, ecc, core);
            status = SEOBInitialConditions_e_anomaly(ICvalues, MfMin, ecc, zeta, xi, core);
        } else
            status = SEOBInitialConditions_Conserve(ICvalues, MfMin, ecc, core);
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    }
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
#if 0
DEBUG_START;
REAL8Vector dbg_xVec, dbg_pVec, dbg_s1Vec, dbg_s2Vec, dbg_sigKerr, dbg_sigStar;
REAL8 dbg_xdata[3], dbg_pdata[3], dbg_s1data[3], dbg_s2data[3], dbg_sigKdata[3], dbg_sigSdata[3];
dbg_xVec.length = dbg_pVec.length = dbg_s1Vec.length = dbg_s2Vec.length = dbg_sigKerr.length = dbg_sigStar.length = 3;
dbg_xVec.data = dbg_xdata;
dbg_pVec.data = dbg_pdata;
dbg_s1Vec.data = dbg_s1data;
dbg_s2Vec.data = dbg_s2data;
dbg_sigKerr.data = dbg_sigKdata;
dbg_sigStar.data = dbg_sigSdata;
// memcpy(dbg_xdata, ICvalues->data, 3*sizeof(REAL8));
// memcpy(dbg_pdata, ICvalues->data+3, 3*sizeof(REAL8));
dbg_xdata[0] = 4.1659938868533658e1;
dbg_xdata[1] = dbg_xdata[2] = 0.0;

dbg_pdata[0] = -1e-2;
dbg_pdata[1] = 1.5936732938778625e-1;
dbg_pdata[2] = -2.6799608562176783e-12;
// memcpy(dbg_s1data, ICvalues->data+6, 3*sizeof(REAL8));
// memcpy(dbg_s2data, ICvalues->data+9, 3*sizeof(REAL8));
dbg_s1data[0] = 0.0;
dbg_s1data[1] = -1.3888888888888892e-01;
dbg_s1data[2] = 6.2500000000000011e-01;

dbg_s2data[0] = 2.7777777777777779e-03;
dbg_s2data[1] = 0.0;
dbg_s2data[2] = 2.2222222222222227e-02;
EOBCalculateSigmaStar(&dbg_sigStar, core->m1, core->m2, &dbg_s1Vec, &dbg_s2Vec);
EOBCalculateSigmaKerr(&dbg_sigKerr, &dbg_s1Vec, &dbg_s2Vec);
REAL8 dbg_Ham = XLALSimIMRSpinPrecEOBHamiltonian(core->eta, &dbg_xVec, &dbg_pVec, &dbg_s1Vec, &dbg_s2Vec, &dbg_sigKerr, &dbg_sigStar, 1, core->seobCoeffs, core->hParams);
DEBUG_END;
if(1) {failed = 1; goto QUIT;}
#endif

    /*
    *
    * 
    *       Evolve EOB trajectory with adaptive sampling 
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Evolve EOB trajectory with adaptive sampling.", this_step);
    INT retLenAdaS = 0;
    REAL8 tendAdaS;  
    // tendAdaS = 20. / mTScaled;  
    tendAdaS  = GET_CONSERV_TIME;
    REAL8 tstartAdaS = 0.;
    REAL8 deltaT_min = 8.0e-5;
    REAL8 EPS_ABS = 1.0e-8;
    REAL8 EPS_REL = 1.0e-8;
    if (core->alignedSpins)
    {
        EPS_REL = 1.0e-9;
        EPS_ABS = 1.0e-10;
    }
    if (hparams->inEPS_REL > 0.)
        EPS_REL = hparams->inEPS_REL;
    if (hparams->inEPS_ABS > 0.)
        EPS_ABS = hparams->inEPS_ABS;
    status = SEOBIntegrateDynamics_Conserve(&dynamicsAdaS, &retLenAdaS, 
        ICvalues, EPS_ABS, EPS_REL, 
        deltaT, deltaT_min, tstartAdaS, tendAdaS, core, core->alignedSpins);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    PRINT_LOG_INFO(LOG_DEBUG, "[%f, %f]AdaS data length = %d", tstartAdaS, tendAdaS, retLenAdaS);

    status = SEOBComputeExtendedSEOBdynamics_Conserve(&seobdynamicsAdaS, dynamicsAdaS, retLenAdaS, core);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    tVecPmodes = CreateREAL8Vector(retLenAdaS);
    memcpy(tVecPmodes->data, seobdynamicsAdaS->tVec, retLenAdaS*sizeof(REAL8));

    REAL8 tmpMf, tmpe0;
    // print_debug("estimate ecc\n");
    // EstimateEquatorialEccentricity(seobdynamicsAdaS->polarrVec[0], seobdynamicsAdaS->polarpphiVec[0], &tmpMf, &tmpe0, MfMin, ecc, core);
    // print_debug("e0 = %.16e, forb = %.16e\n", tmpe0, tmpMf);

    // SEOBComputeExactEquatorialInitialCondition(ICvalues, MfMin, ecc, core);
    /*
    *
    * 
    *       Get final J/L/spins from HiS dynamics at peak of Omega, compute constant angles EulerI2J
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Get final J/L/spins from HiS dynamics at peak of Omega, compute constant angles EulerI2J.", this_step);
    REAL8 tfinal = seobdynamicsAdaS->tVec[retLenAdaS-1];
    status = SEOBInterpolateDynamicsAtTime(&seobvalues_final, tfinal, seobdynamicsAdaS);

    /* Compute final J from dynamics quantities */
    SEOBJfromDynamics(&Jfinal, seobvalues_final, core);
    /*Compute the L-hat vector. Note that it has unit norm */
    SEOBLhatfromDynamics(&Lhatfinal, seobvalues_final, core);
    REAL8 Jmag = sqrt(inner_product3d(Jfinal->data, Jfinal->data));
    /* Cosine of the angle between L-hat and J. Needed to determine
    * the correct sign of the final spin
    */
    REAL8 cos_angle = inner_product3d(Jfinal->data, Lhatfinal->data) / Jmag;
    /* Compute final-J-frame unit vectors e1J, e2J, e3J=Jfinalhat */
    /* Convention: if (ex, ey, ez) is the initial I-frame, e1J chosen such that ex
    * is in the plane (e1J, e3J) and ex.e1J>0 */
    REAL8Vector e1J, e2J, e3J;
    e1J.length = e2J.length = e3J.length = 3;
    REAL8 e1Jdata[3] = {0.};
    REAL8 e2Jdata[3] = {0.};
    REAL8 e3Jdata[3] = {0.};
    e1J.data = e1Jdata;
    e2J.data = e2Jdata;
    e3J.data = e3Jdata;
    SEOBBuildJframeVectors(&e1J, &e2J, &e3J, Jfinal);
    PRINT_LOG_INFO(LOG_DEBUG, "e1J = (%.16e, %.16e, %.16e)\n", e1Jdata[0], e1Jdata[1], e1Jdata[2]);
    PRINT_LOG_INFO(LOG_DEBUG, "e2J = (%.16e, %.16e, %.16e)\n", e2Jdata[0], e2Jdata[1], e2Jdata[2]);
    PRINT_LOG_INFO(LOG_DEBUG, "e3J = (%.16e, %.16e, %.16e)\n", e3Jdata[0], e3Jdata[1], e3Jdata[2]);

    /* Compute Euler angles from initial I-frame to final-J-frame */
    /* Note: if spins are aligned, the function SEOBEulerI2JFromJframeVectors */
    /* becomes ill-defined - just keep these Euler angles to zero then */
    REAL8 alphaI2J = 0., betaI2J = 0., gammaI2J = 0.;
    if (!core->alignedSpins) 
        SEOBEulerI2JFromJframeVectors(&alphaI2J, &betaI2J, &gammaI2J, &e1J, &e2J, &e3J);
    PRINT_LOG_INFO(LOG_DEBUG, "alphaI2J = %.16e, betaI2J = %.16e, gammaI2J = %.16e\n", alphaI2J, betaI2J, gammaI2J);
    /*
    *
    * 
    *       Compute P-frame amp/phase for all modes
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Build the joined dynamics AdaS+HiS up to attachment, joined P-modes AdaS+HiS+RDpatch", this_step);
    SEOBCalculateSphHarmListhlmAmpPhase_noNQC(&listhPlm_AdaS, modes, nmodes,
                                    seobdynamicsAdaS, core);
    /*
    *
    * 
    *       Compute Euler angles J2P from AdaS and HiS dynamics up to attachment
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Compute Euler angles J2P from AdaS and HiS dynamics up to attachment", this_step);

    /* Compute Euler angles J2P from the dynamics before attachment point */
    /* If SpinsAlmostAligned, all Euler angles are set to 0 */
    status = SEOBEulerJ2PFromDynamics(&alphaJ2P, &betaJ2P, &gammaJ2P, 
                &e1J, &e2J, &e3J,
                retLenAdaS, retLenAdaS,
                seobdynamicsAdaS, core);
    if (status == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "SEOBEulerJ2PFromDynamics failed.");
        failed = 1;
        goto QUIT;
    }

    /*
    *
    * 
    *       Compute modes hJlm on the output time series by rotating and interpolating the modes hPlm
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Compute modes hJlm on the output time series by rotating and interpolating the modes hPlm", this_step);
    // status = SEOBRotateInterpolatehJlmReImFromSphHarmListhPlmAmpPhase(
    //     &hJlm, modes, nmodes, modes_lmax, deltaT, retLenAdaS, tVecPmodes,
    //     listhPlm_AdaS, alphaJ2P, betaJ2P, gammaJ2P);
    status = SEOBPlmListToSphHarmTimeSeries(&hJlm, modes, nmodes, modes_lmax, 
        deltaT, retLenAdaS, listhPlm_AdaS, alphaJ2P, betaJ2P, gammaJ2P);
    if ( status == CEV_FAILURE)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in SEOBRotateInterpolatehJlmReImFromSphHarmListhPlmAmpPhase for mode (l,m) = (2,2).");
        failed = 1;
        goto QUIT;
    }

    /*
    *
    * 
    *       Rotate waveform from J-frame to the output I-frame on timeseries-sampling (constant Wigner coeffs)
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Rotate waveform from J-frame to the output I-frame on timeseries-sampling (constant Wigner coeffs)", this_step);

    /* Rotate waveform from J-frame to I-frame */
    status = SEOBRotatehIlmFromhJlm(&hIlm, hJlm, modes_lmax, alphaI2J, betaI2J, gammaI2J, deltaT);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}
    hIlm->tAttach = 0.0;

    /*
    * 
    * 
    *       Output and cleanup
    * 
    * 
    */

    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Output and cleanup", this_step);
    /* Point the output pointers to the relevant time series and return */

    seobdynamicsAdaS->th22Peak = 0;
    all->dyn = seobdynamicsAdaS;
    all->hLM = hIlm;

QUIT:
    PRINT_LOG_INFO(LOG_INFO, "Program END\n\n");
    STRUCTFREE(core, SpinEOBParams);

    STRUCTFREE(dynamicsAdaS, REAL8Array);
    STRUCTFREE(dynamicsInverse, REAL8Array);
    STRUCTFREE(ICvalues, REAL8Vector);
    // STRUCTFREE(seobdynamicsAdaS, SEOBdynamics);
    STRUCTFREE(seobvalues_final, REAL8Vector);

    STRUCTFREE(Jfinal, REAL8Vector);
    STRUCTFREE(Lhatfinal, REAL8Vector);
    STRUCTFREE(sigmaQNMlm0, COMPLEX16Vector);

    STRUCTFREE(listhPlm_AdaS, SphHarmListCAmpPhaseSequence);

    STRUCTFREE(tVecPmodes, REAL8Vector);
    // STRUCTFREE(listhPlm, SphHarmListCAmpPhaseSequence);
    // STRUCTFREE(seobdynamicsAdaSHiS, SEOBdynamics);

    STRUCTFREE(alphaJ2P, REAL8Vector);
    STRUCTFREE(betaJ2P, REAL8Vector);
    STRUCTFREE(gammaJ2P, REAL8Vector);

    // STRUCTFREE(hIlm, SphHarmTimeSeries);
    STRUCTFREE(hJlm, SphHarmTimeSeries);
    if (failed)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Program abort at step %d\n", this_step);
        return CEV_FAILURE;
    }
    return CEV_SUCCESS;
}


INT evolve_adaptive(REAL8 m1,  REAL8 m2, 
           REAL8 s1x, REAL8 s1y, REAL8 s1z, 
           REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 phi0, REAL8 distance,
           REAL8 ecc, REAL8 zeta, REAL8 xi, REAL8 f_min, REAL8 Mf_ref, REAL8 INdeltaT, REAL8 inc,
           INT is_only22,
           HyperParams *hparams, 
           REAL8TimeSeries **hPlusOut,
           REAL8TimeSeries **hCrossOut,
           SEOBCoreOutputs *all)
{
    register INT i;
    INT failed = 0, this_step = 0, status = CEV_FAILURE;
    SpinEOBParams *core = NULL;
    REAL8Vector *ICvalues = NULL;
    REAL8Vector *ICvaluesHiS = NULL;
    REAL8Vector *seobvalues_tstartHiS = NULL;
    REAL8Vector *seobvalues_tPeakOmega = NULL;
    REAL8Vector *seobvalues_test = NULL;
    REAL8Vector* m1rVec = NULL;
    REAL8Array *dynamicsAdaS = NULL;
    REAL8Array *dynamicsHiS = NULL;
    REAL8Vector *Jfinal = NULL;
    REAL8Vector *Lhatfinal = NULL;
    COMPLEX16Vector *sigmaQNMlm0 = NULL;
    SphHarmListEOBNonQCCoeffs *nqcCoeffsList = NULL;
    SphHarmListCAmpPhaseSequence *listhPlm_HiS = NULL;
    SphHarmListCAmpPhaseSequence *listhPlm_HiSRDpatch = NULL;
    SphHarmListCAmpPhaseSequence *listhPlm_AdaS = NULL;
    SphHarmListCAmpPhaseSequence *listhClm = NULL;
    SEOBdynamics *seobdynamicsAdaS = NULL;
    SEOBdynamics *seobdynamicsHiS = NULL;
    
    REAL8Vector *chi1L_tPeakOmega = NULL;
    REAL8Vector *chi2L_tPeakOmega = NULL;

    // Output
    REAL8Vector *tVecPmodes = NULL;
    SphHarmListCAmpPhaseSequence *listhPlm = NULL;
    SEOBdynamics *seobdynamicsAdaSHiS = NULL;
    REAL8Vector *alphaJ2P = NULL;
    REAL8Vector *betaJ2P = NULL;
    REAL8Vector *gammaJ2P = NULL;
    SphHarmTimeSeries *hIlm = NULL;
    SphHarmTimeSeries *hJlm = NULL;
    REAL8TimeSeries *hplusTS = NULL;
    REAL8TimeSeries *hcrossTS = NULL;

    UINT nmodes = 5;
    INT modes_lmax = 5;
    INT modes[5][2] = {{2,2}, {2,1}, {3,3}, {4,4}, {5,5}};
    // memset(modes, 0, 2 * nmodes * sizeof(INT));
    REAL8 mTotal = m1 + m2;
    REAL8 mTScaled = mTotal * CST_MTSUN_SI;
    REAL8 EPS_ALIGN = 1.0e-4;
    REAL8 deltaT = INdeltaT / mTScaled;
    REAL8 amp0 = mTotal * CST_MRSUN_SI / distance / 1e6 / CST_PC_SI;
    REAL8 tStepBack = 200.;
    INT flagZframe = FLAG_SEOBNRv4P_ZFRAME_L;
    hparams->flagZframe = flagZframe;
    /*
    *
    *       Check params
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Check params.", this_step);
    if (m1 < 0 || m2 < 0) {failed = 1; goto QUIT;}
    if (m2 > m1)
        SWAP(m1, m2);
    if (m1/m2 > 100) {failed = 1; goto QUIT;}
    if (sqrt(s1x*s1x+s1y*s1y+s1z*s1z) > 1. ||
        sqrt(s2x*s2x+s2y*s2y+s2z*s2z) > 1)
    {failed = 1; goto QUIT;}

    REAL8 freqMinRad = 0;
    /*Compute the highest initial frequency of the 22 mode */
    XLALEOBHighestInitialFreq(&freqMinRad, mTotal);
    if (f_min > freqMinRad)
    {
        PRINT_LOG_INFO(LOG_WARNING, "Initial frequency is too high, the limit is %4.10f", freqMinRad);
        // failed = 1;
        // goto QUIT;
    }
    /* Check NyquistFrequency */
    UINT ell_max = 4;
    REAL8 INchi1[3] = {s1x, s1y, s1z};
    REAL8 INchi2[3] = {s2x, s2y, s2z};
    if (hparams->Mf_min > hparams->Mf_max && hparams->t_max < 0) {
        status = XLALEOBCheckNyquistFrequency(m1, m2, INchi1, INchi2, INdeltaT, ell_max);
        if (status != CEV_SUCCESS)
        {failed = 1; goto QUIT;}
    }
    /*
    *
    * 
    *       Initialization
    *
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Initialization.", this_step);
    core = CreateSpinEOBParams(m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, ecc, hparams);
    if (!core) {failed = 1; goto QUIT;}
    // if (1) {failed = 1; goto QUIT;}
    tStepBack = GET_MAX(tStepBack, core->hParams->tStepBack);
    /*
    *
    * 
    *       Solve Initial Conditions
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Solve Initial Conditions.", this_step);
    ICvalues = CreateREAL8Vector(14);
    REAL8 MfMin = mTScaled * f_min;
    if (core->hParams && core->hParams->d_ini > 0.0)
    {
        REAL8 d_ini = core->hParams->d_ini; //initial seperation
        REAL8 pr_ini = core->hParams->pr_ini;
        REAL8 pphi_ini = core->hParams->pphi_ini;
        REAL8 ptheta_ini = core->hParams->ptheta_ini;
        REAL8 xSph[3] = {d_ini, 0., 0.};
        REAL8 pSph[3] = {pr_ini, ptheta_ini, pphi_ini};
        REAL8 xCart[3] = {0,0,0};
        REAL8 pCart[3] = {0,0,0};
        if (ptheta_ini > 0)
            core->alignedSpins = TRUE;
        SphericalToCartesian(xCart, pCart, xSph, pSph);
        memset(ICvalues->data, 0, ICvalues->length*sizeof(REAL8));
        memcpy(ICvalues->data, xCart, 3*sizeof(REAL8));
        memcpy(ICvalues->data+3, pCart, 3*sizeof(REAL8));
        memcpy(ICvalues->data+6, core->s1Vec->data, 3*sizeof(REAL8));
        memcpy(ICvalues->data+9, core->s2Vec->data, 3*sizeof(REAL8));
    } else {
        status = SEOBInitialConditions(ICvalues, MfMin, ecc, core);
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    }
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

    /*
    *
    * 
    *       Evolve EOB trajectory with adaptive sampling 
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Evolve EOB trajectory with adaptive sampling.", this_step);
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
    if (hparams->inEPS_REL > 0.)
        EPS_REL = hparams->inEPS_REL;
    if (hparams->inEPS_ABS > 0.)
        EPS_ABS = hparams->inEPS_ABS;
    status = SEOBIntegrateDynamics_adaptive(&dynamicsAdaS, &retLenAdaS, 
        ICvalues, EPS_ABS, EPS_REL, 
        deltaT, deltaT_min, tstartAdaS, tendAdaS, core, core->alignedSpins);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    PRINT_LOG_INFO(LOG_DEBUG, "AdaS data length = %d", retLenAdaS);
    status = SEOBComputeExtendedSEOBdynamics(&seobdynamicsAdaS, dynamicsAdaS, retLenAdaS, core);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    m1rVec = CreateREAL8Vector(retLenAdaS);
    for (i = 0; i < retLenAdaS; i++) 
    {
        m1rVec->data[i] = -1* seobdynamicsAdaS->polarrVec[i];
    }
    UINT index_6M = FindClosestIndex(m1rVec, -6.0);
    REAL8 time_6M = seobdynamicsAdaS->tVec[index_6M];
    tStepBack = GET_MAX(tStepBack, seobdynamicsAdaS->tVec[retLenAdaS-1] - time_6M + 10*deltaT);

    /*
    *
    * 
    *       (DEBUG MODE) Compare new waveform and old waveform
    * 
    * 
    */
    if (INPUT_DEBUG_FLAG == 1)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "DEBUG_FLAG[%d]:Compare new waveform and old waveform", INPUT_DEBUG_FLAG);
        INT dbg_LL;
        INT dbg_MM;
        COMPLEX16TimeSeries *dbg_hLM_new = NULL;
        COMPLEX16TimeSeries *dbg_hLM_old = NULL;
        CHAR dbg_file_dump_waveform[MAX_STR_LEN];
        FILE *dbg_fout;
        for (dbg_LL = 2; dbg_LL <= 4; dbg_LL++)
        {
            for (dbg_MM = 1; dbg_MM <= dbg_LL; dbg_MM++)
            {
                print_debug("Calculate Mode %d %d\n", dbg_LL, dbg_MM);
                // status = eobcore_object_CalculateWaveform_debug(core, dy, &hLM_new, &hLM_old, LL, MM);
                status = dbg_CalculateWaveformFromDynamicsAdaS(seobdynamicsAdaS, core, dbg_LL, dbg_MM, &dbg_hLM_new, &dbg_hLM_old);
                REAL8 dbg_dt = dbg_hLM_new->deltaT;
                // print_debug("deltaT = %f\n", dbg_dt);
                sprintf(dbg_file_dump_waveform, "debug_h%d%d.dat", dbg_LL, dbg_MM);
                dbg_fout = fopen(dbg_file_dump_waveform, "w");
                for(i=0; i < dbg_hLM_new->data->length; i++)
                {
                    fprintf(dbg_fout, "%.16e %.16e %.16e %.16e %.16e\n", i*dbg_dt, 
                        creal(dbg_hLM_new->data->data[i]), cimag(dbg_hLM_new->data->data[i]),
                        creal(dbg_hLM_old->data->data[i]), cimag(dbg_hLM_old->data->data[i]));
                }
                fclose(dbg_fout);
                DestroyCOMPLEX16TimeSeries(dbg_hLM_new);
                dbg_hLM_new = NULL;
                DestroyCOMPLEX16TimeSeries(dbg_hLM_old);
                dbg_hLM_old = NULL;
            }
        }
        failed = 1;
        goto QUIT;
    }
    /*
    *
    * 
    *       Step back and evolve EOB trajectory at high sampling rate (HiS)
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Step back and evolve EOB trajectory at high sampling rate (HiS).", this_step);
    if (!core->alignedSpins)
    {
        EPS_ABS = 1.0e-8;
        EPS_REL = 1.0e-8;
    } else {
        EPS_REL = 1.0e-9;
        EPS_ABS = 1.0e-10;
    }
    REAL8 deltaTHiS = 1. / 50; /* Fixed at 1/50M */
    REAL8 tstartHiSTarget = seobdynamicsAdaS->tVec[retLenAdaS - 1] - tStepBack;
    INT4 indexstartHiS = retLenAdaS - 1; /* index for the AdaS dynamics */
    while ((indexstartHiS > 0) && (seobdynamicsAdaS->tVec[indexstartHiS] > tstartHiSTarget))
        indexstartHiS--;
    REAL8 tstartHiS = seobdynamicsAdaS->tVec[indexstartHiS];
    PRINT_LOG_INFO(LOG_DEBUG, "Interpolate Dynamics At Time %f (%f)/%f", tstartHiS, tstartHiSTarget, seobdynamicsAdaS->tVec[retLenAdaS - 1]);
    status = SEOBInterpolateDynamicsAtTime(&seobvalues_tstartHiS, tstartHiS, seobdynamicsAdaS);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    ICvaluesHiS = CreateREAL8Vector(14);
    if (!ICvaluesHiS) {failed = 1; goto QUIT;}
    memcpy(ICvaluesHiS->data, &(seobvalues_tstartHiS->data[1]), 14 * sizeof(REAL8));
    /* Integrate again the dynamics with a constant high sampling rate */
    core->prev_dr = 0.; // This is used to check whether dr/dt is increasing
                            // in the stopping condition
    core->termination_reason =
        -999; // This is going to store the termination condition for the
                // high-sampling integration
    INT is_re_evolved = FALSE; 
    INT retLenHiS = 0;
    REAL8 tendHiS = tstartHiS + tStepBack; /* Seems it should be ignored anyway because of
                                    integrator->stopontestonly */
    REAL8 rISCO;
HISR:
    rISCO = get_h_rISCO();
    PRINT_LOG_INFO(LOG_DEBUG, "tStepBack = %g", tStepBack);
    PRINT_LOG_INFO(LOG_DEBUG, "Integrate Dynamics");
    status = SEOBIntegrateDynamics(&dynamicsHiS, &retLenHiS, 
        ICvaluesHiS, EPS_ABS, EPS_REL, 
        deltaTHiS, 0., tstartHiS, tendHiS, core, 1);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    PRINT_LOG_INFO(LOG_DEBUG, "HiSR data length = %d", retLenHiS);
    /* Compute derived quantities for the high-sampling dynamics */
    status = SEOBComputeExtendedSEOBdynamics(&seobdynamicsHiS, dynamicsHiS, retLenHiS, core);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}

    /* Find time of peak of omega for the High-sampling dynamics */
    INT foundPeakOmega = 0;
    REAL8 tPeakOmega = 0.;
    SEOBLocateTimePeakOmega(&tPeakOmega, &foundPeakOmega, dynamicsHiS,
                            seobdynamicsHiS, retLenHiS, core);
    PRINT_LOG_INFO(LOG_DEBUG, "Locate timePeak of omega = %f", tPeakOmega);
    /*
    *
    * 
    *       Get final J/L/spins from HiS dynamics at peak of Omega, compute constant angles EulerI2J
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Get final J/L/spins from HiS dynamics at peak of Omega, compute constant angles EulerI2J.", this_step);
    status = SEOBInterpolateDynamicsAtTime(&seobvalues_tPeakOmega, tPeakOmega, seobdynamicsHiS);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}

    /* Interpolate the dynamics to r=10M. This is used as a
    feducial point to evaluate the final mass and spin fits */
    REAL8Vector timeVec;
    timeVec.length = retLenAdaS;
    timeVec.data = seobdynamicsAdaS->tVec;

    // FindClosestIndex requires an *increasing* function of time, so we use -r
    // instead of r
    UINT index_10M = FindClosestIndex(m1rVec, -10.0);
    REAL8 time_10M = timeVec.data[index_10M];
    status = SEOBInterpolateDynamicsAtTime(&seobvalues_test, time_10M, seobdynamicsAdaS);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}

    /* Compute the timeshift to get the attachment point */
    SEOBLFrameVectors(&chi1L_tPeakOmega, &chi2L_tPeakOmega, 
        seobvalues_tPeakOmega, m1, m2, core->hParams->flagZframe);

    PRINT_LOG_INFO(LOG_DEBUG, "chi1L_tPeakOmega = (%.16e, %.16e, %.16e)\n", 
        chi1L_tPeakOmega->data[0], chi1L_tPeakOmega->data[1], chi1L_tPeakOmega->data[2]);
    PRINT_LOG_INFO(LOG_DEBUG, "chi2L_tPeakOmega = (%.16e, %.16e, %.16e)\n",
        chi2L_tPeakOmega->data[0], chi2L_tPeakOmega->data[1], chi2L_tPeakOmega->data[2]);

    REAL8 Deltat22 = XLALSimIMREOBGetNRSpinPeakDeltaTv4(
        2, 2, m1, m2, chi1L_tPeakOmega->data[2], chi2L_tPeakOmega->data[2], core->hParams);

    /* Determine the time of attachment */
    REAL8 tAttach = tPeakOmega - Deltat22;
    PRINT_LOG_INFO(LOG_INFO, "The location of attachment is %f, Deltat22 = %f", tAttach, Deltat22);
    //Compute NQC window factors
    if (ecc != 0.0)
    {
        status = SEOBCalculateNQCWindowFactorsFromDyn(seobdynamicsAdaS, tAttach, time_6M, tstartHiS, 6., 0, core);
        // core->wWind = 1.;
    }
    
    if (is_re_evolved)
    {
        PRINT_LOG_INFO(LOG_INFO, "re-calculating window parameters; rISCO = %f, rmin = %f", rISCO, seobdynamicsHiS->polarrVec[retLenHiS-1]);
        status = SEOBCalculateNQCWindowFactorsFromDyn(seobdynamicsHiS, tAttach, time_6M, tstartHiS, rISCO, 1, core);
    }

    /* Compute final J from dynamics quantities */
    SEOBJfromDynamics(&Jfinal, seobvalues_tPeakOmega, core);
    /*Compute the L-hat vector. Note that it has unit norm */
    SEOBLhatfromDynamics(&Lhatfinal, seobvalues_tPeakOmega, core);
    REAL8 Jmag = sqrt(inner_product3d(Jfinal->data, Jfinal->data));
    /* Cosine of the angle between L-hat and J. Needed to determine
    * the correct sign of the final spin
    */
    REAL8 cos_angle = inner_product3d(Jfinal->data, Lhatfinal->data) / Jmag;
    /* Compute final-J-frame unit vectors e1J, e2J, e3J=Jfinalhat */
    /* Convention: if (ex, ey, ez) is the initial I-frame, e1J chosen such that ex
    * is in the plane (e1J, e3J) and ex.e1J>0 */
    REAL8Vector e1J, e2J, e3J;
    e1J.length = e2J.length = e3J.length = 3;
    REAL8 e1Jdata[3] = {0.};
    REAL8 e2Jdata[3] = {0.};
    REAL8 e3Jdata[3] = {0.};
    e1J.data = e1Jdata;
    e2J.data = e2Jdata;
    e3J.data = e3Jdata;
    SEOBBuildJframeVectors(&e1J, &e2J, &e3J, Jfinal);

    /* Compute Euler angles from initial I-frame to final-J-frame */
    /* Note: if spins are aligned, the function SEOBEulerI2JFromJframeVectors */
    /* becomes ill-defined - just keep these Euler angles to zero then */
    REAL8 alphaI2J = 0., betaI2J = 0., gammaI2J = 0.;
    if (!core->alignedSpins) 
        SEOBEulerI2JFromJframeVectors(&alphaI2J, &betaI2J, &gammaI2J, &e1J, &e2J, &e3J);

    /*
    *
    * 
    *       Compute P-frame amp/phase for all modes on HiS and compute NQC Coeffs
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Compute P-frame amp/phase for all modes on HiS and compute NQC Coeffs.", this_step);
    // REAL8 wWind_old = core->wWind;
    status = SEOBCalculateSphHarmListNQCCoefficientsV4(
            &nqcCoeffsList, modes, nmodes, tPeakOmega, seobdynamicsHiS,
            core, chi1L_tPeakOmega, chi2L_tPeakOmega);
    if ( status != CEV_SUCCESS) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "NQC computation failed.");
        failed = 1;
        goto QUIT;
    }
    if (!is_re_evolved)
    {
        core->wWind = 1.;
    }
    // Check should we need to re-evolve HiSR EOB
// #if 1
    if (ecc != 0.0 && rISCO > 3 && CheckStopCondition(seobdynamicsHiS, core, nqcCoeffsList, tAttach - seobdynamicsHiS->tVec[0]) != CEV_SUCCESS)
    {
#if 1
        // re-evolve EOB
        SET_RISCO(rISCO-1.);
        // PRINT_LOG_INFO(LOG_INFO, "Re-evolve HiSR EOB...");
        // print_debug("Re-evolve EOB...");
        is_re_evolved = TRUE;
        STRUCTFREE(seobdynamicsHiS, SEOBdynamics);
        STRUCTFREE(Jfinal, REAL8Vector);
        STRUCTFREE(Lhatfinal, REAL8Vector);
        STRUCTFREE(seobvalues_tPeakOmega, REAL8Vector);
        STRUCTFREE(seobvalues_test, REAL8Vector);
        STRUCTFREE(chi1L_tPeakOmega, REAL8Vector);
        STRUCTFREE(chi2L_tPeakOmega, REAL8Vector);
        STRUCTFREE(dynamicsHiS, REAL8Array);
        // STRUCTFREE(seobvalues_tstartHiS, REAL8Vector);
        STRUCTFREE(nqcCoeffsList, SphHarmListEOBNonQCCoeffs);
        // STRUCTFREE(ICvaluesHiS, REAL8Vector);
        goto HISR;
#endif
    }
    /*
    *
    *       Compute P-frame amp/phase for all modes on HiS, now including NQC
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Compute P-frame amp/phase for all modes on HiS, now including NQC.", this_step);
    // /* We now include the NQC when generating modes */
    // flagNQC = 1;

    /* Compute amplitude and phase of the P-frame modes hPlm on high sampling,
    * with NQC */
    SEOBCalculateSphHarmListhlmAmpPhase(&listhPlm_HiS, modes, nmodes,
                                        seobdynamicsHiS, nqcCoeffsList,
                                        core, 1);

    /*
    *
    * 
    *       (DEBUG MODE) Compare new waveform and old waveform
    * 
    * 
    */
    if (INPUT_DEBUG_FLAG == 2)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "DEBUG_FLAG[%d]:Compare new waveform and old waveform with correction", INPUT_DEBUG_FLAG);
        print_debug("tAttach = %f\n", tAttach);
        INT dbg_LL;
        INT dbg_MM;
        COMPLEX16TimeSeries *dbg_hLM_new = NULL;
        COMPLEX16TimeSeries *dbg_hLM_old = NULL;
        CHAR dbg_file_dump_waveform[MAX_STR_LEN];
        FILE *dbg_fout;
        for (dbg_LL = 2; dbg_LL <= 4; dbg_LL++)
        {
            for (dbg_MM = 1; dbg_MM <= dbg_LL; dbg_MM++)
            {
                print_debug("Calculate Mode %d %d\n", dbg_LL, dbg_MM);
                // status = eobcore_object_CalculateWaveform_debug(core, dy, &hLM_new, &hLM_old, LL, MM);
                status = dbg_CalculateWaveformFromDynamicsAdaS(seobdynamicsAdaS, core, dbg_LL, dbg_MM, &dbg_hLM_new, &dbg_hLM_old);
                REAL8 dbg_dt = dbg_hLM_new->deltaT;
                sprintf(dbg_file_dump_waveform, "debug_h%d%d.dat", dbg_LL, dbg_MM);
                dbg_fout = fopen(dbg_file_dump_waveform, "w");
                for(i=0; i < dbg_hLM_new->data->length; i++)
                {
                    fprintf(dbg_fout, "%.16e %.16e %.16e %.16e %.16e\n", i*deltaT - tAttach, 
                        creal(dbg_hLM_new->data->data[i]), cimag(dbg_hLM_new->data->data[i]),
                        creal(dbg_hLM_old->data->data[i]), cimag(dbg_hLM_old->data->data[i]));
                }
                fclose(dbg_fout);
                DestroyCOMPLEX16TimeSeries(dbg_hLM_new);
                dbg_hLM_new = NULL;
                DestroyCOMPLEX16TimeSeries(dbg_hLM_old);
                dbg_hLM_old = NULL;
            }
        }
        failed = 1;
        goto QUIT;
    }

    if (INPUT_DEBUG_FLAG == 3)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "DEBUG_FLAG[%d]:Calculate NQC correction series", INPUT_DEBUG_FLAG);
        INT dbg_LL = 2;
        INT dbg_MM = 2;
        COMPLEX16TimeSeries *dbg_nqc = NULL;
        CHAR dbg_file_dump_nqc[MAX_STR_LEN];
        FILE *dbg_fout;
        sprintf(dbg_file_dump_nqc, "debug_nqc_h%d%d.dat", dbg_LL, dbg_MM);
        dbg_fout = fopen(dbg_file_dump_nqc, "w");
        print_debug("Calculate Mode %d %d\n", dbg_LL, dbg_MM);
        status = dbg_CalculateNQCTimeSeries(seobdynamicsHiS, core, nqcCoeffsList, dbg_LL, dbg_MM, &dbg_nqc);
        REAL8 dbg_dt = dbg_nqc->deltaT;
        REAL8 nWind;
        print_debug("tWind = %f, wWind = %f", core->tWind, core->wWind);
        for(i=0; i < dbg_nqc->data->length; i++)
        {
            nWind = NQCWindow(seobdynamicsHiS->tVec[i], core->tWind, core->wWind);
            // print_debug("nWind = %f\n", nWind);
            fprintf(dbg_fout, "%.16e %.16e %.16e %.16e\n", i*dbg_dt + seobdynamicsHiS->tVec[0] - tAttach, 
                creal(dbg_nqc->data->data[i]), cimag(dbg_nqc->data->data[i]), nWind);
        }
        fclose(dbg_fout);
        STRUCTFREE(dbg_nqc, COMPLEX16TimeSeries);
        // failed = 1;
        // goto QUIT;
    }
    /*
    *
    * 
    *       Attach RD to the P-frame waveform
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Attach RD to the P-frame waveform.", this_step);
    REAL8 finalMass = 0., finalSpin = 0.;
    status = SEOBGetFinalSpinMass(&finalMass, &finalSpin, seobvalues_test, core);
    if (status != CEV_SUCCESS) 
    {failed = 1; goto QUIT;}

    /* The function above returns only the magnitude of the spin.
    *  We pick the direction based on whether Lhat \cdot J is positive
    *  or negative */
    if (cos_angle < 0)
    {
        finalSpin *= -1;
    }
    /* finalSpin interpolation is available only between -0.9996 and 0.9996 */
    /* Set finalSpin to +/- 0.9996 if it is out of this range */
    if (finalSpin < -0.9996)
        finalSpin = -0.9996;
    if (finalSpin > 0.9996)
        finalSpin = 0.9996;

    PRINT_LOG_INFO(LOG_INFO, "final mass = %e, final spin = %e, cos_angle = %e", finalMass, finalSpin, cos_angle);

    /* Estimate leading QNM , to determine how long the ringdown
    * patch will be */
    /* NOTE: XLALSimIMREOBGenerateQNMFreqV2Prec returns the complex frequency in
    * physical units... */
    COMPLEX16Vector sigmaQNM220estimatephysicalVec;
    COMPLEX16 sigmaQNM220estimatephysical = 0.;
    sigmaQNM220estimatephysicalVec.length = 1;
    sigmaQNM220estimatephysicalVec.data = &sigmaQNM220estimatephysical;
    status = XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec(&sigmaQNM220estimatephysicalVec, m1,
                                                m2, finalMass, finalSpin, 2,
                                                2, 1);
    if (status == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec for mode (l,m) = (2,2).");
        failed = 1;
        goto QUIT;
    }
    COMPLEX16 sigmaQNM220estimate = mTScaled * sigmaQNM220estimatephysical;

    /* Length of RD patch, 40 e-folds of decay of the estimated QNM220 */
    UINT retLenRDPatch =
    (UINT)ceil(EFOLDS / (cimag(sigmaQNM220estimate) * deltaTHiS));


    /* Attach RD to the P-frame modes */
    /* Vector holding the values of the 0-th overtone QNM complex frequencies for
    * the modes (l,m) */
    // NOTE: the QNM complex frequencies are computed inside
    // SEOBAttachRDToSphHarmListhPlm
    status = SEOBAttachRDToSphHarmListhPlm(
        &listhPlm_HiSRDpatch, &sigmaQNMlm0, modes, nmodes, finalMass, finalSpin,
        listhPlm_HiS, deltaTHiS, retLenHiS, retLenRDPatch, tAttach,
        seobvalues_tPeakOmega, seobdynamicsHiS, core);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}

    /*
    *
    * 
    *       Build the joined dynamics AdaS+HiS up to attachment, joined P-modes AdaS+HiS+RDpatch 
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Build the joined dynamics AdaS+HiS up to attachment, joined P-modes AdaS+HiS+RDpatch", this_step);

    /* Compute amplitude and phase of the P-frame modes hPlm on adaptive sampling,
    * with NQC */
    // flagNQC = 1;
    SEOBCalculateSphHarmListhlmAmpPhase(&listhPlm_AdaS, modes, nmodes,
                                    seobdynamicsAdaS, nqcCoeffsList,
                                    core, 1);

    /* Vector of times for the P-modes, that will be used for interpolation:
    * joining AdaS and HiS+RDpatch */
    UINT retLenPmodes = 0;
    /* First junction at indexAdaSHiS, tAdaSHiS */
    UINT indexJoinHiS = 0;
    REAL8 tJoinHiS = 0.;
    /* Second junction at seobdynamicsAdaS, tJoinAttach */
    UINT indexJoinAttach = 0;
    REAL8 tJoinAttach = 0.;
    /* Construct the joined vector of times (AdaS+HiS+RDpatch) and keep the
    * jonction indices and times */
    status = SEOBJoinTimeVector(&tVecPmodes, &retLenPmodes, &tJoinHiS, &indexJoinHiS,
                    &tJoinAttach, &indexJoinAttach, retLenRDPatch, deltaTHiS,
                    tstartHiS, tAttach, seobdynamicsAdaS, seobdynamicsHiS);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}
    PRINT_LOG_INFO(LOG_DEBUG, "tJoinAttach = %.16e, indexJoinAttach = %u", tJoinAttach, indexJoinAttach);
    /* Copy dynamics from AdaS<HiS and HiS<tAttach to form joined dynamics, ending
    * at the last time sample <tAttach */
    // NOTE: we cut the dynamics at tAttach, as we will extend the Euler
    // angles for t>=tAttach -- but we could also choose to finish at tPeakOmega
    // which is used for final-J and for the final mass/spin fit

    status = SEOBJoinDynamics(&seobdynamicsAdaSHiS, seobdynamicsAdaS, seobdynamicsHiS,
                indexJoinHiS, indexJoinAttach);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}

    /* Copy waveform modes from AdaS and HiS+RDpatch - adjusting 2pi-phase shift
    * at the junction point AdaS/HiS */
    status = SEOBJoinSphHarmListhlm(&listhPlm, listhPlm_AdaS, listhPlm_HiSRDpatch, modes,
                        nmodes, indexstartHiS);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}

    /* Get the time of the frame-invariant amplitude peak */
    REAL8 tPeak = 0;
    UINT indexPeak = 0;
    // NOTE: peak amplitude using l=2 only: h22 required, and h21 used if present
    status = SEOBAmplitudePeakFromAmp22Amp21(&tPeak, &indexPeak, listhPlm, modes, nmodes,
                                tVecPmodes);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}
    PRINT_LOG_INFO(LOG_DEBUG, "tPeak = %.16e, indexPeak = %u", tPeak, indexPeak);
    /*
    *
    * 
    *       Compute Euler angles J2P from AdaS and HiS dynamics up to attachment
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Compute Euler angles J2P from AdaS and HiS dynamics up to attachment", this_step);
    /* Compute Euler angles J2P from the dynamics before attachment point */
    /* If SpinsAlmostAligned, all Euler angles are set to 0 */
    status = SEOBEulerJ2PFromDynamics(&alphaJ2P, &betaJ2P, &gammaJ2P, 
                &e1J, &e2J, &e3J,
                retLenPmodes, indexJoinAttach,
                seobdynamicsAdaSHiS, core);
    if (status == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "SEOBEulerJ2PFromDynamics failed.");
        failed = 1;
        goto QUIT;
    }

    /*
    *
    * 
    *       Compute Euler angles J2P extension after attachment
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Compute Euler angles J2P extension after attachment", this_step);

    /* Compute Euler angles J2P according to the prescription flagEulerextension
    * after attachment point */
    /* If SpinsAlmostAligned, all Euler angles are set to 0 */
    /* NOTE: Regardless of the mode content of hPlm, the frame extension at the
    * moment is based on sigmaQNM22, sigmaQNM21 */
    COMPLEX16 sigmaQNM220 = 0., sigmaQNM210 = 0.;
    COMPLEX16Vector sigmaQNM220physicalVec, sigmaQNM210physicalVec;
    sigmaQNM220physicalVec.length = 1;
    sigmaQNM210physicalVec.length = 1;
    COMPLEX16 sigmaQNM220physicalval = 0.;
    COMPLEX16 sigmaQNM210physicalval = 0.;
    sigmaQNM220physicalVec.data = &sigmaQNM220physicalval;
    sigmaQNM210physicalVec.data = &sigmaQNM210physicalval;
    /* NOTE: XLALSimIMREOBGenerateQNMFreqV2Prec returns the complex frequency in
    * physical units... */
    status = XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec(&sigmaQNM220physicalVec, m1,
                                                m2, finalMass, finalSpin, 2,
                                                2, 1);
    if (status == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec for mode (l,m) = (2,2).");
        failed = 1;
        goto QUIT;
    }
    status = XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec(&sigmaQNM210physicalVec, m1,
                                                m2, finalMass, finalSpin, 2,
                                                1, 1);
    if (status == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec for mode (l,m) = (2,1).");
        failed = 1;
        goto QUIT;
    }
    sigmaQNM220 = mTScaled * sigmaQNM220physicalVec.data[0];
    sigmaQNM210 = mTScaled * sigmaQNM210physicalVec.data[0];
    INT flip = 1;
    if (cos_angle < 0)
        flip = -1;
    // flagEulerextension = 0
    status = SEOBEulerJ2PPostMergerExtension(
        alphaJ2P, betaJ2P, gammaJ2P, sigmaQNM220, sigmaQNM210, tVecPmodes,
        retLenPmodes, indexJoinAttach, core, flip);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}

    /*
    *
    * 
    *       Compute modes hJlm on the output time series by rotating and interpolating the modes hPlm
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Compute modes hJlm on the output time series by rotating and interpolating the modes hPlm", this_step);

    /* Determine the length of the fixed-sampling output time series */
    UINT retLenTS = floor(((tVecPmodes)->data[retLenPmodes - 1] - (tVecPmodes)->data[0]) / deltaT);

    /* Rotate waveform from P-frame to J-frame */
    // flagSymmetrizehPlminusm = 1
    status = SEOBRotateInterpolatehJlmReImFromSphHarmListhPlmAmpPhase(
        &hJlm, &listhClm, modes, nmodes, modes_lmax, deltaT, retLenTS, tVecPmodes,
        listhPlm, alphaJ2P, betaJ2P, gammaJ2P);
    if ( status == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in SEOBRotateInterpolatehJlmReImFromSphHarmListhPlmAmpPhase for mode (l,m) = (2,2).");
        failed = 1;
        goto QUIT;
    }

    /*
    *
    * 
    *       Rotate waveform from J-frame to the output I-frame on timeseries-sampling (constant Wigner coeffs)
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Rotate waveform from J-frame to the output I-frame on timeseries-sampling (constant Wigner coeffs)", this_step);

    /* Rotate waveform from J-frame to I-frame */
    status = SEOBRotatehIlmFromhJlm(&hIlm, hJlm, modes_lmax, alphaI2J, betaI2J, gammaI2J, deltaT);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}
    hIlm->tAttach = tAttach;
    PRINT_LOG_INFO(LOG_DEBUG, "alphaI2J, betaI2J, gammaI2J = %.16e, %.16e, %.16e", alphaI2J, betaI2J, gammaI2J);

    /*
    *
    * 
    *       Compute h20 from I-frame waveform on timeseries sampling
    * 
    * 
    */
    // this_step++;
    // PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Compute h20 from I-frame waveform on timeseries sampling", this_step);

    // status = SEOBCalculateh20(hIlm);
    // {failed = 1; goto QUIT;}

    /*
    *
    * 
    *       Compute hplus, hcross from I-frame waveform on timeseries sampling
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Compute hplus, hcross from I-frame waveform on timeseries sampling", this_step);

    /* GPS time for output time series and modes */
    // LIGOTimeGPS tGPS = LIGOTIMEGPSZERO;
    // XLALGPSAdd(&tGPS, -mTScaled *
    //                     tPeak); /* tPeak converted back to dimensionfull time */
    REAL8 tGPS = -mTScaled * tPeak;
    /* Create output timeseries for hplus, hcross */
    /* Use the dimensionfull INdeltaT (s) as time step */
    // hplusTS = CreateREAL8TimeSeries("H_PLUS", &tGPS, 0.0, INdeltaT,
    //                                 &lalStrainUnit, retLenTS);
    // hcrossTS = CreateREAL8TimeSeries("H_CROSS", &tGPS, 0.0, INdeltaT,
    //                                 &lalStrainUnit, retLenTS);
    hplusTS = CreateREAL8TimeSeries(tGPS, INdeltaT, retLenTS);
    hcrossTS = CreateREAL8TimeSeries(tGPS, INdeltaT, retLenTS);
    /* Compute hplus, hcross from hIlm */
    // NOTE: azimuthal angle of the observer entering the -2Ylm is pi/2-phi
    // according to LAL conventions
    status = SEOBComputehplushcrossFromhIlm(hplusTS, hcrossTS, modes_lmax, hIlm, amp0,
        inc, phi0, is_only22);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}

    /*
    *
    * 
    *       Output and cleanup
    * 
    * 
    */

    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Output and cleanup", this_step);
    /* Point the output pointers to the relevant time series and return */
    (*hPlusOut) = hplusTS;
    (*hCrossOut) = hcrossTS;
    seobdynamicsAdaSHiS->th22Peak = tAttach;
    all->dyn = seobdynamicsAdaSHiS;
    all->hLM = hIlm;
    all->Plm = listhClm;
#if 0
    /* Output vector gathering quantities related to merger (similar to previous
    * AttachParams) */
    /* Format: tPeakOmega tAttach tPeak Jfinalx Jfinaly Jfinalz finalMassfit
    * finalSpinfit termination_reason [sigmaQNMlm0Re sigmaQNMlm0Im for lm in
    * modes] */
    /* NOTE: the size of this output vector depends on the number of modes, due to
    * the sigmaQNM */
    if (!((*mergerParams) = CreateREAL8Vector(9 + 2 * nmodes))) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failed to allocate REAL8Vector mergerParams.");
        failed = 1;
        goto QUIT;
    }
    (*mergerParams)->data[0] = tPeakOmega;
    (*mergerParams)->data[1] = tAttach;
    (*mergerParams)->data[2] = tPeak;
    (*mergerParams)->data[3] = Jfinal->data[0];
    (*mergerParams)->data[4] = Jfinal->data[1];
    (*mergerParams)->data[5] = Jfinal->data[2];
    (*mergerParams)->data[6] = finalMass;
    (*mergerParams)->data[7] = finalSpin;
    (*mergerParams)->data[8] = seobParams.termination_reason;
    for (UINT4 nmode = 0; nmode < nmodes; nmode++) {
    (*mergerParams)->data[9 + 2 * nmode] = creal(sigmaQNMlm0->data[nmode]);
    (*mergerParams)->data[9 + 2 * nmode + 1] = cimag(sigmaQNMlm0->data[nmode]);
    }


    /* Additional outputs */

    /* Dynamics */
    // NOTE: casting to REAL8Vector due to the SWIG wrapping
    *seobdynamicsAdaSVector =
    XLALCreateREAL8Vector(v4PdynamicsVariables * seobdynamicsAdaS->length);
    *seobdynamicsHiSVector =
    XLALCreateREAL8Vector(v4PdynamicsVariables * seobdynamicsHiS->length);
    *seobdynamicsAdaSHiSVector =
    XLALCreateREAL8Vector(v4PdynamicsVariables * seobdynamicsAdaSHiS->length);
    memcpy((*seobdynamicsAdaSVector)->data, seobdynamicsAdaS->array->data,
    (v4PdynamicsVariables * seobdynamicsAdaS->length) * sizeof(REAL8));
    memcpy((*seobdynamicsHiSVector)->data, seobdynamicsHiS->array->data,
    (v4PdynamicsVariables * seobdynamicsHiS->length) * sizeof(REAL8));
    memcpy((*seobdynamicsAdaSHiSVector)->data, seobdynamicsAdaSHiS->array->data,
    (v4PdynamicsVariables * seobdynamicsAdaSHiS->length) * sizeof(REAL8));

    /* Modes in the P-frame */
    // NOTE: casting to REAL8Vector due to the SWIG wrapping
    // NOTE: in the output, real amplitude instead of complex envelope
    /* (2,2) */
    *hP22_amp = XLALCreateREAL8Vector(retLenPmodes);
    *hP22_phase = XLALCreateREAL8Vector(retLenPmodes);
    SphHarmListCAmpPhaseSequence *hP22 =
    SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 2, 2);
    if (hP22 == NULL) {
    memset((*hP22_amp)->data, 0, retLenPmodes * sizeof(REAL8));
    memset((*hP22_phase)->data, 0, retLenPmodes * sizeof(REAL8));
    } else {
    memcpy((*hP22_amp)->data, hP22->campphase->camp_real->data,
        retLenPmodes * sizeof(REAL8));
    memcpy((*hP22_phase)->data, hP22->campphase->phase->data,
        retLenPmodes * sizeof(REAL8));
    }
    /* (2,1) */
    *hP21_amp = XLALCreateREAL8Vector(retLenPmodes);
    *hP21_phase = XLALCreateREAL8Vector(retLenPmodes);
    SphHarmListCAmpPhaseSequence *hP21 =
    SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 2, 1);
    if (hP21 == NULL) {
    memset((*hP21_amp)->data, 0, retLenPmodes * sizeof(REAL8));
    memset((*hP21_phase)->data, 0, retLenPmodes * sizeof(REAL8));
    } else {
    memcpy((*hP21_amp)->data, hP21->campphase->camp_real->data,
        retLenPmodes * sizeof(REAL8));
    memcpy((*hP21_phase)->data, hP21->campphase->phase->data,
        retLenPmodes * sizeof(REAL8));
    }
    /* (3,3) */
    *hP33_amp = XLALCreateREAL8Vector(retLenPmodes);
    *hP33_phase = XLALCreateREAL8Vector(retLenPmodes);
    SphHarmListCAmpPhaseSequence *hP33 =
    SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 3, 3);
    if (hP33 == NULL) {
    memset((*hP33_amp)->data, 0, retLenPmodes * sizeof(REAL8));
    memset((*hP33_phase)->data, 0, retLenPmodes * sizeof(REAL8));
    } else {
    memcpy((*hP33_amp)->data, hP33->campphase->camp_real->data,
        retLenPmodes * sizeof(REAL8));
    memcpy((*hP33_phase)->data, hP33->campphase->phase->data,
        retLenPmodes * sizeof(REAL8));
    }
    /* (4,4) */
    *hP44_amp = XLALCreateREAL8Vector(retLenPmodes);
    *hP44_phase = XLALCreateREAL8Vector(retLenPmodes);
    SphHarmListCAmpPhaseSequence *hP44 =
    SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 4, 4);
    if (hP44 == NULL) {
    memset((*hP44_amp)->data, 0, retLenPmodes * sizeof(REAL8));
    memset((*hP44_phase)->data, 0, retLenPmodes * sizeof(REAL8));
    } else {
    memcpy((*hP44_amp)->data, hP44->campphase->camp_real->data,
        retLenPmodes * sizeof(REAL8));
    memcpy((*hP44_phase)->data, hP44->campphase->phase->data,
        retLenPmodes * sizeof(REAL8));
    }
    /* (5,5) */
    *hP55_amp = XLALCreateREAL8Vector(retLenPmodes);
    *hP55_phase = XLALCreateREAL8Vector(retLenPmodes);
    SphHarmListCAmpPhaseSequence *hP55 =
    SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 5, 5);
    if (hP55 == NULL) {
    memset((*hP55_amp)->data, 0, retLenPmodes * sizeof(REAL8));
    memset((*hP55_phase)->data, 0, retLenPmodes * sizeof(REAL8));
    } else {
    memcpy((*hP55_amp)->data, hP55->campphase->camp_real->data,
        retLenPmodes * sizeof(REAL8));
    memcpy((*hP55_phase)->data, hP55->campphase->phase->data,
        retLenPmodes * sizeof(REAL8));
    }
#endif

QUIT:
    PRINT_LOG_INFO(LOG_INFO, "Program END\n\n");
    STRUCTFREE(core, SpinEOBParams);

    STRUCTFREE(dynamicsAdaS, REAL8Array);
    STRUCTFREE(ICvalues, REAL8Vector);
    STRUCTFREE(seobdynamicsAdaS, SEOBdynamics);

    STRUCTFREE(dynamicsHiS, REAL8Array);
    STRUCTFREE(ICvaluesHiS, REAL8Vector);
    STRUCTFREE(seobvalues_tstartHiS, REAL8Vector);
    STRUCTFREE(seobdynamicsHiS, SEOBdynamics);

    STRUCTFREE(seobvalues_tPeakOmega, REAL8Vector);
    STRUCTFREE(seobvalues_test, REAL8Vector);
    STRUCTFREE(m1rVec,REAL8Vector);
    STRUCTFREE(Jfinal, REAL8Vector);
    STRUCTFREE(Lhatfinal, REAL8Vector);
    STRUCTFREE(sigmaQNMlm0, COMPLEX16Vector);

    STRUCTFREE(chi2L_tPeakOmega, REAL8Vector);
    STRUCTFREE(chi1L_tPeakOmega, REAL8Vector);
    STRUCTFREE(nqcCoeffsList, SphHarmListEOBNonQCCoeffs);

    STRUCTFREE(listhPlm_HiS, SphHarmListCAmpPhaseSequence);
    STRUCTFREE(listhPlm_HiSRDpatch, SphHarmListCAmpPhaseSequence);
    STRUCTFREE(listhPlm_AdaS, SphHarmListCAmpPhaseSequence);

    STRUCTFREE(tVecPmodes, REAL8Vector);
    STRUCTFREE(listhPlm, SphHarmListCAmpPhaseSequence);
    // STRUCTFREE(seobdynamicsAdaSHiS, SEOBdynamics);

    STRUCTFREE(alphaJ2P, REAL8Vector);
    STRUCTFREE(betaJ2P, REAL8Vector);
    STRUCTFREE(gammaJ2P, REAL8Vector);

    // STRUCTFREE(hIlm, SphHarmTimeSeries);
    STRUCTFREE(hJlm, SphHarmTimeSeries);
    if (failed)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Program abort at step %d\n", this_step);
        return CEV_FAILURE;
    }
    return CEV_SUCCESS;
}






/* ------------------------------------------------------------------------------------------------ */
/*                                                                                                  */
/*                                                                                                  */
/*                                                                                                  */
/*                                                                                                  */
/*                                                                                                  */
/*                                                                                                  */
/*                                                                                                  */
/*                                    Spin-Aligned Faster version                                   */
/*                                                                                                  */
/*                                                                                                  */
/*                                                                                                  */
/*                                                                                                  */
/*                                                                                                  */
/*                                                                                                  */
/*                                                                                                  */
/* ------------------------------------------------------------------------------------------------ */

INT evolve_SA(REAL8 m1,  REAL8 m2, 
           REAL8 s1z, 
           REAL8 s2z,
           REAL8 ecc, REAL8 zeta, REAL8 f_min, REAL8 Mf_ref, REAL8 INdeltaT,
           HyperParams *hparams, 
           REAL8Vector **tRetVec,
           SphHarmListCAmpPhaseSequence **hlm,
           int is_noringdown, int is_dyn_debug,
           SEOBSAdynamics **dyn_debug)
{
    register INT i;
    INT failed = 0, this_step = 0, status = CEV_FAILURE;
    SpinEOBParams *core = NULL;
    REAL8Vector *ICvalues = NULL;
    REAL8Vector *ICvaluesHiS = NULL;
    REAL8Vector *ICvalues_SA = NULL;
    REAL8Vector *ICvaluesHiS_SA = NULL;
    REAL8Vector *seobvalues_tstartHiS = NULL;
    REAL8Vector *seobvalues_tPeakOmega = NULL;
    REAL8Vector *seobvalues_test = NULL;
    REAL8Vector* m1rVec = NULL;
    REAL8Array *dynamicsAdaS = NULL;
    REAL8Array *dynamicsInverse = NULL;
    REAL8Array *dynamicsHiS = NULL;
    REAL8Vector *Jfinal = NULL;
    REAL8Vector *Lhatfinal = NULL;
    COMPLEX16Vector *sigmaQNMlm0 = NULL;
    SphHarmListEOBNonQCCoeffs *nqcCoeffsList = NULL;
    SphHarmListCAmpPhaseSequence *listhPlm_HiS = NULL;
    SphHarmListCAmpPhaseSequence *listhPlm_HiSRDpatch = NULL;
    SphHarmListCAmpPhaseSequence *listhPlm_AdaS = NULL;
    SEOBSAdynamics *seobdynamicsAdaS = NULL;
    SEOBSAdynamics *seobdynamicsHiS = NULL;
    
    REAL8Vector *chi1L_tPeakOmega = NULL;
    REAL8Vector *chi2L_tPeakOmega = NULL;

    // Output
    REAL8Vector *tVecPmodes = NULL;
    SphHarmListCAmpPhaseSequence *listhPlm = NULL;
    SEOBSAdynamics *seobdynamicsAdaSHiS = NULL;
    // REAL8Vector *alphaJ2P = NULL;
    // REAL8Vector *betaJ2P = NULL;
    // REAL8Vector *gammaJ2P = NULL;
    // SphHarmTimeSeries *hIlm = NULL;
    // SphHarmTimeSeries *hJlm = NULL;
    // REAL8TimeSeries *hplusTS = NULL;
    // REAL8TimeSeries *hcrossTS = NULL;

    UINT nmodes = 5;
    INT modes_lmax = 5;
    INT modes[5][2] = {{2,2}, {2,1}, {3,3}, {4,4}, {5,5}};
    // memset(modes, 0, 2 * nmodes * sizeof(INT));
    REAL8 mTotal = m1 + m2;
    REAL8 mTScaled = mTotal * CST_MTSUN_SI;
    REAL8 EPS_ALIGN = 1.0e-4;
    REAL8 deltaT = INdeltaT / mTScaled;
    // REAL8 amp0 = mTotal * CST_MRSUN_SI / distance / 1e6 / CST_PC_SI;
    REAL8 tStepBack = 200.;
    INT flagZframe = FLAG_SEOBNRv4P_ZFRAME_L;
    hparams->flagZframe = flagZframe;
    /*
    *
    *       Check params
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Check params.", this_step);
    if (m1 < 0 || m2 < 0) {failed = 1; goto QUIT;}
    if (m2 > m1)
        SWAP(m1, m2);
    // if (m1/m2 > 100) {failed = 1; goto QUIT;}
    if (fabs(s1z) > 1. ||
        fabs(s2z) > 1)
    {failed = 1; goto QUIT;}

    REAL8 freqMinRad = 0;
    /*Compute the highest initial frequency of the 22 mode */
    XLALEOBHighestInitialFreq(&freqMinRad, mTotal);
    if (f_min > freqMinRad)
    {
        PRINT_LOG_INFO(LOG_WARNING, "Initial frequency is too high, the limit is %4.10f", freqMinRad);
        // failed = 1;
        // goto QUIT;
    }
    /* Check NyquistFrequency */
    UINT ell_max = 4;
    REAL8 INchi1[3] = {0, 0, s1z};
    REAL8 INchi2[3] = {0, 0, s2z};
    if (hparams->Mf_min > hparams->Mf_max && hparams->t_max < 0) {
        status = XLALEOBCheckNyquistFrequency(m1, m2, INchi1, INchi2, INdeltaT, ell_max);
        if (status != CEV_SUCCESS)
        {failed = 1; goto QUIT;}
    }
    /*
    *
    * 
    *       Initialization
    *
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ S-A Initialization.", this_step);
    core = CreateSpinEOBParams(m1, m2, 0, 0, s1z, 0, 0, s2z, ecc, hparams);
    if (!core) {failed = 1; goto QUIT;}
    tStepBack = GET_MAX(tStepBack, core->hParams->tStepBack);
    /*
    *
    * 
    *       Solve Initial Conditions
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Solve Initial Conditions.", this_step);
    ICvalues = CreateREAL8Vector(14);
    REAL8 MfMin = mTScaled * f_min;
    if (core->hParams && core->hParams->d_ini > 0.0)
    {
        REAL8 d_ini = core->hParams->d_ini; //initial seperation
        REAL8 pr_ini = core->hParams->pr_ini;
        REAL8 pphi_ini = core->hParams->pphi_ini;
        REAL8 ptheta_ini = core->hParams->ptheta_ini;
        REAL8 xSph[3] = {d_ini, 0., 0.};
        REAL8 pSph[3] = {pr_ini, 0.0, pphi_ini};
        REAL8 xCart[3] = {0,0,0};
        REAL8 pCart[3] = {0,0,0};
        SphericalToCartesian(xCart, pCart, xSph, pSph);
        memset(ICvalues->data, 0, ICvalues->length*sizeof(REAL8));
        memcpy(ICvalues->data, xCart, 3*sizeof(REAL8));
        memcpy(ICvalues->data+3, pCart, 3*sizeof(REAL8));
        memcpy(ICvalues->data+6, core->s1Vec->data, 3*sizeof(REAL8));
        memcpy(ICvalues->data+9, core->s2Vec->data, 3*sizeof(REAL8));
    } else {
        if (ecc > 0.0 && get_egw_flag()) {
            //status = SEOBInitialConditions_egw(ICvalues, MfMin, ecc, core);
            status = SEOBInitialConditions_e_anomaly(ICvalues, Mf_ref, ecc, zeta, 0, core);
            // if (status != CEV_SUCCESS && ecc < 0.05)
            // {
            //     status = SEOBInitialConditions(ICvalues, Mf_ref, ecc, core);
            // }
        } else 
            status = SEOBInitialConditions(ICvalues, Mf_ref, ecc, core);
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    }
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
    ICvalues_SA = CreateREAL8Vector(4);
    {
        REAL8 temp_r = sqrt(ICvalues->data[0] * ICvalues->data[0] +
                            ICvalues->data[1] * ICvalues->data[1] +
                            ICvalues->data[2] * ICvalues->data[2]);
        REAL8 temp_phi = ICvalues->data[12];

        ICvalues_SA->data[0] = temp_r;   // General form of r
        ICvalues_SA->data[1] = temp_phi; // phi
        ICvalues_SA->data[2] = 
            ICvalues->data[3] * cos(temp_phi) +
            ICvalues->data[4] * sin(temp_phi); // p_r^*
        ICvalues_SA->data[3] =
            temp_r * (ICvalues->data[4] * cos(temp_phi) -
                    ICvalues->data[3] * sin(temp_phi)); // p_phi
    }
#if 0
    {
        REAL8 tmp_omega0, tmp_e0;
        status = SEOBConvertInitConditionFromGeneralToV1(ICvalues, core, &tmp_omega0, &tmp_e0);
        print_debug("omega0, e0 = %.16e, %.16e\n", tmp_omega0, tmp_e0);
    }
#endif

#if 0
DEBUG_START;
REAL8Vector dbg_xVec, dbg_pVec, dbg_s1Vec, dbg_s2Vec, dbg_sigKerr, dbg_sigStar;
REAL8 dbg_xdata[3], dbg_pdata[3], dbg_s1data[3], dbg_s2data[3], dbg_sigKdata[3], dbg_sigSdata[3];
dbg_xVec.length = dbg_pVec.length = dbg_s1Vec.length = dbg_s2Vec.length = dbg_sigKerr.length = dbg_sigStar.length = 3;
dbg_xVec.data = dbg_xdata;
dbg_pVec.data = dbg_pdata;
dbg_s1Vec.data = dbg_s1data;
dbg_s2Vec.data = dbg_s2data;
dbg_sigKerr.data = dbg_sigKdata;
dbg_sigStar.data = dbg_sigSdata;
// memcpy(dbg_xdata, ICvalues->data, 3*sizeof(REAL8));
// memcpy(dbg_pdata, ICvalues->data+3, 3*sizeof(REAL8));
// dbg_xdata[0] = 4.1659938868533658e1;
dbg_xdata[0] = 2.2786576810580382e+01;
dbg_xdata[1] = dbg_xdata[2] = 0.0;

// dbg_pdata[0] = -1e-2;
// dbg_pdata[1] = 1.5936732938778625e-1;
// dbg_pdata[2] = -2.6799608562176783e-12;
dbg_pdata[0] = -2.0223298742475641e-04;
dbg_pdata[1] = 2.2970285283699207e-01;
dbg_pdata[2] = 0.0;
// memcpy(dbg_s1data, ICvalues->data+6, 3*sizeof(REAL8));
// memcpy(dbg_s2data, ICvalues->data+9, 3*sizeof(REAL8));
// dbg_s1data[0] = 0.0;
// dbg_s1data[1] = -1.3888888888888892e-01;
// dbg_s1data[2] = 6.2500000000000011e-01;
dbg_s1data[0] = dbg_s1data[1] = 0.0;
dbg_s1data[2] = -8.100000000000000533e-01 ;

// dbg_s2data[0] = 2.7777777777777779e-03;
// dbg_s2data[1] = 0.0;
// dbg_s2data[2] = 2.2222222222222227e-02;
dbg_s2data[0] = dbg_s2data[1] = 0.0;
dbg_s2data[2] = -8.100000000000000533e-01 ;

REAL8 dbg_m1, dbg_m2, dbg_eta, dbg_a, dbg_S_con;
dbg_m1 = 48.3871;
dbg_m2 = 11.6129;
dbg_eta = dbg_m1 * dbg_m2 / (dbg_m1 + dbg_m2) / (dbg_m1 + dbg_m2);
// EOBCalculateSigmaStar(&dbg_sigStar, dbg_m1, dbg_m2, &dbg_s1Vec, &dbg_s2Vec);
// EOBCalculateSigmaKerr(&dbg_sigKerr, &dbg_s1Vec, &dbg_s2Vec);

SpinEOBHSACoeffs dbg_seobCoeffs;
memset(&dbg_seobCoeffs, 0, sizeof(dbg_seobCoeffs));

CalculateSpinEOBHSACoeffs(dbg_m1, dbg_m2, dbg_s1data[2], dbg_s2data[2], &dbg_seobCoeffs);
REAL8 dbg_invcsi;
REAL8 dbg_Ham = EOBSAHamiltonian(dbg_xdata[0], dbg_pdata[0], dbg_pdata[1]*dbg_xdata[0], &dbg_seobCoeffs, &dbg_invcsi);
DEBUG_END;
if(1) {failed = 1; goto QUIT;}
#endif

    /*
    *
    * 
    *       Evolve EOB trajectory with adaptive sampling 
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Evolve EOB trajectory with adaptive sampling.", this_step);
    INT retLenAdaS = 0, retLenInverse = 0;
    REAL8 tendAdaS = 20. / mTScaled;    
    REAL8 tstartAdaS = 0.;
    REAL8 deltaT_min = 8.0e-5;
    REAL8 EPS_REL = 1.0e-10;
    REAL8 EPS_ABS = 1.0e-10;
    if (hparams->inEPS_REL > 0.)
        EPS_REL = hparams->inEPS_REL;
    if (hparams->inEPS_ABS > 0.)
        EPS_ABS = hparams->inEPS_ABS;
    if (MfMin < Mf_ref)
    {
        status = SEOBIntegrateDynamics_SA_inverse(&dynamicsInverse, &retLenInverse, ICvalues_SA, EPS_ABS, EPS_REL,
           deltaT, deltaT_min, tstartAdaS, tendAdaS, core);
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    }
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    if (MfMin < hparams->Mf_max || hparams->t_max > 0.0)
    {
        /**
         * @brief 20240304 by X.L.:
         *          if Mf_max > MfMin,
         *              let integrate stops at omega = pi*Mf_max
         */
        status = SEOBIntegrateDynamics_SA_withFMax(&dynamicsAdaS, &retLenAdaS, 
            ICvalues_SA, EPS_ABS, EPS_REL, 
            deltaT, deltaT_min, tstartAdaS, tendAdaS, core);
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
        if (dynamicsInverse)
            SEOBConcactInverseDynToAdaSDyn_SA(&dynamicsAdaS, dynamicsInverse, &retLenAdaS, retLenInverse);
        status = SEOBComputeExtendedSEOBSAdynamics(&seobdynamicsAdaSHiS, dynamicsAdaS, retLenAdaS, core);    
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
        if (MfMin > Mf_ref)
        {
            status = CutSEOBSAdynamics(&seobdynamicsAdaSHiS, MfMin);
            retLenAdaS = seobdynamicsAdaSHiS->length;
            if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
        }
        PRINT_LOG_INFO(LOG_DEBUG, "Get retLenAdaS = %d", retLenAdaS);
        tVecPmodes = CreateREAL8Vector(retLenAdaS);
        memcpy(tVecPmodes->data, seobdynamicsAdaSHiS->tVec, retLenAdaS*sizeof(REAL8));

        SEOBSACalculateSphHarmListhlmAmpPhase(&listhPlm, modes, nmodes,
                                            seobdynamicsAdaSHiS, NULL,
                                            core, 0);
        
        *tRetVec = tVecPmodes;
        *hlm = listhPlm;
        if (is_dyn_debug)
        {
            PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Concact Dynamics", this_step);
            *dyn_debug = seobdynamicsAdaSHiS;
            // *dyn_debug = CreateSEOBSAdynamics(seobdynamicsAdaSHiS->length);
            // for (int j = 0; j < v4SAdynamicsVariables; j++) 
            // {
            //     // print_debug("Copy truncated dynamics 1 - v4PdynamicsVariables data fields: %d\n", j);
            //     memcpy(&((*dyn_debug)->array->data[j * seobdynamicsAdaSHiS->length]),
            //         &(seobdynamicsAdaSHiS->array->data[j * seobdynamicsAdaSHiS->length]),
            //         seobdynamicsAdaSHiS->length * sizeof(REAL8));
            // }
            // status = SEOBSAJoinDynamics(dyn_debug, seobdynamicsAdaS, seobdynamicsHiS,
            //             indexJoinHiS, indexJoinAttach);
        }

        PRINT_LOG_INFO(LOG_DEBUG, "Finished");
        goto QUIT;
    }
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    status = SEOBIntegrateDynamics_SA(&dynamicsAdaS, &retLenAdaS, 
        ICvalues_SA, EPS_ABS, EPS_REL, 
        deltaT, deltaT_min, tstartAdaS, tendAdaS, core);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    // print_debug("retLenAdaS = %d\n", retLenAdaS);
    if (dynamicsInverse)
        SEOBConcactInverseDynToAdaSDyn_SA(&dynamicsAdaS, dynamicsInverse, &retLenAdaS, retLenInverse);
    // print_debug("retLenAdaS = %d\n", retLenAdaS);

    // {
    //     print_debug("test dynSA\n");
    //     INT jj = 0;
    //     print_debug("%d, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e\n", jj, 
    //         seobdynamicsAdaS->tVec[jj], seobdynamicsAdaS->rVec[jj], seobdynamicsAdaS->phiVec[jj],
    //         seobdynamicsAdaS->prTVec[jj], seobdynamicsAdaS->pphiVec[jj], seobdynamicsAdaS->drVec[jj],
    //         seobdynamicsAdaS->dphiVec[jj], seobdynamicsAdaS->dprTVec[jj], seobdynamicsAdaS->dpphiVec[jj],
    //         seobdynamicsAdaS->HVec[jj]);
    //     print_debug("test dynSA done\n");
    // }

    // PRINT_LOG_INFO(LOG_INFO, "END, AdaS data length = %d", retLenAdaS);
// failed = 1; goto QUIT;
    status = SEOBComputeExtendedSEOBSAdynamics(&seobdynamicsAdaS, dynamicsAdaS, retLenAdaS, core);
    // {
    //     print_debug("test dynSA\n");
    //     INT jj = 0;
    //     print_debug("%d, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e\n", jj, 
    //         seobdynamicsAdaS->tVec[jj], seobdynamicsAdaS->rVec[jj], seobdynamicsAdaS->phiVec[jj],
    //         seobdynamicsAdaS->prTVec[jj], seobdynamicsAdaS->pphiVec[jj], seobdynamicsAdaS->drVec[jj],
    //         seobdynamicsAdaS->dphiVec[jj], seobdynamicsAdaS->dprTVec[jj], seobdynamicsAdaS->dpphiVec[jj],
    //         seobdynamicsAdaS->HVec[jj]);
    //     print_debug("test dynSA done\n");
    // }
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    if (MfMin > Mf_ref)
    {
        status = CutSEOBSAdynamics(&seobdynamicsAdaS, MfMin);
        retLenAdaS = seobdynamicsAdaS->length;
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    }
    m1rVec = CreateREAL8Vector(retLenAdaS);
    for (i = 0; i < retLenAdaS; i++) 
    {
        m1rVec->data[i] = -1* seobdynamicsAdaS->rVec[i];
    }
    UINT index_6M = FindClosestIndex(m1rVec, -6.0);
    REAL8 time_6M = seobdynamicsAdaS->tVec[index_6M];
    tStepBack = GET_MAX(tStepBack, seobdynamicsAdaS->tVec[retLenAdaS-1] - time_6M + 10*deltaT);
#if 0
    FILE *out = fopen( "debug_dynamics_AdaS.dat","w");
    for (int ii=0; ii<retLenAdaS; ii++)
    {
        // t, x, y, z, px, py, pz
    fprintf(out, "%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n",
      seobdynamicsAdaS->tVec[ii],
      seobdynamicsAdaS->posVecx[ii],
      seobdynamicsAdaS->posVecy[ii],
      seobdynamicsAdaS->posVecz[ii],
      seobdynamicsAdaS->momVecx[ii],
      seobdynamicsAdaS->momVecy[ii],
      seobdynamicsAdaS->momVecz[ii],
      seobdynamicsAdaS->s1Vecx[ii],
      seobdynamicsAdaS->s1Vecy[ii],
      seobdynamicsAdaS->s1Vecz[ii],
      seobdynamicsAdaS->s2Vecx[ii],
      seobdynamicsAdaS->s2Vecy[ii],
      seobdynamicsAdaS->s2Vecz[ii]);
    }
    fclose(out);
#endif

    /*
    *
    * 
    *       Step back and evolve EOB trajectory at high sampling rate (HiS)
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Step back and evolve EOB trajectory at high sampling rate (HiS).", this_step);
    REAL8 deltaTHiS = 1. / 50; /* Fixed at 1/50M */
    REAL8 tstartHiSTarget = seobdynamicsAdaS->tVec[retLenAdaS - 1] - tStepBack;
    INT4 indexstartHiS = retLenAdaS - 1; /* index for the AdaS dynamics */
    while ((indexstartHiS > 0) && (seobdynamicsAdaS->tVec[indexstartHiS] > tstartHiSTarget))
        indexstartHiS--;
    REAL8 tstartHiS = seobdynamicsAdaS->tVec[indexstartHiS];
    PRINT_LOG_INFO(LOG_DEBUG, "Interpolate Dynamics At Time %f (%f)/%f", tstartHiS, tstartHiSTarget, seobdynamicsAdaS->tVec[retLenAdaS - 1]);
    status = SEOBInterpolateSADynamicsAtTime(&seobvalues_tstartHiS, tstartHiS, seobdynamicsAdaS);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    ICvaluesHiS_SA = CreateREAL8Vector(4);
    if (!ICvaluesHiS_SA) {failed = 1; goto QUIT;}
    memcpy(ICvaluesHiS_SA->data, &(seobvalues_tstartHiS->data[1]), 4 * sizeof(REAL8));
    PRINT_LOG_INFO(LOG_DEBUG, "initConds = (%.16e, %.16e, %.16e, %.16e)\n",
        ICvaluesHiS_SA->data[0], ICvaluesHiS_SA->data[1], ICvaluesHiS_SA->data[2], ICvaluesHiS_SA->data[3]);
    /* Integrate again the dynamics with a constant high sampling rate */
    core->prev_dr = 0.; // This is used to check whether dr/dt is increasing
                            // in the stopping condition
    core->termination_reason =
        -999; // This is going to store the termination condition for the
                // high-sampling integration
    INT is_re_evolved = FALSE; 
    INT retLenHiS = 0;
    REAL8 tendHiS = tstartHiS + tStepBack; /* Seems it should be ignored anyway because of
                                    integrator->stopontestonly */
    REAL8 rISCO;
HISR:
    rISCO = get_h_rISCO();
    PRINT_LOG_INFO(LOG_DEBUG, "tStepBack = %g", tStepBack);
    status = SEOBIntegrateDynamics_SA(&dynamicsHiS, &retLenHiS, 
        ICvaluesHiS_SA, EPS_ABS, EPS_REL, 
        deltaTHiS, 0., tstartHiS, tendHiS, core);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    // PRINT_LOG_INFO(LOG_INFO, "END, HiSR data length = %d", retLenHiS);
    /* Compute derived quantities for the high-sampling dynamics */
    status = SEOBComputeExtendedSEOBSAdynamics(&seobdynamicsHiS, dynamicsHiS, retLenHiS, core);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
// print_debug("sizeof dynamicsAdaS = (%d, %d)\n", dynamicsAdaS->dimLength->data[0], dynamicsAdaS->dimLength->data[1]);
// print_debug("sizeof dynamicsHiS = (%d, %d)\n", dynamicsAdaS->dimLength->data[0], dynamicsAdaS->dimLength->data[1]);

    /* Find time of peak of omega for the High-sampling dynamics */
    INT foundPeakOmega = 0;
    REAL8 tPeakOmega = 0.;
    SEOBSALocateTimePeakOmega(&tPeakOmega, &foundPeakOmega,
                            seobdynamicsHiS, retLenHiS, core);
    // {
    //     int j;
    //     REAL8 *vector;
    //     for (i=0; i<seobdynamicsHiS->length; i++)
    //     {
    //         print_out("%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n", 
    //             seobdynamicsHiS->tVec[i], 
    //             seobdynamicsHiS->rVec[i], seobdynamicsHiS->phiVec[i],
    //             seobdynamicsHiS->prTVec[i], seobdynamicsHiS->pphiVec[i],
    //             seobdynamicsHiS->drVec[i], seobdynamicsHiS->dphiVec[i],
    //             seobdynamicsHiS->dprTVec[i], seobdynamicsHiS->dpphiVec[i],
    //             seobdynamicsHiS->HVec[i]);
    //     }
    // }
    // failed = 1; goto QUIT;

    PRINT_LOG_INFO(LOG_DEBUG, "Locate timePeak of omega = %f", tPeakOmega);
    /*
    *
    * 
    *       Get final J/L/spins from HiS dynamics at peak of Omega, compute constant angles EulerI2J
    * 
    * 
    */

    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Get final J/L/spins from HiS dynamics at peak of Omega.", this_step);
    status = SEOBInterpolateSADynamicsAtTime(&seobvalues_tPeakOmega, tPeakOmega, seobdynamicsHiS);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}

    /* Interpolate the dynamics to r=10M. This is used as a
    feducial point to evaluate the final mass and spin fits */
    REAL8Vector timeVec;
    timeVec.length = retLenAdaS;
    timeVec.data = seobdynamicsAdaS->tVec;

    // FindClosestIndex requires an *increasing* function of time, so we use -r
    // instead of r
    UINT index_10M = FindClosestIndex(m1rVec, -10.0);
    REAL8 time_10M = timeVec.data[index_10M];
    status = SEOBInterpolateSADynamicsAtTime(&seobvalues_test, time_10M, seobdynamicsAdaS);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}

    /* Compute the timeshift to get the attachment point */
    // SEOBLFrameVectors(&chi1L_tPeakOmega, &chi2L_tPeakOmega, 
    //     seobvalues_tPeakOmega, m1, m2, core->hParams->flagZframe);
    REAL8 Deltat22 = XLALSimIMREOBGetNRSpinPeakDeltaTv4(
        2, 2, m1, m2, core->chi1, core->chi2, core->hParams);

    /* Determine the time of attachment */
    REAL8 tAttach = tPeakOmega - Deltat22;
    PRINT_LOG_INFO(LOG_INFO, "The location of attachment is %f, Deltat22 = %f", tAttach, Deltat22);
    if (core->hParams->zero_dyncoaphase)
    {
        REAL8 de_phi;
        status = SetZeroPhaseAtTimeSA(seobdynamicsHiS, tAttach, &de_phi);
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
        OrbitalPhaseReduceSA(seobdynamicsAdaS, de_phi);
    }

    //Compute NQC window factors
    if (ecc != 0.0)
    {
        status = SEOBSACalculateNQCWindowFactorsFromDyn(seobdynamicsAdaS, tAttach, time_6M, tstartHiS, 6., 0, core);
        // core->wWind = 1.;
    }
    
    if (is_re_evolved)
    {
        PRINT_LOG_INFO(LOG_INFO, "re-calculating window parameters; rISCO = %f, rmin = %f", rISCO, seobdynamicsHiS->rVec[retLenHiS-1]);
        status = SEOBSACalculateNQCWindowFactorsFromDyn(seobdynamicsAdaS, tAttach, time_6M, tstartHiS, rISCO, 1, core);
        // status = SEOBCalculateNQCWindowFactorsFromDyn(seobdynamicsHiS, tAttach, time_6M, tstartHiS, rISCO, 1, core);
    }

#if 0
    /* Compute final J from dynamics quantities */
    SEOBJfromDynamics(&Jfinal, seobvalues_tPeakOmega, core);
    /*Compute the L-hat vector. Note that it has unit norm */
    SEOBLhatfromDynamics(&Lhatfinal, seobvalues_tPeakOmega, core);
    REAL8 Jmag = sqrt(inner_product3d(Jfinal->data, Jfinal->data));
    /* Cosine of the angle between L-hat and J. Needed to determine
    * the correct sign of the final spin
    */
    REAL8 cos_angle = inner_product3d(Jfinal->data, Lhatfinal->data) / Jmag;
    /* Compute final-J-frame unit vectors e1J, e2J, e3J=Jfinalhat */
    /* Convention: if (ex, ey, ez) is the initial I-frame, e1J chosen such that ex
    * is in the plane (e1J, e3J) and ex.e1J>0 */
    REAL8Vector e1J, e2J, e3J;
    e1J.length = e2J.length = e3J.length = 3;
    REAL8 e1Jdata[3] = {0.};
    REAL8 e2Jdata[3] = {0.};
    REAL8 e3Jdata[3] = {0.};
    e1J.data = e1Jdata;
    e2J.data = e2Jdata;
    e3J.data = e3Jdata;
    SEOBBuildJframeVectors(&e1J, &e2J, &e3J, Jfinal);

    /* Compute Euler angles from initial I-frame to final-J-frame */
    /* Note: if spins are aligned, the function SEOBEulerI2JFromJframeVectors */
    /* becomes ill-defined - just keep these Euler angles to zero then */
    REAL8 alphaI2J = 0., betaI2J = 0., gammaI2J = 0.;
    if (!core->alignedSpins) 
        SEOBEulerI2JFromJframeVectors(&alphaI2J, &betaI2J, &gammaI2J, &e1J, &e2J, &e3J);
#endif

    /*
    *
    * 
    *       Compute P-frame amp/phase for all modes on HiS and compute NQC Coeffs
    * 
    * 
    */

    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Compute P-frame amp/phase for all modes on HiS and compute NQC Coeffs.", this_step);
    // REAL8 wWind_old = core->wWind;
    status = SEOBSACalculateSphHarmListNQCCoefficientsV4(
            &nqcCoeffsList, modes, nmodes, tPeakOmega, seobdynamicsHiS,
            core);

    if ( status != CEV_SUCCESS) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "NQC computation failed.");
        failed = 1;
        goto QUIT;
    }
    if (!is_re_evolved)
    {
        core->wWind = 1.;
    }
    // Check should we need to re-evolve HiSR EOB
// #if 1
    if (ecc != 0.0 && rISCO > 3 && CheckStopConditionSA(seobdynamicsHiS, core, nqcCoeffsList, tAttach - seobdynamicsHiS->tVec[0]) != CEV_SUCCESS)
    {
#if 1
        // re-evolve EOB
        SET_RISCO(rISCO-1.);
        // PRINT_LOG_INFO(LOG_INFO, "Re-evolve HiSR EOB...");
        // print_debug("Re-evolve EOB...");
        is_re_evolved = TRUE;
        STRUCTFREE(seobdynamicsHiS, SEOBSAdynamics);
        // STRUCTFREE(Jfinal, REAL8Vector);
        STRUCTFREE(Lhatfinal, REAL8Vector);
        STRUCTFREE(seobvalues_tPeakOmega, REAL8Vector);
        STRUCTFREE(seobvalues_test, REAL8Vector);
        STRUCTFREE(chi1L_tPeakOmega, REAL8Vector);
        STRUCTFREE(chi2L_tPeakOmega, REAL8Vector);
        // STRUCTFREE(dynamicsHiS, REAL8Array);
        // STRUCTFREE(seobvalues_tstartHiS, REAL8Vector);
        STRUCTFREE(nqcCoeffsList, SphHarmListEOBNonQCCoeffs);
        // STRUCTFREE(ICvaluesHiS, REAL8Vector);
        goto HISR;
#endif
    }

    /*
    *
    *       Compute P-frame amp/phase for all modes on HiS, now including NQC
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Compute P-frame amp/phase for all modes on HiS, now including NQC.", this_step);
    // /* We now include the NQC when generating modes */
    // flagNQC = 1;

    /* Compute amplitude and phase of the P-frame modes hPlm on high sampling,
    * with NQC */
    SEOBSACalculateSphHarmListhlmAmpPhase(&listhPlm_HiS, modes, nmodes,
                                        seobdynamicsHiS, nqcCoeffsList,
                                        core, 1);
    if (is_noringdown)
    {
        // Attach AdaS and HiS waveform and quit
        PRINT_LOG_INFO(LOG_INFO, "Attach AdaS and HiS waveform and quit.");
        SEOBSACalculateSphHarmListhlmAmpPhase(&listhPlm_AdaS, modes, nmodes,
                                        seobdynamicsAdaS, nqcCoeffsList,
                                        core, 1);
        REAL8 finalSpinSign = 1.;
        REAL8 tmp = core->eta*seobvalues_tPeakOmega->data[4] + 
            (core->m1/mTotal)*(core->m1/mTotal)*core->chi1 + 
            (core->m2/mTotal)*(core->m2/mTotal)*core->chi2;
        // print_debug("seobvalues_tPeakOmega = (%.16e, %.16e, %.16e, %.16e)\n", 
        //     seobvalues_tPeakOmega->data[0], seobvalues_tPeakOmega->data[1], seobvalues_tPeakOmega->data[2], seobvalues_tPeakOmega->data[3]);
        if (tmp < 0)
        {
            finalSpinSign = -1.;
        }
        PRINT_LOG_INFO(LOG_INFO, "Attach AdaS and HiS time vector.");
        status = SEOBSAAttachAdaSandHiSRTimeVector(&tVecPmodes, 
            seobdynamicsAdaS, seobdynamicsHiS, 
            tstartHiS);

        if (status != CEV_SUCCESS)
        {failed = 1; goto QUIT;}
        for (i=0; i<tVecPmodes->length; i++)
            tVecPmodes->data[i] += finalSpinSign*tAttach;

        PRINT_LOG_INFO(LOG_INFO, "Attach AdaS and HiS waveform at %.16e\n", tAttach);
        status = SEOBJoinSphHarmListhlm(&listhPlm, listhPlm_AdaS, listhPlm_HiS, modes,
                            nmodes, indexstartHiS);
        if (status != CEV_SUCCESS)
        {failed = 1; goto QUIT;}
        *tRetVec = tVecPmodes;
        *hlm = listhPlm;
        goto QUIT;
    }
    /*
    *
    * 
    *       Attach RD to the P-frame waveform
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Attach RD to the P-frame waveform.", this_step);
    REAL8 finalMass = 0., finalSpin = 0.;
    status = SEOBSAGetFinalSpinMass(&finalMass, &finalSpin, core);
    if (status != CEV_SUCCESS) 
    {failed = 1; goto QUIT;}
    PRINT_LOG_INFO(LOG_DEBUG, "finalMass = %.16e, finalSpin = %.16e\n", finalMass, finalSpin);
    /* The function above returns only the magnitude of the spin.
    *  We pick the direction based on whether Lhat \cdot J is positive
    *  or negative */
    REAL8 cos_angle = core->eta*seobvalues_tPeakOmega->data[4] + 
        (core->m1/mTotal)*(core->m1/mTotal)*core->chi1 + 
        (core->m2/mTotal)*(core->m2/mTotal)*core->chi2;
    // print_debug("seobvalues_tPeakOmega = (%.16e, %.16e, %.16e, %.16e)\n", 
    //     seobvalues_tPeakOmega->data[0], seobvalues_tPeakOmega->data[1], seobvalues_tPeakOmega->data[2], seobvalues_tPeakOmega->data[3]);
    if (cos_angle < 0)
    {
        finalSpin *= -1;
    }
    /* finalSpin interpolation is available only between -0.9996 and 0.9996 */
    /* Set finalSpin to +/- 0.9996 if it is out of this range */
    if (finalSpin < -0.9996)
        finalSpin = -0.9996;
    if (finalSpin > 0.9996)
        finalSpin = 0.9996;

    // PRINT_LOG_INFO(LOG_INFO, "final mass = %e, final spin = %e, cos_angle = %e", finalMass, finalSpin, cos_angle);

    /* Estimate leading QNM , to determine how long the ringdown
    * patch will be */
    /* NOTE: XLALSimIMREOBGenerateQNMFreqV2Prec returns the complex frequency in
    * physical units... */
    COMPLEX16Vector sigmaQNM220estimatephysicalVec;
    COMPLEX16 sigmaQNM220estimatephysical = 0.;
    sigmaQNM220estimatephysicalVec.length = 1;
    sigmaQNM220estimatephysicalVec.data = &sigmaQNM220estimatephysical;
    status = XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec(&sigmaQNM220estimatephysicalVec, m1,
                                                m2, finalMass, finalSpin, 2,
                                                2, 1);
    if (status == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec for mode (l,m) = (2,2).");
        failed = 1;
        goto QUIT;
    }
    COMPLEX16 sigmaQNM220estimate = mTScaled * sigmaQNM220estimatephysical;
    PRINT_LOG_INFO(LOG_DEBUG, "sigmaQNM220 = %.16e + i%.16e\n", creal(sigmaQNM220estimate), cimag(sigmaQNM220estimate));

    /* Length of RD patch, 40 e-folds of decay of the estimated QNM220 */
    UINT retLenRDPatch =
    (UINT)ceil(EFOLDS / (cimag(sigmaQNM220estimate) * deltaTHiS));


    /* Attach RD to the P-frame modes */
    /* Vector holding the values of the 0-th overtone QNM complex frequencies for
    * the modes (l,m) */
    // NOTE: the QNM complex frequencies are computed inside
    // SEOBAttachRDToSphHarmListhPlm
    status = SEOBSAAttachRDToSphHarmListhPlm(
        &listhPlm_HiSRDpatch, &sigmaQNMlm0, modes, nmodes, finalMass, finalSpin,
        listhPlm_HiS, deltaTHiS, retLenHiS, retLenRDPatch, tAttach,
        seobvalues_tPeakOmega, seobdynamicsHiS, core);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}

    /*
    *
    * 
    *       Build the joined dynamics AdaS+HiS up to attachment, joined P-modes AdaS+HiS+RDpatch 
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Build the joined dynamics AdaS+HiS up to attachment, joined P-modes AdaS+HiS+RDpatch", this_step);

    /* Compute amplitude and phase of the P-frame modes hPlm on adaptive sampling,
    * with NQC */
    // flagNQC = 1;

    SEOBSACalculateSphHarmListhlmAmpPhase(&listhPlm_AdaS, modes, nmodes,
                                    seobdynamicsAdaS, nqcCoeffsList,
                                    core, 1);

    /* Vector of times for the P-modes, that will be used for interpolation:
    * joining AdaS and HiS+RDpatch */
    UINT retLenPmodes = 0;
    /* First junction at indexAdaSHiS, tAdaSHiS */
    UINT indexJoinHiS = 0;
    REAL8 tJoinHiS = 0.;
    /* Second junction at seobdynamicsAdaS, tJoinAttach */
    UINT indexJoinAttach = 0;
    REAL8 tJoinAttach = 0.;
    /* Construct the joined vector of times (AdaS+HiS+RDpatch) and keep the
    * jonction indices and times */
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Concact Time Vector", this_step);
    // {
    //     INT jj = 0;
    //     print_debug("%d, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e\n", jj, 
    //         seobdynamicsAdaS->tVec[jj], seobdynamicsAdaS->rVec[jj], seobdynamicsAdaS->phiVec[jj],
    //         seobdynamicsAdaS->prTVec[jj], seobdynamicsAdaS->pphiVec[jj], seobdynamicsAdaS->drVec[jj],
    //         seobdynamicsAdaS->dphiVec[jj], seobdynamicsAdaS->dprTVec[jj], seobdynamicsAdaS->dpphiVec[jj],
    //         seobdynamicsAdaS->HVec[jj]);
    // }
    status = SEOBSAJoinTimeVector(&tVecPmodes, &retLenPmodes, &tJoinHiS, &indexJoinHiS,
                    &tJoinAttach, &indexJoinAttach, retLenRDPatch, deltaTHiS,
                    tstartHiS, tAttach, seobdynamicsAdaS, seobdynamicsHiS);
 
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}

    /* Copy dynamics from AdaS<HiS and HiS<tAttach to form joined dynamics, ending
    * at the last time sample <tAttach */
    // NOTE: we cut the dynamics at tAttach, as we will extend the Euler
    // angles for t>=tAttach -- but we could also choose to finish at tPeakOmega
    // which is used for final-J and for the final mass/spin fit

    // status = SEOBJoinDynamics(&seobdynamicsAdaSHiS, seobdynamicsAdaS, seobdynamicsHiS,
    //             indexJoinHiS, indexJoinAttach);
    // if (status != CEV_SUCCESS)
    // {failed = 1; goto QUIT;}

    /* Copy waveform modes from AdaS and HiS+RDpatch - adjusting 2pi-phase shift
    * at the junction point AdaS/HiS */
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Concact Spherical Harmonic list hlm", this_step);

    status = SEOBJoinSphHarmListhlm(&listhPlm, listhPlm_AdaS, listhPlm_HiSRDpatch, modes,
                        nmodes, indexstartHiS);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}
    *tRetVec = tVecPmodes;
    *hlm = listhPlm;
    if (is_dyn_debug)
    {
        PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Concact Dynamics", this_step);
        status = SEOBSAJoinDynamics(dyn_debug, seobdynamicsAdaS, seobdynamicsHiS,
                    indexJoinHiS, indexJoinAttach);
    }
    // PRINT_LOG_INFO(LOG_DEBUG, "len_timeP = %d, len_Plm = %d\n", tVecPmodes->length, listhPlm->campphase->camp_imag->length);

    // status = SEOBRotateInterpolatehJlmReImFromSphHarmListhPlmAmpPhase(
    //     &hJlm, modes, nmodes, modes_lmax, deltaT, retLenTS, tVecPmodes,
    //     listhPlm, alphaJ2P, betaJ2P, gammaJ2P);
    // if ( status == CEV_FAILURE) 
    // {
    //     PRINT_LOG_INFO(LOG_CRITICAL, "failure in SEOBRotateInterpolatehJlmReImFromSphHarmListhPlmAmpPhase for mode (l,m) = (2,2).");
    //     failed = 1;
    //     goto QUIT;
    // }

#if 0

    /* Get the time of the frame-invariant amplitude peak */
    REAL8 tPeak = 0;
    UINT indexPeak = 0;
    // NOTE: peak amplitude using l=2 only: h22 required, and h21 used if present
    status = SEOBAmplitudePeakFromAmp22Amp21(&tPeak, &indexPeak, listhPlm, modes, nmodes,
                                tVecPmodes);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}


    /*
    *
    * 
    *       Compute Euler angles J2P from AdaS and HiS dynamics up to attachment
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Compute Euler angles J2P from AdaS and HiS dynamics up to attachment", this_step);
    /* Compute Euler angles J2P from the dynamics before attachment point */
    /* If SpinsAlmostAligned, all Euler angles are set to 0 */
    status = SEOBEulerJ2PFromDynamics(&alphaJ2P, &betaJ2P, &gammaJ2P, 
                &e1J, &e2J, &e3J,
                retLenPmodes, indexJoinAttach,
                seobdynamicsAdaSHiS, core);
    if (status == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "SEOBEulerJ2PFromDynamics failed.");
        failed = 1;
        goto QUIT;
    }

    /*
    *
    * 
    *       Compute Euler angles J2P extension after attachment
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Compute Euler angles J2P extension after attachment", this_step);

    /* Compute Euler angles J2P according to the prescription flagEulerextension
    * after attachment point */
    /* If SpinsAlmostAligned, all Euler angles are set to 0 */
    /* NOTE: Regardless of the mode content of hPlm, the frame extension at the
    * moment is based on sigmaQNM22, sigmaQNM21 */
    COMPLEX16 sigmaQNM220 = 0., sigmaQNM210 = 0.;
    COMPLEX16Vector sigmaQNM220physicalVec, sigmaQNM210physicalVec;
    sigmaQNM220physicalVec.length = 1;
    sigmaQNM210physicalVec.length = 1;
    COMPLEX16 sigmaQNM220physicalval = 0.;
    COMPLEX16 sigmaQNM210physicalval = 0.;
    sigmaQNM220physicalVec.data = &sigmaQNM220physicalval;
    sigmaQNM210physicalVec.data = &sigmaQNM210physicalval;
    /* NOTE: XLALSimIMREOBGenerateQNMFreqV2Prec returns the complex frequency in
    * physical units... */
    status = XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec(&sigmaQNM220physicalVec, m1,
                                                m2, finalMass, finalSpin, 2,
                                                2, 1);
    if (status == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec for mode (l,m) = (2,2).");
        failed = 1;
        goto QUIT;
    }
    status = XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec(&sigmaQNM210physicalVec, m1,
                                                m2, finalMass, finalSpin, 2,
                                                1, 1);
    if (status == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec for mode (l,m) = (2,1).");
        failed = 1;
        goto QUIT;
    }
    sigmaQNM220 = mTScaled * sigmaQNM220physicalVec.data[0];
    sigmaQNM210 = mTScaled * sigmaQNM210physicalVec.data[0];
    INT flip = 1;
    if (cos_angle < 0)
        flip = -1;
    // flagEulerextension = 0
    status = SEOBEulerJ2PPostMergerExtension(
        alphaJ2P, betaJ2P, gammaJ2P, sigmaQNM220, sigmaQNM210, tVecPmodes,
        retLenPmodes, indexJoinAttach, core, flip);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}

    /*
    *
    * 
    *       Compute modes hJlm on the output time series by rotating and interpolating the modes hPlm
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Compute modes hJlm on the output time series by rotating and interpolating the modes hPlm", this_step);

    /* Determine the length of the fixed-sampling output time series */
    UINT retLenTS = floor(((tVecPmodes)->data[retLenPmodes - 1] - (tVecPmodes)->data[0]) / deltaT);

    /* Rotate waveform from P-frame to J-frame */
    // flagSymmetrizehPlminusm = 1
    status = SEOBRotateInterpolatehJlmReImFromSphHarmListhPlmAmpPhase(
        &hJlm, modes, nmodes, modes_lmax, deltaT, retLenTS, tVecPmodes,
        listhPlm, alphaJ2P, betaJ2P, gammaJ2P);
    if ( status == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in SEOBRotateInterpolatehJlmReImFromSphHarmListhPlmAmpPhase for mode (l,m) = (2,2).");
        failed = 1;
        goto QUIT;
    }
    /*
    *
    * 
    *       Rotate waveform from J-frame to the output I-frame on timeseries-sampling (constant Wigner coeffs)
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Rotate waveform from J-frame to the output I-frame on timeseries-sampling (constant Wigner coeffs)", this_step);

    /* Rotate waveform from J-frame to I-frame */
    status = SEOBRotatehIlmFromhJlm(&hIlm, hJlm, modes_lmax, alphaI2J, betaI2J, gammaI2J, deltaT);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}
    hIlm->tAttach = tAttach;

    /*
    *
    * 
    *       Compute h20 from I-frame waveform on timeseries sampling
    * 
    * 
    */
    // this_step++;
    // PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Compute h20 from I-frame waveform on timeseries sampling", this_step);

    // status = SEOBCalculateh20(hIlm);
    // {failed = 1; goto QUIT;}

    /*
    *
    * 
    *       Compute hplus, hcross from I-frame waveform on timeseries sampling
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Compute hplus, hcross from I-frame waveform on timeseries sampling", this_step);

    /* GPS time for output time series and modes */
    // LIGOTimeGPS tGPS = LIGOTIMEGPSZERO;
    // XLALGPSAdd(&tGPS, -mTScaled *
    //                     tPeak); /* tPeak converted back to dimensionfull time */
    REAL8 tGPS = -mTScaled * tPeak;
    /* Create output timeseries for hplus, hcross */
    /* Use the dimensionfull INdeltaT (s) as time step */
    // hplusTS = CreateREAL8TimeSeries("H_PLUS", &tGPS, 0.0, INdeltaT,
    //                                 &lalStrainUnit, retLenTS);
    // hcrossTS = CreateREAL8TimeSeries("H_CROSS", &tGPS, 0.0, INdeltaT,
    //                                 &lalStrainUnit, retLenTS);
    hplusTS = CreateREAL8TimeSeries(tGPS, INdeltaT, retLenTS);
    hcrossTS = CreateREAL8TimeSeries(tGPS, INdeltaT, retLenTS);
    /* Compute hplus, hcross from hIlm */
    // NOTE: azimuthal angle of the observer entering the -2Ylm is pi/2-phi
    // according to LAL conventions
    status = SEOBComputehplushcrossFromhIlm(hplusTS, hcrossTS, modes_lmax, hIlm, amp0,
        inc, phi0, is_only22);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}

    /*
    *
    * 
    *       Output and cleanup
    * 
    * 
    */

    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Output and cleanup", this_step);
    /* Point the output pointers to the relevant time series and return */
    (*hPlusOut) = hplusTS;
    (*hCrossOut) = hcrossTS;
    seobdynamicsAdaSHiS->th22Peak = tAttach;
    all->dyn = seobdynamicsAdaSHiS;
    all->hLM = hIlm;
#if 0
    /* Output vector gathering quantities related to merger (similar to previous
    * AttachParams) */
    /* Format: tPeakOmega tAttach tPeak Jfinalx Jfinaly Jfinalz finalMassfit
    * finalSpinfit termination_reason [sigmaQNMlm0Re sigmaQNMlm0Im for lm in
    * modes] */
    /* NOTE: the size of this output vector depends on the number of modes, due to
    * the sigmaQNM */
    if (!((*mergerParams) = CreateREAL8Vector(9 + 2 * nmodes))) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failed to allocate REAL8Vector mergerParams.");
        failed = 1;
        goto QUIT;
    }
    (*mergerParams)->data[0] = tPeakOmega;
    (*mergerParams)->data[1] = tAttach;
    (*mergerParams)->data[2] = tPeak;
    (*mergerParams)->data[3] = Jfinal->data[0];
    (*mergerParams)->data[4] = Jfinal->data[1];
    (*mergerParams)->data[5] = Jfinal->data[2];
    (*mergerParams)->data[6] = finalMass;
    (*mergerParams)->data[7] = finalSpin;
    (*mergerParams)->data[8] = seobParams.termination_reason;
    for (UINT4 nmode = 0; nmode < nmodes; nmode++) {
    (*mergerParams)->data[9 + 2 * nmode] = creal(sigmaQNMlm0->data[nmode]);
    (*mergerParams)->data[9 + 2 * nmode + 1] = cimag(sigmaQNMlm0->data[nmode]);
    }


    /* Additional outputs */

    /* Dynamics */
    // NOTE: casting to REAL8Vector due to the SWIG wrapping
    *seobdynamicsAdaSVector =
    XLALCreateREAL8Vector(v4PdynamicsVariables * seobdynamicsAdaS->length);
    *seobdynamicsHiSVector =
    XLALCreateREAL8Vector(v4PdynamicsVariables * seobdynamicsHiS->length);
    *seobdynamicsAdaSHiSVector =
    XLALCreateREAL8Vector(v4PdynamicsVariables * seobdynamicsAdaSHiS->length);
    memcpy((*seobdynamicsAdaSVector)->data, seobdynamicsAdaS->array->data,
    (v4PdynamicsVariables * seobdynamicsAdaS->length) * sizeof(REAL8));
    memcpy((*seobdynamicsHiSVector)->data, seobdynamicsHiS->array->data,
    (v4PdynamicsVariables * seobdynamicsHiS->length) * sizeof(REAL8));
    memcpy((*seobdynamicsAdaSHiSVector)->data, seobdynamicsAdaSHiS->array->data,
    (v4PdynamicsVariables * seobdynamicsAdaSHiS->length) * sizeof(REAL8));

    /* Modes in the P-frame */
    // NOTE: casting to REAL8Vector due to the SWIG wrapping
    // NOTE: in the output, real amplitude instead of complex envelope
    /* (2,2) */
    *hP22_amp = XLALCreateREAL8Vector(retLenPmodes);
    *hP22_phase = XLALCreateREAL8Vector(retLenPmodes);
    SphHarmListCAmpPhaseSequence *hP22 =
    SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 2, 2);
    if (hP22 == NULL) {
    memset((*hP22_amp)->data, 0, retLenPmodes * sizeof(REAL8));
    memset((*hP22_phase)->data, 0, retLenPmodes * sizeof(REAL8));
    } else {
    memcpy((*hP22_amp)->data, hP22->campphase->camp_real->data,
        retLenPmodes * sizeof(REAL8));
    memcpy((*hP22_phase)->data, hP22->campphase->phase->data,
        retLenPmodes * sizeof(REAL8));
    }
    /* (2,1) */
    *hP21_amp = XLALCreateREAL8Vector(retLenPmodes);
    *hP21_phase = XLALCreateREAL8Vector(retLenPmodes);
    SphHarmListCAmpPhaseSequence *hP21 =
    SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 2, 1);
    if (hP21 == NULL) {
    memset((*hP21_amp)->data, 0, retLenPmodes * sizeof(REAL8));
    memset((*hP21_phase)->data, 0, retLenPmodes * sizeof(REAL8));
    } else {
    memcpy((*hP21_amp)->data, hP21->campphase->camp_real->data,
        retLenPmodes * sizeof(REAL8));
    memcpy((*hP21_phase)->data, hP21->campphase->phase->data,
        retLenPmodes * sizeof(REAL8));
    }
    /* (3,3) */
    *hP33_amp = XLALCreateREAL8Vector(retLenPmodes);
    *hP33_phase = XLALCreateREAL8Vector(retLenPmodes);
    SphHarmListCAmpPhaseSequence *hP33 =
    SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 3, 3);
    if (hP33 == NULL) {
    memset((*hP33_amp)->data, 0, retLenPmodes * sizeof(REAL8));
    memset((*hP33_phase)->data, 0, retLenPmodes * sizeof(REAL8));
    } else {
    memcpy((*hP33_amp)->data, hP33->campphase->camp_real->data,
        retLenPmodes * sizeof(REAL8));
    memcpy((*hP33_phase)->data, hP33->campphase->phase->data,
        retLenPmodes * sizeof(REAL8));
    }
    /* (4,4) */
    *hP44_amp = XLALCreateREAL8Vector(retLenPmodes);
    *hP44_phase = XLALCreateREAL8Vector(retLenPmodes);
    SphHarmListCAmpPhaseSequence *hP44 =
    SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 4, 4);
    if (hP44 == NULL) {
    memset((*hP44_amp)->data, 0, retLenPmodes * sizeof(REAL8));
    memset((*hP44_phase)->data, 0, retLenPmodes * sizeof(REAL8));
    } else {
    memcpy((*hP44_amp)->data, hP44->campphase->camp_real->data,
        retLenPmodes * sizeof(REAL8));
    memcpy((*hP44_phase)->data, hP44->campphase->phase->data,
        retLenPmodes * sizeof(REAL8));
    }
    /* (5,5) */
    *hP55_amp = XLALCreateREAL8Vector(retLenPmodes);
    *hP55_phase = XLALCreateREAL8Vector(retLenPmodes);
    SphHarmListCAmpPhaseSequence *hP55 =
    SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 5, 5);
    if (hP55 == NULL) {
    memset((*hP55_amp)->data, 0, retLenPmodes * sizeof(REAL8));
    memset((*hP55_phase)->data, 0, retLenPmodes * sizeof(REAL8));
    } else {
    memcpy((*hP55_amp)->data, hP55->campphase->camp_real->data,
        retLenPmodes * sizeof(REAL8));
    memcpy((*hP55_phase)->data, hP55->campphase->phase->data,
        retLenPmodes * sizeof(REAL8));
    }
#endif
#endif
QUIT:
    PRINT_LOG_INFO(LOG_INFO, "Program END\n\n");
    STRUCTFREE(core, SpinEOBParams);

    // STRUCTFREE(dynamicsAdaS, REAL8Array);
    STRUCTFREE(dynamicsInverse, REAL8Array);
    STRUCTFREE(ICvalues, REAL8Vector);
    STRUCTFREE(ICvalues_SA, REAL8Vector);
    STRUCTFREE(seobdynamicsAdaS, SEOBSAdynamics);

    // STRUCTFREE(dynamicsHiS, REAL8Array);
    STRUCTFREE(ICvaluesHiS, REAL8Vector);
    STRUCTFREE(ICvaluesHiS_SA, REAL8Vector);
    STRUCTFREE(seobvalues_tstartHiS, REAL8Vector);
    STRUCTFREE(seobdynamicsHiS, SEOBSAdynamics);

    STRUCTFREE(seobvalues_tPeakOmega, REAL8Vector);
    STRUCTFREE(seobvalues_test, REAL8Vector);
    STRUCTFREE(m1rVec,REAL8Vector);
    STRUCTFREE(Jfinal, REAL8Vector);
    STRUCTFREE(Lhatfinal, REAL8Vector);
    STRUCTFREE(sigmaQNMlm0, COMPLEX16Vector);

    STRUCTFREE(chi2L_tPeakOmega, REAL8Vector);
    STRUCTFREE(chi1L_tPeakOmega, REAL8Vector);
    STRUCTFREE(nqcCoeffsList, SphHarmListEOBNonQCCoeffs);

    STRUCTFREE(listhPlm_HiS, SphHarmListCAmpPhaseSequence);
    STRUCTFREE(listhPlm_HiSRDpatch, SphHarmListCAmpPhaseSequence);
    STRUCTFREE(listhPlm_AdaS, SphHarmListCAmpPhaseSequence);

    // STRUCTFREE(tVecPmodes, REAL8Vector);
    // STRUCTFREE(listhPlm, SphHarmListCAmpPhaseSequence);
    // STRUCTFREE(seobdynamicsAdaSHiS, SEOBSAdynamics);

    // STRUCTFREE(alphaJ2P, REAL8Vector);
    // STRUCTFREE(betaJ2P, REAL8Vector);
    // STRUCTFREE(gammaJ2P, REAL8Vector);

    // STRUCTFREE(hIlm, SphHarmTimeSeries);
    // STRUCTFREE(hJlm, SphHarmTimeSeries);
    if (failed)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Program abort at step %d\n", this_step);
        return CEV_FAILURE;
    }
    return CEV_SUCCESS;
}



/* ------------------------------------------------------------------------------------------------ */
/*                                                                                                  */
/*                                                                                                  */
/*                                                                                                  */
/*                                                                                                  */
/*                                                                                                  */
/*                                                                                                  */
/*                                                                                                  */
/*                                    Elliptical Precession Codee                                   */
/*                                                                                                  */
/*                                                                                                  */
/*                                                                                                  */
/*                                                                                                  */
/*                                                                                                  */
/*                                                                                                  */
/*                                                                                                  */
/* ------------------------------------------------------------------------------------------------ */

INT evolve_prec(REAL8 m1,  REAL8 m2, 
           REAL8 s1x, REAL8 s1y, REAL8 s1z, 
           REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 phi0, REAL8 distance,
           REAL8 ecc, REAL8 zeta, REAL8 xi, REAL8 f_min, REAL8 Mf_ref, REAL8 INdeltaT, REAL8 inc,
           HyperParams *hparams, 
           SEOBPrecCoreOutputs *all)
{
    register INT i;
    INT failed = 0, this_step = 0, status = CEV_FAILURE;
    SpinEOBParams *core = NULL;
    REAL8Vector *ICvalues = NULL;
    REAL8Vector *ICvaluesHiS = NULL;
    REAL8Vector *seobvalues_tstartHiS = NULL;
    REAL8Vector *seobvalues_tPeakOmega = NULL;
    REAL8Vector *seobvalues_test = NULL;
    REAL8Vector* m1rVec = NULL;
    REAL8Array *dynamicsAdaS = NULL;
    REAL8Array *dynamicsInverse = NULL;
    REAL8Array *dynamicsHiS = NULL;
    REAL8Vector *Jfinal = NULL;
    REAL8Vector *Lhatfinal = NULL;
    COMPLEX16Vector *sigmaQNMlm0 = NULL;
    SphHarmListEOBNonQCCoeffs *nqcCoeffsList = NULL;
    SphHarmListCAmpPhaseSequence *listhPlm_HiS = NULL;
    SphHarmListCAmpPhaseSequence *listhPlm_HiSRDpatch = NULL;
    SphHarmListCAmpPhaseSequence *listhPlm_AdaS = NULL;
    SphHarmListCAmpPhaseSequence *listhClm = NULL;
    SEOBPrecdynamics *seobdynamicsAdaS = NULL;
    SEOBPrecdynamics *seobdynamicsHiS = NULL;
    
    REAL8Vector *chi1L_tPeakOmega = NULL;
    REAL8Vector *chi2L_tPeakOmega = NULL;

    // Output
    REAL8Vector *tVecPmodes = NULL;
    SphHarmListCAmpPhaseSequence *listhPlm = NULL;
    SEOBPrecdynamics *seobdynamicsAdaSHiS = NULL;
    REAL8Vector *alphaJ2P = NULL;
    REAL8Vector *betaJ2P = NULL;
    REAL8Vector *gammaJ2P = NULL;
    SphHarmTimeSeries *hIlm = NULL;
    SphHarmTimeSeries *hJlm = NULL;
    REAL8TimeSeries *hplusTS = NULL;
    REAL8TimeSeries *hcrossTS = NULL;

    UINT nmodes = 5;
    INT modes_lmax = 5;
    INT modes[5][2] = {{2,2}, {2,1}, {3,3}, {4,4}, {5,5}};
    // memset(modes, 0, 2 * nmodes * sizeof(INT));
    REAL8 mTotal = m1 + m2;
    REAL8 mTScaled = mTotal * CST_MTSUN_SI;
    REAL8 EPS_ALIGN = 1.0e-4;
    REAL8 deltaT = INdeltaT / mTScaled;
    REAL8 amp0 = mTotal * CST_MRSUN_SI / distance / 1e6 / CST_PC_SI;
    REAL8 tStepBack = 200.;
    INT flagZframe = FLAG_SEOBNRv4P_ZFRAME_L;
    hparams->flagZframe = flagZframe;
    /*
    *
    *       Check params
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Check params.", this_step);
    if (m1 < 0 || m2 < 0) {failed = 1; goto QUIT;}
    if (m2 > m1)
        SWAP(m1, m2);
    if (m1/m2 > 100) {failed = 1; goto QUIT;}
    if (sqrt(s1x*s1x+s1y*s1y+s1z*s1z) > 1. ||
        sqrt(s2x*s2x+s2y*s2y+s2z*s2z) > 1)
    {failed = 1; goto QUIT;}

    REAL8 freqMinRad = 0;
    /*Compute the highest initial frequency of the 22 mode */
    XLALEOBHighestInitialFreq(&freqMinRad, mTotal);
    if (f_min > freqMinRad)
    {
        PRINT_LOG_INFO(LOG_WARNING, "Initial frequency = %.16e is too high, the limit is %4.10f", f_min, freqMinRad);
        // failed = 1;
        // goto QUIT;
    }
    // print_debug("check NyquistFrequency\n");
    /* Check NyquistFrequency */
    UINT ell_max = 4;
    REAL8 INchi1[3] = {s1x, s1y, s1z};
    REAL8 INchi2[3] = {s2x, s2y, s2z};
    // print_debug("m1 = %.16e, m2 = %.16e\nINchi1 = (%.16e, %.16e, %.16e)\nINchi2 = (%.16e, %.16e, %.16e)\nINdeltaT = %.16e\n",
    //     m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, INdeltaT);
    if (hparams->Mf_min > hparams->Mf_max && hparams->t_max < 0) {
        status = XLALEOBCheckNyquistFrequency(m1, m2, INchi1, INchi2, INdeltaT, ell_max);
        if (status != CEV_SUCCESS)
        {failed = 1; goto QUIT;}
    }
    // print_debug("here\n");
    /*
    *
    * 
    *       Initialization
    *
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Initialization.", this_step);
    // print_debug("here\n");
    core = CreateSpinEOBParams(m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, ecc, hparams);
    if (!core) {failed = 1; goto QUIT;}
    tStepBack = GET_MAX(tStepBack, core->hParams->tStepBack);
    /*
    *
    * 
    *       Solve Initial Conditions
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Solve Initial Conditions.", this_step);
    ICvalues = CreateREAL8Vector(14);
    REAL8 MfMin = mTScaled * f_min;
    if (core->hParams && core->hParams->d_ini > 0.0)
    {
        REAL8 d_ini = core->hParams->d_ini; //initial seperation
        REAL8 pr_ini = core->hParams->pr_ini;
        REAL8 pphi_ini = core->hParams->pphi_ini;
        REAL8 ptheta_ini = core->hParams->ptheta_ini;
        REAL8 xSph[3] = {d_ini, 0., 0.};
        REAL8 pSph[3] = {pr_ini, ptheta_ini, pphi_ini};
        REAL8 xCart[3] = {0,0,0};
        REAL8 pCart[3] = {0,0,0};
        SphericalToCartesian(xCart, pCart, xSph, pSph);
        memset(ICvalues->data, 0, ICvalues->length*sizeof(REAL8));
        memcpy(ICvalues->data, xCart, 3*sizeof(REAL8));
        memcpy(ICvalues->data+3, pCart, 3*sizeof(REAL8));
        memcpy(ICvalues->data+6, core->s1Vec->data, 3*sizeof(REAL8));
        memcpy(ICvalues->data+9, core->s2Vec->data, 3*sizeof(REAL8));
    } else {
        if (ecc > 0.0 && get_egw_flag()) {
            //status = SEOBInitialConditions_egw(ICvalues, MfMin, ecc, core);
            status = SEOBInitialConditions_e_anomaly(ICvalues, Mf_ref, ecc, zeta, xi, core);
        } else
            status = SEOBInitialConditions(ICvalues, Mf_ref, ecc, core);
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    }
    PRINT_LOG_INFO(LOG_DEBUG, "initial conditions:");
    PRINT_LOG_INFO(LOG_DEBUG, "(x,y,z) = (%.16e %.16e %.16e)", 
            ICvalues->data[0], ICvalues->data[1], ICvalues->data[2]);
    PRINT_LOG_INFO(LOG_DEBUG, "(px,py,pz) = (%.16e %.16e %.16e)",
            ICvalues->data[3], ICvalues->data[4], ICvalues->data[5]);
    PRINT_LOG_INFO(LOG_DEBUG, "(S1x,S1y,S1z) = (%.16e %.16e %.16e)",
            ICvalues->data[6], ICvalues->data[7], ICvalues->data[8]);
    PRINT_LOG_INFO(LOG_DEBUG, "(S2x,S2y,S2z) = (%.16e %.16e %.16e)",
            ICvalues->data[9], ICvalues->data[10], ICvalues->data[11]);
    PRINT_LOG_INFO(LOG_DEBUG, "MfMin = %f, MfMax = %f, deltaT = %f\n", MfMin, core->hParams->Mf_max, deltaT);
    PRINT_LOG_INFO(LOG_DEBUG, "e0 = %f\n", ecc);
    REAL8 L0Vec[3] = {0,0,0};
    cross_product3d(ICvalues->data, ICvalues->data+3, L0Vec);
    core->J0Vec->data[0] += L0Vec[0];
    core->J0Vec->data[1] += L0Vec[1];
    core->J0Vec->data[2] += L0Vec[2];

#if 0
    // here we convert general Initial Condition to V1
    {
        REAL8 tmp_omega0, tmp_e0;
        status = SEOBConvertInitConditionFromGeneralToV1(ICvalues, core, &tmp_omega0, &tmp_e0);
        print_debug("omega0, e0 = %.16e, %.16e\n", tmp_omega0, tmp_e0);
    }
#endif
    // print_debug("g_ulPrintDebugLogFlag = %zu", g_ulPrintDebugLogFlag);
    // if (1) {failed = 1; goto QUIT;}

    /*
    *
    * 
    *       Evolve EOB trajectory with adaptive sampling 
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Evolve EOB trajectory with adaptive sampling.", this_step);
    INT retLenAdaS = 0, retLenInverse = 0;
    REAL8 tendAdaS = 20. / mTScaled;    
    REAL8 tstartAdaS = 0.;
    REAL8 deltaT_min = 8.0e-5;
    REAL8 EPS_ABS = 1.0e-8;
    REAL8 EPS_REL = 1.0e-8;
    if (hparams->inEPS_REL > 0.)
        EPS_REL = hparams->inEPS_REL;
    if (hparams->inEPS_ABS > 0.)
        EPS_ABS = hparams->inEPS_ABS;
    if (MfMin < Mf_ref)
    {
        // print_debug("MfMin, Mfref = %.16e, %.16e\n", MfMin, Mf_ref);
        status = SEOBIntegrateDynamics_prec_inverse(&dynamicsInverse, &retLenInverse, ICvalues, EPS_ABS, EPS_REL,
           deltaT, deltaT_min, tstartAdaS, tendAdaS, core, 0);
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    }
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    if (MfMin < hparams->Mf_max || hparams->t_max > 0.0)
    {
        /**
         * @brief 20240304 by X.L.:
         *          if Mf_max > MfMin,
         *              let integrate stops at omega = pi*Mf_max
         */
        status = SEOBIntegrateDynamics_prec_withFMax(&dynamicsAdaS, &retLenAdaS, 
            ICvalues, EPS_ABS, EPS_REL, 
            deltaT, deltaT_min, tstartAdaS, tendAdaS, core, core->alignedSpins);
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
        if (dynamicsInverse)
            SEOBConcactInverseDynToAdaSDynPrec(&dynamicsAdaS, dynamicsInverse, &retLenAdaS, retLenInverse);
        status = SEOBComputeExtendedSEOBPrecdynamics(&seobdynamicsAdaSHiS, dynamicsAdaS, retLenAdaS, core);    
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
        if (MfMin > Mf_ref)
        {
            status = CutSEOBPrecdynamics(&seobdynamicsAdaSHiS, MfMin);
            retLenAdaS = seobdynamicsAdaSHiS->length;
            if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
        }
        PRINT_LOG_INFO(LOG_DEBUG, "Get retLenAdaS = %d", retLenAdaS);
        tVecPmodes = CreateREAL8Vector(retLenAdaS);
        memcpy(tVecPmodes->data, seobdynamicsAdaSHiS->tVec, retLenAdaS*sizeof(REAL8));
        REAL8 tEndAtFMax = seobdynamicsAdaSHiS->tVec[retLenAdaS-1];
        PRINT_LOG_INFO(LOG_INFO, "Step %d_ Get final J/L/spins from HiS dynamics at peak of Omega, compute constant angles EulerI2J.", this_step);
        status = SEOBInterpolatePrecDynamicsAtTime(&seobvalues_tPeakOmega, tEndAtFMax, seobdynamicsAdaSHiS);
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
        SEOBLFrameVectors(&chi1L_tPeakOmega, &chi2L_tPeakOmega, 
            seobvalues_tPeakOmega, m1, m2, core->hParams->flagZframe);
        PRINT_LOG_INFO(LOG_DEBUG, "chi1L_tPeakOmega = (%.16e, %.16e, %.16e)\n", 
            chi1L_tPeakOmega->data[0], chi1L_tPeakOmega->data[1], chi1L_tPeakOmega->data[2]);
        PRINT_LOG_INFO(LOG_DEBUG, "chi2L_tPeakOmega = (%.16e, %.16e, %.16e)\n",
            chi2L_tPeakOmega->data[0], chi2L_tPeakOmega->data[1], chi2L_tPeakOmega->data[2]);
        /* Compute final J from dynamics quantities */
        SEOBJfromDynamics(&Jfinal, seobvalues_tPeakOmega, core);
        /*Compute the L-hat vector. Note that it has unit norm */
        // SEOBLhatfromDynamics(&Lhatfinal, seobvalues_tPeakOmega, core);
        PRINT_LOG_INFO(LOG_DEBUG, "Jfinal = (%.16e, %.16e, %.16e)", Jfinal->data[0], Jfinal->data[1], Jfinal->data[2]);
        // PRINT_LOG_INFO(LOG_DEBUG, "Lhatfinal = (%.16e, %.16e, %.16e)", Lhatfinal->data[0], Lhatfinal->data[1], Lhatfinal->data[2]);

        REAL8Vector e1J_fmax, e2J_fmax, e3J_fmax;
        e1J_fmax.length = e2J_fmax.length = e3J_fmax.length = 3;
        REAL8 e1Jdata_fmax[3] = {0.};
        REAL8 e2Jdata_fmax[3] = {0.};
        REAL8 e3Jdata_fmax[3] = {0.};
        e1J_fmax.data = e1Jdata_fmax;
        e2J_fmax.data = e2Jdata_fmax;
        e3J_fmax.data = e3Jdata_fmax;
        SEOBBuildJframeVectors(&e1J_fmax, &e2J_fmax, &e3J_fmax, Jfinal);
        PRINT_LOG_INFO(LOG_DEBUG, "e1J = (%.16e, %.16e, %.16e)", e1J_fmax.data[0], e1J_fmax.data[1], e1J_fmax.data[2]);
        PRINT_LOG_INFO(LOG_DEBUG, "e2J = (%.16e, %.16e, %.16e)", e2J_fmax.data[0], e2J_fmax.data[1], e2J_fmax.data[2]);
        PRINT_LOG_INFO(LOG_DEBUG, "e2J = (%.16e, %.16e, %.16e)", e3J_fmax.data[0], e3J_fmax.data[1], e3J_fmax.data[2]);

        /* Compute Euler angles from initial I-frame to final-J-frame */
        /* Note: if spins are aligned, the function SEOBEulerI2JFromJframeVectors */
        /* becomes ill-defined - just keep these Euler angles to zero then */
        REAL8 alphaI2J_fmax = 0., betaI2J_fmax = 0., gammaI2J_fmax = 0.;
        SEOBEulerI2JFromJframeVectors(&alphaI2J_fmax, &betaI2J_fmax, &gammaI2J_fmax, &e1J_fmax, &e2J_fmax, &e3J_fmax);

        SEOBPrecCalculateSphHarmListhlmAmpPhase(&listhPlm, modes, nmodes,
                                            seobdynamicsAdaSHiS, NULL,
                                            core, 0);
        PRINT_LOG_INFO(LOG_DEBUG, "Calculate hPlms done");
        status = SEOBPrecEulerJ2PFromDynamics(&alphaJ2P, &betaJ2P, &gammaJ2P, 
                    &e1J_fmax, &e2J_fmax, &e3J_fmax,
                    retLenAdaS, retLenAdaS-1,
                    seobdynamicsAdaSHiS, core);
        PRINT_LOG_INFO(LOG_DEBUG, "Calculate Euler J2P done");
        UINT retLenTS_fmax = floor(((tVecPmodes)->data[retLenAdaS - 1] - (tVecPmodes)->data[0]) / deltaT);
        status = SEOBRotateInterpolatehJlmReImFromSphHarmListhPlmAmpPhase(
            &hJlm, &listhClm, modes, nmodes, modes_lmax, deltaT, retLenTS_fmax, tVecPmodes,
            listhPlm, alphaJ2P, betaJ2P, gammaJ2P);
        if ( status == CEV_FAILURE) 
        {
            PRINT_LOG_INFO(LOG_CRITICAL, "failure in SEOBRotateInterpolatehJlmReImFromSphHarmListhPlmAmpPhase for mode (l,m) = (2,2).");
            failed = 1;
            goto QUIT;
        }
        PRINT_LOG_INFO(LOG_DEBUG, "Calculate hJlm done");
        status = SEOBRotatehIlmFromhJlm(&hIlm, hJlm, modes_lmax, alphaI2J_fmax, betaI2J_fmax, gammaI2J_fmax, deltaT);
        if (status != CEV_SUCCESS)
        {failed = 1; goto QUIT;}
        PRINT_LOG_INFO(LOG_DEBUG, "Calculate hIlm done");
        seobdynamicsAdaSHiS->th22Peak = tEndAtFMax;
        all->dyn = seobdynamicsAdaSHiS;
        all->hLM = hIlm;
        all->tVec = tVecPmodes;
        all->Plm = listhClm;
        PRINT_LOG_INFO(LOG_DEBUG, "Finished");
        goto QUIT;
    }
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    /******************************************/
    status = SEOBIntegrateDynamics_prec(&dynamicsAdaS, &retLenAdaS, 
        ICvalues, EPS_ABS, EPS_REL, 
        deltaT, deltaT_min, tstartAdaS, tendAdaS, core, 0);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    PRINT_LOG_INFO(LOG_DEBUG, "AdaS data length = %d", retLenAdaS);
    if (dynamicsInverse)
        SEOBConcactInverseDynToAdaSDynPrec(&dynamicsAdaS, dynamicsInverse, &retLenAdaS, retLenInverse);
    status = SEOBComputeExtendedSEOBPrecdynamics(&seobdynamicsAdaS, dynamicsAdaS, retLenAdaS, core);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    if (MfMin > Mf_ref)
    {
        status = CutSEOBPrecdynamics(&seobdynamicsAdaS, MfMin);
        retLenAdaS = seobdynamicsAdaS->length;
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    }
    m1rVec = CreateREAL8Vector(retLenAdaS);
    for (i = 0; i < retLenAdaS; i++) 
    {
        m1rVec->data[i] = -1* seobdynamicsAdaS->polarrVec[i];
    }
    UINT index_6M = FindClosestIndex(m1rVec, -6.0);
    REAL8 time_6M = seobdynamicsAdaS->tVec[index_6M];
    tStepBack = GET_MAX(tStepBack, seobdynamicsAdaS->tVec[retLenAdaS-1] - time_6M + 10*deltaT);
    
    /*
    *
    * 
    *       Step back and evolve EOB trajectory at high sampling rate (HiS)
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Step back and evolve EOB trajectory at high sampling rate (HiS).", this_step);
    REAL8 deltaTHiS = 1. / 50; /* Fixed at 1/50M */
    REAL8 tstartHiSTarget = seobdynamicsAdaS->tVec[retLenAdaS - 1] - tStepBack;
    INT4 indexstartHiS = retLenAdaS - 1; /* index for the AdaS dynamics */
    while ((indexstartHiS > 0) && (seobdynamicsAdaS->tVec[indexstartHiS] > tstartHiSTarget))
        indexstartHiS--;
    REAL8 tstartHiS = seobdynamicsAdaS->tVec[indexstartHiS];
    PRINT_LOG_INFO(LOG_DEBUG, "Interpolate Dynamics At Time %f (%f)/%f", tstartHiS, tstartHiSTarget, seobdynamicsAdaS->tVec[retLenAdaS - 1]);
    status = SEOBInterpolatePrecDynamicsAtTime(&seobvalues_tstartHiS, tstartHiS, seobdynamicsAdaS);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    ICvaluesHiS = CreateREAL8Vector(14);
    if (!ICvaluesHiS) {failed = 1; goto QUIT;}
    memcpy(ICvaluesHiS->data, &(seobvalues_tstartHiS->data[1]), 14 * sizeof(REAL8));
    /* Integrate again the dynamics with a constant high sampling rate */
    core->prev_dr = 0.; // This is used to check whether dr/dt is increasing
                            // in the stopping condition
    core->termination_reason =
        -999; // This is going to store the termination condition for the
                // high-sampling integration
    INT is_re_evolved = FALSE; 
    INT retLenHiS = 0;
    REAL8 tendHiS = tstartHiS + tStepBack; /* Seems it should be ignored anyway because of
                                    integrator->stopontestonly */
    REAL8 rISCO;
HISR:
    rISCO = get_h_rISCO();
    PRINT_LOG_INFO(LOG_DEBUG, "tStepBack = %g", tStepBack);
    PRINT_LOG_INFO(LOG_DEBUG, "Integrate Dynamics");
    status = SEOBIntegrateDynamics_prec(&dynamicsHiS, &retLenHiS, 
        ICvaluesHiS, EPS_ABS, EPS_REL, 
        deltaTHiS, 0., tstartHiS, tendHiS, core, 1);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    PRINT_LOG_INFO(LOG_DEBUG, "HiSR data length = %d", retLenHiS);
    status = SEOBComputeExtendedSEOBPrecdynamics(&seobdynamicsHiS, dynamicsHiS, retLenHiS, core);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
    /* Find time of peak of omega for the High-sampling dynamics */
    INT foundPeakOmega = 0;
    REAL8 tPeakOmega = 0.;
    SEOBPrecLocateTimePeakOmega(&tPeakOmega, &foundPeakOmega, dynamicsHiS,
                            seobdynamicsHiS, retLenHiS, core);
    PRINT_LOG_INFO(LOG_DEBUG, "Locate timePeak of omega = %f", tPeakOmega);
    /*
    *
    * 
    *       Get final J/L/spins from HiS dynamics at peak of Omega, compute constant angles EulerI2J
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Get final J/L/spins from HiS dynamics at peak of Omega, compute constant angles EulerI2J.", this_step);
    status = SEOBInterpolatePrecDynamicsAtTime(&seobvalues_tPeakOmega, tPeakOmega, seobdynamicsHiS);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}

    /* Interpolate the dynamics to r=10M. This is used as a
    feducial point to evaluate the final mass and spin fits */
    REAL8Vector timeVec;
    timeVec.length = retLenAdaS;
    timeVec.data = seobdynamicsAdaS->tVec;

    // FindClosestIndex requires an *increasing* function of time, so we use -r
    // instead of r
    UINT index_10M = FindClosestIndex(m1rVec, -10.0);
    // UINT index_10M = 1;
    REAL8 time_10M = timeVec.data[index_10M];
    PRINT_LOG_INFO(LOG_DEBUG, "Get r = 10M position at %d in (%.5f, %.5f)", index_10M,
        m1rVec->data[0], m1rVec->data[retLenAdaS-1]);
    status = SEOBInterpolatePrecDynamicsAtTime(&seobvalues_test, time_10M, seobdynamicsAdaS);
    if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}

    /* Compute the timeshift to get the attachment point */
    SEOBLFrameVectors(&chi1L_tPeakOmega, &chi2L_tPeakOmega, 
        seobvalues_tPeakOmega, m1, m2, core->hParams->flagZframe);
    PRINT_LOG_INFO(LOG_DEBUG, "chi1L_tPeakOmega = (%.16e, %.16e, %.16e)\n", 
        chi1L_tPeakOmega->data[0], chi1L_tPeakOmega->data[1], chi1L_tPeakOmega->data[2]);
    PRINT_LOG_INFO(LOG_DEBUG, "chi2L_tPeakOmega = (%.16e, %.16e, %.16e)\n",
        chi2L_tPeakOmega->data[0], chi2L_tPeakOmega->data[1], chi2L_tPeakOmega->data[2]);
    REAL8 Deltat22 = XLALSimIMREOBGetNRSpinPeakDeltaTv4(
        2, 2, m1, m2, chi1L_tPeakOmega->data[2], chi2L_tPeakOmega->data[2], core->hParams);

    /* Determine the time of attachment */
    REAL8 tAttach = tPeakOmega - Deltat22;
    PRINT_LOG_INFO(LOG_INFO, "The location of attachment is %f, Deltat22 = %f", tAttach, Deltat22);
    if (core->hParams->zero_dyncoaphase)
    {
        REAL8 de_phiD, de_phiM, de_phi;
        status = SetZeroPhaseAtTimePrec(seobdynamicsHiS, tAttach, &de_phiD, &de_phiM, &de_phi);
        if (status != CEV_SUCCESS) {failed = 1; goto QUIT;}
        OrbitalPhaseReducePrec(seobdynamicsAdaS, de_phiD, de_phiM, de_phi);
    }
    //Compute NQC window factors
    if (ecc != 0.0)
    {
        status = SEOBPrecCalculateNQCWindowFactorsFromDyn(seobdynamicsAdaS, tAttach, time_6M, tstartHiS, 6., 0, core);
        // core->wWind = 1.;
    }

    if (is_re_evolved)
    {
        PRINT_LOG_INFO(LOG_INFO, "re-calculating window parameters; rISCO = %f, rmin = %f", rISCO, seobdynamicsHiS->polarrVec[retLenHiS-1]);
        status = SEOBPrecCalculateNQCWindowFactorsFromDyn(seobdynamicsHiS, tAttach, time_6M, tstartHiS, rISCO, 1, core);
    }

    /* Compute final J from dynamics quantities */
    SEOBJfromDynamics(&Jfinal, seobvalues_tPeakOmega, core);

    /*Compute the L-hat vector. Note that it has unit norm */
    SEOBLhatfromDynamics(&Lhatfinal, seobvalues_tPeakOmega, core);

    REAL8 Jmag = sqrt(inner_product3d(Jfinal->data, Jfinal->data));
    /* Cosine of the angle between L-hat and J. Needed to determine
    * the correct sign of the final spin
    */

    REAL8 cos_angle = inner_product3d(Jfinal->data, Lhatfinal->data) / Jmag;
    REAL8Vector e1J, e2J, e3J;
    e1J.length = e2J.length = e3J.length = 3;
    REAL8 e1Jdata[3] = {0.};
    REAL8 e2Jdata[3] = {0.};
    REAL8 e3Jdata[3] = {0.};
    e1J.data = e1Jdata;
    e2J.data = e2Jdata;
    e3J.data = e3Jdata;
    SEOBBuildJframeVectors(&e1J, &e2J, &e3J, Jfinal);
    /* Compute Euler angles from initial I-frame to final-J-frame */
    /* Note: if spins are aligned, the function SEOBEulerI2JFromJframeVectors */
    /* becomes ill-defined - just keep these Euler angles to zero then */
    REAL8 alphaI2J = 0., betaI2J = 0., gammaI2J = 0.;
    SEOBEulerI2JFromJframeVectors(&alphaI2J, &betaI2J, &gammaI2J, &e1J, &e2J, &e3J);
    /*
    *
    * 
    *       Compute P-frame amp/phase for all modes on HiS and compute NQC Coeffs
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Compute P-frame amp/phase for all modes on HiS and compute NQC Coeffs.", this_step);
    // REAL8 wWind_old = core->wWind;
    status = SEOBPrecCalculateSphHarmListNQCCoefficientsV4(
            &nqcCoeffsList, modes, nmodes, tPeakOmega, seobdynamicsHiS,
            core, chi1L_tPeakOmega, chi2L_tPeakOmega);
    if ( status != CEV_SUCCESS) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "NQC computation failed.");
        failed = 1;
        goto QUIT;
    }
    core->wWind = 1.;
    /*
    *
    *       Compute P-frame amp/phase for all modes on HiS, now including NQC
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Compute P-frame amp/phase for all modes on HiS, now including NQC.", this_step);
    SEOBPrecCalculateSphHarmListhlmAmpPhase(&listhPlm_HiS, modes, nmodes,
                                        seobdynamicsHiS, nqcCoeffsList,
                                        core, 1);
    /*
    *
    * 
    *       Attach RD to the P-frame waveform
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Attach RD to the P-frame waveform.", this_step);
    REAL8 finalMass = 0., finalSpin = 0.;
    status = SEOBGetFinalSpinMass(&finalMass, &finalSpin, seobvalues_test, core);
    if (status != CEV_SUCCESS) 
    {failed = 1; goto QUIT;}

    /* The function above returns only the magnitude of the spin.
    *  We pick the direction based on whether Lhat \cdot J is positive
    *  or negative */
    if (cos_angle < 0)
    {
        finalSpin *= -1;
    }
    /* finalSpin interpolation is available only between -0.9996 and 0.9996 */
    /* Set finalSpin to +/- 0.9996 if it is out of this range */
    if (finalSpin < -0.9996)
        finalSpin = -0.9996;
    if (finalSpin > 0.9996)
        finalSpin = 0.9996;

    PRINT_LOG_INFO(LOG_INFO, "final mass = %e, final spin = %e", finalMass, finalSpin);

    /* Estimate leading QNM , to determine how long the ringdown
    * patch will be */
    /* NOTE: XLALSimIMREOBGenerateQNMFreqV2Prec returns the complex frequency in
    * physical units... */
    COMPLEX16Vector sigmaQNM220estimatephysicalVec;
    COMPLEX16 sigmaQNM220estimatephysical = 0.;
    sigmaQNM220estimatephysicalVec.length = 1;
    sigmaQNM220estimatephysicalVec.data = &sigmaQNM220estimatephysical;
    status = XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec(&sigmaQNM220estimatephysicalVec, m1,
                                                m2, finalMass, finalSpin, 2,
                                                2, 1);
    if (status == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec for mode (l,m) = (2,2).");
        failed = 1;
        goto QUIT;
    }
    COMPLEX16 sigmaQNM220estimate = mTScaled * sigmaQNM220estimatephysical;

    /* Length of RD patch, 40 e-folds of decay of the estimated QNM220 */
    UINT retLenRDPatch =
    (UINT)ceil(EFOLDS / (cimag(sigmaQNM220estimate) * deltaTHiS));

    /* Attach RD to the P-frame modes */
    /* Vector holding the values of the 0-th overtone QNM complex frequencies for
    * the modes (l,m) */
    // NOTE: the QNM complex frequencies are computed inside
    // SEOBAttachRDToSphHarmListhPlm
    status = SEOBPrecAttachRDToSphHarmListhPlm(
        &listhPlm_HiSRDpatch, &sigmaQNMlm0, modes, nmodes, finalMass, finalSpin,
        listhPlm_HiS, deltaTHiS, retLenHiS, retLenRDPatch, tAttach,
        seobvalues_tPeakOmega, seobdynamicsHiS, core);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}

    /*
    *
    * 
    *       Build the joined dynamics AdaS+HiS up to attachment, joined P-modes AdaS+HiS+RDpatch 
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Build the joined dynamics AdaS+HiS up to attachment, joined P-modes AdaS+HiS+RDpatch", this_step);
    /* Compute amplitude and phase of the P-frame modes hPlm on adaptive sampling,
    * with NQC */
    // flagNQC = 1;
    SEOBPrecCalculateSphHarmListhlmAmpPhase(&listhPlm_AdaS, modes, nmodes,
                                    seobdynamicsAdaS, nqcCoeffsList,
                                    core, 1);
    /* Vector of times for the P-modes, that will be used for interpolation:
    * joining AdaS and HiS+RDpatch */
    UINT retLenPmodes = 0;
    /* First junction at indexAdaSHiS, tAdaSHiS */
    UINT indexJoinHiS = 0;
    REAL8 tJoinHiS = 0.;
    /* Second junction at seobdynamicsAdaS, tJoinAttach */
    UINT indexJoinAttach = 0;
    REAL8 tJoinAttach = 0.;
    /* Construct the joined vector of times (AdaS+HiS+RDpatch) and keep the
    * jonction indices and times */
    status = SEOBPrecJoinTimeVector(&tVecPmodes, &retLenPmodes, &tJoinHiS, &indexJoinHiS,
                    &tJoinAttach, &indexJoinAttach, retLenRDPatch, deltaTHiS,
                    tstartHiS, tAttach, seobdynamicsAdaS, seobdynamicsHiS);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}
    PRINT_LOG_INFO(LOG_DEBUG, "tJoinAttach = %.16e, indexJoinAttach = %u", tJoinAttach, indexJoinAttach);
    /* Copy dynamics from AdaS<HiS and HiS<tAttach to form joined dynamics, ending
    * at the last time sample <tAttach */
    // NOTE: we cut the dynamics at tAttach, as we will extend the Euler
    // angles for t>=tAttach -- but we could also choose to finish at tPeakOmega
    // which is used for final-J and for the final mass/spin fit

    status = SEOBPrecJoinDynamics(&seobdynamicsAdaSHiS, seobdynamicsAdaS, seobdynamicsHiS,
                indexJoinHiS, indexJoinAttach);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}

    /* Copy waveform modes from AdaS and HiS+RDpatch - adjusting 2pi-phase shift
    * at the junction point AdaS/HiS */
    status = SEOBJoinSphHarmListhlm(&listhPlm, listhPlm_AdaS, listhPlm_HiSRDpatch, modes,
                        nmodes, indexstartHiS);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}

    /* Get the time of the frame-invariant amplitude peak */
    REAL8 tPeak = 0;
    UINT indexPeak = 0;
    // NOTE: peak amplitude using l=2 only: h22 required, and h21 used if present
    status = SEOBAmplitudePeakFromAmp22Amp21(&tPeak, &indexPeak, listhPlm, modes, nmodes,
                                tVecPmodes);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}
    PRINT_LOG_INFO(LOG_DEBUG, "tPeak = %.16e, indexPeak = %u", tPeak, indexPeak);
    /*
    *
    * 
    *       Compute Euler angles J2P from AdaS and HiS dynamics up to attachment
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Compute Euler angles J2P from AdaS and HiS dynamics up to attachment", this_step);
    /* Compute Euler angles J2P from the dynamics before attachment point */
    /* If SpinsAlmostAligned, all Euler angles are set to 0 */
    status = SEOBPrecEulerJ2PFromDynamics(&alphaJ2P, &betaJ2P, &gammaJ2P, 
                &e1J, &e2J, &e3J,
                retLenPmodes, indexJoinAttach,
                seobdynamicsAdaSHiS, core);
    if (status == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "SEOBEulerJ2PFromDynamics failed.");
        failed = 1;
        goto QUIT;
    }
    /*
    * 
    * 
    *       Compute Euler angles J2P extension after attachment
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Compute Euler angles J2P extension after attachment", this_step);

    /* Compute Euler angles J2P according to the prescription flagEulerextension
    * after attachment point */
    /* If SpinsAlmostAligned, all Euler angles are set to 0 */
    /* NOTE: Regardless of the mode content of hPlm, the frame extension at the
    * moment is based on sigmaQNM22, sigmaQNM21 */
    COMPLEX16 sigmaQNM220 = 0., sigmaQNM210 = 0.;
    COMPLEX16Vector sigmaQNM220physicalVec, sigmaQNM210physicalVec;
    sigmaQNM220physicalVec.length = 1;
    sigmaQNM210physicalVec.length = 1;
    COMPLEX16 sigmaQNM220physicalval = 0.;
    COMPLEX16 sigmaQNM210physicalval = 0.;
    sigmaQNM220physicalVec.data = &sigmaQNM220physicalval;
    sigmaQNM210physicalVec.data = &sigmaQNM210physicalval;
    /* NOTE: XLALSimIMREOBGenerateQNMFreqV2Prec returns the complex frequency in
    * physical units... */
    status = XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec(&sigmaQNM220physicalVec, m1,
                                                m2, finalMass, finalSpin, 2,
                                                2, 1);
    if (status == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec for mode (l,m) = (2,2).");
        failed = 1;
        goto QUIT;
    }
    status = XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec(&sigmaQNM210physicalVec, m1,
                                                m2, finalMass, finalSpin, 2,
                                                1, 1);
    if (status == CEV_FAILURE) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec for mode (l,m) = (2,1).");
        failed = 1;
        goto QUIT;
    }
    sigmaQNM220 = mTScaled * sigmaQNM220physicalVec.data[0];
    sigmaQNM210 = mTScaled * sigmaQNM210physicalVec.data[0];
    INT flip = 1;
    if (cos_angle < 0)
        flip = -1;
    // flagEulerextension = 0
    status = SEOBEulerJ2PPostMergerExtension(
        alphaJ2P, betaJ2P, gammaJ2P, sigmaQNM220, sigmaQNM210, tVecPmodes,
        retLenPmodes, indexJoinAttach, core, flip);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}

    /*
    *
    * 
    *       Compute modes hJlm on the output time series by rotating and interpolating the modes hPlm
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Compute modes hJlm on the output time series by rotating and interpolating the modes hPlm", this_step);

    /* Determine the length of the fixed-sampling output time series */
    UINT retLenTS = floor(((tVecPmodes)->data[retLenPmodes - 1] - (tVecPmodes)->data[0]) / deltaT);

    /* Rotate waveform from P-frame to J-frame */
    // flagSymmetrizehPlminusm = 1
    status = SEOBRotateInterpolatehJlmReImFromSphHarmListhPlmAmpPhase(
        &hJlm, &listhClm, modes, nmodes, modes_lmax, deltaT, retLenTS, tVecPmodes,
        listhPlm, alphaJ2P, betaJ2P, gammaJ2P);
    if ( status == CEV_FAILURE )
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failure in SEOBRotateInterpolatehJlmReImFromSphHarmListhPlmAmpPhase for mode (l,m) = (2,2).");
        failed = 1;
        goto QUIT;
    }

    /*
    * 
    * 
    *       Rotate waveform from J-frame to the output I-frame on timeseries-sampling (constant Wigner coeffs)
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "STEP %d_ Rotate waveform from J-frame to the output I-frame on timeseries-sampling (constant Wigner coeffs)", this_step);

    /* Rotate waveform from J-frame to I-frame */
    status = SEOBRotatehIlmFromhJlm(&hIlm, hJlm, modes_lmax, alphaI2J, betaI2J, gammaI2J, deltaT);
    if (status != CEV_SUCCESS)
    {failed = 1; goto QUIT;}
    hIlm->tAttach = tAttach;
    PRINT_LOG_INFO(LOG_DEBUG, "alphaI2J, betaI2J, gammaI2J = %.16e, %.16e, %.16e", alphaI2J, betaI2J, gammaI2J);

    /*
    *
    * 
    *       Output 
    * 
    * 
    */
    this_step++;
    PRINT_LOG_INFO(LOG_INFO, "Step %d_ Set Output.", this_step);
    all->dyn = seobdynamicsAdaSHiS;
    all->hLM = hIlm;
    all->tVec = tVecPmodes;
    all->Plm = listhClm;

QUIT:

    PRINT_LOG_INFO(LOG_INFO, "Program END\n\n");
    STRUCTFREE(core, SpinEOBParams);

    STRUCTFREE(dynamicsAdaS, REAL8Array);
    STRUCTFREE(dynamicsInverse, REAL8Array);
    STRUCTFREE(ICvalues, REAL8Vector);

    STRUCTFREE(dynamicsHiS, REAL8Array);
    STRUCTFREE(ICvaluesHiS, REAL8Vector);
    STRUCTFREE(seobvalues_tstartHiS, REAL8Vector);
    STRUCTFREE(seobdynamicsAdaS, SEOBPrecdynamics);
    STRUCTFREE(seobdynamicsHiS, SEOBPrecdynamics);

    STRUCTFREE(seobvalues_tPeakOmega, REAL8Vector);
    STRUCTFREE(seobvalues_test, REAL8Vector);
    STRUCTFREE(m1rVec,REAL8Vector);
    STRUCTFREE(Jfinal, REAL8Vector);
    STRUCTFREE(Lhatfinal, REAL8Vector);
    STRUCTFREE(sigmaQNMlm0, COMPLEX16Vector);

    STRUCTFREE(chi2L_tPeakOmega, REAL8Vector);
    STRUCTFREE(chi1L_tPeakOmega, REAL8Vector);
    STRUCTFREE(nqcCoeffsList, SphHarmListEOBNonQCCoeffs);

    STRUCTFREE(listhPlm_HiS, SphHarmListCAmpPhaseSequence);
    STRUCTFREE(listhPlm_HiSRDpatch, SphHarmListCAmpPhaseSequence);
    STRUCTFREE(listhPlm_AdaS, SphHarmListCAmpPhaseSequence);

    // STRUCTFREE(tVecPmodes, REAL8Vector);
    STRUCTFREE(listhPlm, SphHarmListCAmpPhaseSequence);
    // STRUCTFREE(seobdynamicsAdaSHiS, SEOBPrecdynamics);

    STRUCTFREE(alphaJ2P, REAL8Vector);
    STRUCTFREE(betaJ2P, REAL8Vector);
    STRUCTFREE(gammaJ2P, REAL8Vector);

    // STRUCTFREE(hIlm, SphHarmTimeSeries);
    STRUCTFREE(hJlm, SphHarmTimeSeries);
    if (failed)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Program abort at step %d\n", this_step);
        return CEV_FAILURE;
    }
    return CEV_SUCCESS;
}
