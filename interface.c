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
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include "myLog.h"
#include "myFileIO.h"

int debug();

INT get_PrecFlag();
void set_PrecFlag(INT flag);

INT set_egw_flag(INT flag);

INT evolve(REAL8 m1,  REAL8 m2, 
           REAL8 s1x, REAL8 s1y, REAL8 s1z, 
           REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 phi0, REAL8 distance,
           REAL8 ecc, REAL8 f_min, REAL8 INdeltaT, REAL8 inc,
           INT is_only22,
           HyperParams *hparams, 
           REAL8TimeSeries **hPlusOut,
           REAL8TimeSeries **hCrossOut,
           SEOBCoreOutputs *all);

INT evolve_adaptive(REAL8 m1,  REAL8 m2, 
           REAL8 s1x, REAL8 s1y, REAL8 s1z, 
           REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 phi0, REAL8 distance,
           REAL8 ecc, REAL8 f_min, REAL8 INdeltaT, REAL8 inc,
           INT is_only22,
           HyperParams *hparams, 
           REAL8TimeSeries **hPlusOut,
           REAL8TimeSeries **hCrossOut,
           SEOBCoreOutputs *all);

INT evolve_conserv(REAL8 m1,  REAL8 m2, 
           REAL8 s1x, REAL8 s1y, REAL8 s1z, 
           REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 phi0, REAL8 distance,
           REAL8 ecc, REAL8 f_min, REAL8 INdeltaT, REAL8 inc,
           HyperParams *hparams, 
           SEOBCoreOutputs *all);

INT evolve_SA(REAL8 m1,  REAL8 m2, 
           REAL8 s1z, 
           REAL8 s2z,
           REAL8 ecc, REAL8 f_min, REAL8 INdeltaT,
           HyperParams *hparams, 
           REAL8Vector **tRetVec,
           SphHarmListCAmpPhaseSequence **hlm,
           int is_noringdown, int is_dyn_debug,
           SEOBSAdynamics **dyn_debug);
INT evolve_prec(REAL8 m1,  REAL8 m2, 
           REAL8 s1x, REAL8 s1y, REAL8 s1z, 
           REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 phi0, REAL8 distance,
           REAL8 ecc, REAL8 f_min, REAL8 INdeltaT, REAL8 inc,
           HyperParams *hparams, 
           SEOBPrecCoreOutputs *all);

INT choose_debug(INT debug_id,
    REAL8 m1,  REAL8 m2, 
    REAL8 s1x, REAL8 s1y, REAL8 s1z, 
    REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 phi0, REAL8 distance,
    REAL8 ecc, REAL8 f_min, REAL8 INdeltaT, REAL8 inc, HyperParams *hparams);

typedef struct {
    INT use_geom;
    REAL8 m1;
    REAL8 m2;
    REAL8 s1x;
    REAL8 s1y;
    REAL8 s1z;
    REAL8 s2x;
    REAL8 s2y;
    REAL8 s2z;
    REAL8 ecc;
    REAL8 inc;
    REAL8 phiRef;
    REAL8 beta;
    REAL8 f_min;
    REAL8 distance;
    REAL8 srate;
    REAL8 deltaT;
    INT is_only22;
    INT is_constp;
    INT conserve_flag;
    REAL8 conserve_time;
    INT is_noringdown;
    INT prec_flag;
    INT debug_id;
    REAL8 d_ini;
    REAL8 pr_ini;
    REAL8 pphi_ini;
    REAL8 ptheta_ini;
    INT flagTuning;
    REAL8 tStepBack;
    REAL8 sl_p;
    REAL8 x0;
    INT code_version;
    INT log_level;
    REAL8 risco;
    REAL8 KK;
    REAL8 dSS;
    REAL8 dSO;
    REAL8 dtPeak;
    INT egw_flag;
    INT ret_dyn;
    REAL8 inEPS_REL;
    REAL8 inEPS_ABS;
    INT is_coframe;
}pyInputParams_t;

typedef struct {
    INT length;
    REAL8Vector *time;
    REAL8Vector *timeM;
    REAL8Vector *hplus;
    REAL8Vector *hcross;
    REAL8Vector *h22_real;
    REAL8Vector *h22_imag;

    REAL8Vector *h21_real;
    REAL8Vector *h21_imag;

    REAL8Vector *h33_real;
    REAL8Vector *h33_imag;

    REAL8Vector *h44_real;
    REAL8Vector *h44_imag;

    REAL8Vector *h55_real;
    REAL8Vector *h55_imag;
}pyOutputStruct_t;

typedef struct {
    INT length;
    REAL8Vector *timeM;
    REAL8Vector *xVec;
    REAL8Vector *yVec;
    REAL8Vector *zVec;
    REAL8Vector *pTxVec;
    REAL8Vector *pTyVec;
    REAL8Vector *pTzVec;
    REAL8Vector *vxVec;
    REAL8Vector *vyVec;
    REAL8Vector *vzVec;
    REAL8Vector *s1xVec;
    REAL8Vector *s1yVec;
    REAL8Vector *s1zVec;
    REAL8Vector *s2xVec;
    REAL8Vector *s2yVec;
    REAL8Vector *s2zVec;
    REAL8Vector *prTDotVec;
    REAL8Vector *hamVec;
}pyDynOutputStruct_t;

pyOutputStruct_t *CreatePyOutputStruct_t(INT length)
{
    pyOutputStruct_t *output = (pyOutputStruct_t *) MYMalloc(sizeof(pyOutputStruct_t));
    output->length = length;
    output->time = CreateREAL8Vector(length);
    output->timeM = CreateREAL8Vector(length);
    output->hplus = CreateREAL8Vector(length);
    output->hcross = CreateREAL8Vector(length);
    output->h22_real = CreateREAL8Vector(length);
    output->h22_imag = CreateREAL8Vector(length);

    output->h21_real = CreateREAL8Vector(length);
    output->h21_imag = CreateREAL8Vector(length);

    output->h33_real = CreateREAL8Vector(length);
    output->h33_imag = CreateREAL8Vector(length);

    output->h44_real = CreateREAL8Vector(length);
    output->h44_imag = CreateREAL8Vector(length);

    output->h55_real = CreateREAL8Vector(length);
    output->h55_imag = CreateREAL8Vector(length);
    return output;
}

pyDynOutputStruct_t *CreatepyDynOutputStruct_t(INT length)
{
    pyDynOutputStruct_t *output = (pyDynOutputStruct_t *) MYMalloc(sizeof(pyDynOutputStruct_t));
    output->length = length;
    output->timeM = CreateREAL8Vector(length);

    output->xVec = CreateREAL8Vector(length);
    output->yVec = CreateREAL8Vector(length);
    output->zVec = CreateREAL8Vector(length);

    output->vxVec = CreateREAL8Vector(length);
    output->vyVec = CreateREAL8Vector(length);
    output->vzVec = CreateREAL8Vector(length);

    output->pTxVec = CreateREAL8Vector(length);
    output->pTyVec = CreateREAL8Vector(length);
    output->pTzVec = CreateREAL8Vector(length);

    output->s1xVec = CreateREAL8Vector(length);
    output->s1yVec = CreateREAL8Vector(length);
    output->s1zVec = CreateREAL8Vector(length);

    output->s2xVec = CreateREAL8Vector(length);
    output->s2yVec = CreateREAL8Vector(length);
    output->s2zVec = CreateREAL8Vector(length);
    output->prTDotVec = CreateREAL8Vector(length);
    output->hamVec = CreateREAL8Vector(length);
    return output;
}

void DestroypyOutputStruct_t(pyOutputStruct_t *out)
{
    if (!out)
        return;
    STRUCTFREE(out->time, REAL8Vector);
    STRUCTFREE(out->timeM, REAL8Vector);
    STRUCTFREE(out->hplus, REAL8Vector);
    STRUCTFREE(out->hcross, REAL8Vector);

    STRUCTFREE(out->h22_real, REAL8Vector);
    STRUCTFREE(out->h22_imag, REAL8Vector);

    STRUCTFREE(out->h21_real, REAL8Vector);
    STRUCTFREE(out->h21_imag, REAL8Vector);

    STRUCTFREE(out->h33_real, REAL8Vector);
    STRUCTFREE(out->h33_imag, REAL8Vector);

    STRUCTFREE(out->h44_real, REAL8Vector);
    STRUCTFREE(out->h44_imag, REAL8Vector);

    STRUCTFREE(out->h55_real, REAL8Vector);
    STRUCTFREE(out->h55_imag, REAL8Vector);
    MYFree(out);
    out = NULL;
    return;
}

void DestroypyDynOutputStruct_t(pyDynOutputStruct_t *out)
{
    if (!out)
        return;
    STRUCTFREE(out->timeM, REAL8Vector);
    STRUCTFREE(out->xVec, REAL8Vector);
    STRUCTFREE(out->yVec, REAL8Vector);
    STRUCTFREE(out->zVec, REAL8Vector);
    STRUCTFREE(out->vxVec, REAL8Vector);
    STRUCTFREE(out->vyVec, REAL8Vector);
    STRUCTFREE(out->vzVec, REAL8Vector);
    STRUCTFREE(out->pTxVec, REAL8Vector);
    STRUCTFREE(out->pTyVec, REAL8Vector);
    STRUCTFREE(out->pTzVec, REAL8Vector);
    STRUCTFREE(out->s1xVec, REAL8Vector);
    STRUCTFREE(out->s1yVec, REAL8Vector);
    STRUCTFREE(out->s1zVec, REAL8Vector);
    STRUCTFREE(out->s2xVec, REAL8Vector);
    STRUCTFREE(out->s2yVec, REAL8Vector);
    STRUCTFREE(out->s2zVec, REAL8Vector);
    STRUCTFREE(out->prTDotVec, REAL8Vector);
    STRUCTFREE(out->hamVec, REAL8Vector);
    MYFree(out);
    return;
}

void convert_SEOBPrecCoreOutputs_to_pyOutputStruct_t(INT is_only22, REAL8 mtot, REAL8 dL, REAL8 inc, REAL8 phic, REAL8 beta, SEOBPrecCoreOutputs *All_prec, pyOutputStruct_t **ret);
void convert_SEOBCoreOutputs_to_pyOutputStruct_t(INT is_only22, REAL8 mtot, REAL8 dL, REAL8 inc, REAL8 phic, REAL8 beta, SEOBCoreOutputs *All_prec, pyOutputStruct_t **ret);
void convert_SphHarmListCAmpPhaseSequence_to_pyOutputStruct_t(INT is_only22, REAL8 mtot, REAL8 dL, REAL8 inc, REAL8 phic, REAL8 beta, REAL8Vector *tVec, SphHarmListCAmpPhaseSequence *hLM, SEOBSAdynamics *dyn_debug, pyOutputStruct_t **ret);
void convert_PrecSphHarmListCAmpPhaseSequence_to_pyOutputStruct_t(INT is_only22, REAL8 mtot, REAL8 dL, REAL8 inc, REAL8 phic, REAL8 beta, REAL8Vector *tVec, SphHarmListCAmpPhaseSequence *PLM, pyOutputStruct_t **ret);

void convert_SEOBPrecCoreOutputs_to_pyDynOutputStruct_t(SEOBPrecCoreOutputs *All_prec, pyDynOutputStruct_t **ret);
void convert_SEOBCoreOutputs_to_pyDynOutputStruct_t(SEOBCoreOutputs *All, pyDynOutputStruct_t **ret);
void convert_SEOBSAdynamics_to_pyDynOutputStruct_t(SEOBSAdynamics *dyn_debug, REAL8 m1, REAL8 m2, REAL8 chi1, REAL8 chi2, pyDynOutputStruct_t **ret);

INT generate_waveform(pyInputParams_t *params, pyOutputStruct_t **output, pyDynOutputStruct_t **dynoutput)
{
    clock_t tstart, tend;
    tstart = clock();
    HyperParams hparams;
    CtrlParams ctpms;
    ctpms.level = params->log_level;
    strncpy(ctpms.flog, "None", STR_COMM_SIZE);
    strncpy(ctpms.fdebug, "conserve.h5", STR_COMM_SIZE);
    INT use_geom;
    use_geom = params->use_geom;
    hparams.d_ini = params->d_ini;
    hparams.pr_ini = params->pr_ini;
    hparams.pphi_ini = params->pphi_ini;
    hparams.ptheta_ini = params->ptheta_ini;
    hparams.flagTuning = params->flagTuning;
    hparams.tStepBack = params->tStepBack;
    hparams.sl_p = params->sl_p;
    hparams.x0 = params->x0;
    hparams.KK = params->KK;
    hparams.dSS = params->dSS;
    hparams.dSO = params->dSO;
    hparams.dtPeak = params->dtPeak;
    hparams.flagZframe = FLAG_SEOBNRv4P_ZFRAME_L;
    hparams.inEPS_REL = params->inEPS_REL;
    hparams.inEPS_ABS = params->inEPS_ABS;
    hparams.is_coframe = params->is_coframe;
    set_egw_flag(params->egw_flag);
    SET_CODE_VERSION(params->code_version);
    set_PrecFlag(params->prec_flag);
    if (params->risco > 0.0)
        SET_RISCO(params->risco);
    if (params->conserve_time > 0.0)
        SET_CONSERV(params->conserve_flag, params->conserve_time);
    if (use_geom)
    {
        REAL8 mTScaled = (params->m1+params->m2) * CST_MTSUN_SI;
        params->f_min = params->f_min / mTScaled;
    }
    if (ctpms.level >= 10)
    {
        SET_INPUT_DEBUG_FLAG(ctpms.level - 10);
    }
    if (strcmp(ctpms.flog, "None")!=0)
    {
        LOG_SetPrintLogPlaceFlag(1);
    }
    LOG_SetPrintDebugLogFlag(ctpms.level);
    LOG_Init(ctpms.flog, 40960);
    PRINT_LOG_INFO(LOG_DEBUG, "params: m1, m2 = (%g, %g)\n\tchi1Vec = (%g, %g, %g), chi2Vec = (%g, %g, %g)\n\teccentricity = %g, deltaT = %g,\n\tf-min=%g, inclination=%g, phiRef=%g\n", 
        params->m1, params->m2, params->s1x, params->s1y, params->s1z, params->s2x, params->s2y, params->s2z, params->ecc, params->deltaT, params->f_min, params->inc * 180/CST_PI, params->beta*180/CST_PI);
    // print_debug("params: m1, m2 = (%g, %g)\n\tchi1Vec = (%g, %g, %g), chi2Vec = (%g, %g, %g)\n\teccentricity = %g, deltaT = %g,\n\tf-min=%g, inclination=%g, phiRef=%g\n", 
    //     params->m1, params->m2, params->s1x, params->s1y, params->s1z, params->s2x, params->s2y, params->s2z, params->ecc, params->deltaT, params->f_min, params->inc * 180/CST_PI, params->phiRef*180/CST_PI);
    // return 0;
    if (params->debug_id)
    {
        choose_debug(params->debug_id, params->m1, params->m2,
            params->s1x, params->s1y, params->s1z,
            params->s2x, params->s2y, params->s2z,
            params->beta, params->distance, params->ecc, 
            params->f_min, params->deltaT, params->inc,
            &hparams);
        return CEV_SUCCESS;
    }
    // call waveform generator
    INT status, is_failed;
    is_failed = 0;
    if (get_PrecFlag() != 0) 
    {
        SEOBPrecCoreOutputs *All_prec = (SEOBPrecCoreOutputs *)MYCalloc(1, sizeof(SEOBPrecCoreOutputs));
        PRINT_LOG_INFO(LOG_DEBUG, "prec code");
        status = evolve_prec(params->m1, params->m2, 
            params->s1x, params->s1y, params->s1z, 
            params->s2x, params->s2y, params->s2z, 
            params->beta, params->distance, params->ecc, params->f_min,
            params->deltaT, params->inc, &hparams, All_prec);
        if (status == CEV_SUCCESS)
        {
            if (hparams.is_coframe)
                convert_PrecSphHarmListCAmpPhaseSequence_to_pyOutputStruct_t(params->is_only22, params->m1 + params->m2, params->distance, params->inc, params->phiRef, params->beta, All_prec->tVec, All_prec->Plm, output);
            else
                convert_SEOBPrecCoreOutputs_to_pyOutputStruct_t(params->is_only22, params->m1 + params->m2, params->distance, params->inc, params->phiRef, params->beta, All_prec, output);
            if (params->ret_dyn)
                convert_SEOBPrecCoreOutputs_to_pyDynOutputStruct_t(All_prec, dynoutput);
        }
        else
            is_failed = 1;
        STRUCTFREE(All_prec, SEOBPrecCoreOutputs);
    } else if (CONSERVE_FLAG != 0) {
        PRINT_LOG_INFO(LOG_DEBUG, "conserve code");
        SEOBCoreOutputs *All = (SEOBCoreOutputs *)MYCalloc(1, sizeof(SEOBCoreOutputs));
        status = evolve_conserv(params->m1, params->m2, 
            params->s1x, params->s1y, params->s1z, 
            params->s2x, params->s2y, params->s2z, 
            params->beta, params->distance, params->ecc, params->f_min,
            params->deltaT, params->inc, &hparams, All);
        if (status == CEV_SUCCESS)
        {
            convert_SEOBCoreOutputs_to_pyOutputStruct_t(params->is_only22, params->m1 + params->m2, params->distance, params->inc, params->phiRef, params->beta, All, output);
            if (params->ret_dyn)
                convert_SEOBCoreOutputs_to_pyDynOutputStruct_t(All, dynoutput);
        }
        else
            is_failed = 1;
        STRUCTFREE(All, SEOBCoreOutputs);
    } else if ( (!params->is_constp) && (sqrt(params->s1x*params->s1x + params->s1y*params->s1y) < 1e-5 && sqrt(params->s2x*params->s2x + params->s2y*params->s2y) < 1e-5)) {
        /* default */
        PRINT_LOG_INFO(LOG_DEBUG, "default SA code");
        REAL8Vector *tVec = NULL;
        SphHarmListCAmpPhaseSequence *hLM = NULL;
        SEOBSAdynamics *dyn_debug = NULL;
        status = evolve_SA(params->m1, params->m2, params->s1z, params->s2z, params->ecc, params->f_min, params->deltaT, &hparams, &tVec, &hLM, params->is_noringdown, params->ret_dyn, &dyn_debug);
        if (status == CEV_SUCCESS)
        {
            convert_SphHarmListCAmpPhaseSequence_to_pyOutputStruct_t(params->is_only22, params->m1 + params->m2, params->distance, params->inc, params->phiRef, params->beta, tVec, hLM, dyn_debug, output);
            if (params->ret_dyn)
                convert_SEOBSAdynamics_to_pyDynOutputStruct_t(dyn_debug, params->m1, params->m2, params->s1z, params->s2z, dynoutput);
        }
        else
            is_failed = 1;
        STRUCTFREE(tVec, REAL8Vector);
        STRUCTFREE(hLM, SphHarmListCAmpPhaseSequence);
        STRUCTFREE(dyn_debug, SEOBSAdynamics);
    } else {
        if (params->is_constp)
        {
            REAL8TimeSeries *hplus = NULL;
            REAL8TimeSeries *hcross = NULL;
            SEOBCoreOutputs *All = (SEOBCoreOutputs *)MYCalloc(1, sizeof(SEOBCoreOutputs));
            status = evolve(params->m1, params->m2, 
                params->s1x, params->s1y, params->s1z, 
                params->s2x, params->s2y, params->s2z, 
                params->beta, params->distance, params->ecc, params->f_min, 
                params->deltaT, params->inc, params->is_only22,
                &hparams, &hplus, &hcross, All);
            if (status == CEV_SUCCESS)
            {
                convert_SEOBCoreOutputs_to_pyOutputStruct_t(params->is_only22, params->m1 + params->m2, params->distance, params->inc, params->phiRef, params->beta, All, output);
                if (params->ret_dyn)
                    convert_SEOBCoreOutputs_to_pyDynOutputStruct_t(All, dynoutput);
            }
            else
                is_failed = 1;
            STRUCTFREE(All, SEOBCoreOutputs);
            STRUCTFREE(hplus, REAL8TimeSeries);
            STRUCTFREE(hcross, REAL8TimeSeries);
        } else {
            REAL8TimeSeries *hplus = NULL;
            REAL8TimeSeries *hcross = NULL;
            SEOBCoreOutputs *All = (SEOBCoreOutputs *)MYCalloc(1, sizeof(SEOBCoreOutputs));
            status = evolve_adaptive(params->m1, params->m2, 
                params->s1x, params->s1y, params->s1z, 
                params->s2x, params->s2y, params->s2z, 
                params->beta, params->distance, params->ecc, params->f_min, 
                params->deltaT, params->inc, params->is_only22,
                &hparams, &hplus, &hcross, All);
            if (status == CEV_SUCCESS)
            {
                convert_SEOBCoreOutputs_to_pyOutputStruct_t(params->is_only22, params->m1 + params->m2, params->distance, params->inc, params->phiRef, params->beta, All, output);
                if (params->ret_dyn)
                    convert_SEOBCoreOutputs_to_pyDynOutputStruct_t(All, dynoutput);
            }
            else
                is_failed = 1;
            STRUCTFREE(All, SEOBCoreOutputs);
            STRUCTFREE(hplus, REAL8TimeSeries);
            STRUCTFREE(hcross, REAL8TimeSeries);
        }
    }
    tend = clock();
    PRINT_LOG_INFO(LOG_INFO, "Time Cost: %fs\n", ((REAL8)(tend-tstart)/CLOCKS_PER_SEC));  
    if (is_failed) return CEV_FAILURE;
    return CEV_SUCCESS;
}


static INT find_exact_amp_peak(REAL8Vector *tMVec, REAL8Vector *amp22)
{
    REAL8 tpeak22, amp;
    INT i, ii, length = tMVec->length;
    // gsl_spline *spline = NULL;
    // gsl_interp_accel *acc = NULL;
    // spline = gsl_spline_alloc (gsl_interp_cspline, length);
    // acc = gsl_interp_accel_alloc ();
    // gsl_spline_init(spline, tMVec->data, amp22->data, length);
    amp = amp22->data[0];
    ii = 0;
    for(i=1; i<length; i++)
    {
        if (amp < amp22->data[i])
        {
            ii = i;
            amp = amp22->data[i];
        }
    }
    // gsl_spline_free(spline);
    // gsl_interp_accel_free(acc);
    return ii;
}

static void apply_phic_on_hpc(REAL8Vector *tMVec, REAL8Vector *hplus, REAL8Vector *hcross, INT ipeak22, REAL8 phic)
{
    REAL8 hp, hc, phi0, dphi;
    hp = hplus->data[ipeak22];
    hc = hcross->data[ipeak22];
    phi0 = atan2(-hc, hp);
    dphi = phic - phi0;
    for(int i=0; i<hplus->length; i++)
    {
        hp = hplus->data[i];
        hc = hcross->data[i];
        hplus->data[i] = hp * cos(dphi) + hc * sin(dphi);
        hcross->data[i] = hc * cos(dphi) - hp * sin(dphi);
    }
    return;
}

void convert_SEOBPrecCoreOutputs_to_pyDynOutputStruct_t(SEOBPrecCoreOutputs *All_prec, pyDynOutputStruct_t **ret)
{
    INT i, length;
    length = All_prec->dyn->length;
    pyDynOutputStruct_t *output = CreatepyDynOutputStruct_t(length);
    memcpy(output->timeM->data, All_prec->dyn->tVec, length*sizeof(REAL8));
    memcpy(output->xVec->data, All_prec->dyn->posVecx, length*sizeof(REAL8));
    memcpy(output->yVec->data, All_prec->dyn->posVecy, length*sizeof(REAL8));
    memcpy(output->zVec->data, All_prec->dyn->posVecz, length*sizeof(REAL8));
    memcpy(output->pTxVec->data, All_prec->dyn->momTVecx, length*sizeof(REAL8));
    memcpy(output->pTyVec->data, All_prec->dyn->momTVecy, length*sizeof(REAL8));
    memcpy(output->pTzVec->data, All_prec->dyn->momTVecz, length*sizeof(REAL8));
    memcpy(output->vxVec->data, All_prec->dyn->velVecx, length*sizeof(REAL8));
    memcpy(output->vyVec->data, All_prec->dyn->velVecy, length*sizeof(REAL8));
    memcpy(output->vzVec->data, All_prec->dyn->velVecz, length*sizeof(REAL8));
    memcpy(output->s1xVec->data, All_prec->dyn->s1Vecx, length*sizeof(REAL8));
    memcpy(output->s1yVec->data, All_prec->dyn->s1Vecy, length*sizeof(REAL8));
    memcpy(output->s1zVec->data, All_prec->dyn->s1Vecz, length*sizeof(REAL8));
    memcpy(output->s2xVec->data, All_prec->dyn->s2Vecx, length*sizeof(REAL8));
    memcpy(output->s2yVec->data, All_prec->dyn->s2Vecy, length*sizeof(REAL8));
    memcpy(output->s2zVec->data, All_prec->dyn->s2Vecz, length*sizeof(REAL8));
    memcpy(output->prTDotVec->data, All_prec->dyn->prTDotVec, length*sizeof(REAL8));
    memcpy(output->hamVec->data, All_prec->dyn->HamVec, length*sizeof(REAL8));
    *ret = output;
    return;
}

void convert_SEOBPrecCoreOutputs_to_pyOutputStruct_t(INT is_only22, REAL8 mtot, REAL8 dL, REAL8 inc, REAL8 phic, REAL8 beta, SEOBPrecCoreOutputs *All_prec, pyOutputStruct_t **ret)
{
    INT i;
    REAL8 mT, amp0;
    mT = mtot * CST_MTSUN_SI;
    amp0 = mtot * CST_MRSUN_SI / dL / 1e6 / CST_PC_SI;
    COMPLEX16TimeSeries *h22 = XLALSphHarmTimeSeriesGetMode(All_prec->hLM, 2, 2);
    COMPLEX16TimeSeries *h21 = XLALSphHarmTimeSeriesGetMode(All_prec->hLM, 2, 1);
    COMPLEX16TimeSeries *h33 = XLALSphHarmTimeSeriesGetMode(All_prec->hLM, 3, 3);
    COMPLEX16TimeSeries *h44 = XLALSphHarmTimeSeriesGetMode(All_prec->hLM, 4, 4);
    COMPLEX16TimeSeries *h55 = XLALSphHarmTimeSeriesGetMode(All_prec->hLM, 5, 5);
    REAL8 deltaT = h22->deltaT;
    INT length = h22->data->length;
    pyOutputStruct_t *output = CreatePyOutputStruct_t(length);
    REAL8Vector *amp22 = CreateREAL8Vector(length);
    for(i=0; i<length; i++)
    {
        output->timeM->data[i] = deltaT*i;
        output->time->data[i] = output->timeM->data[i] * mT;
        output->h22_real->data[i] = creal(h22->data->data[i]);
        output->h22_imag->data[i] = cimag(h22->data->data[i]);
        output->h21_real->data[i] = creal(h21->data->data[i]);
        output->h21_imag->data[i] = cimag(h21->data->data[i]);
        output->h33_real->data[i] = creal(h33->data->data[i]);
        output->h33_imag->data[i] = cimag(h33->data->data[i]),
        output->h44_real->data[i] = creal(h44->data->data[i]);
        output->h44_imag->data[i] = cimag(h44->data->data[i]);
        output->h55_real->data[i] = creal(h55->data->data[i]);
        output->h55_imag->data[i] = cimag(h55->data->data[i]);
        amp22->data[i] = cabs(h22->data->data[i]);
        output->hplus->data[i] = 0.0;
        output->hcross->data[i] = 0.0;
#if 0
        COMPLEX16 sYlm, hpc_contrib;
        // 2, 2
        SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 2, 2, &sYlm);
        hpc_contrib = sYlm * (output->h22_real->data[i] + I*output->h22_imag->data[i]);
        output->hplus->data[i] += amp0 * creal(hpc_contrib);
        output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
        // 2, -2
        SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 2, -2, &sYlm);
        hpc_contrib = sYlm * (output->h22_real->data[i] - I*output->h22_imag->data[i]);
        output->hplus->data[i] += amp0 * creal(hpc_contrib);
        output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
        if (!is_only22)
        {
            // 2, 1
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 2, 1, &sYlm);
            hpc_contrib = sYlm * (output->h21_real->data[i] + I*output->h21_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            // 2, -1
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 2, -1, &sYlm);
            hpc_contrib = sYlm * (output->h21_real->data[i] - I*output->h21_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            // 3, 3
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 3, 3, &sYlm);
            hpc_contrib = sYlm * (output->h33_real->data[i] + I*output->h33_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            // 3, -3
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 3, -3, &sYlm);
            hpc_contrib = -sYlm * (output->h33_real->data[i] - I*output->h33_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            // 4, 4
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 4, 4, &sYlm);
            hpc_contrib = sYlm * (output->h44_real->data[i] + I*output->h44_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            // 4, -4
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 4, -4, &sYlm);
            hpc_contrib = sYlm * (output->h44_real->data[i] - I*output->h44_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            // 5, 5
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 5, 5, &sYlm);
            hpc_contrib = sYlm * (output->h55_real->data[i] + I*output->h55_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            // 5, -5
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 5, -5, &sYlm);
            hpc_contrib = -sYlm * (output->h55_real->data[i] - I*output->h55_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
        }
#endif
    }
    INT l, m;
    COMPLEX16 hpc_contrib, sYlm;
    for (l = 2; l <= 5; l++) 
    {
        for (m = -l; m <= l; m++) 
        {
            if (is_only22 && (l != 2 || abs(m) != 2))
                continue;
            SpinWeightedSphericalHarmonic(inc, CST_PI / 2. - beta, -2, l, m, &sYlm);
            COMPLEX16TimeSeries *hIlm = XLALSphHarmTimeSeriesGetMode(All_prec->hLM, l, m);
            for (i = 0; i < length; i++) 
            {
                hpc_contrib = sYlm * hIlm->data->data[i];
                output->hplus->data[i] += amp0 * creal(hpc_contrib);
                output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            }
        }
    }
    INT ipeak22 = find_exact_amp_peak(output->timeM, amp22);
    apply_phic_on_hpc(output->timeM, output->hplus, output->hcross, ipeak22, phic);
    STRUCTFREE(amp22, REAL8Vector);
    *ret = output;
    return;
}

void convert_SEOBCoreOutputs_to_pyDynOutputStruct_t(SEOBCoreOutputs *All, pyDynOutputStruct_t **ret)
{
    INT i, length;
    length = All->dyn->length;
    pyDynOutputStruct_t *output = CreatepyDynOutputStruct_t(length);
    memcpy(output->timeM->data, All->dyn->tVec, length*sizeof(REAL8));
    memcpy(output->xVec->data, All->dyn->posVecx, length*sizeof(REAL8));
    memcpy(output->yVec->data, All->dyn->posVecy, length*sizeof(REAL8));
    memcpy(output->zVec->data, All->dyn->posVecz, length*sizeof(REAL8));
    memcpy(output->pTxVec->data, All->dyn->momVecx, length*sizeof(REAL8));
    memcpy(output->pTyVec->data, All->dyn->momVecy, length*sizeof(REAL8));
    memcpy(output->pTzVec->data, All->dyn->momVecz, length*sizeof(REAL8));
    memcpy(output->vxVec->data, All->dyn->velVecx, length*sizeof(REAL8));
    memcpy(output->vyVec->data, All->dyn->velVecy, length*sizeof(REAL8));
    memcpy(output->vzVec->data, All->dyn->velVecz, length*sizeof(REAL8));
    memcpy(output->s1xVec->data, All->dyn->s1Vecx, length*sizeof(REAL8));
    memcpy(output->s1yVec->data, All->dyn->s1Vecy, length*sizeof(REAL8));
    memcpy(output->s1zVec->data, All->dyn->s1Vecz, length*sizeof(REAL8));
    memcpy(output->s2xVec->data, All->dyn->s2Vecx, length*sizeof(REAL8));
    memcpy(output->s2yVec->data, All->dyn->s2Vecy, length*sizeof(REAL8));
    memcpy(output->s2zVec->data, All->dyn->s2Vecz, length*sizeof(REAL8));
    memcpy(output->prTDotVec->data, All->dyn->polarprDotVec, length*sizeof(REAL8));
    memcpy(output->hamVec->data, All->dyn->hamVec, length*sizeof(REAL8));
    *ret = output;
    return;
}

void convert_SEOBCoreOutputs_to_pyOutputStruct_t(INT is_only22, REAL8 mtot, REAL8 dL, REAL8 inc, REAL8 phic, REAL8 beta, SEOBCoreOutputs *All, pyOutputStruct_t **ret)
{
    INT i;
    REAL8 mT, amp0;
    mT = mtot * CST_MTSUN_SI;
    amp0 = mtot * CST_MRSUN_SI / dL / 1e6 / CST_PC_SI;
    COMPLEX16TimeSeries *h22 = XLALSphHarmTimeSeriesGetMode(All->hLM, 2, 2);
    COMPLEX16TimeSeries *h21 = XLALSphHarmTimeSeriesGetMode(All->hLM, 2, 1);
    COMPLEX16TimeSeries *h33 = XLALSphHarmTimeSeriesGetMode(All->hLM, 3, 3);
    COMPLEX16TimeSeries *h44 = XLALSphHarmTimeSeriesGetMode(All->hLM, 4, 4);
    COMPLEX16TimeSeries *h55 = XLALSphHarmTimeSeriesGetMode(All->hLM, 5, 5);

    REAL8 deltaT = h22->deltaT;
    INT length = h22->data->length;
    REAL8Vector *amp22 = CreateREAL8Vector(length);
    pyOutputStruct_t *output = CreatePyOutputStruct_t(length);
    // print_debug("inc = %.16e, phi = %.16e\n", inc, beta);
    for(i=0;i<length;i++)
    {
        output->timeM->data[i] = deltaT*i;
        output->time->data[i] = output->timeM->data[i] * mT;
        output->h22_real->data[i] = creal(h22->data->data[i]);
        output->h22_imag->data[i] = cimag(h22->data->data[i]);
        output->h21_real->data[i] = creal(h21->data->data[i]);
        output->h21_imag->data[i] = cimag(h21->data->data[i]);
        output->h33_real->data[i] = creal(h33->data->data[i]);
        output->h33_imag->data[i] = cimag(h33->data->data[i]),
        output->h44_real->data[i] = creal(h44->data->data[i]);
        output->h44_imag->data[i] = cimag(h44->data->data[i]);
        output->h55_real->data[i] = creal(h55->data->data[i]);
        output->h55_imag->data[i] = cimag(h55->data->data[i]);
        amp22->data[i] = cabs(h22->data->data[i]);
        output->hplus->data[i] = 0.0;
        output->hcross->data[i] = 0.0;
#if 0
        COMPLEX16 sYlm, hpc_contrib;
        // 2, 2
        SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 2, 2, &sYlm);
        hpc_contrib = sYlm * (output->h22_real->data[i] + I*output->h22_imag->data[i]);
        output->hplus->data[i] += amp0 * creal(hpc_contrib);
        output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
        // 2, -2
        SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 2, -2, &sYlm);
        hpc_contrib = sYlm * (output->h22_real->data[i] - I*output->h22_imag->data[i]);
        output->hplus->data[i] += amp0 * creal(hpc_contrib);
        output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
        if (!is_only22) 
        {
            // 2, 1
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 2, 1, &sYlm);
            hpc_contrib = sYlm * (output->h21_real->data[i] + I*output->h21_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            // 2, -1
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 2, -1, &sYlm);
            hpc_contrib = sYlm * (output->h21_real->data[i] - I*output->h21_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            // 3, 3
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 3, 3, &sYlm);
            hpc_contrib = sYlm * (output->h33_real->data[i] + I*output->h33_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            // 3, -3
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 3, -3, &sYlm);
            hpc_contrib = -sYlm * (output->h33_real->data[i] - I*output->h33_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            // 4, 4
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 4, 4, &sYlm);
            hpc_contrib = sYlm * (output->h44_real->data[i] + I*output->h44_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            // 4, -4
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 4, -4, &sYlm);
            hpc_contrib = sYlm * (output->h44_real->data[i] - I*output->h44_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            // 5, 5
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 5, 5, &sYlm);
            hpc_contrib = sYlm * (output->h55_real->data[i] + I*output->h55_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            // 5, -5
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 5, -5, &sYlm);
            hpc_contrib = -sYlm * (output->h55_real->data[i] - I*output->h55_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
        }
#endif
    }
    
    INT l, m;
    COMPLEX16 hpc_contrib, sYlm;
    for (l = 2; l <= 5; l++) 
    {
        for (m = -l; m <= l; m++) 
        {
            if (is_only22 && (l != 2 || abs(m) != 2))
                continue;
            SpinWeightedSphericalHarmonic(inc, CST_PI / 2. - beta, -2, l, m, &sYlm);
            COMPLEX16TimeSeries *hIlm = XLALSphHarmTimeSeriesGetMode(All->hLM, l, m);
            for (i = 0; i < length; i++) 
            {
                hpc_contrib = sYlm * hIlm->data->data[i];
                output->hplus->data[i] += amp0 * creal(hpc_contrib);
                output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            }
        }
    }
    INT ipeak22 = find_exact_amp_peak(output->timeM, amp22);
    apply_phic_on_hpc(output->timeM, output->hplus, output->hcross, ipeak22, phic);
    STRUCTFREE(amp22, REAL8Vector);
    *ret = output;
    return;
}

void convert_SEOBSAdynamics_to_pyDynOutputStruct_t(SEOBSAdynamics *dyn_debug, REAL8 m1, REAL8 m2, REAL8 chi1, REAL8 chi2, pyDynOutputStruct_t **ret)
{
    INT i, length;
    length = dyn_debug->length;
    REAL8 r, phi, vr, vf, prT, pphi;
    pyDynOutputStruct_t *output = CreatepyDynOutputStruct_t(length);
    REAL8 mtot = m1 + m2;
    REAL8 s1NormFac = m1*m1/mtot/mtot;
    REAL8 s2NormFac = m2*m2/mtot/mtot;
    memset(output->zVec->data, 0, length*sizeof(REAL8));
    memset(output->pTzVec->data, 0, length*sizeof(REAL8));
    memset(output->vzVec->data, 0, length*sizeof(REAL8));
    memset(output->s1xVec->data, 0, length*sizeof(REAL8));
    memset(output->s1yVec->data, 0, length*sizeof(REAL8));
    memset(output->s2xVec->data, 0, length*sizeof(REAL8));
    memset(output->s2yVec->data, 0, length*sizeof(REAL8));
    for (i=0; i<length; i++)
    {
        r = dyn_debug->rVec[i];
        phi = dyn_debug->phiVec[i];
        vr = dyn_debug->drVec[i];
        vf = r * dyn_debug->dphiVec[i];
        prT = dyn_debug->prTVec[i];
        pphi = dyn_debug->pphiVec[i];
        // print_debug("[%d](r, phi, vr, vf, prT, pphi) = (%f, %f, %f, %f, %f)\n",
        //     i, r, phi, vr, vf, prT, pphi);
        output->timeM->data[i] = dyn_debug->tVec[i];
        output->xVec->data[i] = r * cos(phi);
        output->yVec->data[i] = r * sin(phi);
        // output->zVec->data[i] = 0.0;
        output->pTxVec->data[i] = prT * cos(phi) - pphi * sin(phi)/r;
        output->pTyVec->data[i] = prT * sin(phi) + pphi * cos(phi)/r;
        // output->pTzVec->data[i] = 0.0;
        output->vxVec->data[i] = vr * cos(phi) - vf * sin(phi);
        output->vyVec->data[i] = vr * sin(phi) + vf * cos(phi);
        // output->vzVec->data[i] = 0.0;
        output->s1zVec->data[i] = chi1 * s1NormFac;
        output->s2zVec->data[i] = chi2 * s2NormFac;
        output->prTDotVec->data[i] = dyn_debug->dprTVec[i];
        output->hamVec->data[i] = dyn_debug->HVec[i];
    }
    *ret = output;
    return;
}

void convert_PrecSphHarmListCAmpPhaseSequence_to_pyOutputStruct_t(INT is_only22, REAL8 mtot, REAL8 dL, REAL8 inc, REAL8 phic, REAL8 beta, REAL8Vector *tVec, SphHarmListCAmpPhaseSequence *PLM, pyOutputStruct_t **ret)
{
    INT i;
    REAL8 mT, amp0;
    mT = mtot * CST_MTSUN_SI;
    amp0 = mtot * CST_MRSUN_SI / dL / 1e6 / CST_PC_SI;
    CAmpPhaseSequence *h22 = SphHarmListCAmpPhaseSequence_GetMode(PLM, 2, 2)->campphase;
    CAmpPhaseSequence *h21 = SphHarmListCAmpPhaseSequence_GetMode(PLM, 2, 1)->campphase;
    CAmpPhaseSequence *h33 = SphHarmListCAmpPhaseSequence_GetMode(PLM, 3, 3)->campphase;
    CAmpPhaseSequence *h44 = SphHarmListCAmpPhaseSequence_GetMode(PLM, 4, 4)->campphase;
    CAmpPhaseSequence *h55 = SphHarmListCAmpPhaseSequence_GetMode(PLM, 5, 5)->campphase;
    // REAL8Vector *tVec = h22->xdata;
    INT length = tVec->length;
    pyOutputStruct_t *output = CreatePyOutputStruct_t(length);
    REAL8Vector *amp22 = CreateREAL8Vector(length);
    for(i=0;i<length;i++)
    {
        output->timeM->data[i] = tVec->data[i];
        output->time->data[i] = output->timeM->data[i] * mT;
        output->h22_real->data[i] = h22->camp_real->data[i]*cos(h22->phase->data[i]);
        output->h22_imag->data[i] = h22->camp_real->data[i]*sin(h22->phase->data[i]);
        output->h21_real->data[i] = h21->camp_real->data[i]*cos(h21->phase->data[i]);
        output->h21_imag->data[i] = h21->camp_real->data[i]*sin(h21->phase->data[i]);
        output->h33_real->data[i] = h33->camp_real->data[i]*cos(h33->phase->data[i]);
        output->h33_imag->data[i] = h33->camp_real->data[i]*sin(h33->phase->data[i]);
        output->h44_real->data[i] = h44->camp_real->data[i]*cos(h44->phase->data[i]);
        output->h44_imag->data[i] = h44->camp_real->data[i]*sin(h44->phase->data[i]);
        output->h55_real->data[i] = h55->camp_real->data[i]*cos(h55->phase->data[i]);
        output->h55_imag->data[i] = h55->camp_real->data[i]*sin(h55->phase->data[i]);
        output->hplus->data[i] = 0.0;
        output->hcross->data[i] = 0.0;
        amp22->data[i] = h22->camp_real->data[i];
    }

    INT l, m;
    COMPLEX16 hpc_contrib, sYlm;
    for (l = 2; l <= 5; l++) 
    {
        for (m = -l; m <= l; m++) 
        {
            if (is_only22 && (l != 2 || abs(m) != 2))
                continue;
            SpinWeightedSphericalHarmonic(inc, CST_PI / 2. - beta, -2, l, m, &sYlm);
            // COMPLEX16TimeSeries *hIlm = XLALSphHarmTimeSeriesGetMode(All->hLM, l, m);
            CAmpPhaseSequence *hIlm = SphHarmListCAmpPhaseSequence_GetMode(PLM, l, m)->campphase;
            for (i = 0; i < length; i++) 
            {
                hpc_contrib = sYlm * (hIlm->camp_real->data[i]*cos(hIlm->phase->data[i]) + I*hIlm->camp_real->data[i]*sin(hIlm->phase->data[i]));
                output->hplus->data[i] += amp0 * creal(hpc_contrib);
                output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            }
        }
    }
    INT ipeak22 = find_exact_amp_peak(output->timeM, amp22);
    apply_phic_on_hpc(output->timeM, output->hplus, output->hcross, ipeak22, phic);
    STRUCTFREE(amp22, REAL8Vector);
    *ret = output;
    return;
}

void convert_SphHarmListCAmpPhaseSequence_to_pyOutputStruct_t(INT is_only22, REAL8 mtot, REAL8 dL, REAL8 inc, REAL8 phic, REAL8 beta, REAL8Vector *tVec, SphHarmListCAmpPhaseSequence *hLM, SEOBSAdynamics *dyn_debug, pyOutputStruct_t **ret)
{
    INT i;
    REAL8 mT, amp0;
    mT = mtot * CST_MTSUN_SI;
    amp0 = mtot * CST_MRSUN_SI / dL / 1e6 / CST_PC_SI;
    CAmpPhaseSequence *h22 = SphHarmListCAmpPhaseSequence_GetMode(hLM, 2, 2)->campphase;
    CAmpPhaseSequence *h21 = SphHarmListCAmpPhaseSequence_GetMode(hLM, 2, 1)->campphase;
    CAmpPhaseSequence *h33 = SphHarmListCAmpPhaseSequence_GetMode(hLM, 3, 3)->campphase;
    CAmpPhaseSequence *h44 = SphHarmListCAmpPhaseSequence_GetMode(hLM, 4, 4)->campphase;
    CAmpPhaseSequence *h55 = SphHarmListCAmpPhaseSequence_GetMode(hLM, 5, 5)->campphase;
    INT length = tVec->length;
    pyOutputStruct_t *output = CreatePyOutputStruct_t(length);
    REAL8Vector *amp22 = CreateREAL8Vector(length);
    for(i=0;i<length;i++)
    {
        output->timeM->data[i] = tVec->data[i];
        output->time->data[i] = output->timeM->data[i] * mT;
        output->h22_real->data[i] = h22->camp_real->data[i]*cos(h22->phase->data[i]);
        output->h22_imag->data[i] = h22->camp_real->data[i]*sin(h22->phase->data[i]);
        output->h21_real->data[i] = h21->camp_real->data[i]*cos(h21->phase->data[i]);
        output->h21_imag->data[i] = h21->camp_real->data[i]*sin(h21->phase->data[i]);
        output->h33_real->data[i] = h33->camp_real->data[i]*cos(h33->phase->data[i]);
        output->h33_imag->data[i] = h33->camp_real->data[i]*sin(h33->phase->data[i]);
        output->h44_real->data[i] = h44->camp_real->data[i]*cos(h44->phase->data[i]);
        output->h44_imag->data[i] = h44->camp_real->data[i]*sin(h44->phase->data[i]);
        output->h55_real->data[i] = h55->camp_real->data[i]*cos(h55->phase->data[i]);
        output->h55_imag->data[i] = h55->camp_real->data[i]*sin(h55->phase->data[i]);
        output->hplus->data[i] = 0.0;
        output->hcross->data[i] = 0.0;
        amp22->data[i] = h22->camp_real->data[i];
#if 0
        COMPLEX16 sYlm, hpc_contrib;
        // 2, 2
        SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 2, 2, &sYlm);
        hpc_contrib = sYlm * (output->h22_real->data[i] + I*output->h22_imag->data[i]);
        output->hplus->data[i] += amp0 * creal(hpc_contrib);
        output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
        // 2, -2
        SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 2, -2, &sYlm);
        hpc_contrib = sYlm * (output->h22_real->data[i] - I*output->h22_imag->data[i]);
        output->hplus->data[i] += amp0 * creal(hpc_contrib);
        output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
        if (!is_only22) 
        {
            // 2, 1
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 2, 1, &sYlm);
            hpc_contrib = sYlm * (output->h21_real->data[i] + I*output->h21_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            // 2, -1
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 2, -1, &sYlm);
            hpc_contrib = sYlm * (output->h21_real->data[i] - I*output->h21_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            // 3, 3
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 3, 3, &sYlm);
            hpc_contrib = sYlm * (output->h33_real->data[i] + I*output->h33_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            // 3, -3
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 3, -3, &sYlm);
            hpc_contrib = -sYlm * (output->h33_real->data[i] - I*output->h33_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            // 4, 4
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 4, 4, &sYlm);
            hpc_contrib = sYlm * (output->h44_real->data[i] + I*output->h44_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            // 4, -4
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 4, -4, &sYlm);
            hpc_contrib = sYlm * (output->h44_real->data[i] - I*output->h44_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            // 5, 5
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 5, 5, &sYlm);
            hpc_contrib = sYlm * (output->h55_real->data[i] + I*output->h55_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            // 5, -5
            SpinWeightedSphericalHarmonic(inc, CST_PI/2. - beta, -2, 5, -5, &sYlm);
            hpc_contrib = -sYlm * (output->h55_real->data[i] - I*output->h55_imag->data[i]);
            output->hplus->data[i] += amp0 * creal(hpc_contrib);
            output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
        }
#endif
    }

    INT l, m;
    COMPLEX16 hpc_contrib, sYlm;
    for (l = 2; l <= 5; l++) 
    {
        for (m = -l; m <= l; m++) 
        {
            if (is_only22 && (l != 2 || abs(m) != 2))
                continue;
            SpinWeightedSphericalHarmonic(inc, CST_PI / 2. - beta, -2, l, m, &sYlm);
            // COMPLEX16TimeSeries *hIlm = XLALSphHarmTimeSeriesGetMode(All->hLM, l, m);
            CAmpPhaseSequence *hIlm = SphHarmListCAmpPhaseSequence_GetMode(hLM, l, m)->campphase;
            for (i = 0; i < length; i++) 
            {
                hpc_contrib = sYlm * (hIlm->camp_real->data[i]*cos(hIlm->phase->data[i]) + I*hIlm->camp_real->data[i]*sin(hIlm->phase->data[i]));
                output->hplus->data[i] += amp0 * creal(hpc_contrib);
                output->hcross->data[i] += -amp0 * cimag(hpc_contrib);
            }
        }
    }
    INT ipeak22 = find_exact_amp_peak(output->timeM, amp22);
    apply_phic_on_hpc(output->timeM, output->hplus, output->hcross, ipeak22, phic);
    STRUCTFREE(amp22, REAL8Vector);
    *ret = output;
    return;
}


// Create Spin Params
NewtonMultipolePrefixes *test_interface(REAL8 m1, REAL8 m2)
{
    NewtonMultipolePrefixes *prefixes = (NewtonMultipolePrefixes *) MYCalloc(1, sizeof(NewtonMultipolePrefixes));
    if (XLALSimIMREOBComputeNewtonMultipolePrefixes(prefixes, m1, m2) != CEV_SUCCESS)
    {
        return NULL;
    }
    INT l, m;
    for(l=2; l<8; l++)
        for (m=1; m<=l; m++)
            print_debug("prefixes[%d][%d] = %.5e + I %.5e\n", 
                l, m,
                creal(prefixes->values[l][m]), 
                cimag(prefixes->values[l][m]));
    prefixes->values[0][0] = 1.0;
    return prefixes;
}