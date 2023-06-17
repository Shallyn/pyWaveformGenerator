/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "pUtils.h"
#include "myLog.h"
#include <gsl/gsl_spline.h>
#define GSL_START \
          gsl_error_handler_t *saveGSLErrorHandler_; \
          saveGSLErrorHandler_ = gsl_set_error_handler_off();

#define GSL_END \
          gsl_set_error_handler( saveGSLErrorHandler_ );

/* 
 * 
 * whether use egw 
 * 0: not use
 * 1: egw
 * 
 */
INT g_egw_flag = 0;
INT get_egw_flag()
{
    return g_egw_flag;
}
void set_egw_flag(INT flag)
{
    g_egw_flag = flag;
    return;
}

/* 
 * 
 * version of code 
 * 0: SEOBNRv4PHM
 * 1: SEOBNRP
 * 
 */
INT g_DebugFlag = 0;
INT get_DebugFlag()
{
    return g_DebugFlag;
}
void set_DebugFlag(INT flag)
{
    g_DebugFlag = flag;
    return;
}

INT g_InputDebugFlag = 0;
INT get_InputDebugFlag()
{
    return g_InputDebugFlag;
}
void set_InputDebugFlag(INT flag)
{
    g_InputDebugFlag = flag;
    return;
}


INT g_VersionFlag = 1;
INT get_VersionFlag()
{
    return g_VersionFlag;
}
void set_VersionFlag(INT flag)
{
    g_VersionFlag = flag;
    return;
}

INT g_ConserveFlag = 0;
REAL8 g_ConserveTime = 1000.;
INT get_ConserveFlag()
{
    return g_ConserveFlag;
}
REAL8 get_ConserveTime()
{
    return g_ConserveTime;
}
void set_ConserveFlag(INT flag, REAL8 ctime)
{
    g_ConserveFlag = flag;
    g_ConserveTime = ctime;
    return;
}


REAL8 *rotate_vector(const REAL8 v[], const REAL8 k[], const REAL8 theta,REAL8 result[] )
{
    // Rodrigues rotation formula, see https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    UINT ii;
    REAL8 kcrossv[3]={0,0,0};
    cross_product3d(k,v,kcrossv);
    REAL8 kdotv = inner_product3d(k,v);
    for ( ii=0;ii<3;ii++)
    {
        result[ii] = v[ii]*cos(theta)+kcrossv[ii]*sin(theta)+k[ii]*kdotv*(1-cos(theta));
    }
    return result;
}


REAL8 inner_product3d(const REAL8 values1[], const REAL8 values2[])
{
    REAL8 result = 0;
    INT i;
    for (i=0; i<3; i++)
        result += values1[i] * values2[i];
    return result;
}

REAL8 *cross_product3d(const REAL8 values1[],
                        const REAL8 values2[],
                        REAL8 result[] )
{
    result[0] = values1[1]*values2[2] - values1[2]*values2[1];
    result[1] = values1[2]*values2[0] - values1[0]*values2[2];
    result[2] = values1[0]*values2[1] - values1[1]*values2[0];

    return result;
}

REAL8 SEOBCalculateChiS(REAL8 chi1dotZ, REAL8 chi2dotZ) {
  return 0.5 * (chi1dotZ + chi2dotZ);
}
REAL8 SEOBCalculateChiA(REAL8 chi1dotZ, REAL8 chi2dotZ) {
  return 0.5 * (chi1dotZ - chi2dotZ);
}

INT EOBCalculateSigmaStar(REAL8Vector *sigmaStar,
                          REAL8 mass1,
                          REAL8 mass2,
                          REAL8Vector *s1norm,
                          REAL8Vector *s2norm)

{
    UINT i;
    for (i = 0; i < 3; i++)
    {
        sigmaStar->data[i] = (mass2/mass1 * s1norm->data[i] + mass1/mass2 * s2norm->data[i]);
    }
    return CEV_SUCCESS;
}

INT EOBCalculateSigmaKerr(REAL8Vector *sigmaKerr,
                          REAL8Vector *s1norm,
                          REAL8Vector *s2norm)
{
    UINT i;
    
    for (i=0; i<3; i++)
    {
        sigmaKerr->data[i] = s1norm->data[i] + s2norm->data[i];
    }
    return CEV_SUCCESS;
}

INT XLALSimIMRSpinEOBCalculateSigmaStar(REAL8Vector *sigmaStar,
                          REAL8 mass1,
                          REAL8 mass2,
                          REAL8Vector *s1,
                          REAL8Vector *s2)

{
    UINT i;
    REAL8 mTotal = mass1 + mass2, mT2;
    mT2 = mTotal * mTotal;
    for (i = 0; i < 3; i++)
    {
        sigmaStar->data[i] = (mass2/mass1 * s1->data[i] + mass1/mass2 * s2->data[i])/mT2;
    }
    return CEV_SUCCESS;
}

INT XLALSimIMRSpinEOBCalculateSigmaKerr(REAL8Vector *sigmaKerr,
                                    REAL8 mass1,
                                    REAL8 mass2,
                                    REAL8Vector *s1,
                                    REAL8Vector *s2)
{
    UINT i;
    REAL8 mTotal = mass1 + mass2, mT2;
    mT2 = mTotal * mTotal;
    for (i=0; i<3; i++)
    {
        sigmaKerr->data[i] = (s1->data[i] + s2->data[i])/mT2;
    }
    return CEV_SUCCESS;
}

SpinEOBDynamics *CreateSpinEOBDynamics(UINT length)
{
    SpinEOBDynamics *dy = NULL;
    dy = (SpinEOBDynamics*)MYMalloc(sizeof(*dy));
    dy->t0 = 0.0;
    dy->length = length;
    dy->tVec = 
        dy->xVec = dy->yVec = dy->zVec = 
        dy->pxVec = dy->pyVec = dy->pzVec = NULL;
    dy->tVec = CreateREAL8Vector(length);
    if (!dy->tVec) {STRUCTFREE(dy, SpinEOBDynamics); return NULL;}
    dy->xVec = CreateREAL8Vector(length);
    if (!dy->xVec) {STRUCTFREE(dy, SpinEOBDynamics); return NULL;}
    dy->yVec = CreateREAL8Vector(length);
    if (!dy->yVec) {STRUCTFREE(dy, SpinEOBDynamics); return NULL;}
    dy->zVec = CreateREAL8Vector(length);
    if (!dy->zVec) {STRUCTFREE(dy, SpinEOBDynamics); return NULL;}
    dy->pxVec = CreateREAL8Vector(length);
    if (!dy->pxVec) {STRUCTFREE(dy, SpinEOBDynamics); return NULL;}
    dy->pyVec = CreateREAL8Vector(length);
    if (!dy->pyVec) {STRUCTFREE(dy, SpinEOBDynamics); return NULL;}
    dy->pzVec = CreateREAL8Vector(length);
    if (!dy->pzVec) {STRUCTFREE(dy, SpinEOBDynamics); return NULL;}
    return dy;
}

void DestroySpinEOBDynamics(SpinEOBDynamics *dy)
{
    if (!dy)
        return;
    
    STRUCTFREE(dy->tVec, REAL8Vector);
    STRUCTFREE(dy->xVec, REAL8Vector);
    STRUCTFREE(dy->yVec, REAL8Vector);
    STRUCTFREE(dy->zVec, REAL8Vector);
    STRUCTFREE(dy->pxVec, REAL8Vector);
    STRUCTFREE(dy->pyVec, REAL8Vector);
    STRUCTFREE(dy->pzVec, REAL8Vector);
    dy->length = 0;
    dy->t0 = 0.0;
    MYFree(dy);
    return;
}

void DestroySpinEOBParams(SpinEOBParams *params)
{
    if (!params)
        return;
    
    STRUCTFREE(params->s1Vec, REAL8Vector);
    STRUCTFREE(params->s2Vec, REAL8Vector);
    STRUCTFREE(params->J0Vec, REAL8Vector);
    STRUCTFREE(params->sigmaStar, REAL8Vector);
    STRUCTFREE(params->sigmaKerr, REAL8Vector);
    if (params->seobCoeffs) 
    {MYFree(params->seobCoeffs); params->seobCoeffs = NULL;}
    if (params->seobCoeffConsts) 
    {MYFree(params->seobCoeffConsts); params->seobCoeffConsts = NULL;}
    if (params->hParams) 
    {MYFree(params->hParams); params->hParams = NULL;}
    if (params->prefixes)
    {MYFree(params->prefixes); params->prefixes = NULL;}
    if (params->hCoeffs)
    {MYFree(params->hCoeffs); params->hCoeffs = NULL;}
    if (params->nqcCoeffs)
    {MYFree(params->nqcCoeffs); params->nqcCoeffs = NULL;}
    if (params->eccCoeffs)
    {MYFree(params->eccCoeffs); params->eccCoeffs = NULL;}
    if (params->saCoeffs)
    {MYFree(params->saCoeffs); params->saCoeffs = NULL;}
    MYFree(params);
    return;
}

SEOBdynamics *CreateSEOBdynamics(INT length)
{
    SEOBdynamics *ret = NULL;
    ret = (SEOBdynamics *) MYMalloc(sizeof(SEOBdynamics));
    if (!ret)
        return NULL;
    ret->length = length;
    ret->array = (REAL8Array *) CreateREAL8Array(2, v4PdynamicsVariables, length);
    if (!ret->array)
        return NULL;
    /* Set array pointers corresponding to the data vectors (successive vectors of
    * length retLen in the 1D array data) */
    ret->tVec = ret->array->data;
    ret->posVecx = ret->array->data + 1 * length;
    ret->posVecy = ret->array->data + 2 * length;
    ret->posVecz = ret->array->data + 3 * length;
    ret->momVecx = ret->array->data + 4 * length;
    ret->momVecy = ret->array->data + 5 * length;
    ret->momVecz = ret->array->data + 6 * length;
    ret->s1Vecx = ret->array->data + 7 * length;
    ret->s1Vecy = ret->array->data + 8 * length;
    ret->s1Vecz = ret->array->data + 9 * length;
    ret->s2Vecx = ret->array->data + 10 * length;
    ret->s2Vecy = ret->array->data + 11 * length;
    ret->s2Vecz = ret->array->data + 12 * length;
    ret->phiDMod = ret->array->data + 13 * length;
    ret->phiMod = ret->array->data + 14 * length;
    ret->velVecx = ret->array->data + 15 * length;
    ret->velVecy = ret->array->data + 16 * length;
    ret->velVecz = ret->array->data + 17 * length;
    ret->polarrVec = ret->array->data + 18 * length;
    ret->polarphiVec = ret->array->data + 19 * length;
    ret->polarprVec = ret->array->data + 20 * length;
    ret->polarpphiVec = ret->array->data + 21 * length;
    ret->omegaVec = ret->array->data + 22 * length;
    ret->s1dotZVec = ret->array->data + 23 * length;
    ret->s2dotZVec = ret->array->data + 24 * length;
    ret->hamVec = ret->array->data + 25 * length;

    ret->chiAxVec = ret->array->data + 26 * length;
    ret->chiAyVec = ret->array->data + 27 * length;
    ret->chiAzVec = ret->array->data + 28 * length;
    ret->chiSxVec = ret->array->data + 29 * length;
    ret->chiSyVec = ret->array->data + 30 * length;
    ret->chiSzVec = ret->array->data + 31 * length;

    ret->fluxVec = ret->array->data + 32 * length;
    ret->FrVec = ret->array->data + 33 * length;
    ret->FfVec = ret->array->data + 34 * length;
    ret->polarprDotVec = ret->array->data + 35 * length;

    return ret;
}

void DestroySEOBdynamics(SEOBdynamics *seobdynamics)
{
    STRUCTFREE(seobdynamics->array, REAL8Array);
    MYFree(seobdynamics);
    return;
}

SEOBSAdynamics *CreateSEOBSAdynamics(INT length)
{
    SEOBSAdynamics *ret = NULL;
    ret = (SEOBSAdynamics *) MYMalloc(sizeof(SEOBSAdynamics));
    if (!ret)
        return NULL;
    ret->length = length;
    ret->array = (REAL8Array *) CreateREAL8Array(2, v4SAdynamicsVariables, length);
    if (!ret->array)
        return NULL;
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
    return ret;
}

void DestroySEOBSAdynamics(SEOBSAdynamics *seobdynamics)
{
    STRUCTFREE(seobdynamics->array, REAL8Array);
    MYFree(seobdynamics);
    return;
}

SEOBPrecdynamics *CreateSEOBPrecdynamics(INT length)
{
    SEOBPrecdynamics *ret = NULL;
    ret = (SEOBPrecdynamics *) MYMalloc(sizeof(SEOBPrecdynamics));
    if (!ret)
        return NULL;
    ret->length = length;
    ret->array = (REAL8Array *) CreateREAL8Array(2, v4PrecdynamicsVariables, length);
    if (!ret->array)
        return NULL;
    ret->tVec = ret->array->data;
    ret->posVecx = ret->array->data + 1 * length;
    ret->posVecy = ret->array->data + 2 * length;
    ret->posVecz = ret->array->data + 3 * length;
    ret->momTVecx = ret->array->data + 4 * length;
    ret->momTVecy = ret->array->data + 5 * length;
    ret->momTVecz = ret->array->data + 6 * length;
    ret->s1Vecx = ret->array->data + 7 * length;
    ret->s1Vecy = ret->array->data + 8 * length;
    ret->s1Vecz = ret->array->data + 9 * length;
    ret->s2Vecx = ret->array->data + 10 * length;
    ret->s2Vecy = ret->array->data + 11 * length;
    ret->s2Vecz = ret->array->data + 12 * length;
    ret->phiDMod = ret->array->data + 13 * length;
    ret->phiMod = ret->array->data + 14 * length;
    ret->velVecx = ret->array->data + 15 * length;
    ret->velVecy = ret->array->data + 16 * length;
    ret->velVecz = ret->array->data + 17 * length;
    ret->HamVec = ret->array->data + 18 * length;
    ret->fluxVec = ret->array->data + 19 * length;
    ret->prTDotVec = ret->array->data + 20 * length;

    ret->polarrVec = ret->array->data + 21 * length;
    ret->polarphiVec = ret->array->data + 22 * length;
    ret->polarprTVec = ret->array->data + 23 * length;
    ret->polarpphiVec = ret->array->data + 24 * length;
    ret->omegaVec = ret->array->data + 25 * length;

    ret->nchiaVec = ret->array->data + 26 * length;
    ret->nchisVec = ret->array->data + 27 * length;
    ret->lchiaVec = ret->array->data + 28 * length;
    ret->lchisVec = ret->array->data + 29 * length;
    ret->echiaVec = ret->array->data + 30 * length;
    ret->echisVec = ret->array->data + 31 * length;

    ret->s1dotZVec = ret->array->data + 32 * length;
    ret->s2dotZVec = ret->array->data + 33 * length;

    ret->chi1chi1 = ret->array->data + 34 * length;
    ret->chi1chi2 = ret->array->data + 35 * length;
    ret->chi2chi2 = ret->array->data + 36 * length;

    ret->JnVec = ret->array->data + 37 * length;
    ret->JlVec = ret->array->data + 38 * length;
    ret->JeVec = ret->array->data + 39 * length;

    return ret;
}

void DestroySEOBPrecdynamics(SEOBPrecdynamics *seobdynamics)
{
    STRUCTFREE(seobdynamics->array, REAL8Array);
    MYFree(seobdynamics);
    return;
}


int SphHarmListEOBNonQCCoeffs_AddMode(
    SphHarmListEOBNonQCCoeffs **list_prepended,       /* List structure to prepend to */
    EOBNonQCCoeffs *nqcCoeffs, /* Mode data to be added */
    UINT l,                   /*< Mode number l */
    INT m                     /*< Mode number m */
) 
{
    SphHarmListEOBNonQCCoeffs *list;
    /* Check if the node with this mode already exists */
    list = *list_prepended;
    while (list) 
    {
        if (l == list->l && m == list->m)
            break;
        list = list->next;
    }
    if (list) 
    { /* We don't allow for the case where the mode already exists in
                    the list*/
        PRINT_LOG_INFO(LOG_CRITICAL, "tried to add an already existing mode to a SphHarmListCAmpPhaseSequence.");
        return CEV_FAILURE;
    } 
    else 
    {
        list = MYMalloc(sizeof(SphHarmListEOBNonQCCoeffs));
    }
    list->l = l;
    list->m = m;
    if (nqcCoeffs) 
    {
        list->nqcCoeffs = nqcCoeffs;
    } else 
    {
        list->nqcCoeffs = NULL;
    }
    if (*list_prepended) 
    {
        list->next = *list_prepended;
    } else 
    {
        list->next = NULL;
    }
    *list_prepended = list;

    return CEV_SUCCESS;
}

SphHarmListEOBNonQCCoeffs *SphHarmListEOBNonQCCoeffs_GetMode(
    SphHarmListEOBNonQCCoeffs
        *list, /* List structure to get a particular mode from */
    UINT l,   /*< Mode number l */
    INT m     /*< Mode number m */
) 
{
    if (!list)
        return NULL;

    SphHarmListEOBNonQCCoeffs *itr = list;
    while (itr->l != l || itr->m != m) 
    {
        itr = itr->next;
        if (!itr)
        return NULL;
    }
    /* Return a pointer to a SphHarmListCAmpPhaseSequence */
    return itr;
}

int DestroySphHarmListEOBNonQCCoeffs
(
    SphHarmListEOBNonQCCoeffs *list /* Pointer to list to be destroyed */
)
{
    SphHarmListEOBNonQCCoeffs *pop;
    while ((pop = list))
    {
        if (pop->nqcCoeffs)
            MYFree(pop->nqcCoeffs);
        /* Go to next element */
        list = pop->next;
        /* Free structure itself */
        MYFree(pop);
    }
    return CEV_SUCCESS;
}

int CAmpPhaseSequence_Init(
    CAmpPhaseSequence **campphase, /* Double pointer to campphase sequence */
    int size                       /* Size of data */
) 
{
    /* Check input pointer */
    if (!campphase) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "input double pointer is NULL.");
        return CEV_FAILURE;
    }
    if (*campphase) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "input pointer is not NULL.");
        return CEV_FAILURE;
    }

    /* Allocate structure */
    *campphase = MYMalloc(sizeof(CAmpPhaseSequence));

    /* Allocate storage for data and initialize to 0 */
    if (!((*campphase)->xdata = CreateREAL8Vector(size))) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failed to create REAL8Vector xdata.");
        return CEV_FAILURE;
    }
    if (!((*campphase)->camp_real = CreateREAL8Vector(size))) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failed to create REAL8Vector camp_real.");
        return CEV_FAILURE;
    }
    if (!((*campphase)->camp_imag = CreateREAL8Vector(size))) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failed to create REAL8Vector camp_imag.");
        return CEV_FAILURE;
    }
    if (!((*campphase)->phase = CreateREAL8Vector(size))) 
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "failed to create REAL8Vector phase.");
        return CEV_FAILURE;
    }
    memset((*campphase)->xdata->data, 0, size * sizeof(REAL8));
    memset((*campphase)->camp_real->data, 0, size * sizeof(REAL8));
    memset((*campphase)->camp_imag->data, 0, size * sizeof(REAL8));
    memset((*campphase)->phase->data, 0, size * sizeof(REAL8));

    return CEV_SUCCESS;
}


void DestroyCAmpPhaseSequence(
    CAmpPhaseSequence *campphase /* Pointer to structure to be destroyed */
) 
{
    /* Raise an error if NULL pointer */
    if (!campphase)
        return;
    STRUCTFREE(campphase->xdata, REAL8Vector);
    STRUCTFREE(campphase->camp_real, REAL8Vector);
    STRUCTFREE(campphase->camp_imag, REAL8Vector);
    STRUCTFREE(campphase->phase, REAL8Vector);
    MYFree(campphase);

    return;
}

int SphHarmListCAmpPhaseSequence_AddMode(
    SphHarmListCAmpPhaseSequence *
        *list_prepended,          /* List structure to prepend to */
    CAmpPhaseSequence *campphase, /* Mode data to be added */
    UINT l,                      /*< Mode number l */
    INT  m                        /*< Mode number m */
) 
{
    SphHarmListCAmpPhaseSequence *list;
    /* Check if the node with this mode already exists */
    list = *list_prepended;
    while (list) 
    {
        if (l == list->l && m == list->m)
            break;
        list = list->next;
    }
    if (list) 
    { /* We don't allow for the case where the mode already exists in
                    the list*/
        PRINT_LOG_INFO(LOG_CRITICAL, "tried to add an already existing mode to a SphHarmListCAmpPhaseSequence.");
        return CEV_FAILURE;
    } else {
        list = MYMalloc(sizeof(SphHarmListCAmpPhaseSequence));
    }
    list->l = l;
    list->m = m;
    if (campphase) {
        list->campphase = campphase;
    } else {
        list->campphase = NULL;
    }
    if (*list_prepended) {
        list->next = *list_prepended;
    } else {
        list->next = NULL;
    }
    *list_prepended = list;

    return CEV_SUCCESS;
}

SphHarmListCAmpPhaseSequence *SphHarmListCAmpPhaseSequence_GetMode(
        SphHarmListCAmpPhaseSequence *list, /* List structure to get a particular mode from */
        UINT l,   /*< Mode number l */
        INT m     /*< Mode number m */
) 
{
    if (!list)
        return NULL;

    SphHarmListCAmpPhaseSequence *itr = list;
    while (itr->l != l || itr->m != m) 
    {
        itr = itr->next;
        if (!itr)
        return NULL;
    }
    /* Return a pointer to a SphHarmListCAmpPhaseSequence */
    return itr;
}


int DestroySphHarmListCAmpPhaseSequence(
    SphHarmListCAmpPhaseSequence *list /* Pointer to list to be destroyed */
) 
{
    SphHarmListCAmpPhaseSequence *pop;
    while ((pop = list)) 
    {
        if (pop->campphase) 
        { /* Internal CAmpPhaseSequence is freed */
            STRUCTFREE(pop->campphase, CAmpPhaseSequence);
        }
        /* Go to nsext element */
        list = pop->next;
        /* Free structure itself */
        MYFree(pop);
    }

    return CEV_SUCCESS;
}

/**
 * Set the tdata member for *all* nodes in the list.
 */
void XLALSphHarmTimeSeriesSetTData(
            SphHarmTimeSeries *ts, /**< Linked list to be prepended */
            REAL8Vector* tdata /**< series of time data*/
)
{
	while( ts )
    {
		ts->tdata = tdata;
		ts = ts->next;
	}
}

/**
 * Prepend a node to a linked list of SphHarmTimeSeries, or create a new head
 */
SphHarmTimeSeries* XLALSphHarmTimeSeriesAddMode(
            SphHarmTimeSeries *appended, /**< Linked list to be prepended */
            const COMPLEX16TimeSeries* inmode, /**< Time series of h_lm mode being prepended */
            UINT l, /**< l index of h_lm mode being prepended */
            INT m /**< m index of h_lm mode being prepended */
)
{
    SphHarmTimeSeries* ts;

    // Check if the node with this l, m already exists
    ts = appended;
    while( ts )
    {
        if( l == ts->l && m == ts->m )
        {
            break;
        }
        ts = ts->next;
    }

    if( ts )
    {
        STRUCTFREE( ts->mode , COMPLEX16TimeSeries);
        ts->mode = CutCOMPLEX16TimeSeries( inmode, 0, inmode->data->length);
        return appended;
    } else {
        ts = MYMalloc( sizeof(SphHarmTimeSeries) );
    }

    ts->l = l;
    ts->m = m;
    // Cut returns a new series using a slice of the original. I ask it to
    // return a new one for the full data length --- essentially a duplication
    if( inmode )
    {
        ts->mode = CutCOMPLEX16TimeSeries( inmode, 0, inmode->data->length);
    } else 
    {
        ts->mode = NULL;
    }

    if( appended )
    {
        ts->next = appended;
        ts->tdata = appended->tdata;
    } else {
        ts->next = NULL;
        ts->tdata = NULL;
    }

    return ts;
}

/**
 * Get the time series of a waveform's (l,m) spherical harmonic mode from a
 * SphHarmTimeSeries linked list. Returns a pointer to its COMPLEX16TimeSeries
 */
COMPLEX16TimeSeries* XLALSphHarmTimeSeriesGetMode(
            SphHarmTimeSeries *ts, /**< linked list to extract mode from */
            UINT l, /**< l index of h_lm mode to get */
            INT m /**< m index of h_lm mode to get */
)
{
    if( !ts ) return NULL;

    SphHarmTimeSeries *itr = ts;
    while( itr->l != l || itr->m != m )
    {
        itr = itr->next;
        if( !itr ) return NULL;
    }
    return itr->mode;
}

/**
 * Convention for Wigner matrices (mp stands for m', * for conjugation):
 * for the active rotation from the I-frame to the P-frame, parameterized by
 * Euler angles (alpha, beta, gamma) in the ZYZ convention \f[ h^{P}_{lm} =
 * \sum_{m'} D^{l}_{m m'} h^{I}_{lm'}\f] \f[ h^{I}_{lm} = \sum_{m'} D^{l}_{m
 * m'}* h^{P}_{lm'}\f] \f[ D^{l}_{m m'} = d^{l}_{m m'}(\beta) \exp{i m \alpha}
 * \exp{i m' \gamma}\f] with the notation \f$ c,s = \cos, \sin (\beta/2)\f$, \f$
 * k_{min} = \max(0,m-m'), k_{max} = \min(l+m, l-m')\f$:
 *
 * \f[d^{l}_{m m'}(\beta) = \sum_{k=k_{min}}^{k_{max}} \frac{(-1)^{k+m'-m}}{k!}
 * \frac{\sqrt{(l+m)!(l-m)!(l+m')!(l-m')!}}{(l+m-k)!(l-m'-k)!(k-m+m')!}
 * c^{2l+m-m'-2k} s^{2k-m+m'}\f]
 */
REAL8 SEOBWignerDAmp(UINT l, INT m, INT mp, REAL8 beta) 
{
  return WignerdMatrix(l, m, mp, beta);
}

REAL8 SEOBWignerDPhase(INT m, INT mp, REAL8 alpha, REAL8 gamma) 
{
  return m * alpha + mp * gamma;
}

void DestroySphHarmTimeSeries(SphHarmTimeSeries *ts)
{
    SphHarmTimeSeries* pop;
    while( (pop = ts) )
    {
        STRUCTFREE( pop->mode , COMPLEX16TimeSeries);
        // The tdata pointer is shared so we delete on the last node
        if( pop->next == NULL)
        {
            STRUCTFREE( pop->tdata , REAL8Vector);
        }
        ts = pop->next;
        MYFree( pop );
    }
}

void DestroySEOBCoreOutputs(SEOBCoreOutputs *all)
{
    if (all->hLM)
    {
        STRUCTFREE(all->hLM, SphHarmTimeSeries);
    }
    if (all->dyn)
    {
        STRUCTFREE(all->dyn, SEOBdynamics);
    }
    if (all->flux)
    {
        STRUCTFREE(all->flux, REAL8Vector);
    }
    if (all->Plm)
    {
        STRUCTFREE(all->Plm, SphHarmListCAmpPhaseSequence);
    }
    MYFree(all);
    return;
}

void DestroySEOBPrecCoreOutputs(SEOBPrecCoreOutputs *all)
{
    if (all->tVec)
    {
        STRUCTFREE(all->tVec, REAL8Vector);
    }
    if (all->hLM)
    {
        STRUCTFREE(all->hLM, SphHarmTimeSeries);
    }
    if (all->dyn)
    {
        STRUCTFREE(all->dyn, SEOBPrecdynamics);
    }
    if (all->Plm)
    {
        STRUCTFREE(all->Plm, SphHarmListCAmpPhaseSequence);
    }
    MYFree(all);
    return;
}


int XLALEOBFindRobustPeak(REAL8 *tPeakQuant, REAL8Vector *tVec,
                                REAL8Vector *quantVec,
                                UINT window_width) {
    // We begin at the end and go backwards
    UINT vlen = tVec->length, kk;
    UINT local_argmax = 0;
    UINT idx_global = 0;
    UINT lo, hi; // Bounds of local window
    // Global argmax over the entire array
    UINT global_arg_max = argmax(quantVec);
    REAL8 global_max = quantVec->data[global_arg_max];
    REAL8Vector *sl = NULL;
    REAL8 curr_max = 0;
    for ( kk = vlen - window_width - 1; kk > window_width; kk--) 
    {
        lo = kk - window_width;
        hi = kk + window_width +
            1; // Slice function does not return the value at the end
        sl = get_slice(quantVec, lo, hi);
        local_argmax = argmax(sl);
        if (sl->data[local_argmax] > sl->data[0] &&
            sl->data[local_argmax] > sl->data[sl->length - 1]) 
        {
            // We have *a* local max, figure out it's global index
            // Is the local argmax the largest local argmax so far?
            if (sl->data[local_argmax] > curr_max) 
            {
                curr_max = sl->data[local_argmax];
                idx_global = lo + local_argmax;
            }
        }
        STRUCTFREE(sl, REAL8Vector);
    }
    *tPeakQuant = 0;
    // Conditions under which we pick the last point of the dynamics:
    // i) we found no local maxes at all
    // ii) the global arg max is larger than the largest of the local maxes
    // by more than 2 % of the largest maxes value (ideally they should be equal)
    // iii) the  peak is  so close to end that we can't interpolate below.

    if (idx_global == 0 ||
        ((quantVec->data[global_arg_max] - quantVec->data[idx_global]) /
            quantVec->data[idx_global] >
        0.1) ||
        (idx_global > tVec->length - 4)) 
    {
        PRINT_LOG_INFO(LOG_WARNING, "Warning no local max found, using last point");
        *tPeakQuant = tVec->data[tVec->length - 1];
        return CEV_SUCCESS;
    }


    // We have a well-behaved local max. Get the time more accurately.
    // Now we need to interpolate and then set the derivative of the interpolant
    // to 0. We solve this via bisection in an interval of 3 points to the left
    // and right of the argmax.
    GSL_START;
    gsl_spline *spline = NULL;
    gsl_interp_accel *acc = NULL;
    spline = gsl_spline_alloc(gsl_interp_cspline, quantVec->length);
    acc = gsl_interp_accel_alloc();

    REAL8 time1 = tVec->data[idx_global - 3];
    REAL8 time2 = tVec->data[idx_global + 3];
    REAL8 timePeak = 0;
    REAL8 omegaDerivMid = 0;
    INT status, is_failed = 0;
    status = gsl_spline_init(spline, tVec->data, quantVec->data, quantVec->length);
    if (status != GSL_SUCCESS)
    {
        is_failed = 1;
        goto QUIT;
    }
    REAL8 omegaDeriv1 = gsl_spline_eval_deriv(spline, time1, acc);
    do 
    {
        timePeak = (time1 + time2) / 2.;
        omegaDerivMid = gsl_spline_eval_deriv(spline, timePeak, acc);

        if (omegaDerivMid * omegaDeriv1 < 0.0) 
        {
            time2 = timePeak;
        } 
        else 
        {
            omegaDeriv1 = omegaDerivMid;
            time1 = timePeak;
        }
    } while (time2 - time1 > 1.0e-8);
    *tPeakQuant = timePeak;
QUIT:
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    GSL_END;
    if (is_failed)
        return CEV_FAILURE;
    return CEV_SUCCESS;
}

