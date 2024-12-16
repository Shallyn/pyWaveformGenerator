/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_PUTILS__
#define __INCLUDE_PUTILS__
#include "pStructs.h"
#include "mySpherical.h"
INT get_egw_flag();

INT get_DebugFlag();
void set_DebugFlag(INT flag);
#define IS_DEBUG (get_DebugFlag())
#define DEBUG_START \
do{set_DebugFlag(1);}while(0)
#define DEBUG_END \
do{set_DebugFlag(0);}while(0)

INT get_InputDebugFlag();
void set_InputDebugFlag(INT flag);
#define INPUT_DEBUG_FLAG (get_InputDebugFlag())
#define SET_INPUT_DEBUG_FLAG(flag) \
do{set_InputDebugFlag(flag);}while(0)

INT get_VersionFlag();
void set_VersionFlag(INT flag);
#define CODE_VERSION (get_VersionFlag())
#define SET_CODE_VERSION(ver) \
do{set_VersionFlag(ver);}while(0)

INT get_ConserveFlag();
REAL8 get_ConserveTime();
void set_ConserveFlag(INT flag, REAL8 ctime);
#define CONSERVE_FLAG (get_ConserveFlag())
#define GET_CONSERV_TIME (get_ConserveTime())
#define SET_CONSERV(flag, ctime) \
do{set_ConserveFlag(flag, ctime);}while(0)

REAL8 *rotate_vector(const REAL8 v[], const REAL8 k[], const REAL8 theta,REAL8 result[] );
REAL8 inner_product3d(const REAL8 values1[], const REAL8 values2[]);
REAL8 *cross_product3d(const REAL8 values1[],
                        const REAL8 values2[],
                        REAL8 result[] );
REAL8 SEOBWignerDAmp(UINT l, INT m, INT mp, REAL8 beta);
REAL8 SEOBWignerDPhase(INT m, INT mp, REAL8 alpha, REAL8 gamma);

REAL8 SEOBCalculateChiS(REAL8 chi1dotZ, REAL8 chi2dotZ);
REAL8 SEOBCalculateChiA(REAL8 chi1dotZ, REAL8 chi2dotZ);
                                   
SpinEOBDynamics *CreateSpinEOBDynamics(UINT length);
void DestroySpinEOBDynamics(SpinEOBDynamics *dy);
void DestroySpinEOBParams(SpinEOBParams *params);

INT EOBCalculateSigmaStar(REAL8Vector *sigmaStar,
                          REAL8 mass1,
                          REAL8 mass2,
                          REAL8Vector *s1norm,
                          REAL8Vector *s2norm);

INT EOBCalculateSigmaKerr(REAL8Vector *sigmaKerr,
                          REAL8Vector *s1norm,
                          REAL8Vector *s2norm);
INT XLALSimIMRSpinEOBCalculateSigmaStar(REAL8Vector *sigmaStar,
                          REAL8 mass1,
                          REAL8 mass2,
                          REAL8Vector *s1,
                          REAL8Vector *s2);
INT XLALSimIMRSpinEOBCalculateSigmaKerr(REAL8Vector *sigmaKerr,
                                    REAL8 mass1,
                                    REAL8 mass2,
                                    REAL8Vector *s1,
                                    REAL8Vector *s2);
SEOBdynamics *CreateSEOBdynamics(INT length);
void DestroySEOBdynamics(SEOBdynamics *seobdynamics);

SEOBSAdynamics *CreateSEOBSAdynamics(INT length);
void DestroySEOBSAdynamics(SEOBSAdynamics *seobdynamics);

SEOBPrecdynamics *CreateSEOBPrecdynamics(INT length);
void DestroySEOBPrecdynamics(SEOBPrecdynamics *seobdynamics);

int SphHarmListEOBNonQCCoeffs_AddMode(
    SphHarmListEOBNonQCCoeffs **list_prepended,       /* List structure to prepend to */
    EOBNonQCCoeffs *nqcCoeffs, /* Mode data to be added */
    UINT l,                   /*< Mode number l */
    INT m                     /*< Mode number m */
);
SphHarmListEOBNonQCCoeffs *SphHarmListEOBNonQCCoeffs_GetMode(
    SphHarmListEOBNonQCCoeffs
        *list, /* List structure to get a particular mode from */
    UINT l,   /*< Mode number l */
    INT m     /*< Mode number m */
);

int DestroySphHarmListEOBNonQCCoeffs
(
    SphHarmListEOBNonQCCoeffs *list /* Pointer to list to be destroyed */
);

int CAmpPhaseSequence_Init(
    CAmpPhaseSequence **campphase, /* Double pointer to campphase sequence */
    int size                       /* Size of data */
);

void DestroyCAmpPhaseSequence(
    CAmpPhaseSequence *campphase /* Pointer to structure to be destroyed */
);


int SphHarmListCAmpPhaseSequence_AddMode(
    SphHarmListCAmpPhaseSequence *
        *list_prepended,          /* List structure to prepend to */
    CAmpPhaseSequence *campphase, /* Mode data to be added */
    UINT l,                      /*< Mode number l */
    INT  m                        /*< Mode number m */
);

SphHarmListCAmpPhaseSequence *SphHarmListCAmpPhaseSequence_GetMode(
        SphHarmListCAmpPhaseSequence *list, /* List structure to get a particular mode from */
        UINT l,   /*< Mode number l */
        INT m     /*< Mode number m */
);

int DestroySphHarmListCAmpPhaseSequence(
    SphHarmListCAmpPhaseSequence *list /* Pointer to list to be destroyed */
);

void XLALSphHarmTimeSeriesSetTData(
            SphHarmTimeSeries *ts, /**< Linked list to be prepended */
            REAL8Vector* tdata /**< series of time data*/
);

SphHarmTimeSeries* XLALSphHarmTimeSeriesAddMode(
            SphHarmTimeSeries *appended, /**< Linked list to be prepended */
            const COMPLEX16TimeSeries* inmode, /**< Time series of h_lm mode being prepended */
            UINT l, /**< l index of h_lm mode being prepended */
            INT m /**< m index of h_lm mode being prepended */
);

COMPLEX16TimeSeries* XLALSphHarmTimeSeriesGetMode(
            SphHarmTimeSeries *ts, /**< linked list to extract mode from */
            UINT l, /**< l index of h_lm mode to get */
            INT m /**< m index of h_lm mode to get */
);
void DestroySphHarmTimeSeries(SphHarmTimeSeries *ts);
void DestroySEOBCoreOutputs(SEOBCoreOutputs *all);
void DestroySEOBPrecCoreOutputs(SEOBPrecCoreOutputs *all);

int XLALEOBFindRobustPeak(REAL8 *tPeakQuant, REAL8Vector *tVec,
                                REAL8Vector *quantVec,
                                UINT window_width);
#endif

