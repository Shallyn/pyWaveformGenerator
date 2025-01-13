/**
 * Writer: Xiaolin.liu
 * xiaolin.liu@mail.bnu.edu.cn
 *
 * This module contains basic functions for  calculation.
 * Functions list:
 * Kernel:
 * 20xx.xx.xx, LOC
 **/

#ifndef __INCLUDE_PCORE__
#define __INCLUDE_PCORE__
#include "newFactorizedWaveform.h"
#include "newFactorizedWaveformPrec.h"
#include "pBHRingdown.h"
#include "pEnergyFlux.h"
#include "pFactorizedWaveform.h"
#include "pHamiltonian.h"
#include "pInitialCondition.h"
#include "pInitialConditionExact.h"
#include "pMemory.h"
#include "pNQCcorrection.h"
#include "pRK4pdeIntegrator.h"
#include "pUtils.h"

void set_h_rISCO(REAL8 rISCO);
REAL8 get_h_rISCO();
#define SET_RISCO(RISCO)                                                                                               \
    do                                                                                                                 \
    {                                                                                                                  \
        set_h_rISCO(RISCO);                                                                                            \
    } while (0)

SpinEOBParams *CreateSpinEOBParams(REAL8 m1, REAL8 m2, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z,
                                   REAL8 e0, HyperParams *params);

INT SEOBInitialConditions(REAL8Vector *ICvalues, REAL8 MfMin, REAL8 ecc, SpinEOBParams *seobParams);

INT SEOBInitialConditions_Conserve(REAL8Vector *ICvalues, REAL8 MfMin, REAL8 ecc, SpinEOBParams *seobParams);

INT SEOBInitialConditions_e_anomaly(REAL8Vector *ICvalues, REAL8 MfMin, REAL8 ecc, REAL8 zeta, REAL8 xi,
                                    SpinEOBParams *seobParams);

INT CalculateAOmegaFromrpphi(REAL8 r, REAL8 pphi, SpinEOBParams *core, REAL8 *omegaOut);

REAL8 SEOBCalculatetplspin(REAL8 m1, REAL8 m2, REAL8 eta, REAL8 chi1dotZ, REAL8 chi2dotZ);
REAL8 FindClosestValueInIncreasingVector(REAL8Vector *vec, /**<< Input: monotonically increasing vector */
                                         REAL8 value       /**<< Input: value to look for */
);

INT SEOBIntegrateDynamics(REAL8Array **dynamics, INT *retLenOut, REAL8Vector *ICvalues, REAL8 EPS_ABS, REAL8 EPS_REL,
                          REAL8 deltaT, REAL8 deltaT_min, REAL8 tstart, REAL8 tend, SpinEOBParams *seobParams,
                          INT flagConstantSampling);

INT SEOBIntegrateDynamics_withfMax(REAL8Array **dynamics, INT *retLenOut, REAL8Vector *ICvalues, REAL8 EPS_ABS,
                                   REAL8 EPS_REL, REAL8 deltaT, REAL8 deltaT_min, REAL8 tstart, REAL8 tend,
                                   SpinEOBParams *seobParams, INT flagConstantSampling);

INT SEOBIntegrateDynamics_inverse(REAL8Array **dynamics, INT *retLenOut, REAL8Vector *ICvalues, REAL8 EPS_ABS,
                                  REAL8 EPS_REL, REAL8 deltaT, REAL8 deltaT_min, REAL8 tstart, REAL8 tend,
                                  SpinEOBParams *seobParams, INT flagConstantSampling);
INT SEOBIntegrateDynamics_SA_inverse(REAL8Array **dynamics, INT *retLenOut, REAL8Vector *ICvalues, REAL8 EPS_ABS,
                                     REAL8 EPS_REL, REAL8 deltaT, REAL8 deltaT_min, REAL8 tstart, REAL8 tend,
                                     SpinEOBParams *seobParams);
void SEOBConcactInverseDynToAdaSDyn_SA(REAL8Array **dyn_out, REAL8Array *dyn_inv, INT *retLen_out, INT retLen_inv);
void SEOBConcactInverseDynToAdaSDyn(REAL8Array **dyn_out, REAL8Array *dyn_inv, INT *retLen_out, INT retLen_inv);

INT SEOBIntegrateDynamics_adaptive(REAL8Array **dynamics, INT *retLenOut, REAL8Vector *ICvalues, REAL8 EPS_ABS,
                                   REAL8 EPS_REL, REAL8 deltaT, REAL8 deltaT_min, REAL8 tstart, REAL8 tend,
                                   SpinEOBParams *seobParams, INT flagConstantSampling);

INT SEOBIntegrateDynamics_Conserve(REAL8Array **dynamics, INT *retLenOut, REAL8Vector *ICvalues, REAL8 EPS_ABS,
                                   REAL8 EPS_REL, REAL8 deltaT, REAL8 deltaT_min, REAL8 tstart, REAL8 tend,
                                   SpinEOBParams *seobParams, INT flagConstantSampling);

INT SEOBIntegrateDynamics_SA(REAL8Array **dynamics, INT *retLenOut, REAL8Vector *ICvalues, REAL8 EPS_ABS, REAL8 EPS_REL,
                             REAL8 deltaT, REAL8 deltaT_min, REAL8 tstart, REAL8 tend, SpinEOBParams *seobParams);

INT SEOBIntegrateDynamics_SA_withFMax(REAL8Array **dynamics, INT *retLenOut, REAL8Vector *ICvalues, REAL8 EPS_ABS,
                                      REAL8 EPS_REL, REAL8 deltaT, REAL8 deltaT_min, REAL8 tstart, REAL8 tend,
                                      SpinEOBParams *seobParams);

INT CutSEOBdynamics(SEOBdynamics **eobdyn, REAL8 MfMin);
INT CutSEOBSAdynamics(SEOBSAdynamics **eobdyn, REAL8 MfMin);
void OrbitalPhaseReduce(SEOBdynamics *dyn, REAL8 phiD, REAL8 phiM, REAL8 phi);
INT SetZeroPhaseAtTime(SEOBdynamics *dyn, REAL8 t, REAL8 *ret_dphiD, REAL8 *ret_dphiM, REAL8 *ret_dphi);
void OrbitalPhaseReduceSA(SEOBSAdynamics *dyn, REAL8 phi);
INT SetZeroPhaseAtTimeSA(SEOBSAdynamics *dyn, REAL8 t, REAL8 *ret_dphi);

int XLALEOBCheckNyquistFrequency(REAL8 m1, REAL8 m2, REAL8 spin1[3], REAL8 spin2[3], REAL8 deltaT, UINT ell_max);
INT SEOBComputeExtendedSEOBdynamics(SEOBdynamics **seobdynamics, REAL8Array *dynamics, INT retLen,
                                    SpinEOBParams *seobParams);
INT SEOBComputeExtendedSEOBSAdynamics(SEOBSAdynamics **seobsadynamics, REAL8Array *dynamics, INT length,
                                      SpinEOBParams *seobParams);

INT SEOBComputeExtendedSEOBdynamics_Conserve(SEOBdynamics **seobdynamics, REAL8Array *dynamics, INT retLen,
                                             SpinEOBParams *seobParams);

int SEOBInterpolateDynamicsAtTime(REAL8Vector **seobdynamics_values, /**<< Output: pointer to vector for
                                                                        seobdynamics interpolated values */
                                  REAL8 t,                           /**<< Input: time at which to evaluate */
                                  SEOBdynamics *seobdynamics         /**<< Input: SEOB dynamics */
);

int SEOBInterpolateSADynamicsAtTime(REAL8Vector **seobdynamics_values, /**<< Output: pointer to vector for
                                                                          seobdynamics interpolated values */
                                    REAL8 t,                           /**<< Input: time at which to evaluate */
                                    SEOBSAdynamics *seobdynamics       /**<< Input: SEOB dynamics */
);

int SEOBLocateTimePeakOmega(REAL8 *tPeakOmega,          /**<< Output: time of peak of Omega if found (see inside
                                         XLALSimLocateOmegaTime for what is returned otherwise) */
                            INT *foundPeakOmega,        /**<< Output: flag indicating wether tPeakOmega has
                                                           been found */
                            REAL8Array *dynamics,       /**<< Input: array for dynamics */
                            SEOBdynamics *seobdynamics, /**<< Input: SEOB dynamics object */
                            UINT retLen,                /**<< Input: length of dynamics */
                            SpinEOBParams *seobParams   /**<< SEOB params */
);

int SEOBSALocateTimePeakOmega(REAL8 *tPeakOmega,            /**<< Output: time of peak of Omega if found (see inside
                                             XLALSimLocateOmegaTime for what is returned otherwise) */
                              INT *foundPeakOmega,          /**<< Output: flag indicating wether tPeakOmega has
                                                               been found */
                              SEOBSAdynamics *seobdynamics, /**<< Input: SEOB dynamics object */
                              UINT retLen,                  /**<< Input: length of dynamics */
                              SpinEOBParams *seobParams     /**<< SEOB params */
);

int SEOBLFrameVectors(REAL8Vector **S1,        /**<<Output: S1 in L-n frame */
                      REAL8Vector **S2,        /**<<Output: S2 in L-n frame */
                      REAL8Vector *seobvalues, /**<<Input: vector of extended dynamics */
                      REAL8 m1,                /**<<Input: mass of the first object in solar masses */
                      REAL8 m2,                /**<<Input: mass of the second object in solar masses */
                      INT flagZframe);

int SEOBJfromDynamics(REAL8Vector **J,          /**<< Output: pointer to vector J */
                      REAL8Vector *seobvalues,  /**<< Input: vector for
                                                   extended dynamics values */
                      SpinEOBParams *seobParams /**<< SEOB params */
);

int SEOBLhatfromDynamics(REAL8Vector **L,          /**<< Output: pointer to vector L */
                         REAL8Vector *seobvalues,  /**<< Input: vector for extended dynamics values */
                         SpinEOBParams *seobParams /**<< SEOB params */
);

int SEOBBuildJframeVectors(REAL8Vector *e1J, /**<< Output: vector for e1J, already allocated */
                           REAL8Vector *e2J, /**<< Output: vector for e2J, already allocated */
                           REAL8Vector *e3J, /**<< Output: vector for e3J, already allocated */
                           REAL8Vector *JVec /**<< Input: vector J */
);

int SEOBEulerI2JFromJframeVectors(REAL8 *alphaI2J,  /**<< Output: Euler angle alpha I2J */
                                  REAL8 *betaI2J,   /**<< Output: Euler angle beta I2J */
                                  REAL8 *gammaI2J,  /**<< Output: Euler angle gamma I2J */
                                  REAL8Vector *e1J, /**<< Input: unit Jframe vector e1J */
                                  REAL8Vector *e2J, /**<< Input: unit Jframe vector e2J */
                                  REAL8Vector *e3J  /**<< Input: unit Jframe vector e3J */
);

int SEOBCalculateSphHarmListNQCCoefficientsV4(
    SphHarmListEOBNonQCCoeffs **nqcCoeffsList, /**<< Output: non-quasi-circular coefficients as a list
                                                  for each mode */
    INT modes[][2],                            /**<< Input: array of modes (l,m) */
    UINT nmodes,                               /**<< Input: number of modes (l,m) */
    REAL8 tPeakOmega,                          /**<< Input: time of peak of Omega */
    SEOBdynamics *seobdynamics,                /**<< Input: SEOB dynamics */
    SpinEOBParams *seobParams,                 /**<< Input: SEOB params */
    REAL8Vector *chi1_omegaPeak,               /**<< Input: dimensionless spin 1 at peak of
                                                  omega in L_N frame */
    REAL8Vector *chi2_omegaPeak                /**<< Input: dimensionless spin 2 at peak of
                                                   omega in L_N frame */
);

int SEOBSACalculateSphHarmListNQCCoefficientsV4(
    SphHarmListEOBNonQCCoeffs **nqcCoeffsList, /**<< Output: non-quasi-circular coefficients as a list
                                                  for each mode */
    INT modes[][2],                            /**<< Input: array of modes (l,m) */
    UINT nmodes,                               /**<< Input: number of modes (l,m) */
    REAL8 tPeakOmega,                          /**<< Input: time of peak of Omega */
    SEOBSAdynamics *seobdynamics,              /**<< Input: SEOB dynamics */
    SpinEOBParams *seobParams);

int SEOBCalculateSphHarmListhlmAmpPhase(SphHarmListCAmpPhaseSequence **listhlm,   /**<< Output: list of modes for hlm */
                                        INT modes[][2],                           /**<< Input: array of modes (l,m) */
                                        UINT nmodes,                              /**<< Input: number of modes (l,m) */
                                        SEOBdynamics *seobdynamics,               /**<< Input: SEOB dynamics */
                                        SphHarmListEOBNonQCCoeffs *listnqcCoeffs, /**<< Input: list of NQCs */
                                        SpinEOBParams *seobParams,                /**<< SEOB params */
                                        UINT flagNQC /**<< flag to choose wether or not to include NQC */
);

int SEOBSACalculateSphHarmListhlmAmpPhase(SphHarmListCAmpPhaseSequence **listhlm, /**<< Output: list of modes for hlm */
                                          INT modes[][2],                         /**<< Input: array of modes (l,m) */
                                          UINT nmodes,                            /**<< Input: number of modes (l,m) */
                                          SEOBSAdynamics *seobdynamics,           /**<< Input: SEOB dynamics */
                                          SphHarmListEOBNonQCCoeffs *listnqcCoeffs, /**<< Input: list of NQCs */
                                          SpinEOBParams *seobParams,                /**<< SEOB params */
                                          UINT flagNQC /**<< flag to choose wether or not to include NQC */
);

int SEOBCalculateSphHarmListhlmAmpPhase_noNQC(
    SphHarmListCAmpPhaseSequence **listhlm, /**<< Output: list of modes for hlm */
    INT modes[][2],                         /**<< Input: array of modes (l,m) */
    UINT nmodes,                            /**<< Input: number of modes (l,m) */
    SEOBdynamics *seobdynamics,             /**<< Input: SEOB dynamics */
    SpinEOBParams *seobParams               /**<< SEOB params */
);

int SEOBGetFinalSpinMass(REAL8 *finalMass,         /**<< Output: final mass computed from fit (scaled by M) */
                         REAL8 *finalSpin,         /**<< Output: final spin computed from fit
                                                      (dimensionless) */
                         REAL8Vector *seobvalues,  /**<< Input: vector for dynamics values at time of
                                                      peak of omega */
                         SpinEOBParams *seobParams /**<< Input: SEOB params */
);

int SEOBSAGetFinalSpinMass(REAL8 *finalMass,         /**<< Output: final mass computed
                                                        from fit (scaled by M) */
                           REAL8 *finalSpin,         /**<< Output: final spin computed
                                                        from fit (dimensionless) */
                           SpinEOBParams *seobParams /**<< Input: SEOB params */
);

INT SEOBAttachRDToSphHarmListhPlm(SphHarmListCAmpPhaseSequence **listhPlm_RDattached, /**<< Output: list of
                                                               extended modes hlm with RD attached */
                                  COMPLEX16Vector **sigmaQNMlm0, /**<< Output: list of QNM complex frequency
                                                       for modes lm, 0th overtone (dimensionless) */
                                  INT4 modes[][2],               /**<< Input: array of modes (l,m) */
                                  UINT nmodes,                   /**<< Input: number of modes (l,m) */
                                  REAL8 finalMass, /**<< Input: final mass computed from fit (scaled by M) */
                                  REAL8 finalSpin, /**<< Input: final spin computed from fit (dimensionless) */
                                  SphHarmListCAmpPhaseSequence *listhPlm, /**<< Input: list of modes hlm */
                                  REAL8 deltaT,                           /**<< Input: time step */
                                  UINT retLen,        /**<< Input: length of the input modes and dynamics */
                                  UINT retLenRDPatch, /**<< Input: length of the ringdown patch */
                                  REAL8 tAttach,      /**<< Input: time of attachment */
                                  // REAL8 tStart, /**<< Input: starting time (of the HiS part) */
                                  REAL8Vector *seobvalues,    /**<< Input: vector for dynamics values at time of
                                                                 peak of omega */
                                  SEOBdynamics *seobdynamics, /**<< Input: SEOB dynamics */
                                  SpinEOBParams *seobParams   /**<< SEOB params */
);

INT SEOBSAAttachRDToSphHarmListhPlm(SphHarmListCAmpPhaseSequence **listhPlm_RDattached, /**<< Output: list of
                                                                 extended modes hlm with RD attached */
                                    COMPLEX16Vector **sigmaQNMlm0, /**<< Output: list of QNM complex frequency
                                                         for modes lm, 0th overtone (dimensionless) */
                                    INT4 modes[][2],               /**<< Input: array of modes (l,m) */
                                    UINT nmodes,                   /**<< Input: number of modes (l,m) */
                                    REAL8 finalMass, /**<< Input: final mass computed from fit (scaled by M) */
                                    REAL8 finalSpin, /**<< Input: final spin computed from fit (dimensionless) */
                                    SphHarmListCAmpPhaseSequence *listhPlm, /**<< Input: list of modes hlm */
                                    REAL8 deltaT,                           /**<< Input: time step */
                                    UINT retLen,        /**<< Input: length of the input modes and dynamics */
                                    UINT retLenRDPatch, /**<< Input: length of the ringdown patch */
                                    REAL8 tAttach,      /**<< Input: time of attachment */
                                    // REAL8 tStart, /**<< Input: starting time (of the HiS part) */
                                    REAL8Vector *seobvalues,      /**<< Input: vector for dynamics values at time of
                                                                     peak of omega */
                                    SEOBSAdynamics *seobdynamics, /**<< Input: SEOB dynamics */
                                    SpinEOBParams *seobParams     /**<< SEOB params */
);

int SEOBJoinTimeVector(REAL8Vector **tVecPmodes,       /**<< Output: vector of times for P-modes
                                                          (AdaS+HiS+RDpatch) */
                       UINT *retLenPmodes,             /**<< Output: length of output vector of times for
                                                           P-modes */
                       REAL8 *tJoinHiS,                /**<< Output: first time >= tstartHiS */
                       UINT *indexJoinHiS,             /**<< Output: first index >= tstartHiS */
                       REAL8 *tJoinAttach,             /**<< Output: first time >= tAttach */
                       UINT *indexJoinAttach,          /**<< Output: first index >= tAttach */
                       UINT retLenHiSRDpatch,          /**<< Input: length of RD patch to be added at the
                                                           end of HiS with the same constant sampling */
                       REAL8 deltaTHiS,                /**<< Input: time step for the high sampling */
                       REAL8 tstartHiS,                /**<< Input: time of start of HiS */
                       REAL8 tAttach,                  /**<< Input: time of attachment */
                       SEOBdynamics *seobdynamicsAdaS, /**<< Input: SEOB dynamics with adaptive-sampling */
                       SEOBdynamics *seobdynamicsHiS   /**<< Input: SEOB dynamics with high-sampling */
);

int SEOBSAJoinTimeVector(REAL8Vector **tVecPmodes,         /**<< Output: vector of times for P-modes
                                                              (AdaS+HiS+RDpatch) */
                         UINT *retLenPmodes,               /**<< Output: length of output vector of times for
                                                               P-modes */
                         REAL8 *tJoinHiS,                  /**<< Output: first time >= tstartHiS */
                         UINT *indexJoinHiS,               /**<< Output: first index >= tstartHiS */
                         REAL8 *tJoinAttach,               /**<< Output: first time >= tAttach */
                         UINT *indexJoinAttach,            /**<< Output: first index >= tAttach */
                         UINT retLenHiSRDpatch,            /**<< Input: length of RD patch to be added at the
                                                               end of HiS with the same constant sampling */
                         REAL8 deltaTHiS,                  /**<< Input: time step for the high sampling */
                         REAL8 tstartHiS,                  /**<< Input: time of start of HiS */
                         REAL8 tAttach,                    /**<< Input: time of attachment */
                         SEOBSAdynamics *seobdynamicsAdaS, /**<< Input: SEOB dynamics with adaptive-sampling */
                         SEOBSAdynamics *seobdynamicsHiS   /**<< Input: SEOB dynamics with high-sampling */
);

int SEOBSAAttachAdaSandHiSRTimeVector(REAL8Vector **tVecPmodes, SEOBSAdynamics *seobdynamicsAdaS,
                                      SEOBSAdynamics *seobdynamicsHiS, REAL8 tstartHiS);

int SEOBJoinDynamics(SEOBdynamics **seobdynamicsJoined, /**<< Output: pointer to joined dynamics */
                     SEOBdynamics *seobdynamics1,       /**<< Input: first dynamics */
                     SEOBdynamics *seobdynamics2,       /**<< Input: second dynamics */
                     UINT indexJoin12,                  /**<< Input: index where to join the two dynamics */
                     UINT indexEnd2                     /**<< Input: index of the joined dynamics where to stop
                                                            dynamics 2 (excluded) */
);

int SEOBSAJoinDynamics(SEOBSAdynamics **seobdynamicsJoined, /**<< Output: pointer to joined dynamics */
                       SEOBSAdynamics *seobdynamics1,       /**<< Input: first dynamics */
                       SEOBSAdynamics *seobdynamics2,       /**<< Input: second dynamics */
                       UINT indexJoin12,                    /**<< Input: index where to join the two dynamics */
                       UINT indexEnd2                       /**<< Input: index of the joined dynamics where to stop
                                                                dynamics 2 (excluded) */
);

int SEOBAmplitudePeakFromAmp22Amp21(REAL8 *tPeak,                           /**<< Output: time of peak */
                                    UINT *indexPeak,                        /**<< Output: index of peak */
                                    SphHarmListCAmpPhaseSequence *listhPlm, /**<< Input: list of modes hlm */
                                    INT4 modes[][2],                        /**<< Input: array of modes (l,m) */
                                    UINT nmodes,                            /**<< Input: number of modes (l,m) */
                                    REAL8Vector *tVec                       /**<< Input: vector of times */
);

int SEOBJoinSphHarmListhlm(SphHarmListCAmpPhaseSequence **listhlm_joined, /**<< Output: list of joined modes */
                           SphHarmListCAmpPhaseSequence *listhlm_1,       /**<< Input: list of modes 1 */
                           SphHarmListCAmpPhaseSequence *listhlm_2,       /**<< Input: list of modes 2 */
                           INT4 modes[][2],                               /**<< Input: array of modes (l,m) */
                           UINT nmodes,                                   /**<< Input: number of modes (l,m) */
                           UINT indexJoin12 /**<< Input: index where to join the two dynamics */
);

int SEOBEulerJ2PFromDynamics(REAL8Vector **alphaJ2P,     /**<< Output: pointer to vector for alpha J2P */
                             REAL8Vector **betaJ2P,      /**<< Output: pointer to vector for beta J2P */
                             REAL8Vector **gammaJ2P,     /**<< Output: pointer to vector for gamma J2P */
                             REAL8Vector *e1J,           /**<< Input: unit Jframe vector e1J */
                             REAL8Vector *e2J,           /**<< Input: unit Jframe vector e2J */
                             REAL8Vector *e3J,           /**<< Input: unit Jframe vector e3J */
                             UINT retLen,                /**<< Input: total length of Euler angles data to be allocated
                                                             (length of P-modes) */
                             UINT indexStop,             /**<< Input: index where we stop the computation (excluded,
                                                             index of time of attachment) */
                             SEOBdynamics *seobdynamics, /**<<Input: SEOB dynamics (joined AdaS+HiS, up
                                                            to tAttach) */
                             SpinEOBParams *seobParams   /**<< SEOB params */
);

int SEOBEulerJ2PPostMergerExtension(REAL8Vector *alphaJ2P, /**<< Output: vector for alpha J2P, already allocated */
                                    REAL8Vector *betaJ2P,  /**<< Output: vector for beta J2P, already allocated */
                                    REAL8Vector *gammaJ2P, /**<< Output: vector for gamma J2P, already allocated */
                                    COMPLEX16 sigmaQNM220, /**<< Input: complex frequency for QNM 22, 0th overtone */
                                    COMPLEX16 sigmaQNM210, /**<< Input: complex frequency for QNM 21, 0th overtone */
                                    REAL8Vector *tVec,     /**<< Input: time vector for Euler angles data (length
                                                              of P-modes) */
                                    UINT retLen,           /**<< Input: total length of Euler angles data (length of
                                                               P-modes) */
                                    UINT indexStart,       /**<< Input: index where we start the extension (included,
                                                               index of time of attachment) */
                                    SpinEOBParams *seobParams, /**<< SEOB params */
                                    INT flip /** << a flag of whether to flip the sign of the precession
                                              * frequency
                                              */
);

int SEOBRotateInterpolatehJlmReImFromSphHarmListhPlmAmpPhase(

    SphHarmTimeSeries **hJlm,                /**<< Output: hJlm time series, will contain
                                                complex values on fixed sampling */
    SphHarmListCAmpPhaseSequence **listhClm, /**<< Output: list of C-frame modes hClm */
    INT modes[][2],                          /**<< Input: array of modes (l,m) */
    UINT nmodes,                             /**<< Input: number of modes (l,m) */
    INT modes_lmax,                          /**<< Input: maximum value of l in modes (l,m) */
    REAL8 deltaT,                            /**<< Input: time step for the hJlm timeseries */
    UINT retLenTS,                           /**<< Input: number of samples for the hJlm timeseries */
    REAL8Vector *tVecPmodes,                 /**<< Input: irregular time vector on which the
                                                hPlm and Euler angles are given */
    SphHarmListCAmpPhaseSequence *listhPlm,  /**<< Input: list of P-frame modes hPlm */
    REAL8Vector *alphaJ2P,                   /**<< Input: vector for Euler angle alpha J2P */
    REAL8Vector *betaJ2P,                    /**<< Input: vector for Euler angle beta J2P */
    REAL8Vector *gammaJ2P                    /**<< Input: vector for Euler angle gamma J2P */
);

int SEOBPlmListToSphHarmTimeSeries(

    SphHarmTimeSeries **hJlm, /**<< Output: hJlm time series, will contain
                                 complex values on fixed sampling */
    INT modes[][2],           /**<< Input: array of modes (l,m) */
    UINT nmodes,              /**<< Input: number of modes (l,m) */
    INT modes_lmax,           /**<< Input: maximum value of l in modes (l,m) */
    REAL8 deltaT,             /**<< Input: time step for the hJlm timeseries */
    UINT retLenTS,            /**<< Input: number of samples for the hJlm timeseries */
    SphHarmListCAmpPhaseSequence *listhPlm, REAL8Vector *alphaJ2P, /**<< Input: vector for Euler angle alpha J2P */
    REAL8Vector *betaJ2P,                                          /**<< Input: vector for Euler angle beta J2P */
    REAL8Vector *gammaJ2P                                          /**<< Input: vector for Euler angle gamma J2P */
);

int SEOBRotatehIlmFromhJlm(SphHarmTimeSeries **hIlm, /**<< Output: hIlm time series, complex values on
                                                        fixed sampling */
                           SphHarmTimeSeries *hJlm,  /**<< Output: hJlm time series, complex values on
                                                        fixed sampling */
                           INT modes_lmax,           /**<< Input: maximum value of l in modes (l,m) */
                           REAL8 alphaI2J,           /**<< Input: Euler angle alpha I->J */
                           REAL8 betaI2J,            /**<< Input: Euler angle beta I->J */
                           REAL8 gammaI2J,           /**<< Input: Euler angle gamma I->J */
                           REAL8 deltaT              /**<< Input: time step, necessary to initialize new timeseries
                                                      */
);

int SEOBComputehplushcrossFromhIlm(REAL8TimeSeries *hplusTS,  /**<< Output: time series for hplus, already created */
                                   REAL8TimeSeries *hcrossTS, /**<< Output: time series for hplus, already created */
                                   INT modes_lmax,            /**<< Input: maximum value of l */
                                   SphHarmTimeSeries *hIlm,   /**<< Input: list with time series for each mode hIlm */
                                   REAL8 amp0,                /**<< Input: amplitude prefactor */
                                   REAL8 inc,                 /**<< Input: inclination */
                                   REAL8 phi,                 /**<< Input: phase */
                                   INT is_only22);

INT SEOBCalculateNQCWindowFactorsFromDyn(SEOBdynamics *dyn, REAL8 thPeak, REAL8 tr6M, REAL8 tHiStart, REAL8 tThresh,
                                         INT is_first, SpinEOBParams *seobParams);

INT SEOBSACalculateNQCWindowFactorsFromDyn(SEOBSAdynamics *dyn, REAL8 thPeak, REAL8 tr6M, REAL8 tHiStart, REAL8 tThresh,
                                           INT is_first, SpinEOBParams *seobParams);

INT dbg_CalculateWaveformFromDynamicsAdaS(SEOBdynamics *seobdynamics, SpinEOBParams *seobParams, INT ModeL, INT ModeM,
                                          COMPLEX16TimeSeries **out_new, COMPLEX16TimeSeries **out_old);
INT dbg_CalculateNQCTimeSeries(SEOBdynamics *seobdynamics, SpinEOBParams *seobParams,
                               SphHarmListEOBNonQCCoeffs *listnqcCoeffs, INT ModeL, INT ModeM,
                               COMPLEX16TimeSeries **ret);

INT CheckStopCondition(SEOBdynamics *seobdynamics, SpinEOBParams *seobParams, SphHarmListEOBNonQCCoeffs *listnqcCoeffs,
                       REAL8 tPeak);

INT CheckStopConditionSA(SEOBSAdynamics *seobdynamics, SpinEOBParams *seobParams,
                         SphHarmListEOBNonQCCoeffs *listnqcCoeffs, REAL8 tPeak);

INT SEOBInitialConditions_egw(REAL8Vector *ICvalues, REAL8 MfMin, REAL8 ecc, SpinEOBParams *seobParams);

INT SEOBIntegrateDynamics_egw(REAL8Array **dynamics, INT *retLenOut, REAL8Vector *ICvalues, REAL8 EPS_ABS,
                              REAL8 EPS_REL, REAL8 deltaT, REAL8 deltaT_min, REAL8 tstart, REAL8 tend,
                              SpinEOBParams *seobParams, INT flagConstantSampling);
INT CalculateSAh22SeriesFromrpphi(REAL8 r, REAL8 pphi, SpinEOBParams *core, REAL8 *omega22, REAL8 *egw);
INT find_SApphi_from_rpm(REAL8 rp, REAL8 rm, SpinEOBParams *core, REAL8 *pphi);
INT EvaluateOmega22SA_form_rpphi(REAL8 r, REAL8 pphi, SpinEOBParams *seobParams, REAL8 *omega22);
REAL8 calc_egw_from_e22(REAL8 e22);
REAL8 calc_e22_from_egw(REAL8 egw);
INT find_SACircpphi_from_r(REAL8 r, SpinEOBParams *core, REAL8 *pphi);

#endif
