/**
 * Writer: Xiaolin.liu
 * xiaolin.liu@mail.bnu.edu.cn
 *
 * This module contains basic functions for  calculation.
 * Functions list:
 * Kernel:
 * 20xx.xx.xx, LOC
 **/

#ifndef __INCLUDE_PBHRINGDOWN__
#define __INCLUDE_PBHRINGDOWN__
#include "pUtils.h"

INT4 XLALSimIMREOBGenerateQNMFreqV2Prec(
    COMPLEX16Vector *modefreqs, /**<< OUTPUT, complex freqs of overtones in unit of Hz */
    const REAL8 mass1,          /**<< The mass of the 1st component (in Solar masses) */
    const REAL8 mass2,          /**<< The mass of the 2nd component (in Solar masses) */
    const REAL8 spin1[3],       /**<< The spin of the 1st object; only needed for
                                   spin waveforms */
    const REAL8 spin2[3],       /**<< The spin of the 2nd object; only needed for
                                   spin waveforms */
    UINT l,                     /**<< The l value of the mode in question */
    INT m,                      /**<< The m value of the mode in question */
    UINT nmodes                 /**<< The number of overtones that should be included (max 8) */
);

INT XLALSimIMREOBFinalMassSpinPrec(REAL8 *finalMass,     /**<< OUTPUT, the final mass (scaled by original total
                                                            mass) */
                                   REAL8 *finalSpin,     /**<< OUTPUT, the final spin (scaled by final mass) */
                                   const REAL8 mass1,    /**<< The mass of the 1st component of the system */
                                   const REAL8 mass2,    /**<< The mass of the 2nd component of the system */
                                   const REAL8 spin1[3], /**<< The spin of the 1st object; only needed for
                                                            spin waveforms */
                                   const REAL8 spin2[3]  /**<< The spin of the 2nd object; only needed for spin
                                                            waveforms */
);

INT XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec(
    COMPLEX16Vector *modefreqs, /**<< OUTPUT, complex freqs of overtones in unit of Hz */
    const REAL8 mass1,          /**<< The mass of the 1st component (in Solar masses) */
    const REAL8 mass2,          /**<< The mass of the 2nd component (in Solar masses) */
    const REAL8 finalMass,      /**<< The mass of the final BH (scaled by original
                                   total mass) */
    REAL8 finalSpin,            /**<< The dimensionless spin of the final BH */
    UINT l,                     /**<< The l value of the mode in question */
    INT m,                      /**<< The m value of the mode in question */
    UINT nmodes                 /**<< The number of overtones that should be included (max 8) */
);

INT XLALSimIMREOBAttachFitRingdown(REAL8Vector *signal1, /**<< OUTPUT, Real of inspiral waveform to which we
                                                            attach ringdown */
                                   REAL8Vector *signal2, /**<< OUTPUT, Imag of inspiral waveform to which we
                                                            attach ringdown */
                                   const INT l,          /**<< Current mode l */
                                   const INT m,          /**<< Current mode m */
                                   const REAL8 dt,       /**<< Sample time step (in seconds) */
                                   const REAL8 mass1,    /**<< First component mass (in Solar masses) */
                                   const REAL8 mass2,    /**<< Second component mass (in Solar masses) */
                                   const REAL8 spin1x,   /**<<The spin of the first object;  */
                                   const REAL8 spin1y,   /**<<The spin of the first object;  */
                                   const REAL8 spin1z,   /**<<The spin of the first object;  */
                                   const REAL8 spin2x,   /**<<The spin of the second object; */
                                   const REAL8 spin2y,   /**<<The spin of the second object; */
                                   const REAL8 spin2z,   /**<<The spin of the second object; */
                                   const REAL8 finalM,   /**<<The spin of the final BH;      */
                                   const REAL8 finalS,   /**<< The magnitude of the final spin */
                                   REAL8Vector *timeVec, /**<< Vector containing the time values */
                                   REAL8Vector *matchrange,
                                   /**<< Time values chosen as points for performing comb matching */
                                   UINT *indAmpMax /**<<This is used only for SEOBNRv4HM, it is needed to
                                                      remember the attaching point of the 22 mode */
);

#endif
