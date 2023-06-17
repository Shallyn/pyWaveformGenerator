/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "myOptparser.h"
#include "pCore.h"
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include "myLog.h"
#include "myFileIO.h"

int debug();

INT get_PrecFlag();
void set_PrecFlag(INT flag);

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

void print_waveform(REAL8TimeSeries *hp, REAL8TimeSeries *hc);
void print_SphHarmTimeSeries(SphHarmTimeSeries *hlm, INT mode);
void print_SEOBdynamics(SEOBdynamics *dyn);
void print_SEOBPrecdynamics(SEOBPrecdynamics *dyn);
void print_SphHarmListCAmpPhaseSequence(REAL8Vector *tVec, SphHarmListCAmpPhaseSequence *hLM);
void print_SphHarmListCAmpPhaseSequenceAndDynamics(REAL8Vector *tVec, SphHarmListCAmpPhaseSequence *hLM, SEOBSAdynamics *dyn_debug);
void print_SEOBPrecWaveform(SphHarmTimeSeries *hLM);
void print_SphHarmListCAmpPhaseSequence_debug(SphHarmListCAmpPhaseSequence *hLM);


typedef struct tagGSParams {
    REAL8 m1;                 /**< mass of companion 1 */
    REAL8 m2;
    REAL8 s1x;
    REAL8 s1y;
    REAL8 s1z;
    REAL8 s2x;
    REAL8 s2y;
    REAL8 s2z;
    REAL8 ecc;                 /**< eccentricity at start frequency */
    REAL8 f_min;              /**< start frequency */
    REAL8 p;                  /**< initial semi-latus of the orbit */
    REAL8 distance;
    REAL8 inc;
    INT is_only22;
    INT mode;
    REAL8 phiRef;
    REAL8 deltaT;             /**< sampling interval */
    REAL8 tStepBack;
    INT is_constsp;
    INT is_noringdown;
    INT is_prec;
    INT debug_id;
} InputParams;

void dump_All_to_hdf5(SEOBCoreOutputs* All, InputParams *p, CHAR *fname);

INT usage(const CHAR *program)
{
    INT a,c;
    print_err("usage: %s [options]\n", program);
    print_err("\t-h, --help\tprint this message and exit\n");
    print_err("\t-g, --geometric\tuse geometric units.\n");
    print_err("\t-A --constant-sampling\n\t\tuse constant sampling.\n");
    print_err("\t-F --no-ringdown\n\t\tremove ringdown part.\n");
    print_err("\t-j, --only22\tonly use l,|m| = (2, 2) mode\n");
    print_err("\t-v, --version\tVersion of code[1]\n\t\t0: SEOBNRv4PHM\n\t\t1: SEOBNREPHM\n");
    print_err("\t-P --prec\n\t\tuse prec code(in test)[0].\n");
    print_err("\t-6=DEBUG --debug\n\t\trun debug[0].\n");

    print_err("\t-m M1, --m1=M1\n\t\tcomponent mass 1  [%g]\n", DEFAULT_m1);
    print_err("\t-M M2, --m2=M2\n\t\tcomponent mass 2  [%g]\n", DEFAULT_m2);
    print_err("\t-x S1X, --spin1x=S1X\n\t\tspin chi z for component 1 [%g]\n", DEFAULT_sx);
    print_err("\t-y S1Y, --spin1y=S1Y\n\t\tspin chi z for component 1 [%g]\n", DEFAULT_sy);
    print_err("\t-z S1Z, --spin1z=S1Z\n\t\tspin chi z for component 1 [%g]\n", DEFAULT_sz);
    print_err("\t-X S2X, --spin2x=S2X\n\t\tspin chi z for component 2 [%g]\n", DEFAULT_sx);
    print_err("\t-Y S2Y, --spin2y=S2Y\n\t\tspin chi z for component 2 [%g]\n", DEFAULT_sy);
    print_err("\t-Z S2Z, --spin2z=S2Z\n\t\tspin chi z for component 2 [%g]\n", DEFAULT_sz);
    print_err("\t-e ECC, --eccentricity=ECC\n\t\tinitial eccentricity of orbit [%g]\n", DEFAULT_ecc);
    print_err("\t-s SEMI-LATUS, --semi-latus=SEMI-LATUS\n\t\tinitial semi-latus of orbit [%g]\n", DEFAULT_semilatus);

    print_err("\t-R SRATE, --sample-rate=SRATE\n\t\t sample rate [%g]\n", DEFAULT_deltat);
    print_err("\t-f FMIN, --f-min=FMIN           \n\t\tfrequency to start waveform [%g]\n", DEFAULT_f_min);
    print_err("\t-i INCLINATION, --inclination=INCLINATION           \n\t\t initial inclination angle to orbit plane in degrees[%g]\n", DEFAULT_inclination);
    print_err("\t-u PHIREF, --phiRef=PHIREF  \n\t\ttreference phase in degrees [%g]", DEFAULT_PHIREF);
    print_err("\t-D DISTANCE, --distance=DISTANCE  \n\t\tLum distance in Mpc of the source, if > 0[100]\n");

    print_err("\t-d D, --d=D  \n\t\tInitial seperation, if > 0, use given initial condition[%g]\n", -1);
    print_err("\t-r PR, --pr=PR  \n\t\tInitial radius momentum [%g]\n", 0);
    print_err("\t-p PPHI, --pphi=PPHI  \n\t\tinitial phi augular momentum, if > 0, use given initial condition[%g]\n", -1);
    print_err("\t-t PTHETA, --ptheta=PTHETA  \n\t\tInitial theta augular momentum, if > 0, use given initial condition[%g]\n", -1);

    print_err("\t-I RISCO, --risco=RISCO  \n\t\tA hyper parameter that controls stop condition when (e0>0)[%g]\n", 6.);

    print_err("\t-k KK, --KK=KK                 \n\t\tAdjustable Coefficient K.\n");
    print_err("\t-< DSS, --dss=dSS                 \n\t\tAdjustable Coefficient dSS.\n");
    print_err("\t-> DSO, --dso=dSO                \n\t\tAdjustable Coefficient dSO.\n");
    print_err("\t-T DTPEAK --dtPeak=DTPEAK             \n\t\tAdjustable Coefficient dtPeak.\n");

    print_err("\t-B TSTEPBACK --tstepback=TSTEPBACK             \n\t\tStep back time for re-evolving[200].\n");

    print_err("\t-C CONSERVETIME --conserve=CONSERVETIME\n\t\tConserve dynamics with time limit.\n");

    print_err("\t-E MODE --mode=MODE        \n\t\tThe output mode.\n");
    print_err("\t\t0:\t\tdefault output, print time[s] hPlus hCross\n");
    print_err("\t\t10l+|m|:\t\toutput (l,|m|) mode, print time[M], hlm_real, hlm_imag\n");
    print_err("\t\t-1:\t\toutput all modes, print time[M], 22, 21, 33, 44, 55\n");
    print_err("\t\t-2:\t\toutput dynamics, print time[M], x, y, z, vx, vy, vz, px, py, pz\n");
    print_err("\t-l LOG_LEVEL, --log-level=LOG_LEVEL\n\t\tDebug level for logging [0]\n");
    print_err("\t-o LOG_FILE, --log-file=LOG_FILE\n\t\tWill dump loggings to here.\n");
    print_err("\t-U DEBUG_FILE, --debug-file=DEBUG_FILE\n\t\tWill dump debug data to here.\n");
    print_err("\n\n");
    return CEV_SUCCESS;
}

InputParams parseargs(INT argc, CHAR **argv, HyperParams *hparams, CtrlParams *ctpms)
{
    InputParams p;
    INT use_geom = 0;
    p.m1 = DEFAULT_m1;
    p.m2 = DEFAULT_m2;
    p.s1x = DEFAULT_sx;
    p.s1y = DEFAULT_sy;
    p.s1z = DEFAULT_sz;
    p.s2x = DEFAULT_sx;
    p.s2y = DEFAULT_sy;
    p.s2z = DEFAULT_sz;
    p.ecc = DEFAULT_ecc;
    p.inc = DEFAULT_inclination * CST_PI / 180.;
    p.phiRef = DEFAULT_PHIREF * CST_PI / 180.;
    p.f_min = DEFAULT_f_min;
    p.distance = 100.;
    p.deltaT = 1./DEFAULT_deltat;
    p.is_only22 = 0;
    p.mode = 0;
    p.is_constsp = 0;
    p.is_noringdown = 0;
    p.is_prec = 0;
    p.debug_id = 0;
    hparams->d_ini = -1;
    hparams->pr_ini = 0;
    hparams->pphi_ini = 0;
    hparams->ptheta_ini = 0;
    hparams->flagTuning = 0;
    hparams->tStepBack = 200.;
    hparams->sl_p = DEFAULT_semilatus;
    hparams->x0 = cos(DEFAULT_inclination * CST_PI / 180.);
    SET_CODE_VERSION(1);
    SET_CONSERV(0, 0.0);
    set_PrecFlag(0);
    extern CHAR *EXT_optarg;
    extern INT EXT_optind;
    ctpms->level = 1;
    strncpy(ctpms->flog, "None", STR_COMM_SIZE);
    strncpy(ctpms->fdebug, "conserve.h5", STR_COMM_SIZE);
    OPTION long_options[] = {
        {"help", opt_no_argument, 0, 'h'},
        {"only22", opt_no_argument, 0, 'j'},
        {"geometric", opt_no_argument, 0, 'g'},
        {"constant-sampling", opt_no_argument, 0, 'A'},
        {"no-ringdown", opt_no_argument, 0, 'F'},
        {"prec", opt_required_argument, 0, 'P'},
        {"debug", opt_required_argument, 0, '6'},
        {"version", opt_required_argument, 0, 'v'},

        {"m1", opt_required_argument, 0, 'm'},
        {"m2", opt_required_argument, 0, 'M'},
        {"spin1x", opt_required_argument, 0, 'x'},
        {"spin1y", opt_required_argument, 0, 'y'},
        {"spin1z", opt_required_argument, 0, 'z'},
        {"spin2x", opt_required_argument, 0, 'X'},
        {"spin2y", opt_required_argument, 0, 'Y'},
        {"spin2z", opt_required_argument, 0, 'Z'},
        {"eccentricity", opt_required_argument, 0, 'e'},
        {"semi-latus", opt_required_argument, 0, 's'},

        {"sample-rate", opt_required_argument, 0, 'R'},
        {"f-min", opt_required_argument, 0, 'f'},
        {"inclination", opt_required_argument, 0, 'i'},
        {"phiRef", opt_required_argument, 0, 'u'},
        {"distance", opt_required_argument, 0, 'D'},

        {"d", opt_required_argument, 0, 'd'},
        {"pr", opt_required_argument, 0, 'r'},
        {"pphi", opt_required_argument, 0, 'p'},
        {"ptheta", opt_required_argument, 0, 't'},

        {"risco", opt_required_argument, 0, 'I'},

        {"KK", opt_required_argument , 0 ,'k'},
        {"dss", opt_required_argument , 0 ,'<'},
        {"dso", opt_required_argument , 0 ,'>'},
        {"dtPeak", opt_required_argument ,0 ,'T'},

        {"tstepback", opt_required_argument ,0 ,'B'},
        {"conserve", opt_required_argument ,0 ,'C'},

        {"mode", opt_required_argument, 0, 'E'},

        {"log-level", opt_required_argument , 0, 'l'},
        {"log-file", opt_required_argument , 0, 'o'},
        {"debug-file", opt_required_argument , 0, 'U'},
        {0, 0, 0, 0}
    };
    CHAR args[] =
    "h:j:g:A:v:P:m:M:x:y:z:X:Y:Z:e:s:R:f:i:u:D:d:r:p:t:I:k:<:>:T:B:C:E:l:o:U";
    while (1)
    {
        INT option_index = 0;
        INT c;
        c = getopt_long_only(argc, argv, args, long_options, &option_index);
        if (c == -1)
            break;
        switch (c)
        {
            case 0:
                if (long_options[option_index].flag)
                    break;
                else
                {
                    print_err("error parsing option %s with argument %s\n", long_options[option_index].name, EXT_optarg);
                    exit(1);
                }
            case 'h':
                usage(argv[0]);
                exit(0);
            case 'j':
                p.is_only22 = 1;
                break;
            case 'g':
                use_geom = 1;
                break;
            case 'A':
                p.is_constsp = 1;
                break;
            case 'F':
                p.is_noringdown = 1;
                break;
            case 'v':
                SET_CODE_VERSION(atoi(EXT_optarg));
                break;
            case 'P':
                // p.is_prec = atoi(EXT_optarg);
                set_PrecFlag(atoi(EXT_optarg));
                // print_debug("prec flag = %d\n", get_PrecFlag());
                break;
            case '6':
                p.debug_id = atoi(EXT_optarg);
                break;
            case 'm':
                p.m1 = atof(EXT_optarg);
                break;
            case 'M':
                p.m2 = atof(EXT_optarg);
                break;
            case 'x':
                p.s1x = atof(EXT_optarg);
                break;
            case 'y':
                p.s1y = atof(EXT_optarg);
                break;
            case 'z':
                p.s1z = atof(EXT_optarg);
                break;
            case 'X':
                p.s2x = atof(EXT_optarg);
                break;
            case 'Y':
                p.s2y = atof(EXT_optarg);
                break;
            case 'Z':
                p.s2z = atof(EXT_optarg);
                break;
            case 'e':
                p.ecc = atof(EXT_optarg);
                break;
            case 's':
                hparams->sl_p = atof(EXT_optarg);
                break;
            case 'R':
                p.deltaT = 1./atof(EXT_optarg);
                break;
            case 'f':
                p.f_min = atof(EXT_optarg);
                break;
            case 'i':
                p.inc = atof(EXT_optarg) * CST_PI/180;
                hparams->x0 = cos(p.inc);
                break;
            case 'u':
                p.phiRef = atof(EXT_optarg) * CST_PI/180;
                break;
            case 'D':
                p.distance = atof(EXT_optarg);
                break;
            case 'd':
                hparams->d_ini = atof(EXT_optarg);
                break;
            case 'r':
                hparams->pr_ini = atof(EXT_optarg);
                break;
            case 'p':
                hparams->pphi_ini = atof(EXT_optarg);
                break;
            case 't':
                hparams->ptheta_ini = atof(EXT_optarg);
                break;
            case 'I':
            {
                REAL8 risco = atof(EXT_optarg);
                SET_RISCO(risco);
                break;
            }
            case 'k':
                hparams->KK = atof(EXT_optarg);
                hparams->flagTuning = 1;
                break;
            case '<':
                hparams->dSS = atof(EXT_optarg);
                hparams->flagTuning = 1;
                break;
            case '>':
                hparams->dSO = atof(EXT_optarg);
                hparams->flagTuning = 1;
                break;
            case 'T':
                hparams->dtPeak = atof(EXT_optarg);
                hparams->flagTuning = 1;
                break;
            case 'B':
                hparams->tStepBack = atof(EXT_optarg);
                break;
            case 'C':
                // print_debug("here\n");
                SET_CONSERV(1, atof(EXT_optarg));
                break;
            case 'E':
                p.mode = atoi(EXT_optarg);
                break;
            case 'l':
                ctpms->level = atoi(EXT_optarg);
                // print_debug("ctpms.level = %d\n", ctpms->level);
                break;
            case 'O':
                strncpy(ctpms->flog, EXT_optarg, STR_COMM_SIZE);
                break;
            case 'U':
                strncpy(ctpms->fdebug, EXT_optarg, STR_COMM_SIZE);
                break;
            default:
                print_err("unknown error while parsing options\n");
                exit(1);
        }
    }
    if (EXT_optind < argc)
    {
        print_err("extraneous command line arguments:\n");
        while (EXT_optind < argc)
            print_err("%s\n", argv[EXT_optind++]);
        exit(1);
    }
    if (use_geom)
    {
        REAL8 mTScaled = (p.m1+p.m2) * CST_MTSUN_SI;
        p.f_min = p.f_min / mTScaled;
    }
    if (ctpms->level >= 10)
    {
        SET_INPUT_DEBUG_FLAG(ctpms->level - 10);
    }
    return p;
}

enum {nprefix = 2};


INT main(INT argc, CHAR **argv)
{
#if 0
    return debug();
#endif
    clock_t tstart, tend;
    tstart = clock();
    InputParams p;
    INT status;
    CtrlParams ctpms;
    HyperParams hparams;
    SphHarmListCAmpPhaseSequence *hLM = NULL;
    REAL8Vector *tVec = NULL;
    REAL8TimeSeries *hplus = NULL;
    REAL8TimeSeries *hcross = NULL;
    SEOBCoreOutputs *All = (SEOBCoreOutputs *)MYCalloc(1, sizeof(SEOBCoreOutputs));
    p = parseargs(argc, argv, &hparams, &ctpms);
    if (strcmp(ctpms.flog, "None")!=0)
    {
        LOG_SetPrintLogPlaceFlag(1);
    }
    LOG_SetPrintDebugLogFlag(ctpms.level);
    LOG_Init(ctpms.flog, 40960);
    // print_debug("g_ulPrintDebugLogFlag = %zu\n", g_ulPrintDebugLogFlag);
    PRINT_LOG_INFO(LOG_DEBUG, "CMD:\n%s --m1=%g --m2=%g --spin1x=%g --spin1y=%g --spin1z=%g --spin2x=%g --spin2y=%g --spin2z=%g --eccentricity=%g --delta-t=%g --f-min=%g --inclination=%g --phiRef=%g\n", 
        argv[0], p.m1, p.m2, p.s1x, p.s1y, p.s1z, p.s2x, p.s2y, p.s2z, p.ecc, p.deltaT, p.f_min, p.inc * 180/CST_PI, p.phiRef*180/CST_PI);
    PRINT_LOG_INFO(LOG_INFO, "CodeVersionFlag = %d\n", CODE_VERSION);
    // print_debug("prec flag = %d\n", ctpms.level);
    if (p.debug_id)
    {
        choose_debug(p.debug_id, p.m1, p.m2, 
            p.s1x, p.s1y, p.s1z, 
            p.s2x, p.s2y, p.s2z, 
            p.phiRef, p.distance, p.ecc, p.f_min,
            p.deltaT, p.inc, &hparams);
        return 0;
    }
    if (get_PrecFlag() != 0) {
        PRINT_LOG_INFO(LOG_DEBUG, "prec code");
        SEOBPrecCoreOutputs *All_prec = (SEOBPrecCoreOutputs *)MYCalloc(1, sizeof(SEOBPrecCoreOutputs));
        status = evolve_prec(p.m1, p.m2, 
            p.s1x, p.s1y, p.s1z, 
            p.s2x, p.s2y, p.s2z, 
            p.phiRef, p.distance, p.ecc, p.f_min,
            p.deltaT, p.inc, &hparams, All_prec);
        if (status != CEV_SUCCESS)
            goto EXIT;
        switch(p.mode)
        {
            case -2:
                PRINT_LOG_INFO(LOG_DEBUG, "Dump dynamics");
                print_SEOBPrecdynamics(All_prec->dyn);
                break;
            case -3:
                PRINT_LOG_INFO(LOG_DEBUG, "Dump Plm");
                print_SphHarmListCAmpPhaseSequence_debug(All_prec->Plm);
                break;
            case -1:
            default:
                PRINT_LOG_INFO(LOG_DEBUG, "Dump waveform");
                print_SEOBPrecWaveform(All_prec->hLM);
        }
        STRUCTFREE(All_prec, SEOBPrecCoreOutputs);
        goto EXIT;
    } else if (CONSERVE_FLAG!=0) {
        print_debug("it is conserved\n");
        status = evolve_conserv(p.m1, p.m2, 
            p.s1x, p.s1y, p.s1z, 
            p.s2x, p.s2y, p.s2z, 
            p.phiRef, p.distance, p.ecc, p.f_min,
            p.deltaT, p.inc, &hparams, All);
        if (status != CEV_SUCCESS)
            goto EXIT;
        dump_All_to_hdf5(All, &p, ctpms.fdebug);
    } else if (sqrt(p.s1x*p.s1x + p.s1y*p.s1y) < 1e-5 && sqrt(p.s2x*p.s2x + p.s2y*p.s2y) < 1e-5) {
        /* default */
        SEOBSAdynamics *dyn_debug = NULL;
        status = evolve_SA(p.m1, p.m2, p.s1z, p.s2z, p.ecc, p.f_min, p.deltaT, &hparams, &tVec, &hLM, p.is_noringdown, p.mode==-2, &dyn_debug);
        if (status != CEV_SUCCESS)
            goto EXIT;
        if (p.mode==-2)
            print_SphHarmListCAmpPhaseSequenceAndDynamics(tVec, hLM, dyn_debug);
        else
            print_SphHarmListCAmpPhaseSequence(tVec, hLM);
        STRUCTFREE(dyn_debug, SEOBSAdynamics);
    } else {
        if (p.is_constsp)
            status = evolve(p.m1, p.m2, 
                p.s1x, p.s1y, p.s1z, 
                p.s2x, p.s2y, p.s2z, 
                p.phiRef, p.distance, p.ecc, p.f_min, 
                p.deltaT, p.inc, p.is_only22,
                &hparams, &hplus, &hcross, All);
        else
            status = evolve_adaptive(p.m1, p.m2, 
                p.s1x, p.s1y, p.s1z, 
                p.s2x, p.s2y, p.s2z, 
                p.phiRef, p.distance, p.ecc, p.f_min, 
                p.deltaT, p.inc, p.is_only22,
                &hparams, &hplus, &hcross, All);
        if (status != CEV_SUCCESS)
            goto EXIT;
        // Output
        switch(p.mode)
        {
            case 0:
                // default
                print_waveform(hplus, hcross);
                break;
            case -2:
                print_SEOBdynamics(All->dyn);
                break;
            case -3:
                PRINT_LOG_INFO(LOG_DEBUG, "Dump Plm");
                print_SphHarmListCAmpPhaseSequence_debug(All->Plm);
                break;
            default:
                print_SphHarmTimeSeries(All->hLM, p.mode);
        }
    }

EXIT:
    STRUCTFREE(hplus, REAL8TimeSeries);
    STRUCTFREE(hcross, REAL8TimeSeries);
    STRUCTFREE(hLM, SphHarmListCAmpPhaseSequence);
    STRUCTFREE(tVec, REAL8Vector);
    STRUCTFREE(All, SEOBCoreOutputs);
    CheckMemoryLeak();
    tend = clock();
    PRINT_LOG_INFO(LOG_INFO, "Time Cost: %fs\n", ((REAL8)(tend-tstart)/CLOCKS_PER_SEC));  
    return 0;
}

void print_SphHarmListCAmpPhaseSequence_debug(SphHarmListCAmpPhaseSequence *hLM)
{
    INT i, length;
    CAmpPhaseSequence *h22 = SphHarmListCAmpPhaseSequence_GetMode(hLM, 2, 2)->campphase;
    CAmpPhaseSequence *h21 = SphHarmListCAmpPhaseSequence_GetMode(hLM, 2, 1)->campphase;
    CAmpPhaseSequence *h33 = SphHarmListCAmpPhaseSequence_GetMode(hLM, 3, 3)->campphase;
    CAmpPhaseSequence *h44 = SphHarmListCAmpPhaseSequence_GetMode(hLM, 4, 4)->campphase;
    CAmpPhaseSequence *h55 = SphHarmListCAmpPhaseSequence_GetMode(hLM, 5, 5)->campphase;
    length = h22->xdata->length;
    print_out("#time\t#h22r\t#h22i\t#h21r\t#h21i\t#h33r\t#h33i\t#h44r\t#h44i\th55r\th55i\n");
    for(i=0;i<length;i++)
    {
        print_out("%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\n", h22->xdata->data[i],
            h22->camp_real->data[i]*cos(h22->phase->data[i]), h22->camp_real->data[i]*sin(h22->phase->data[i]),
            h21->camp_real->data[i]*cos(h21->phase->data[i]), h21->camp_real->data[i]*sin(h21->phase->data[i]),
            h33->camp_real->data[i]*cos(h33->phase->data[i]), h33->camp_real->data[i]*sin(h33->phase->data[i]),
            h44->camp_real->data[i]*cos(h44->phase->data[i]), h44->camp_real->data[i]*sin(h44->phase->data[i]),
            h55->camp_real->data[i]*cos(h55->phase->data[i]), h55->camp_real->data[i]*sin(h55->phase->data[i]));
    }
    return;
}

void print_SphHarmTimeSeries(SphHarmTimeSeries *hlm, INT mode)
{
    INT i;
    if (mode < 0)
    {
        COMPLEX16TimeSeries *h22 = XLALSphHarmTimeSeriesGetMode(hlm, 2, 2);
        COMPLEX16TimeSeries *h21 = XLALSphHarmTimeSeriesGetMode(hlm, 2, 1);
        COMPLEX16TimeSeries *h33 = XLALSphHarmTimeSeriesGetMode(hlm, 3, 3);
        COMPLEX16TimeSeries *h44 = XLALSphHarmTimeSeriesGetMode(hlm, 4, 4);
        COMPLEX16TimeSeries *h55 = XLALSphHarmTimeSeriesGetMode(hlm, 5, 5);
        REAL8 deltaT = h22->deltaT;
        INT length = h22->data->length;
        print_out("#time\t#h22r\t#h22i\t#h21r\t#h21i\t#h33r\t#h33i\t#h44r\t#h44i\t#h55r\t#h55i\n");
        for(i=0;i<length;i++)
        {
            print_out("%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\n", i*deltaT - hlm->tAttach,
                creal(h22->data->data[i]), cimag(h22->data->data[i]),
                creal(h21->data->data[i]), cimag(h21->data->data[i]),
                creal(h33->data->data[i]), cimag(h33->data->data[i]),
                creal(h44->data->data[i]), cimag(h44->data->data[i]),
                creal(h55->data->data[i]), cimag(h55->data->data[i]));

        }
        return;
    }
    UINT l;
    INT m;
    m = mode % 10;
    l = (mode - m)/10;
    COMPLEX16TimeSeries *hIlmTS = XLALSphHarmTimeSeriesGetMode(hlm, l, m);
    if (!hIlmTS)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "Invalid mode tag %d", mode);
        return;
    }
    REAL8 deltaT = hIlmTS->deltaT;
    INT length = hIlmTS->data->length;
    print_out("#time\t#hreal\t#himag\n");
    for(i=0; i<length; i++)
    {
        print_out("%.18e\t%.18e\t%.18e\n", i*deltaT, 
            creal(hIlmTS->data->data[i]), cimag(hIlmTS->data->data[i]));
    }
    return;
}

void print_SEOBdynamics(SEOBdynamics *dyn)
{
    INT i, length;
    length = dyn->length;
    print_out("#time\t#x\t#y\t#z\t#vx\t#vy\t#vz\t#px\t#py\t#pz\t#ham\t#flux\t#Fr\t#Ff\t#polarprDot\n");
    for(i=0;i<length;i++)
    {
        print_out("%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\n",
            dyn->tVec[i], dyn->posVecx[i], dyn->posVecy[i], dyn->posVecz[i],
            dyn->momVecx[i], dyn->momVecy[i], dyn->momVecz[i],
            dyn->s1Vecx[i], dyn->s1Vecy[i], dyn->s1Vecz[i],
            dyn->s2Vecx[i], dyn->s2Vecy[i], dyn->s2Vecz[i],
            dyn->omegaVec[i], dyn->hamVec[i],
            dyn->fluxVec[i], dyn->polarprDotVec[i]);
    }
    return;
}

void print_SEOBPrecdynamics(SEOBPrecdynamics *dyn)
{
    INT i, length;
    length = dyn->length;
    print_out("#time\t#x\t#y\t#z\t#pTx\t#pTy\t#pTz\ts1x\ts1y\ts1z\ts2x\ts2y\ts2z\tomega\tHam\tflux\tprTDot\n");
    for(i=0;i<length;i++)
    {
        print_out("%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\n",
            dyn->tVec[i], dyn->posVecx[i], dyn->posVecy[i], dyn->posVecz[i],
            dyn->momTVecx[i], dyn->momTVecy[i], dyn->momTVecz[i],
            dyn->s1Vecx[i], dyn->s1Vecy[i], dyn->s1Vecz[i],
            dyn->s2Vecx[i], dyn->s2Vecy[i], dyn->s2Vecz[i],
            dyn->omegaVec[i], dyn->HamVec[i], dyn->fluxVec[i],
            dyn->prTDotVec[i]);
    }
    return;
}

void print_SEOBPrecWaveform(SphHarmTimeSeries *hLM)
{
    INT i;
    COMPLEX16TimeSeries *h22 = XLALSphHarmTimeSeriesGetMode(hLM, 2, 2);
    COMPLEX16TimeSeries *h21 = XLALSphHarmTimeSeriesGetMode(hLM, 2, 1);
    COMPLEX16TimeSeries *h33 = XLALSphHarmTimeSeriesGetMode(hLM, 3, 3);
    COMPLEX16TimeSeries *h44 = XLALSphHarmTimeSeriesGetMode(hLM, 4, 4);
    COMPLEX16TimeSeries *h55 = XLALSphHarmTimeSeriesGetMode(hLM, 5, 5);
    REAL8 deltaT = h22->deltaT;
    INT length = h22->data->length;
    print_out("#time\th22r\th22i\th21r\th21i\th33r\th33i\th44r\th44i\th55r\th55i\n");
    for(i=0;i<length;i++)
    {
        print_out("%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\n",
            deltaT*i - hLM->tAttach, 
            creal(h22->data->data[i]), cimag(h22->data->data[i]),
            creal(h21->data->data[i]), cimag(h21->data->data[i]),
            creal(h33->data->data[i]), cimag(h33->data->data[i]),
            creal(h44->data->data[i]), cimag(h44->data->data[i]),
            creal(h55->data->data[i]), cimag(h55->data->data[i]));
    }
    return;
}

void print_waveform(REAL8TimeSeries *hp, REAL8TimeSeries *hc)
{
    if (hp->data->length != hc->data->length || 
        hp->deltaT != hc->deltaT)
    {
        PRINT_LOG_INFO(LOG_CRITICAL, "hPlus and hCorss has different shape.\n");
        return;
    }
    INT i, length = hp->data->length;
    REAL8 deltaT = hp->deltaT, t0 = hp->epoch;
    print_out("#time\t#hplus\t#hcross\n");
    for (i=0; i<length; i++)
    {
        print_out("%.18e\t%.18e\t%.18e\n", 
            t0 + i*deltaT, hp->data->data[i], 
            hc->data->data[i]);
    }
    return;
}


void dump_All_to_hdf5(SEOBCoreOutputs* All, InputParams *p, CHAR *fname)
{
    UINT length = All->dyn->length;
    DumpREAL8ArrayTohdf5(fname, "dynamics", All->dyn->array, 1);
    REAL8Array *hLMarr = NULL;
    REAL8Vector *tVec = NULL;
    REAL8 *hr, *hi;
    INT i;
    hLMarr = CreateREAL8Array(2, 2, length);
    tVec = CreateREAL8Vector(length);
    for (i=0; i<length; i++)
    {
        tVec->data[i] = All->dyn->tVec[i];
    }
    DumpREAL8VectorTohdf5(fname, "time", tVec, 0);
    writeREAL8Tohdf5(fname, "q", p->m1/p->m2, 0);
    REAL8 s1[3] = {p->s1x, p->s1y, p->s1z};
    REAL8 s2[3] = {p->s2x, p->s2y, p->s2z};
    REAL8Vector s1Vec, s2Vec;
    s1Vec.length = s2Vec.length = 3;
    s1Vec.data = s1;
    s2Vec.data = s2;
    DumpREAL8VectorTohdf5(fname, "s1Vec", &s1Vec, 0);
    DumpREAL8VectorTohdf5(fname, "s2Vec", &s2Vec, 0);
    hr = hLMarr->data;
    hi = hr + length;
    COMPLEX16TimeSeries *h22 = XLALSphHarmTimeSeriesGetMode(All->hLM, 2, 2);
    COMPLEX16TimeSeries *h21 = XLALSphHarmTimeSeriesGetMode(All->hLM, 2, 1);
    COMPLEX16TimeSeries *h33 = XLALSphHarmTimeSeriesGetMode(All->hLM, 3, 3);
    COMPLEX16TimeSeries *h44 = XLALSphHarmTimeSeriesGetMode(All->hLM, 4, 4);
    COMPLEX16TimeSeries *h55 = XLALSphHarmTimeSeriesGetMode(All->hLM, 5, 5);
    // 2,2 mode
    for (i=0; i<length; i++)
    {
        hr[i] = creal(h22->data->data[i]);
        hi[i] = cimag(h22->data->data[i]);
    }
    DumpREAL8ArrayTohdf5(fname, "h22", hLMarr, 0);

    // 2,1 mode
    for (i=0; i<length; i++)
    {
        hr[i] = creal(h21->data->data[i]);
        hi[i] = cimag(h21->data->data[i]);
    }
    DumpREAL8ArrayTohdf5(fname, "h21", hLMarr, 0);

    // 3,3 mode
    for (i=0; i<length; i++)
    {
        hr[i] = creal(h33->data->data[i]);
        hi[i] = cimag(h33->data->data[i]);
    }
    DumpREAL8ArrayTohdf5(fname, "h33", hLMarr, 0);

    // 4,4 mode
    for (i=0; i<length; i++)
    {
        hr[i] = creal(h44->data->data[i]);
        hi[i] = cimag(h44->data->data[i]);
    }
    DumpREAL8ArrayTohdf5(fname, "h44", hLMarr, 0);

    STRUCTFREE(hLMarr, REAL8Array);
    STRUCTFREE(tVec, REAL8Vector);
    return;
}

void print_SphHarmListCAmpPhaseSequenceAndDynamics(REAL8Vector *tVec, SphHarmListCAmpPhaseSequence *hLM, SEOBSAdynamics *dyn_debug)
{
    int i;
    size_t length = tVec->length, len_dyn = dyn_debug->length;
    REAL8 r, dr, omega;
    // print_debug("len(hLM) = %d, len(dyn) = %d\n", tVec->length, dyn_debug->length);
    CAmpPhaseSequence *h22 = SphHarmListCAmpPhaseSequence_GetMode(hLM, 2, 2)->campphase;
    CAmpPhaseSequence *h21 = SphHarmListCAmpPhaseSequence_GetMode(hLM, 2, 1)->campphase;
    CAmpPhaseSequence *h33 = SphHarmListCAmpPhaseSequence_GetMode(hLM, 3, 3)->campphase;
    CAmpPhaseSequence *h44 = SphHarmListCAmpPhaseSequence_GetMode(hLM, 4, 4)->campphase;
    CAmpPhaseSequence *h55 = SphHarmListCAmpPhaseSequence_GetMode(hLM, 5, 5)->campphase;
    print_out("#time\t#h22r\t#h22i\t#h21r\t#h21i\t#h33r\t#h33i\t#h44r\t#h44i\th55r\th55i\tr\trDot\tomega\n");
    for(i=0;i<length;i++)
    {
        if (i < len_dyn)
        {
            r = dyn_debug->rVec[i];
            dr = dyn_debug->drVec[i];
            omega = dyn_debug->dphiVec[i];
        } else {
            r = dyn_debug->rVec[len_dyn-1];
            dr = dyn_debug->drVec[len_dyn-1];
            omega = dyn_debug->dphiVec[len_dyn-1];
        }
        print_out("%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\n", tVec->data[i],
            h22->camp_real->data[i]*cos(h22->phase->data[i]), h22->camp_real->data[i]*sin(h22->phase->data[i]),
            h21->camp_real->data[i]*cos(h21->phase->data[i]), h21->camp_real->data[i]*sin(h21->phase->data[i]),
            h33->camp_real->data[i]*cos(h33->phase->data[i]), h33->camp_real->data[i]*sin(h33->phase->data[i]),
            h44->camp_real->data[i]*cos(h44->phase->data[i]), h44->camp_real->data[i]*sin(h44->phase->data[i]),
            h55->camp_real->data[i]*cos(h55->phase->data[i]), h55->camp_real->data[i]*sin(h55->phase->data[i]),
            r, dr, omega);
    }
    return;
}

void print_SphHarmListCAmpPhaseSequence(REAL8Vector *tVec, SphHarmListCAmpPhaseSequence *hLM)
{
    int i;
    size_t length = tVec->length;
    CAmpPhaseSequence *h22 = SphHarmListCAmpPhaseSequence_GetMode(hLM, 2, 2)->campphase;
    CAmpPhaseSequence *h21 = SphHarmListCAmpPhaseSequence_GetMode(hLM, 2, 1)->campphase;
    CAmpPhaseSequence *h33 = SphHarmListCAmpPhaseSequence_GetMode(hLM, 3, 3)->campphase;
    CAmpPhaseSequence *h44 = SphHarmListCAmpPhaseSequence_GetMode(hLM, 4, 4)->campphase;
    CAmpPhaseSequence *h55 = SphHarmListCAmpPhaseSequence_GetMode(hLM, 5, 5)->campphase;
    print_out("#time\t#h22r\t#h22i\t#h21r\t#h21i\t#h33r\t#h33i\t#h44r\t#h44i\th55r\th55i\n");
    for(i=0;i<length;i++)
    {
        print_out("%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\n", tVec->data[i],
            h22->camp_real->data[i]*cos(h22->phase->data[i]), h22->camp_real->data[i]*sin(h22->phase->data[i]),
            h21->camp_real->data[i]*cos(h21->phase->data[i]), h21->camp_real->data[i]*sin(h21->phase->data[i]),
            h33->camp_real->data[i]*cos(h33->phase->data[i]), h33->camp_real->data[i]*sin(h33->phase->data[i]),
            h44->camp_real->data[i]*cos(h44->phase->data[i]), h44->camp_real->data[i]*sin(h44->phase->data[i]),
            h55->camp_real->data[i]*cos(h55->phase->data[i]), h55->camp_real->data[i]*sin(h55->phase->data[i]));

    }
    return;
}








int debug()
{
#if 1
    {
        clock_t t0, te;
        REAL8 m1 = 5.;
        REAL8 m2 = 1.;
        REAL8 chi1 = 0.9;
        REAL8 chi2 = -0.8;
        REAL8 r = 41.65993886853366;
        REAL8 prT = -0.01;
        REAL8 pphi = 6.639233199936644;
        SpinEOBHSACoeffs hsaCoeffs;
        REAL8 tmp, H;
        int nmax = 2000000;
        CalculateSpinEOBHSACoeffs(m1, m2, chi1, chi2, &hsaCoeffs);
        t0 = clock();
        for (int i=0; i<nmax; i++)
            H = EOBSAHamiltonian(r, prT, pphi, &hsaCoeffs, &tmp);
        te = clock();
        print_debug("Time Cost: %.16e\n", ((REAL8)(te-t0)/CLOCKS_PER_SEC) );
        print_debug("H_sa = %.16e\n", H);

        REAL8 eta = m1*m2 / (m1+m2) / (m1+m2);
        REAL8Vector xVec, pVec, s1Vec, s2Vec, sigmaKerr, sigmaStar;
        REAL8 xData[3] = {41.659938868533658,0,0};
        REAL8 pData[3] = {-0.01, 1.5936732938778625e-1, -2.6799608562176783e-12};
        REAL8 s1Data[3] = {0,0,0};
        REAL8 s2Data[3] = {0,0,0};
        REAL8 sKData[3] = {0,0,0};
        REAL8 sSData[3] = {0,0,0};
        xVec.length = pVec.length = s1Vec.length = s2Vec.length = sigmaKerr.length = sigmaStar.length;
        xVec.data = xData;
        pVec.data = pData;
        s1Vec.data = s1Data;
        s2Vec.data = s2Data;
        sigmaKerr.data = sKData;
        sigmaStar.data = sSData;
        s1Data[2] = chi1 * m1 * m1 / (m1+m2) / (m1+m2);
        s2Data[2] = chi2 * m2 * m2 / (m1+m2) / (m1+m2);
        sSData[2] = (m2/m1)*s1Data[2] + (m1/m2)*s2Data[2];
        sKData[2] = s1Data[2] + s2Data[2];
        REAL8 a = fabs(sKData[2]);
        REAL8 S_con = sKData[2] / (1. - 2.*eta);
        SpinEOBHCoeffs coeffs;
        HyperParams hParams;
        hParams.flagTuning = 0;
        XLALSimIMRCalculateSpinPrecEOBHCoeffs_v2(&coeffs, eta, a, S_con, 4, &hParams);
        t0 = clock();
        for (int i=0; i<nmax; i++)
            H = EOBHamiltonian(eta, &xVec, &pVec, &s1Vec, &s2Vec, &sigmaKerr, &sigmaStar, 1, &coeffs);
        te = clock();
        print_debug("Time Cost: %.16e\n", ((REAL8)(te-t0)/CLOCKS_PER_SEC) );
        print_debug("H = %.16e\n", H);
    }
    return 0;
#elif 0
{
    clock_t t0, te;
    SpinEOBParams *core = NULL;
    REAL8 m1 = 5.;
    REAL8 m2 = 1.;
    REAL8 chi1 = 0.9;
    REAL8 chi2 = -0.8;
    REAL8 r = 41.65993886853366;
    REAL8 prT = -0.01;
    REAL8 pphi = 6.639233199936644;
    HyperParams hParams;
    hParams.flagTuning = 0;
    core = CreateSpinEOBParams(m1, m2, 0, 0, chi1, 0, 0, chi2, 0.2, &hParams);
    REAL8 values[4], dvalues[4];
    values[0] = r;
    values[1] = 0;
    values[2] = prT;
    values[3] = pphi;
    int nmax = 100000;
    t0 = clock();
    for (int i=0; i<nmax; i++)
        XLALSpinAlignedHcapDerivative_SA(0, values, dvalues, (void*)core);
    te = clock();
    print_debug("Time Cost SA: %.16e\n", ((REAL8)(te-t0)/CLOCKS_PER_SEC) );
    print_debug("dvalues = (%.16e, %.16e, %.16e, %.16e)\n\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3]);

    t0 = clock();
    for (int i=0; i<nmax; i++)
        XLALSpinAlignedHcapDerivative(0, values, dvalues, (void*)core);
    te = clock();
    print_debug("Time Cost Old: %.16e\n", ((REAL8)(te-t0)/CLOCKS_PER_SEC) );
    print_debug("dvalues = (%.16e, %.16e, %.16e, %.16e)\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3]);

}
return 0;
#elif 0
{
    clock_t t0, te;
    REAL8 cFr, cFf, r, dr, prDot, e0;
    e0 = 0.1;
    REAL8 m1 = 5.;
    REAL8 m2 = 1.;
    REAL8 chi1 = 0.9;
    REAL8 chi2 = -0.8;
    REAL8 eta = m1*m2/(m1+m2)/(m1+m2);
SpinEOBParams *core = NULL;
HyperParams hParams;
hParams.flagTuning = 0;
core = CreateSpinEOBParams(m1, m2, 0, 0, chi1, 0, 0, chi2, 0.2, &hParams);
    r = 4.1659938868533658e+01, dr = -9.6550207519474334e-03, prDot = -2.5203299166118250e-06;
    EccCorrectionCoeffs eccCoeffs;
    CalculateEccCorrectionCoeffs(eta, chi1, chi2, &eccCoeffs);
    int nmax = 100000;
    t0 = clock();
    for (int i=0; i<nmax; i++)
        CalculateEccCorrectionToFluxV3X(r, dr, prDot, &cFr, &cFf, e0, core->eccCoeffs);
    te = clock();
    print_debug("Time Cost SA: %.16e\n", ((REAL8)(te-t0)/CLOCKS_PER_SEC) );
    // print_debug("cFr = %.16e, cFf = %.16e\n", cFr, cFf);
print_err("\n");
    r = 4.1659938868533658e+01, dr = -9.6550207519474299e-03, prDot = -2.5203299166118250e-06;
    t0 = clock();
    for (int i=0; i<nmax; i++)
        CalculateEccCorrectionToFluxV3(eta, chi1, chi2, r, dr, prDot, &cFr, &cFf, e0);
    te = clock();
    print_debug("Time Cost Old: %.16e\n", ((REAL8)(te-t0)/CLOCKS_PER_SEC) );
    // print_debug("cFr = %.16e, cFf = %.16e\n", cFr, cFf);
}
return 0;
#endif
}