/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_PSTRUCTS__
#define __INCLUDE_PSTRUCTS__
#include "myUtils.h"
#define VERSION 1

typedef enum tagflagSEOBNRv4P_Zframe {
  FLAG_SEOBNRv4P_ZFRAME_L = 0, /**< set Z axis of the P-frame along L */
  FLAG_SEOBNRv4P_ZFRAME_LN = 1 /**< set Z axis of the P-frame along LN */
} flagSEOBNRv4P_Zframe;


// Other params
typedef struct tagHyperParams
{
    REAL8 sl_p; // initial semi-latus
    REAL8 x0; // initial inclination angle
    REAL8 KK;
    REAL8 dSO;
    REAL8 dSS;
    REAL8 dtPeak;
    INT flagTuning;
    INT flagZframe;
    REAL8 d_ini;
    REAL8 pr_ini;
    REAL8 pphi_ini;
    REAL8 ptheta_ini;
    REAL8 tStepBack;
    REAL8 inEPS_REL;
    REAL8 inEPS_ABS;
    INT is_coframe;
    REAL8  Mf_min;
    REAL8  Mf_ref;
    INT    zero_dyncoaphase;
}HyperParams;

typedef struct {
    REAL8 FrPN1_f21;
    REAL8 FrPN1_f20;
    REAL8 FrPN1_f12;
    REAL8 FrPN1_f11;
    REAL8 FrPN1_f10;
    REAL8 FrPN1_f03;
    REAL8 FrPN1_f02;
    REAL8 FrPN1_f01;

    REAL8 FfPN1_f21;
    REAL8 FfPN1_f20;
    REAL8 FfPN1_f12;
    REAL8 FfPN1_f11;
    REAL8 FfPN1_f10;
    REAL8 FfPN1_f03;
    REAL8 FfPN1_f02;
    REAL8 FfPN1_f01;

    REAL8 FrPN32_f21;
    REAL8 FrPN32_f20;
    REAL8 FrPN32_f11;
    REAL8 FrPN32_f10;
    REAL8 FrPN32_f02;
    REAL8 FrPN32_f01;

    REAL8 FfPN32_f30;
    REAL8 FfPN32_f20;
    REAL8 FfPN32_f11;
    REAL8 FfPN32_f10;
    REAL8 FfPN32_f02;
    REAL8 FfPN32_f01;

    REAL8 FrPN2Topfm2_f20;
    REAL8 FrPN2Topfm2_f10;
    REAL8 FrPN2Topfm2_f00;
    REAL8 FrPN2Topfm1_f20;
    REAL8 FrPN2Topfm1_f10;
    REAL8 FrPN2Topfm1_f00;
    REAL8 FrPN2Topf0_f30;
    REAL8 FrPN2Topf0_f21;
    REAL8 FrPN2Topf0_f20;
    REAL8 FrPN2Topf0_f12;
    REAL8 FrPN2Topf0_f11;
    REAL8 FrPN2Topf0_f10;
    REAL8 FrPN2Topf0_f03;
    REAL8 FrPN2Topf0_f02;
    REAL8 FrPN2Topf0_f01;
    REAL8 FrPN2Topf0_f00;

    REAL8 FfPN2Topfm2_f23;
    REAL8 FfPN2Topfm2_f14;
    REAL8 FfPN2Topfm2_f05;
    REAL8 FfPN2Topfm1_f22;
    REAL8 FfPN2Topfm1_f13;
    REAL8 FfPN2Topfm1_f04;
    REAL8 FfPN2Topf0_f30;
    REAL8 FfPN2Topf0_f21;
    REAL8 FfPN2Topf0_f20;
    REAL8 FfPN2Topf0_f12;
    REAL8 FfPN2Topf0_f11;
    REAL8 FfPN2Topf0_f10;
    REAL8 FfPN2Topf0_f03;
    REAL8 FfPN2Topf0_f02;
    REAL8 FfPN2Topf0_f01;
}EccCorrectionCoeffs;


typedef struct tagFacWaveformCoeffs
{
    REAL8 delta22vh3;
    REAL8 delta22vh6;
    REAL8 delta22vh6S;
    REAL8 delta22v8;
    REAL8 delta22v8S;
    REAL8 delta22vh9;
    REAL8 delta22v5;
    REAL8 delta22v6;
    REAL8 delta22v6S;

    REAL8 rho22v2;
    REAL8 rho22v3;
    REAL8 rho22v3S;
    REAL8 rho22v4;
    REAL8 rho22v4S;
    REAL8 rho22v5;
    REAL8 rho22v5S;
    REAL8 rho22v6;
    REAL8 rho22v6S;
    REAL8 rho22v6l;
    REAL8 rho22v7;
    REAL8 rho22v7S;
    REAL8 rho22v8;
    REAL8 rho22v8S;
    REAL8 rho22v8l;
    REAL8 rho22v10;
    REAL8 rho22v10l;

    REAL8 delta21vh3;
    REAL8 delta21vh6;
    REAL8 delta21vh6S;
    REAL8 delta21vh7;
    REAL8 delta21vh7S;
    REAL8 delta21vh9;
    REAL8 delta21v5;
    REAL8 delta21v7;

    REAL8 rho21v1;
    REAL8 rho21v2;
    REAL8 rho21v2S;
    REAL8 rho21v3;
    REAL8 rho21v3S;
    REAL8 rho21v4;
    REAL8 rho21v4S;
    REAL8 rho21v5;
    REAL8 rho21v5S;
    REAL8 rho21v6;
    REAL8 rho21v6S;
    REAL8 rho21v6l;
    REAL8 rho21v7;
    REAL8 rho21v7S;
    REAL8 rho21v7l;
    REAL8 rho21v7lS;
    REAL8 rho21v8;
    REAL8 rho21v8l;
    REAL8 rho21v10;
    REAL8 rho21v10l;

    REAL8 f21v1;
    REAL8 f21v1S;
    REAL8 f21v3;
    REAL8 f21v3S;
    REAL8 f21v4;
    REAL8 f21v5;
    REAL8 f21v6;
    REAL8 f21v7c;
    REAL8 f21v7cEff00;
    REAL8 f21v7cEff10;
    REAL8 f21v7cEff11;
    REAL8 f21v7cEff01;
    REAL8 f21v7cEff02;

    REAL8 delta33vh3;
    REAL8 delta33vh6;
    REAL8 delta33vh6S;
    REAL8 delta33vh9;
    REAL8 delta33v5;
    REAL8 delta33v7;

    REAL8 rho33v2;
    REAL8 rho33v3;
    REAL8 rho33v4;
    REAL8 rho33v4S;
    REAL8 rho33v5;
    REAL8 rho33v5S;
    REAL8 rho33v6;
    REAL8 rho33v6S;
    REAL8 rho33v6l;
    REAL8 rho33v7;
    REAL8 rho33v7S;
    REAL8 rho33v8;
    REAL8 rho33v8l;
    REAL8 rho33v10;
    REAL8 rho33v10l;

    REAL8 f33v3;
    REAL8 f33v4;
    REAL8 f33v5;
    REAL8 f33v6;
    REAL8 f33v3S;
    REAL8 f33vh6;

    REAL8 delta32vh3;
    REAL8 delta32vh4;
    REAL8 delta32vh4S;
    REAL8 delta32vh6;
    REAL8 delta32vh6S;
    REAL8 delta32vh9;

    REAL8 rho32v;
    REAL8 rho32vS;
    REAL8 rho32v2;
    REAL8 rho32v2S;
    REAL8 rho32v3;
    REAL8 rho32v3S;
    REAL8 rho32v4;
    REAL8 rho32v4S;
    REAL8 rho32v5;
    REAL8 rho32v5S;
    REAL8 rho32v6;
    REAL8 rho32v6S;
    REAL8 rho32v6l;
    REAL8 rho32v8;
    REAL8 rho32v8l;

    REAL8 delta31vh3;
    REAL8 delta31vh6;
    REAL8 delta31vh6S;
    REAL8 delta31vh7;
    REAL8 delta31vh7S;
    REAL8 delta31vh9;
    REAL8 delta31v5;

    REAL8 rho31v2;
    REAL8 rho31v3;
    REAL8 rho31v4;
    REAL8 rho31v4S;
    REAL8 rho31v5;
    REAL8 rho31v5S;
    REAL8 rho31v6;
    REAL8 rho31v6S;
    REAL8 rho31v6l;
    REAL8 rho31v7;
    REAL8 rho31v7S;
    REAL8 rho31v8;
    REAL8 rho31v8l;

    REAL8 f31v3;
    REAL8 f31v3S;

    REAL8 delta44vh3;
    REAL8 delta44vh6;
    REAL8 delta44vh6S;
    REAL8 delta44v5;
    REAL8 delta44vh9;

    REAL8 rho44v2;
    REAL8 rho44v3;
    REAL8 rho44v3S;
    REAL8 rho44v4;
    REAL8 rho44v4S;
    REAL8 rho44v5;
    REAL8 rho44v5S;
    REAL8 rho44v6;
    REAL8 rho44v6S;
    REAL8 rho44v6l;
    REAL8 rho44v8;
    REAL8 rho44v8l;
    REAL8 rho44v10;
    REAL8 rho44v10l;

    REAL8 delta43vh3;
    REAL8 delta43vh4;
    REAL8 delta43vh4S;
    REAL8 delta43vh6;

    REAL8 rho43v;
    REAL8 rho43v2;
    REAL8 rho43v4;
    REAL8 rho43v4S;
    REAL8 rho43v5;
    REAL8 rho43v5S;
    REAL8 rho43v6;
    REAL8 rho43v6l;

    REAL8 f43v;
    REAL8 f43vS;

    REAL8 delta42vh3;
    REAL8 delta42vh6;
    REAL8 delta42vh6S;

    REAL8 rho42v2;
    REAL8 rho42v3;
    REAL8 rho42v3S;
    REAL8 rho42v4;
    REAL8 rho42v4S;
    REAL8 rho42v5;
    REAL8 rho42v5S;
    REAL8 rho42v6;
    REAL8 rho42v6S;
    REAL8 rho42v6l;

    REAL8 delta41vh3;
    REAL8 delta41vh4;
    REAL8 delta41vh4S;
    REAL8 delta41vh6;

    REAL8 rho41v;
    REAL8 rho41v2;
    REAL8 rho41v4;
    REAL8 rho41v4S;
    REAL8 rho41v5;
    REAL8 rho41v5S;
    REAL8 rho41v6;
    REAL8 rho41v6l;

    REAL8 f41v;
    REAL8 f41vS;

    REAL8 delta55vh3;
    REAL8 delta55vh6;
    REAL8 delta55vh9;

    REAL8 delta55v5;
    REAL8 rho55v2;
    REAL8 rho55v3;
    REAL8 rho55v3S;
    REAL8 rho55v4;
    REAL8 rho55v4S;
    REAL8 rho55v5;
    REAL8 rho55v5S;
    REAL8 rho55v6;
    REAL8 rho55v6l;
    REAL8 rho55v8;
    REAL8 rho55v8l;
    REAL8 rho55v10;
    REAL8 rho55v10l;
    REAL8 f55v3;
    REAL8 f55v4;
    REAL8 f55v5c;


    REAL8 delta54vh3;
    REAL8 delta54vh4;
    REAL8 delta54vh4S;
    REAL8 rho54v2;
    REAL8 rho54v3;
    REAL8 rho54v3S;
    REAL8 rho54v4;
    REAL8 rho54v4S;

    REAL8 delta53vh3;
    REAL8 rho53v2;
    REAL8 rho53v3;
    REAL8 rho53v3S;
    REAL8 rho53v4;
    REAL8 rho53v4S;
    REAL8 rho53v5;
    REAL8 rho53v5S;

    REAL8 delta52vh3;
    REAL8 delta52vh4;
    REAL8 delta52vh4S;
    REAL8 rho52v2;
    REAL8 rho52v3;
    REAL8 rho52v3S;
    REAL8 rho52v4;
    REAL8 rho52v4S;

    REAL8 delta51vh3;
    REAL8 rho51v2;
    REAL8 rho51v3;
    REAL8 rho51v3S;
    REAL8 rho51v4;
    REAL8 rho51v4S;
    REAL8 rho51v5;
    REAL8 rho51v5S;

    REAL8 delta66vh3;
    REAL8 rho66v2;
    REAL8 rho66v3;
    REAL8 rho66v3S;
    REAL8 rho66v4;
    REAL8 rho66v4S;

    REAL8 delta65vh3;
    REAL8 rho65v2;
    REAL8 rho65v3;
    REAL8 rho65v3S;

    REAL8 delta64vh3;
    REAL8 rho64v2;
    REAL8 rho64v3;
    REAL8 rho64v3S;
    REAL8 rho64v4;
    REAL8 rho64v4S;

    REAL8 delta63vh3;
    REAL8 rho63v2;
    REAL8 rho63v3;
    REAL8 rho63v3S;

    REAL8 delta62vh3;
    REAL8 rho62v2;
    REAL8 rho62v3;
    REAL8 rho62v3S;
    REAL8 rho62v4;
    REAL8 rho62v4S;

    REAL8 delta61vh3;
    REAL8 rho61v2;
    REAL8 rho61v3;
    REAL8 rho61v3S;

    REAL8 delta77vh3;
    REAL8 rho77v2;
    REAL8 rho77v3;
    REAL8 rho77v3S;

    REAL8 rho76v2;

    REAL8 delta75vh3;
    REAL8 rho75v2;
    REAL8 rho75v3;
    REAL8 rho75v3S;

    REAL8 rho74v2;

    REAL8 delta73vh3;
    REAL8 rho73v2;
    REAL8 rho73v3;
    REAL8 rho73v3S;

    REAL8 rho72v2;

    REAL8 delta71vh3;
    REAL8 rho71v2;
    REAL8 rho71v3;
    REAL8 rho71v3S;

    REAL8 rho88v2;
    REAL8 rho87v2;
    REAL8 rho86v2;
    REAL8 rho85v2;
    REAL8 rho84v2;
    REAL8 rho83v2;
    REAL8 rho82v2;
    REAL8 rho81v2;

    // New factorized waveform Coeffs
    // h22
    REAL8 h22T0ff00;
    REAL8 h22T0ff02;
    REAL8 h22T0ff20;
    REAL8 h22T0ff11;

    COMPLEX16 h22T2ff00;
    REAL8 h22T2ff02;
    REAL8 h22T2ff20;
    REAL8 h22T2ff11;
    REAL8 h22T2ff40;
    REAL8 h22T2ff31;
    REAL8 h22T2ff22;
    REAL8 h22T2ff13;
    REAL8 h22T2ff04;

    REAL8 h22T3ff10;
    REAL8 h22T3ff01;
    REAL8 h22T3ff30;
    REAL8 h22T3ff21;
    REAL8 h22T3ff12;
    REAL8 h22T3ff03;

    COMPLEX16 h22T4ff00;
    COMPLEX16 h22T4ff20;
    COMPLEX16 h22T4ff11;
    COMPLEX16 h22T4ff02;
    REAL8 h22T4ff40;
    REAL8 h22T4ff31;
    REAL8 h22T4ff22;
    REAL8 h22T4ff13;
    REAL8 h22T4ff04;
    REAL8 h22T4ff60;
    REAL8 h22T4ff51;
    REAL8 h22T4ff42;
    REAL8 h22T4ff33;
    REAL8 h22T4ff24;
    REAL8 h22T4ff15;
    REAL8 h22T4ff06;

    // h21
    REAL8 h21T0ff01;
    REAL8 h21T1ff00;

    REAL8 h21T2ff21;
    REAL8 h21T2ff12;
    REAL8 h21T2ff10;
    REAL8 h21T2ff03;
    COMPLEX16 h21T2ff01;

    REAL8 h21T3ff20;
    REAL8 h21T3ff11;
    REAL8 h21T3ff02;
    COMPLEX16 h21T3ff00;

    // h33
    REAL8 h33T0ff30;
    REAL8 h33T0ff21;
    REAL8 h33T0ff12;
    REAL8 h33T0ff10;
    REAL8 h33T0ff03;
    REAL8 h33T0ff01;

    REAL8 h33T2ff50;
    REAL8 h33T2ff41;
    REAL8 h33T2ff32;
    REAL8 h33T2ff30;
    REAL8 h33T2ff23;
    REAL8 h33T2ff21;
    REAL8 h33T2ff14;
    REAL8 h33T2ff12;
    COMPLEX16 h33T2ff10;
    REAL8 h33T2ff05;
    REAL8 h33T2ff03;
    COMPLEX16 h33T2ff01;

    REAL8 h33T3ff20;
    REAL8 h33T3ff11;
    REAL8 h33T3ff02;
    REAL8 h33T3ff00;
    REAL8 h33T3ff40;
    REAL8 h33T3ff31;
    REAL8 h33T3ff22;
    REAL8 h33T3ff13;

    // h32
    REAL8 h32T0ff11;
    REAL8 h32T0ff02;
    REAL8 h32T1ff10;
    REAL8 h32T1ff01;

    REAL8 h32T2ff31;
    REAL8 h32T2ff22;
    REAL8 h32T2ff13;
    REAL8 h32T2ff11;
    REAL8 h32T2ff04;
    REAL8 h32T2ff02;

    // h31
    REAL8 h31T0ff30;
    REAL8 h31T0ff21;
    REAL8 h31T0ff12;
    REAL8 h31T0ff10;
    REAL8 h31T0ff03;
    REAL8 h31T0ff01;

    REAL8 h31T2ff50;
    REAL8 h31T2ff41;
    REAL8 h31T2ff32;
    REAL8 h31T2ff30;
    REAL8 h31T2ff23;
    REAL8 h31T2ff21;
    REAL8 h31T2ff14;
    REAL8 h31T2ff12;
    REAL8 h31T2ff10;
    REAL8 h31T2ff05;
    REAL8 h31T2ff03;
    REAL8 h31T2ff01;

    REAL8 h31T3ff20;
    REAL8 h31T3ff11;
    REAL8 h31T3ff02;
    REAL8 h31T3ff00;
    
    // h44
    REAL8 h44T0ff40;
    REAL8 h44T0ff31;
    REAL8 h44T0ff22;
    REAL8 h44T0ff20;
    REAL8 h44T0ff13;
    REAL8 h44T0ff11;
    REAL8 h44T0ff04;
    REAL8 h44T0ff02;
    REAL8 h44T0ff00;

    REAL8 h44T2ff60;
    REAL8 h44T2ff51;
    REAL8 h44T2ff42;
    REAL8 h44T2ff40;
    REAL8 h44T2ff31;
    REAL8 h44T2ff24;
    REAL8 h44T2ff22;
    COMPLEX16 h44T2ff20;
    REAL8 h44T2ff15;
    REAL8 h44T2ff13;
    COMPLEX16 h44T2ff11;
    REAL8 h44T2ff06;
    REAL8 h44T2ff04;
    COMPLEX16 h44T2ff02;
    COMPLEX16 h44T2ff00;

    // h43
    REAL8 h43T0ff21;
    REAL8 h43T0ff12;
    REAL8 h43T0ff03;
    REAL8 h43T0ff01;

    REAL8 h43T1ff20;
    REAL8 h43T1ff11;
    REAL8 h43T1ff02;
    REAL8 h43T1ff00;

    // h42
    REAL8 h42T0ff40;
    REAL8 h42T0ff31;
    REAL8 h42T0ff20;
    REAL8 h42T0ff13;
    REAL8 h42T0ff11;
    REAL8 h42T0ff04;
    REAL8 h42T0ff02;
    REAL8 h42T0ff00;

    REAL8 h42T2ff60;
    REAL8 h42T2ff51;
    REAL8 h42T2ff42;
    REAL8 h42T2ff40;
    REAL8 h42T2ff33;
    REAL8 h42T2ff31;
    REAL8 h42T2ff24;
    REAL8 h42T2ff22;
    REAL8 h42T2ff20;
    REAL8 h42T2ff15;
    REAL8 h42T2ff13;
    REAL8 h42T2ff11;
    REAL8 h42T2ff06;
    REAL8 h42T2ff04;
    REAL8 h42T2ff02;
    REAL8 h42T2ff00;

    // h41
    REAL8 h41T0ff21;
    REAL8 h41T0ff12;
    REAL8 h41T0ff03;
    REAL8 h41T0ff01;

    REAL8 h41T1ff20;
    REAL8 h41T1ff11;
    REAL8 h41T1ff02;
    REAL8 h41T1ff00;
}
FacWaveformCoeffs;

typedef struct tagEOBNonQCCoeffs
{
  REAL8 a1;
  REAL8 a2;
  REAL8 a3;
  REAL8 a3S;
  REAL8 a4;
  REAL8 a5;
  REAL8 b1;
  REAL8 b2;
  REAL8 b3;
  REAL8 b4;
} EOBNonQCCoeffs;


/* Spin EOBH Coeffs */
typedef struct
tagSpinEOBHCoeffs
{
    double KK;
    double k0;
    double k1;
    double k2;
    double k3;
    double k4;
    double k5;
    double k5l;
    double b3;
    double bb3;
    double d1;
    double d1v2;
    double dheffSS;
    double dheffSSv2;
    int updateHCoeffs;
}
SpinEOBHCoeffs;

typedef struct {
    REAL8 eta;
    REAL8 a;
    REAL8 a2;
    REAL8 chi;
    REAL8 sigmaStar;
    REAL8 sigmaKerr;

    REAL8 invm1PlusEtaKK;
    REAL8 k0;
    REAL8 k1;
    REAL8 k2;
    REAL8 k3;
    REAL8 k4;
    REAL8 k5;
    REAL8 k5l;

    REAL8 D2;
    REAL8 D3;

    REAL8 b3;
    REAL8 bb3;

    REAL8 PQ4;
    REAL8 Pds0u;
    REAL8 Pds0p2;
    REAL8 Pds0pn2;

    REAL8 PdsSu2;
    REAL8 PdsSup2;
    REAL8 PdsSpn4;
    REAL8 PdsSp4;
    REAL8 PdsSupn2;
    REAL8 PdsSp2pn2;

    REAL8 PdsKu2;
    REAL8 PdsKup2;
    REAL8 PdsKpn4;
    REAL8 PdsKp4;
    REAL8 PdsKupn2;
    REAL8 PdsKp2pn2;
    REAL8 PdsKu3;

    REAL8 Peffss;
    REAL8 sign;
}SpinEOBHSACoeffs;

typedef struct
tagSEOBHCoeffConstants
{

    double a0k2; //Coefficient of a^0 in k2
    double a1k2; //Coefficient of a^1 in k2

    double a0k3; //Coefficient of a^0 in k3
    double a1k3; //Coefficient of a^1 in k3

    double a0k4; //Coefficient of a^0 in k4
    double a1k4; //Coefficient of a^1 in k4
    double a2k4; //Coefficient of a^2 in k4

    double a0k5; //Coefficient of a^0 in k5
    double a1k5; //Coefficient of a^1 in k5
    double a2k5; //Coefficient of a^2 in k5
}
SEOBHCoeffConstants;

typedef struct tagNewtonMultipolePrefixes
{
    COMPLEX16 values[9][9];
}
NewtonMultipolePrefixes;

/* Spin EOB Core */
typedef struct tagSpinEOBParams
{
    // Intrinsic
    REAL8                   eta;
    REAL8                   m1;
    REAL8                   m2;
    REAL8                   a;
    REAL8                   chi1;
    REAL8                   chi2;
    REAL8                   p0; // initial semi-latus
    REAL8                   e0;
    REAL8                   x0; // initial inclination
    // Numerical Deriv aux
    REAL8                   omega;
    UINT                    omegaPeaked;
    INT                     termination_reason;
    REAL8                   prev_dr;
    REAL8                   spin1z_omegaPeak;
    REAL8                   spin2z_omegaPeak;
    REAL8                   tPeakOmega;

    // Spin Vectors
    REAL8Vector             *chi1Vec;
    REAL8Vector             *chi2Vec;
    REAL8Vector             *s1Vec;
    REAL8Vector             *s2Vec;
    REAL8Vector             *sigmaStar;
    REAL8Vector             *sigmaKerr;
    REAL8Vector             *J0Vec;
    // Structs
    SpinEOBHCoeffs          *seobCoeffs;
    SpinEOBHSACoeffs        *saCoeffs;
    SEOBHCoeffConstants     *seobCoeffConsts;
    HyperParams             *hParams;
    NewtonMultipolePrefixes *prefixes;
    FacWaveformCoeffs       *hCoeffs;
    // FacWaveformCoeffs       *hCoeffs_flux;
    EccCorrectionCoeffs     *eccCoeffs;
    EOBNonQCCoeffs          *nqcCoeffs;

    // hCoeffs cal
    COMPLEX16                   cal21;
    REAL8                       cal21E;
    REAL8                       cal21E1;
    REAL8                       cal21E2;
    REAL8                       cal21E3;
    REAL8                       cal21E4;
    COMPLEX16                   cal55;

    // NQC window
    REAL8                   tWind;
    REAL8                   wWind;

    // Booleans
    INT                     tortoise;
    BOOLEAN                 alignedSpins;
    BOOLEAN                 ignoreflux;
    BOOLEAN                 use_hm;

    // cache
    REAL8                   cache[3];
}SpinEOBParams;

typedef struct tagSpinEOBDynamics
{
    UINT length;
    REAL8 t0;
    REAL8Vector *tVec;
    REAL8Vector *xVec;
    REAL8Vector *yVec;
    REAL8Vector *zVec;
    REAL8Vector *pxVec;
    REAL8Vector *pyVec;
    REAL8Vector *pzVec;
}SpinEOBDynamics;

/* ------------------------ Deriv & Integral struct------------------ */
/* Hcap */
typedef
struct tagHcapDerivParams
{
    const REAL8   *values;
    SpinEOBParams *params;
    UINT         varyParam;
}
HcapDerivParams;

typedef
struct tagHcapSphDeriv2Params
{
    const REAL8     *sphValues;
    SpinEOBParams   *params;
    UINT           varyParam1;
    UINT           varyParam2;
}
HcapSphDeriv2Params;


/**
 * Structure the EOB dynamics for precessing waveforms.
 */
#define v4PdynamicsVariables 36
typedef struct tagSEOBdynamics 
{
    UINT length;
    REAL8 th22Peak;
    REAL8Array *array;

    REAL8 *tVec;

    REAL8 *posVecx;
    REAL8 *posVecy;
    REAL8 *posVecz;

    REAL8 *momVecx;
    REAL8 *momVecy;
    REAL8 *momVecz;

    REAL8 *s1Vecx;
    REAL8 *s1Vecy;
    REAL8 *s1Vecz;

    REAL8 *s2Vecx;
    REAL8 *s2Vecy;
    REAL8 *s2Vecz;

    REAL8 *phiDMod;
    REAL8 *phiMod;

    REAL8 *velVecx;
    REAL8 *velVecy;
    REAL8 *velVecz;

    REAL8 *polarrVec;
    REAL8 *polarphiVec;
    REAL8 *polarprVec;
    REAL8 *polarpphiVec;

    REAL8 *omegaVec;

    REAL8 *s1dotZVec;
    REAL8 *s2dotZVec;

    // chi expanded by nVec, vVec, LNVec
    REAL8 *chiAxVec;
    REAL8 *chiAyVec;
    REAL8 *chiAzVec;
    REAL8 *chiSxVec;
    REAL8 *chiSyVec;
    REAL8 *chiSzVec;

    REAL8 *hamVec;
    REAL8 *fluxVec;

    REAL8 *FrVec;
    REAL8 *FfVec;
    REAL8 *polarprDotVec;
} SEOBdynamics;

#define v4SAdynamicsVariables 10
typedef struct {
    UINT length;
    REAL8 th22Peak;
    REAL8Array *array;

    REAL8 *tVec;
    REAL8 *rVec;
    REAL8 *phiVec;

    REAL8 *prTVec;
    REAL8 *pphiVec;

    REAL8 *drVec;
    REAL8 *dphiVec;

    REAL8 *HVec;
    
    REAL8 *dprTVec;
    REAL8 *dpphiVec;

}SEOBSAdynamics;

typedef struct tagSphHarmListEOBNonQCCoeffs {
    EOBNonQCCoeffs*                        nqcCoeffs; /**< NQC coefficients for this mode. */
    UINT                                   l; /**< Mode number l. */
    INT                                    m; /**< Mode number m. */
    struct tagSphHarmListEOBNonQCCoeffs   *next; /**< Pointer to next element in the list. */
} SphHarmListEOBNonQCCoeffs;

/**
 * Structure to represent a data piece (e.g. a mode hlm), either in frequency or time
 * with complex amplitude (enveloppe) and phase.
 * The mode values are camp * exp(I*phase)
 */
typedef struct tagCAmpPhaseSequence {
    REAL8Vector*                        xdata; /**< Sequence of times or frequencies on which data is given. */
    REAL8Vector*                        camp_real; /**< Sequence for the real part of the complex amplitude (enveloppe). */
    REAL8Vector*                        camp_imag; /**< Sequence for the imag part of the complex amplitude (enveloppe). */
    REAL8Vector*                        phase; /**< Sequence for the phase. */
} CAmpPhaseSequence;

/**
 * Structure to represent linked list of modes
 * with complex amplitude (enveloppe) and phase.
 */
typedef struct tagSphHarmListCAmpPhaseSequence {
    CAmpPhaseSequence*                        campphase; /**< Data for this mode. */
    UINT                                      l; /**< Mode number l. */
    INT                                       m; /**< Mode number m. */
    struct tagSphHarmListCAmpPhaseSequence*   next; /**< Pointer to next element in the list. */
} SphHarmListCAmpPhaseSequence;

/**
 * Structure to carry a collection of spherical harmonic modes in COMPLEX16 
 * time series. Contains convenience getter and setter functions, as well as
 * a convienence "maximum l mode" function. Implemented as a singly forward
 * linked list.
 */
typedef struct tagSphHarmTimeSeries {
    COMPLEX16TimeSeries*            mode; /**< The sequences of sampled data. */
    UINT                            l; /**< Node mode l  */
    INT                             m; /**< Node submode m  */
    REAL8Vector*                  tdata; /**< Timestamp values */
    REAL8                           tAttach;
    struct tagSphHarmTimeSeries*    next; /**< next pointer */
} SphHarmTimeSeries;

typedef struct tagSEOBCoreOutputs
{
    SphHarmTimeSeries   *hLM;
    SEOBdynamics        *dyn;
    REAL8Vector         *flux;
    SphHarmListCAmpPhaseSequence *Plm;
}SEOBCoreOutputs;



typedef struct tagRRForceCoeffs
{
    // fabc -> p^2a pr^2b / r^c
    REAL8 LOPre;
    REAL8 LOf100;
    REAL8 LOf010;
    REAL8 LOf001;

    REAL8 P1Pre;
    REAL8 P1f110;
    REAL8 P1f200;
    REAL8 P1f020;
    REAL8 P1f011;
    REAL8 P1f101;
    REAL8 P1f002;

    REAL8 TPre;
    REAL8 Tf201;
    REAL8 Tf102;
    REAL8 Tf101;
    REAL8 Tf011;
    REAL8 Tf111;
    REAL8 Tf020;
    REAL8 Tf012;
    REAL8 Tf030;
    REAL8 Tf002;

    REAL8 P2Pre;
    REAL8 P2f300;
    REAL8 P2f030;
    REAL8 P2f003;
    REAL8 P2f210;
    REAL8 P2f201;
    REAL8 P2f120;
    REAL8 P2f021;
    REAL8 P2f102;
    REAL8 P2f012;
    REAL8 P2f111;

    REAL8 SOPre;
    REAL8 SOf200;
    REAL8 SOf020;
    REAL8 SOf011;
    REAL8 SOf101;
    REAL8 SOf110;
    REAL8 SOf001;
    REAL8 SOf100;
    REAL8 SOf010;

    REAL8 SSPre;
    REAL8 SSf100;
    REAL8 SSf010;
    REAL8 SSf001;
}RRForceCoeffs;


/*--------------------------------------------------------------*/
/*                                                              */
/*                                                              */
/*                                                              */
/*                            PREC                              */
/*                                                              */
/*                                                              */
/*                                                              */
/*--------------------------------------------------------------*/
#define v4PrecdynamicsVariables 40
#define v4PrecEvolvedynamicsVariables 21

typedef struct {
    UINT length;
    REAL8 th22Peak;
    REAL8Array *array;

    REAL8 *tVec;

    REAL8 *posVecx;
    REAL8 *posVecy;
    REAL8 *posVecz;

    REAL8 *momTVecx;
    REAL8 *momTVecy;
    REAL8 *momTVecz;

    REAL8 *s1Vecx;
    REAL8 *s1Vecy;
    REAL8 *s1Vecz;

    REAL8 *s2Vecx;
    REAL8 *s2Vecy;
    REAL8 *s2Vecz;

    REAL8 *phiDMod;
    REAL8 *phiMod;

    REAL8 *velVecx;
    REAL8 *velVecy;
    REAL8 *velVecz;

    REAL8 *HamVec;
    REAL8 *fluxVec; // for debug
    REAL8 *prTDotVec; // for debug

    // other params
    REAL8 *polarrVec;
    REAL8 *polarphiVec;
    REAL8 *polarprTVec;
    REAL8 *polarpphiVec;
    REAL8 *omegaVec;

    REAL8 *nchiaVec;
    REAL8 *nchisVec;
    REAL8 *lchiaVec;
    REAL8 *lchisVec;
    REAL8 *echiaVec;
    REAL8 *echisVec;

    REAL8 *chi1chi1;
    REAL8 *chi1chi2;
    REAL8 *chi2chi2;

    REAL8 *JnVec;
    REAL8 *JlVec;
    REAL8 *JeVec;

    REAL8 *s1dotZVec;
    REAL8 *s2dotZVec;

}SEOBPrecdynamics;

typedef struct tagSEOBPrecCoreOutputs
{
    REAL8Vector         *tVec;
    SphHarmTimeSeries   *hLM;
    SEOBPrecdynamics    *dyn;
    SphHarmListCAmpPhaseSequence *Plm;
}SEOBPrecCoreOutputs;

typedef struct {
    REAL8 xVec[3];
    REAL8 pTVec[3];
    REAL8 s1Vec[3];
    REAL8 s2Vec[3];

    REAL8 LVec[3];
    REAL8 r;
    REAL8 prT;
    REAL8 pT2;
    
    REAL8 xS1;
    REAL8 xS2;
    REAL8 pTS1;
    REAL8 pTS2;
    REAL8 xpS1;
    REAL8 xpS2;

    REAL8 S1S1;
    REAL8 S1S2;
    REAL8 S2S2;
}SEOBPrecVariables;

typedef struct {
    REAL8 nchia;
    REAL8 nchis;
    REAL8 lchia;
    REAL8 lchis;
    REAL8 echia;
    REAL8 echis;
    REAL8 chi1chi1;
    REAL8 chi1chi2;
    REAL8 chi2chi2;
    REAL8 Jn;
    REAL8 Jl;
    REAL8 Je;
    REAL8 x;
    REAL8 x2;
    REAL8 x3;
    REAL8 x4;
    REAL8 y;
    REAL8 y2;
    REAL8 y3;
    REAL8 y4;
    REAL8 y5;
    REAL8 y6;
    REAL8 z2;
    REAL8 z3;
    REAL8 z4;
    REAL8 z5;
    REAL8 z6;
    REAL8 z7;
    REAL8 z8;
    REAL8 z9;
    REAL8 z10;
    REAL8 z11;
    REAL8 z12;
    REAL8 z13;
    REAL8 z14;
    REAL8 z15;
    REAL8 z16;
    REAL8 z18;
    REAL8 z20;
    REAL8 z21;
    REAL8 z22;
    REAL8 z24;
    REAL8 z30;
}SEOBPrecWaveformVariables;

#endif

