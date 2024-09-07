#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu, 23 Mar 2023 06:10:15 +0000

@author: Shallyn
"""
import numpy as np
import matplotlib.pyplot as plt
import sys, os
from pathlib import Path
import ctypes
from .ctype_structs import *
from .pyUtils import *
#pwd = Path(__file__).absolute().parent
#pwd = Path(sys.path[0])
# pwd = Path(os.getcwd())

#-----Constants-----#
k_B_SI = 1.3806505e-23 # J/K
h_SI = 6.62606896e-34 # J s
ly_SI = 9.4605284e15 # 1 ly to m
AU_SI = 149597870700 # 1 AU to m
pc_SI = 3.08567758e16 # 1 Parsec to m
c_SI = 299792458 # Speed Of Light [m s^-1]
M_Sun_SI = 1.98892e30 # Mass of Sun [kg]
R_Sun_SI = 6.96342e8 # Radius of Sun [m]
alp_GP_SI = 192.85948 * np.pi / 180 # Direction Of Galactic Polar
det_GP_SI = 27.12825 * np.pi / 180 # Direction Of Galactic Polar
l_CP_SI = 122.932
G_SI = 6.672e-11 # Gravitational Constant [m^3 s^-2 kg^-1]
h_0_SI = 0.678 # Hubble Constant
MRSUN_SI = 1.47662504e3 
MTSUN_SI = 4.92549095e-6 
#---------Comm--------#
def dim_t(M):
    return c_SI**3 / ( M * M_Sun_SI * G_SI)

lib_path = Path(__file__).parent / '.libs/libEOB.so'
my_waveform_lib = ctypes.CDLL(lib_path)
my_waveform_lib.CreateREAL8Vector.restype = ctypes.POINTER(pyREAL8Vector)
my_waveform_lib.CreateSpinEOBParams.restype = ctypes.POINTER(pySpinEOBParams)
my_waveform_lib.calculate_QNMFrequency.restype = ctypes.c_double
my_waveform_lib.calculate_QNMFrequencies.restype = ctypes.c_int
my_waveform_lib.calculate_QNMFrequenciesFromFinal.restype = ctypes.c_int

def calculate_QNMFrequenciesFromFinal(mFinal:float, chiFinal:float,
                        modeL:int=2, modeM:int=2, nmodes:int=1):
    FreqRVec = ctypes.POINTER(pyREAL8Vector)()
    FreqIVec = ctypes.POINTER(pyREAL8Vector)()
    value_list = [ctypes.c_double(mFinal),
                  ctypes.c_double(chiFinal),
                  ctypes.c_uint(modeL),
                  ctypes.c_uint(modeM),
                  ctypes.c_size_t(nmodes),
                  ctypes.byref(FreqRVec), ctypes.byref(FreqIVec)]
    ret = my_waveform_lib.calculate_QNMFrequenciesFromFinal(*value_list)
    if ret != 0:
        return None
    npFreqRVec = convert_REAL8Vector_to_numpy(FreqRVec)
    npFreqIVec = convert_REAL8Vector_to_numpy(FreqIVec)
    my_waveform_lib.DestroyREAL8Vector(FreqRVec)
    my_waveform_lib.DestroyREAL8Vector(FreqIVec)
    return npFreqRVec + 1.j*npFreqIVec

def calculate_QNMFrequencies(m1:float, m2:float, 
                           chi1x:float, chi1y:float, chi1z:float, 
                           chi2x:float, chi2y:float, chi2z:float,
                           modeL:int=2, modeM:int=2, nmodes:int=1):
    FreqRVec = ctypes.POINTER(pyREAL8Vector)()
    FreqIVec = ctypes.POINTER(pyREAL8Vector)()
    value_list = [ctypes.c_double(m1),
                  ctypes.c_double(m2),
                  ctypes.c_double(chi1x),
                  ctypes.c_double(chi1y),
                  ctypes.c_double(chi1z),
                  ctypes.c_double(chi2x),
                  ctypes.c_double(chi2y),
                  ctypes.c_double(chi2z),
                  ctypes.c_uint(modeL),
                  ctypes.c_uint(modeM),
                  ctypes.c_size_t(nmodes),
                  ctypes.byref(FreqRVec), ctypes.byref(FreqIVec)]
    ret = my_waveform_lib.calculate_QNMFrequencies(*value_list)
    if ret != 0:
        return None
    npFreqRVec = convert_REAL8Vector_to_numpy(FreqRVec)
    npFreqIVec = convert_REAL8Vector_to_numpy(FreqIVec)
    my_waveform_lib.DestroyREAL8Vector(FreqRVec)
    my_waveform_lib.DestroyREAL8Vector(FreqIVec)
    return npFreqRVec + 1.j*npFreqIVec


def calculate_QNMFrequency(m1:float, m2:float, 
                           chi1x:float, chi1y:float, chi1z:float, 
                           chi2x:float, chi2y:float, chi2z:float,
                           modeL:int=5, modeM:int=5):
    value_list = [ctypes.c_double(m1),
                  ctypes.c_double(m2),
                  ctypes.c_double(chi1x),
                  ctypes.c_double(chi1y),
                  ctypes.c_double(chi1z),
                  ctypes.c_double(chi2x),
                  ctypes.c_double(chi2y),
                  ctypes.c_double(chi2z),
                  ctypes.c_uint(modeL),
                  ctypes.c_uint(modeM)]
    return my_waveform_lib.calculate_QNMFrequency(*value_list)

def convert_REAL8Vector_to_numpy(vec:ctypes.POINTER(pyREAL8Vector)):
    length = vec.contents.length
    ret = np.zeros(length)
    npdata = ret.ctypes.data
    ctypes.memmove(npdata, vec.contents.data, length*ctypes.sizeof(ctypes.c_double))
    return ret

def convert_ndarray_to_REAL8Vector(vec:np.ndarray):
    if len(vec.shape) > 1:
        raise Exception('the shape of vec should be 1-d')
    length = int(len(vec))
    npvec = np.zeros(length)
    npvec[:] = vec[:]
    ret = my_waveform_lib.CreateREAL8Vector(ctypes.c_uint(length))
    ctypes.memmove(ret.contents.data, npvec.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), length * ctypes.sizeof(ctypes.c_double))
    # print(ret.contents.data)
    # print(vec.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    # ctypes.memmove(ret.contents.data, vec.ctypes.data, vec.nbytes)
    return ret

class npDynamicData(object):
    def __init__(self, m1:float, m2:float, dyndata:ctypes.POINTER(pyDynOutputStruct)):
        self.__m1 = m1
        self.__m2 = m2
        self.__timeM = convert_REAL8Vector_to_numpy(dyndata.contents.timeM)
        self.__xVec = convert_REAL8Vector_to_numpy(dyndata.contents.xVec)
        self.__yVec = convert_REAL8Vector_to_numpy(dyndata.contents.yVec)
        self.__zVec = convert_REAL8Vector_to_numpy(dyndata.contents.zVec)
        self.__pTxVec = convert_REAL8Vector_to_numpy(dyndata.contents.pTxVec)
        self.__pTyVec = convert_REAL8Vector_to_numpy(dyndata.contents.pTyVec)
        self.__pTzVec = convert_REAL8Vector_to_numpy(dyndata.contents.pTzVec)
        self.__vxVec = convert_REAL8Vector_to_numpy(dyndata.contents.vxVec)
        self.__vyVec = convert_REAL8Vector_to_numpy(dyndata.contents.vyVec)
        self.__vzVec = convert_REAL8Vector_to_numpy(dyndata.contents.vzVec)
        self.__s1xVec = convert_REAL8Vector_to_numpy(dyndata.contents.s1xVec)
        self.__s1yVec = convert_REAL8Vector_to_numpy(dyndata.contents.s1yVec)
        self.__s1zVec = convert_REAL8Vector_to_numpy(dyndata.contents.s1zVec)
        self.__s2xVec = convert_REAL8Vector_to_numpy(dyndata.contents.s2xVec)
        self.__s2yVec = convert_REAL8Vector_to_numpy(dyndata.contents.s2yVec)
        self.__s2zVec = convert_REAL8Vector_to_numpy(dyndata.contents.s2zVec)
        self.__phiDModVec = convert_REAL8Vector_to_numpy(dyndata.contents.phiDModVec)
        self.__phiModVec = convert_REAL8Vector_to_numpy(dyndata.contents.phiModVec)
        self.__prTDotVec = convert_REAL8Vector_to_numpy(dyndata.contents.prTDotVec)
        self.__hamVec = convert_REAL8Vector_to_numpy(dyndata.contents.hamVec)
    def __getitem__(self, key:int):
        return np.array([self.__xVec[key], self.__yVec[key], self.__zVec[key], 
            self.__pTxVec[key], self.__pTyVec[key], self.__pTzVec[key], 
            self.__s1xVec[key], self.__s1yVec[key], self.__s1zVec[key],
            self.__s2xVec[key], self.__s2yVec[key], self.__s2zVec[key],
            self.__phiDModVec[key], self.__phiModVec[key]])
    def __len__(self):
        return len(self.__timeM)
    @property
    def phiDMod(self):
        return self.__phiDModVec
    @property
    def phiMod(self):
        return self.__phiModVec
    @property
    def m1(self):
        return self.__m1
    @property
    def m2(self):
        return self.__m2
    @property
    def M(self):
        return self.m1 + self.m2
    @property
    def eta(self):
        return self.m1*self.m2/self.M/self.M
    @property
    def dm(self):
        return np.abs(self.m1 - self.m2) / self.M
    @property
    def Q1(self):
        return self.m1 / self.M
    @property
    def Q2(self):
        return self.m2 / self.M
    @property
    def timeM(self):
        return self.__timeM
    
    @property
    def x(self):
        return self.__xVec
    @property
    def y(self):
        return self.__yVec
    @property
    def z(self):
        return self.__zVec
        
    @property
    def pTx(self):
        return self.__pTxVec
    @property
    def pTy(self):
        return self.__pTyVec
    @property
    def pTz(self):
        return self.__pTzVec
    
    @property
    def prTDot(self):
        return self.__prTDotVec
    @property
    def Hreal(self):
        return self.__hamVec

    @property
    def vx(self):
        return self.__vxVec
    @property
    def vy(self):
        return self.__vyVec
    @property
    def vz(self):
        return self.__vzVec

    @property
    def r(self):
        return np.sqrt(self.__xVec*self.__xVec + self.__yVec*self.__yVec + self.__zVec*self.__zVec)
    @property
    def rDot(self):
        return (self.vx*self.x + self.vy*self.y + self.vz*self.z) / self.r
    @property
    def v2(self):
        return self.vx*self.vx + self.vy*self.vy + self.vz*self.vz
    @property
    def prT(self):
        return (self.pTx*self.x + self.pTy*self.y + self.pTz*self.z) / self.r
    
    @property
    def s1Vec(self):
        return np.stack([self.__s1xVec, self.__s1yVec, self.__s1zVec], axis = 1)
    @property
    def s2Vec(self):
        return np.stack([self.__s2xVec, self.__s2yVec, self.__s2zVec], axis = 1)
    @property
    def chi1Vec(self):
        return self.s1Vec / np.power(self.Q1, 2)
    @property
    def chi2Vec(self):
        return self.s2Vec / np.power(self.Q2, 2)
    @property
    def rVec(self):
        return np.stack([self.__xVec, self.__yVec, self.__zVec], axis = 1)
    @property
    def vVec(self):
        return np.stack([self.__vxVec, self.__vyVec, self.__vzVec], axis = 1)
    @property
    def pTVec(self):
        return np.stack([self.__pTxVec, self.__pTyVec, self.__pTzVec], axis = 1)
    @property
    def LVec(self):
        Lx = self.pTz * self.y - self.pTy * self.z
        Ly = self.pTx * self.z - self.pTz * self.x
        Lz = self.pTy * self.x - self.pTx * self.y
        # shape (N, 3)
        return np.stack([Lx, Ly, Lz], axis = 1)
    @property
    def LNVec(self):
        LNx = self.vz * self.y - self.vy * self.z
        LNy = self.vx * self.z - self.vz * self.x
        LNz = self.vy * self.x - self.vx * self.y
        # shape (N, 3)
        return np.stack([LNx, LNy, LNz], axis = 1)
    @property
    def magL(self):
        return np.linalg.norm(self.LVec, axis = 1)
    @property
    def magLN(self):
        return np.linalg.norm(self.LNVec, axis = 1)
    @property
    def omega(self):
        return self.magLN / np.power(self.r, 2)
    

class npWaveformData(object):
    def __init__(self, hcdata:ctypes.POINTER(pyOutputStruct), t0 = 0.0):
        self.__timeM = convert_REAL8Vector_to_numpy(hcdata.contents.timeM)
        self.__timeM = self.__timeM - self.__timeM[0]
        self.__time = convert_REAL8Vector_to_numpy(hcdata.contents.time)
        self.__time = self.__time - self.__time[0] + t0
        self.__hplus = convert_REAL8Vector_to_numpy(hcdata.contents.hplus)
        self.__hcross = convert_REAL8Vector_to_numpy(hcdata.contents.hcross)
        self.__h22_real = convert_REAL8Vector_to_numpy(hcdata.contents.h22_real)
        self.__h22_imag = convert_REAL8Vector_to_numpy(hcdata.contents.h22_imag)
        self.__h21_real = convert_REAL8Vector_to_numpy(hcdata.contents.h21_real)
        self.__h21_imag = convert_REAL8Vector_to_numpy(hcdata.contents.h21_imag)
        self.__h33_real = convert_REAL8Vector_to_numpy(hcdata.contents.h33_real)
        self.__h33_imag = convert_REAL8Vector_to_numpy(hcdata.contents.h33_imag)
        self.__h44_real = convert_REAL8Vector_to_numpy(hcdata.contents.h44_real)
        self.__h44_imag = convert_REAL8Vector_to_numpy(hcdata.contents.h44_imag)
        self.__h55_real = convert_REAL8Vector_to_numpy(hcdata.contents.h55_real)
        self.__h55_imag = convert_REAL8Vector_to_numpy(hcdata.contents.h55_imag)
    
    def __len__(self):
        return len(self.__timeM)

    @property
    def timeM(self):
        return self.__timeM

    @property
    def time(self):
        return self.__time

    @property
    def hpc(self):
        return cTimeSeries(self.__time, self.__hplus, -self.__hcross)
    
    @property
    def h22(self):
        return cTimeSeries(self.__timeM, self.__h22_real, self.__h22_imag)
    
    @property
    def h21(self):
        return cTimeSeries(self.__timeM, self.__h21_real, self.__h21_imag)

    @property
    def h33(self):
        return cTimeSeries(self.__timeM, self.__h33_real, self.__h33_imag)

    @property
    def h44(self):
        return cTimeSeries(self.__timeM, self.__h44_real, self.__h44_imag)

    @property
    def h55(self):
        return cTimeSeries(self.__timeM, self.__h55_real, self.__h55_imag)



class SEOBNRWaveformCaller(object):
    def __init__(self, **kwargs):
        self.__set_default_params()
        self.__parse_params(**kwargs)
    
    def __set_default_params(self):
        self.use_geom = False
        self.m1 = 10.
        self.m2 = 10.
        self.s1x = 0.0
        self.s1y = 0.0
        self.s1z = 0.0
        self.s2x = 0.0
        self.s2y = 0.0
        self.s2z = 0.0
        self.ecc = 0.0
        self.inc = 0.0
        self.phiRef = 0.0
        self.beta = 0.0 # FIXME
        self.f_min = 40.
        self.distance = 100.
        self.srate = 16384.
        self.deltaT = 1./self.srate
        self.is_only22 = False
        self.is_constp = False
        self.conserve_flag = 1
        self.conserve_time = -1.
        self.is_noringdown = False
        self.prec_flag = 0
        self.debug_id = 0
        # hparams
        self.d_ini = -1
        self.pr_ini = 0.0
        self.pphi_ini = 0.0
        self.ptheta_ini = 0.0
        self.flagTuning = 0
        self.tStepBack = 200.
        self.sl_p = -1
        self.x0 = np.cos(0.0)
        self.code_version = 1
        self.log_level = 1
        self.risco = -1.0
        self.KK = 0.0
        self.dSS = 0.0
        self.dSO = 0.0
        self.dtPeak = 0.0
        self.t0 = 0.0
        self.egw_flag = 0
        self.ret_dyn = 0
        self.EPS_REL = -1.
        self.EPS_ABS = -1.
        self.is_coframe = 0 # only valid when prec_flag > 0
        self.use_coaphase = 0 # default
        self.zeta = 0 # anomaly angle zeta in r = p / (1 + e cos(zeta))
        self.xi = 0
        #self.Mf_ref = 0.0049258454886941605 # recomended value: omega_22 = 0.03095
        self.Mf_ref = 0.002
        self.zero_dyncoaphase = 0
        self.f_max = -1.
        self.t_max = -1.
        self.ICValues = None

    def set_params(self, **kwargs):
        self.__parse_params(**kwargs)
    
    def __parse_params(self, **kwargs):
        if 'use_geom' in kwargs:
            self.use_geom = 0 if kwargs['use_geom'] == False else 1
        self.m1 = self.m1 if 'm1' not in kwargs else kwargs['m1']
        self.m2 = self.m2 if 'm2' not in kwargs else kwargs['m2']
        self.s1x = self.s1x if 's1x' not in kwargs else kwargs['s1x']
        self.s1y = self.s1y if 's1y' not in kwargs else kwargs['s1y']
        self.s1z = self.s1z if 's1z' not in kwargs else kwargs['s1z']
        self.s2x = self.s2x if 's2x' not in kwargs else kwargs['s2x']
        self.s2y = self.s2y if 's2y' not in kwargs else kwargs['s2y']
        self.s2z = self.s2z if 's2z' not in kwargs else kwargs['s2z']
        self.ecc = self.ecc if 'ecc' not in kwargs else kwargs['ecc']
        self.inc = self.inc if 'inc' not in kwargs else kwargs['inc']*np.pi/180.
        if 'inc_rad' in kwargs:
            self.inc = kwargs['inc_rad']
        self.phiRef = self.phiRef if 'phiRef' not in kwargs else kwargs['phiRef']*np.pi/180.
        if 'phiRef_rad' in kwargs:
            self.phiRef = kwargs['phiRef_rad']
        self.beta = self.beta if 'beta' not in kwargs else kwargs['beta']*np.pi/180.
        if 'beta_rad' in kwargs:
            self.beta = kwargs['beta_rad']
        self.f_min = self.f_min if 'f_min' not in kwargs else kwargs['f_min']
        self.Mf_ref = self.Mf_ref if 'Mf_ref' not in kwargs else kwargs['Mf_ref']
        if 'f_ref' in kwargs:
            self.Mf_ref = kwargs['f_ref'] / dim_t(self.m1 + self.m2)
        self.distance = self.distance if 'distance' not in kwargs else kwargs['distance']
        if 'srate' in kwargs:
            self.srate = kwargs['srate']
            self.deltaT = 1./self.srate
        if 'deltaT' in kwargs:
            self.deltaT = kwargs['deltaT']
            self.srate = 1./self.deltaT
        if 'is_only22' in kwargs:
            self.is_only22 = 0 if kwargs['is_only22'] == False else 1
        if 'is_constp' in kwargs:
            self.is_constp = 0 if kwargs['is_constp'] == False else 1
        self.conserve_flag = self.conserve_flag if 'conserve_flag' not in kwargs else kwargs['conserve_flag']
        self.conserve_time = self.conserve_time if 'conserve_time' not in kwargs else kwargs['conserve_time']
        if 'is_noringdown' in kwargs:
            self.is_noringdown = 0 if kwargs['is_noringdown'] == False else 1
        self.prec_flag = self.prec_flag if 'prec_flag' not in kwargs else kwargs['prec_flag']
        self.debug_id = self.debug_id if 'debug_id' not in kwargs else kwargs['debug_id']
        # hparams
        self.d_ini = self.d_ini if 'd_ini' not in kwargs else kwargs['d_ini']
        self.pr_ini = self.pr_ini if 'pr_ini' not in kwargs else kwargs['pr_ini']
        self.pphi_ini = self.pphi_ini if 'pphi_ini' not in kwargs else kwargs['pphi_ini']
        self.ptheta_ini = self.ptheta_ini if 'ptheta_ini' not in kwargs else kwargs['ptheta_ini']
        self.flagTuning = self.flagTuning if 'flagTuning' not in kwargs else kwargs['flagTuning']
        self.tStepBack = self.tStepBack if 'tStepBack' not in kwargs else kwargs['tStepBack']
        self.sl_p = self.sl_p if 'sl_p' not in kwargs else kwargs['sl_p']
        self.x0 = np.cos(self.inc)
        self.code_version = self.code_version if 'code_version' not in kwargs else kwargs['code_version']
        self.log_level = self.log_level if 'log_level' not in kwargs else kwargs['log_level']
        self.risco = self.risco if 'risco' not in kwargs else kwargs['risco']
        self.KK = self.KK if 'KK' not in kwargs else kwargs['KK']
        self.dSS = self.dSS if 'dSS' not in kwargs else kwargs['dSS']
        self.dSO = self.dSO if 'dSO' not in kwargs else kwargs['dSO']
        self.dtPeak = self.dtPeak if 'dtPeak' not in kwargs else kwargs['dtPeak']
        self.t0 = self.t0 if 't0' not in kwargs else kwargs['t0']
        self.egw_flag = self.egw_flag if 'egw_flag' not in kwargs else kwargs['egw_flag']
        self.ret_dyn = self.ret_dyn if 'ret_dyn' not in kwargs else kwargs['ret_dyn']
        self.EPS_REL = self.EPS_REL if 'EPS_REL' not in kwargs else kwargs['EPS_REL']
        self.EPS_ABS = self.EPS_ABS if 'EPS_ABS' not in kwargs else kwargs['EPS_ABS']
        if 'is_coframe' in kwargs:
            self.is_coframe = 0 if kwargs['is_coframe'] == False else 1
        if 'use_coaphase' in kwargs:
            self.use_coaphase = 0 if kwargs['use_coaphase'] == False else 1
        self.zeta = self.zeta if 'zeta' not in kwargs else kwargs['zeta']*np.pi/180
        if 'zeta_rad' in kwargs:
            self.zeta = kwargs['zeta_rad']
        # zeta in [0, 2pi]
        self.zeta = np.pi * ((self.zeta / np.pi) % 2)
        self.xi = self.xi if 'xi' not in kwargs else kwargs['xi']*np.pi/180
        if 'xi_rad' in kwargs:
            self.xi = kwargs['xi_rad']
        if 'zero_dyncoaphase' in kwargs:
            self.zero_dyncoaphase = 0 if kwargs['zero_dyncoaphase'] == False else 1
        self.f_max = self.f_max if 'f_max' not in kwargs else kwargs['f_max']
        self.t_max = self.t_max if 't_max' not in kwargs else kwargs['t_max']
        self.ICValues = self.ICValues if 'initial_condition' not in kwargs else kwargs['initial_condition']
    def calculate_hcorrections(self, l:int, m:int, dyn:npDynamicData):
        hparams = pyHyperParams()
        hparams.d_ini = self.d_ini
        hparams.pr_ini = self.pr_ini
        hparams.ptheta_ini = self.ptheta_ini
        hparams.flagTuning = self.flagTuning
        hparams.tStepBack = self.tStepBack
        hparams.sl_p = self.sl_p
        hparams.x0 = self.x0
        hparams.KK = self.KK
        hparams.dSS = self.dSS
        hparams.dSO = self.dSO
        hparams.dtPeak = self.dtPeak
        hparams.flagZframe = 0
        params = my_waveform_lib.CreateSpinEOBParams(ctypes.c_double(self.m1), ctypes.c_double(self.m2),
                                                     ctypes.c_double(self.s1x), ctypes.c_double(self.s1y), 
                                                     ctypes.c_double(self.s1z), ctypes.c_double(self.s2x),
                                                     ctypes.c_double(self.s2y), ctypes.c_double(self.s2z),
                                                     ctypes.c_double(self.ecc), hparams)
        xVec = convert_ndarray_to_REAL8Vector(dyn.x)
        yVec = convert_ndarray_to_REAL8Vector(dyn.y)
        zVec = convert_ndarray_to_REAL8Vector(dyn.z)

        vxVec = convert_ndarray_to_REAL8Vector(dyn.vx)
        vyVec = convert_ndarray_to_REAL8Vector(dyn.vy)
        vzVec = convert_ndarray_to_REAL8Vector(dyn.vz)

        pTxVec = convert_ndarray_to_REAL8Vector(dyn.pTx)
        pTyVec = convert_ndarray_to_REAL8Vector(dyn.pTy)
        pTzVec = convert_ndarray_to_REAL8Vector(dyn.pTz)

        s1xVec = convert_ndarray_to_REAL8Vector(dyn.s1Vec[:,0])
        s1yVec = convert_ndarray_to_REAL8Vector(dyn.s1Vec[:,1])
        s1zVec = convert_ndarray_to_REAL8Vector(dyn.s1Vec[:,2])

        s2xVec = convert_ndarray_to_REAL8Vector(dyn.s2Vec[:,0])
        s2yVec = convert_ndarray_to_REAL8Vector(dyn.s2Vec[:,1])
        s2zVec = convert_ndarray_to_REAL8Vector(dyn.s2Vec[:,2])

        hamVec = convert_ndarray_to_REAL8Vector(dyn.Hreal)
        prTDotVec = convert_ndarray_to_REAL8Vector(dyn.prTDot)
        rholm_real = ctypes.POINTER(pyREAL8Vector)()
        rholm_imag = ctypes.POINTER(pyREAL8Vector)()
        flm_real = ctypes.POINTER(pyREAL8Vector)()
        flm_imag = ctypes.POINTER(pyREAL8Vector)()
        ret = my_waveform_lib.prec_calculateSEOBFactorizedWaveformCorrectionFromDynVectors(params, ctypes.c_int(l), ctypes.c_int(m),
                                xVec, yVec, zVec, 
                                vxVec, vyVec, vzVec, pTxVec, pTyVec, pTzVec, s1xVec, s1yVec, s1zVec,
                                s2xVec, s2yVec, s2zVec, hamVec, prTDotVec, 
                                ctypes.byref(rholm_real), ctypes.byref(rholm_imag),
                                ctypes.byref(flm_real), ctypes.byref(flm_imag))
        if ret != 0:
            return None
        rholm_real_np = convert_REAL8Vector_to_numpy(rholm_real)
        rholm_imag_np = convert_REAL8Vector_to_numpy(rholm_imag)
        flm_real_np = convert_REAL8Vector_to_numpy(flm_real)
        flm_imag_np = convert_REAL8Vector_to_numpy(flm_imag)
        my_waveform_lib.DestroySpinEOBParams(params)
        my_waveform_lib.DestroyREAL8Vector(xVec)
        my_waveform_lib.DestroyREAL8Vector(yVec)
        my_waveform_lib.DestroyREAL8Vector(zVec)
        my_waveform_lib.DestroyREAL8Vector(vxVec)
        my_waveform_lib.DestroyREAL8Vector(vyVec)
        my_waveform_lib.DestroyREAL8Vector(vzVec)
        my_waveform_lib.DestroyREAL8Vector(pTxVec)
        my_waveform_lib.DestroyREAL8Vector(pTyVec)
        my_waveform_lib.DestroyREAL8Vector(pTzVec)
        my_waveform_lib.DestroyREAL8Vector(s1xVec)
        my_waveform_lib.DestroyREAL8Vector(s1yVec)
        my_waveform_lib.DestroyREAL8Vector(s1zVec)
        my_waveform_lib.DestroyREAL8Vector(s2xVec)
        my_waveform_lib.DestroyREAL8Vector(s2yVec)
        my_waveform_lib.DestroyREAL8Vector(s2zVec)
        my_waveform_lib.DestroyREAL8Vector(hamVec)
        my_waveform_lib.DestroyREAL8Vector(prTDotVec)
        my_waveform_lib.DestroyREAL8Vector(rholm_real)
        my_waveform_lib.DestroyREAL8Vector(rholm_imag)
        my_waveform_lib.DestroyREAL8Vector(flm_real)
        my_waveform_lib.DestroyREAL8Vector(flm_imag)
        my_waveform_lib.CheckMemoryLeak()
        return rholm_real_np + 1.j*rholm_imag_np, flm_real_np + 1.j*flm_imag_np
    
    def run(self):
        if self.ICValues is not None and not isinstance(self.ICValues, np.ndarray):
            raise TypeError('initial condition must be numpy array with size (14,)')
        if self.ICValues is None:
            initValues = ctypes.POINTER(pyREAL8Vector)() # should be a NULL pointer
        else:
            initValues = convert_ndarray_to_REAL8Vector(self.ICValues)
        value_list = [ctypes.c_int(self.use_geom),
                      ctypes.c_double(self.m1),
                      ctypes.c_double(self.m2),
                      ctypes.c_double(self.s1x),
                      ctypes.c_double(self.s1y),
                      ctypes.c_double(self.s1z),
                      ctypes.c_double(self.s2x),
                      ctypes.c_double(self.s2y),
                      ctypes.c_double(self.s2z),
                      ctypes.c_double(self.ecc),
                      ctypes.c_double(self.inc),
                      ctypes.c_double(self.phiRef),
                      ctypes.c_double(self.beta),
                      ctypes.c_double(self.f_min),
                      ctypes.c_double(self.distance),
                      ctypes.c_double(self.srate),
                      ctypes.c_double(self.deltaT),
                      ctypes.c_int(self.is_only22),
                      ctypes.c_int(self.is_constp),
                      ctypes.c_int(self.conserve_flag),
                      ctypes.c_double(self.conserve_time),
                      ctypes.c_int(self.is_noringdown),
                      ctypes.c_int(self.prec_flag),
                      ctypes.c_int(self.debug_id),
                      ctypes.c_double(self.d_ini),
                      ctypes.c_double(self.pr_ini),
                      ctypes.c_double(self.pphi_ini),
                      ctypes.c_double(self.ptheta_ini),
                      ctypes.c_int(self.flagTuning),
                      ctypes.c_double(self.tStepBack),
                      ctypes.c_double(self.sl_p),
                      ctypes.c_double(self.x0),
                      ctypes.c_int(self.code_version),
                      ctypes.c_int(self.log_level),
                      ctypes.c_double(self.risco),
                      ctypes.c_double(self.KK),
                      ctypes.c_double(self.dSS),
                      ctypes.c_double(self.dSO),
                      ctypes.c_double(self.dtPeak),
                      ctypes.c_int(self.egw_flag),
                      ctypes.c_int(self.ret_dyn),
                      ctypes.c_double(self.EPS_REL),
                      ctypes.c_double(self.EPS_ABS),
                      ctypes.c_int(self.is_coframe),
                      ctypes.c_int(self.use_coaphase),
                      ctypes.c_double(self.zeta),
                      ctypes.c_double(self.xi),
                      ctypes.c_double(self.Mf_ref),
                      ctypes.c_int(self.zero_dyncoaphase),
                      ctypes.c_double(self.f_max),
                      ctypes.c_double(self.t_max),
                      initValues
                      ]
        input_pms = pyInputParams(*value_list)
        ret_struct = ctypes.POINTER(pyOutputStruct)()
        ret_dynstruct = ctypes.POINTER(pyDynOutputStruct)()
        ret = my_waveform_lib.generate_waveform(ctypes.byref(input_pms), ctypes.byref(ret_struct), ctypes.byref(ret_dynstruct))
        if self.ICValues is not None:
            my_waveform_lib.DestroyREAL8Vector(initValues)
        npdynstruct = None
        if ret == 0:
            npstruct = npWaveformData(ret_struct, t0 = self.t0)
            if self.ret_dyn:
                npdynstruct = npDynamicData(self.m1, self.m2, ret_dynstruct)
        else:
            npstruct = None
        my_waveform_lib.DestroypyOutputStruct_t(ret_struct)
        my_waveform_lib.DestroypyDynOutputStruct_t(ret_dynstruct)
        my_waveform_lib.CheckMemoryLeak()
        return npstruct, npdynstruct
    
def calculate_rho22_from_npdynstruct(dyn:npDynamicData):
    pass

def DestroySpinEOBParams(params:ctypes.POINTER(pyNewtonMultipolePrefixes)):
    my_waveform_lib.DestroySpinEOBParams(params)
    return

my_waveform_lib.test_interface.restype = ctypes.POINTER(pyNewtonMultipolePrefixes)
def test_interface(m1:float, m2:float):
    ret = my_waveform_lib.test_interface(ctypes.c_double(m1), ctypes.c_double(m2))
    for l in range(2,8):
        for m in range(1,l+1):
            print(ret.contents.values[l][m].real + 1.j*ret.contents.values[l][m].imag)
    my_waveform_lib.myFree(ret)
    my_waveform_lib.CheckMemoryLeak()
    return ret


