#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 10:25:38 2019

@author: drizl
"""

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from pathlib import Path
from .detector import Detector
from ..pyUtils import rTimeSeries, cTimeSeries

__all__ = ['DetectorPSD', 'Detector']

LOC = Path(__file__).parent

#-----switch method-----#
class switch(object):
    def __init__(self, value):
        self.value = value
        self.fall = False

    def __iter__(self):
        """Return the match method once, then stop"""
        yield self.match
        raise StopIteration

    def match(self, *args):
        """Indicate whether or not to enter a case suite"""
        if self.fall or not args:
            return True
        elif self.value in args: # changed for v1.5, see below
            self.fall = True
            return True
        else:
            return False

class DetectorPSD(object):
    def __init__(self, name = None, flow = 0, fhigh = None):
        if isinstance(name, DetectorPSD):
            name = name.name
        self._name = name
        self._choose_psd(flow, fhigh)
                
    def __call__(self, *args, **kwargs):
        return self._psd(*args, **kwargs)

    @property
    def frange(self):
        return self._psd(1, True)
    
    def _choose_psd(self, flow = 0, fhigh = None):
        file = None
        for case in switch(self._name):
            if case('ET'):
                file = LOC / 'LIGO-P1600143-v18-ET_D.txt'
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh)
                break
            if case('ET_fit'):
                file = None
                self._psd = PSD_ET_fit
                break
            if case('ET_D_1'):
                file = LOC / 'ET-D-1-0000A-18.txt'
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh)
                break
            if case('ET_D_2'):
                file = LOC / 'ET-D-2-0000A-18.txt'
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh)
                break
            if case('ET_D_3'):
                file = LOC / 'ET-D-3-0000A-18.txt'
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh)
                break
            if case('CE_Pes'):
                file = LOC / 'LIGO-P1600143-v18-CE_Pessimistic.txt'
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh)
                break
            if case('CE_Wide'):
                file = LOC / 'LIGO-P1600143-v18-CE_Wideband.txt'
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh)
                break
            if case('CE_LAL'):
                file = LOC / 'LIGO-P1600143-v18-CE.txt'
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh)
                break
            if case('CE'):
                file = LOC / 'cosmic_explorer.txt'
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh)
                break
            if case('CE_20km_pm'):
                file = LOC / 'cosmic_explorer_20km_pm.txt'
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh)
                break
            if case('CE_20km'):
                file = LOC / 'cosmic_explorer_20km.txt'
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh)
                break
            if case('CE_40km_lf'):
                file = LOC / 'cosmic_explorer_40km_lf.txt'
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh)
                break
            if case('advLIGO'):
                file = LOC / 'LIGO-P1200087-v18-aLIGO_DESIGN.txt'
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh)
                break
            if case('advLIGO_fit'):
                file = None
                self._psd = PSD_advLIGO_fit
                break
            if case('advLIGO_zerodethp'):
                file = LOC /"ZERO_DET_high_P.txt"
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh)
                break
            if case('advLIGO_O4_high'):
                file = LOC /"aligo_O4high.txt"
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh)
                break
            if case('advLIGO_O4_low'):
                file = LOC /"aligo_O4low.txt"
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh)
                break
            if case('advVirgo_O4_high'):
                file = LOC /"avirgo_O4high_NEW.txt"
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh)
                break
            if case('advVirgo_O5_high'):
                file = LOC /"avirgo_O5high_NEW.txt"
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh)
                break
            if case('advVirgo_O5_low'):
                file = LOC /"avirgo_O5low_NEW.txt"
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh)
                break
            if case('L1'):
                file = LOC / 'LIGOLivingston_O3PSD-1241571618-21600.txt'
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh, exp = False)
                break
            if case('H1'):
                file = LOC / 'LIGOHanford_O3PSD-1241571618-21600.txt'
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh, exp = False)
                break
            if case('V1'):
                file = LOC / 'Virgo_O3PSD-1241571618-21600.txt'
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh, exp = False)
                break
            if case('cut'):
                self._psd = get_lowCutPSD(flow = flow, fhigh = fhigh)
                break
            if case('LISA'):
                self._psd = get_PSD_Space_fit(name = 'LISA', flow = flow, fhigh = fhigh)
                break
            if case('Taiji'):
                self._psd = get_PSD_Space_fit(name = 'Taiji', flow = flow, fhigh = fhigh)
                break
            if case('Tianqin'):
                self._psd = get_PSD_Space_fit(name = 'Tianqin', flow = flow, fhigh = fhigh)
                break
            if case('DECIGO_fit'):
                self._psd = PSD_DECIGO_fit
                break
            if case('DECIGO_fit2'):
                self._psd = PSD_DECIGO_fit2
                break
            if case('DECIGO_fit3'):
                self._psd = PSD_DECIGO_fit3
                break
            if case('KAGRA_25Mpc'):
                file = LOC / 'kagra_25Mpc.txt'
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh, exp = True)
                break
            if case('KAGRA_80Mpc'):
                file = LOC / 'kagra_80Mpc.txt'
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh, exp = True)
                break
            if case('KAGRA_128Mpc'):
                file = LOC / 'kagra_128Mpc.txt'
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh, exp = True)
                break
            if case('H1_O3'):
                file = LOC / 'aligo_O3actual_H1.txt'
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh, exp = True)
                break
            if case('L1_O3'):
                file = LOC / 'aligo_O3actual_L1.txt'
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh, exp = True)
                break
            if case('V1_O3'):
                file = LOC / 'avirgo_O3actual.txt'
                self._psd = loadPSD_from_file(file, flow, fhigh = fhigh, exp = True)
                break
            self._psd = lambda x : 1
        self._file = file
            
    @property
    def name(self):
        return self._name
        
    def get_psd_data(self, exp = True):
        if self._file is None:
            return None
        data = np.loadtxt(self._file)
        freq = data[:,0]
        h = data[:,1]
        valpsd = np.zeros(len(h))
        idxsift = np.where(h > 0)
        if exp:
            valpsd[idxsift] = np.exp(2 * np.log(h[idxsift]))
        else:
            valpsd[idxsift] = h[idxsift]
        return freq, valpsd

    def _sim_noise_seg(self, seg, psd_data, srate):
        segLen = len(seg)
        df = srate / segLen
        sigma = np.sqrt(np.abs(psd_data) / df) / 2
        stilde = sigma * (np.random.randn(len(sigma)) + 1.j*np.random.randn(len(sigma)))
        return np.fft.irfft(stilde) * srate

    def _sim_noise(self, stride, psd_data, seg, srate):
        segLen = len(seg)
        if stride == 0:
            return self._sim_noise_seg(seg, psd_data, srate)
        elif stride == segLen:
            seg = self._sim_noise_seg(seg, psd_data, srate)
            stride = 0
        overlap = seg[stride:]
        lenolp = len(overlap)
        seg = self._sim_noise_seg(seg, psd_data, srate)
        for i in range(len(overlap)):
            x = np.cos(np.pi * i / (2*lenolp))
            y = np.sin(np.pi * i / (2*lenolp))
            seg[i] = x * overlap[i] + y * seg[i]
        return seg

    def generate_noise(self, duration, srate):
        ncount = duration * srate
        segdur = 4
        length = int(segdur * srate)
        df = srate / length
        freqs = np.arange(0.0, (length/2+1)*df, df)
        psd_data = self.__call__(freqs)
        psd_data[np.isinf(psd_data)]=0.0
        psd_data[np.isnan(psd_data)]=0.0
        stride = int(length / 2)
        seg = np.zeros(length)
        ret = []
        seg = self._sim_noise(0, psd_data, seg, srate)
        while(1):
            for j in range(0, stride):
                ncount -= 1
                if ncount == 0:
                    return np.asarray(ret)
                ret.append(seg[j])
            seg = self._sim_noise(stride, psd_data, seg, srate)
            # if ret is None:
            #     ret = seg.copy().tolist()
            # else:
            #     ret = np.append(ret, seg).tolist()


    

def loadPSD_from_file(file, flow = 0, fhigh = None, exp = True):
    if fhigh is None:
        fhigh = np.inf
    data = np.loadtxt(file)
    freq = data[:,0]
    h = data[:,1]
    valpsd = np.zeros(len(h))
    idxsift = np.where(h > 0)
    if exp:
        valpsd[idxsift] = np.exp(2 * np.log(h[idxsift]))
    else:
        valpsd[idxsift] = h[idxsift]
    func = InterpolatedUnivariateSpline(freq, valpsd)
    if flow < freq[0]:
        flow = freq[0]
    if fhigh > freq[-1]:
        fhigh = freq[-1]
    def funcPSD(freq, is_get_frange = False):
        if is_get_frange:
            return (flow, fhigh)
        ret = func(freq)
        if hasattr(freq, '__len__'):
            ret[np.where(freq < flow)] = np.inf
            ret[np.where(freq > fhigh)] = np.inf
        elif freq < flow:
            ret = np.inf
        elif freq > fhigh:
            ret = np.inf
        return ret
    return funcPSD
        
def get_lowCutPSD(flow = 0, fhigh = None):
    if fhigh is None:
        fhigh = np.inf
    def funcPSD(freq, is_get_frange = False):
        if is_get_frange:
            return (flow, fhigh)
        if hasattr(freq, '__len__'):
            ret = np.ones(len(freq))
            ret[np.where(freq < flow)] = np.inf
            ret[np.where(freq > fhigh)] = np.inf
        elif freq < flow:
            ret = np.inf
        elif freq > fhigh:
            ret = np.inf
        return ret
    return funcPSD
    

"""
"""
def PSD_DECIGO_fit(f, is_get_frange = False):
    if is_get_frange:
        return (0, np.inf)
    fac_f_4 = np.power(f, -4)
    fac_f2 = np.power(f/7.36, 2)
    return 6.53e-49 * ( 1+fac_f2 ) + \
        4.45e-51*fac_f_4 / (1 + 1/fac_f2) + \
        4.49e-52*fac_f_4

def PSD_DECIGO_fit2(f, is_get_frange = False):
    if is_get_frange:
        return (0, np.inf)
    fac_f_4 = np.power(f, -4)
    fac_f2 = np.power(f, 2)
    return 1.25e-47 + 4.21e-50 * fac_f_4 + 3.92e-49 * fac_f2

def PSD_DECIGO_fit3(f, is_get_frange = False):
    if is_get_frange:
        return (0, np.inf)
    fac_f_4 = np.power(f, -4)
    fac_f2 = np.power(f, 2)
    return 1.88e-48 + 6.31e-51 * fac_f_4 + 5.88e-50 * fac_f2


"""
	arXiv:0908.0353
"""
def PSD_ET_fit(f, is_get_frange = False):
    if is_get_frange:
        return (1, np.inf)
    x = f / 200
    x2 = np.power(x,2)
    x3 = np.power(x,3)
    x4 = np.power(x,4)
    x5 = np.power(x,5)
    x6 = np.power(x,6)
    p1 = np.power(x, -4.05)
    a1p2 = 185.62 * np.power(x, -0.69)
    ret =  1.449e-52 * (p1 + a1p2 + 232.56 * \
                        (1 + 31.18*x - 46.72*x2 + 52.24*x3 - 42.16*x4 + 10.17 * x5 + 11.53 * x6) / \
                        (1 + 13.58*x - 36.46*x2 + 18.56*x3 + 27.43*x4) )
    if hasattr(f, '__len__'):
        f = np.asarray(f)
        if f.any() < 1:
            idx = np.where(f < 1)[0][0]
            ret[idx] = np.inf
        return ret
    elif f < 1:
        return np.inf
    else:
        return ret
    
    


"""
    PhysRevD.71.084008
"""
def PSD_advLIGO_fit(f, is_get_frange = False):
    if is_get_frange:
        return (10, np.inf)
    fac = f / 215
    fac2 = np.power(fac,2)
    fac4 = np.power(fac,4)
    ret = 1e-49 * ( np.power(fac, -4.14) - 5/fac2 + \
                    111 * ( 1 - fac2 + fac4/2 ) / (1 + fac2/2) )
    if hasattr(ret, '__len__'):
        f = np.asarray(f)
        if f.any() < 10:
            ret[np.where(f < 10)] = np.inf
        return ret
    elif f < 10:
        return np.inf
    else:
        return ret

"""
    CQG. 36.105011
"""
def get_PSD_Space_fit(name = 'LISA', flow = 0, fhigh = None):
    C_SI = 299792458
    name_low = name.lower()
    if fhigh is None:
        fhigh = np.inf
    for case in switch(name_low):
        if case('lisa'):
            P_OMS = np.power(1.5e-11, 2)
            fP_acc = lambda f : np.power(3.e-15, 2) * (1 + np.power(4.e-4/f,2))
            L = 2.5e9
            break
        if case('taiji'):
            P_OMS = np.power(8e-11, 2)
            fP_acc = lambda f : np.power(3.e-15, 2) * (1 + np.power(4.e-4/f,2))
            L = 3.e9
            break
        if case('tianqin'):
            P_OMS = np.power(1.e-12, 2)
            fP_acc = lambda f : np.power(1.e-15, 2) * (1 + np.power(1.e-4/f,2))
            L = np.sqrt(3) * 1e8
            break
        raise Exception(f'Unrecognized PSD name {name}')
    def func_PSD(f, is_get_frange = False):
        if is_get_frange:
            return (flow, fhigh)
        fstar = C_SI/(2*np.pi*L)
        ret = (10. / (3.*L*L)) * \
            (P_OMS + 2*(1+ np.power(np.cos(f/fstar),2))*\
            (fP_acc(f)/np.power(2*np.pi*f,4)) ) * \
            (1 + 6.*np.power(f/fstar,2)/10.)
        if hasattr(f, '__len__'):
            f = np.asarray(f)
            ret[np.where(f < flow)] = np.inf
            ret[np.where(f > fhigh)] = np.inf
        elif f < flow:
            ret = np.inf
        elif f > fhigh:
            ret = np.inf
        return ret
    
    return func_PSD
    

class DetectorInfo(object):
    def __init__(self, latitude, longitude, orientation, aperture):
        self.__latitude = latitude
        self.__longitude = longitude
        self.__orientation = orientation
        self.__aperture = aperture

    @property
    def latitude(self):
        return self.__latitude
    @property
    def lbd(self):
        return self.__latitude

    @property
    def longitude(self):
        return self.__longitude
    @property
    def varphi(self):
        return self.__longitude

    @property
    def orientation(self):
        return self.__orientation
    @property
    def gamma(self):
        return self.__orientation

    @property
    def aperture(self):
        return self.__aperture
    @property
    def zeta(self):
        return self.__aperture
# 2207.02771
Detector_CE1 = DetectorInfo(latitude = 46.5*np.pi/180, longitude = -119.4*np.pi/180, orientation = 171*np.pi/180, aperture = np.pi/2)
Detector_CE2 = DetectorInfo(latitude = 30.6*np.pi/180, longitude = -90.8*np.pi/180, orientation = 242.7*np.pi/180, aperture = np.pi/2)
Detector_ET1 = DetectorInfo(latitude = 40.5*np.pi/180, longitude = 9.4*np.pi/180, orientation = 0*np.pi/180, aperture = np.pi/3)
Detector_ET2 = DetectorInfo(latitude = 40.5*np.pi/180, longitude = 9.4*np.pi/180, orientation = 120*np.pi/180, aperture = np.pi/3)
Detector_ET3 = DetectorInfo(latitude = 40.5*np.pi/180, longitude = 9.4*np.pi/180, orientation = 240*np.pi/180, aperture = np.pi/3)
Detector_H1 = DetectorInfo(latitude = 46.5*np.pi/180, longitude = -119.4*np.pi/180, orientation = 171*np.pi/180, aperture = np.pi/2)
Detector_L1 = DetectorInfo(latitude = 30.6*np.pi/180, longitude = -90.8*np.pi/180, orientation = 242.7*np.pi/180, aperture = np.pi/2)
Detector_Virgo = DetectorInfo(latitude = 43.6*np.pi/180, longitude = 10.5*np.pi/180, orientation = 115.6*np.pi/180, aperture = np.pi/2)
Detector_KAGRA = DetectorInfo(latitude = 36.4*np.pi/180, longitude = 137.3*np.pi/180, orientation = 15.4*np.pi/180, aperture = np.pi/2)
DetectorDict = {'CE1':Detector_CE1,
                'CE2':Detector_CE2,
                'ET1': Detector_ET1,
                'ET2': Detector_ET2,
                'ET3': Detector_ET3,
                'H1': Detector_H1,
                'L1': Detector_L1,
                'Virgo':Detector_Virgo,
                'KAGRA':Detector_KAGRA}

c_SI = 299792458 # Speed Of Light [m s^-1]
class AntennaPatternFCalculator(object):
    def __init__(self, psi, ra, dec, latitude, longitude, gamma, zeta = np.pi/2.):
        self.__psi = psi
        self.__alpha = ra
        self.__delta = dec
        self.__lbd = latitude
        self.__varphi = longitude
        self.__gamma = gamma
        self.__zeta = zeta

        self.__alphi = self.__alpha - self.__varphi
        self.__sin_zeta = np.sin(self.__zeta)
        self.__sin_2psi = np.sin(2.*self.__psi)
        self.__cos_2psi = np.cos(2.*self.__psi)
        self.__sin_2alphi = np.sin(2.*self.__alphi)
        self.__cos_2alphi = np.cos(2.*self.__alphi)
        self.__sin_alphi = np.sin(self.__alphi)
        self.__cos_alphi = np.cos(self.__alphi)
        self.__sin_2gamma = np.sin(2.*self.__gamma)
        self.__cos_2gamma = np.cos(2.*self.__gamma)
        self.__sin_2delta = np.sin(2.*self.__delta)
        self.__cos_2delta = np.cos(2.*self.__delta)
        self.__sin_delta = np.sin(self.__delta)
        self.__cos_delta = np.cos(self.__delta)
        self.__sin_2lambda = np.sin(2.*self.__lbd)
        self.__cos_2lambda = np.cos(2.*self.__lbd)
        self.__sin_lambda = np.sin(self.__lbd)
        self.__cos_lambda = np.cos(self.__lbd)
        self.__3MinusCos2Lambda = 3.-self.__cos_2lambda
        self.__3MinusCos2Delta = 3.-self.__cos_2delta

        self.__Gplus = (1./16.)*self.__sin_2gamma*self.__3MinusCos2Lambda*self.__3MinusCos2Delta*self.__cos_2alphi \
            -0.25*self.__cos_2gamma*self.__sin_lambda*self.__3MinusCos2Delta*self.__sin_2alphi \
            +0.25*self.__sin_2gamma*self.__sin_2lambda*self.__sin_2delta*self.__cos_alphi \
            -0.5*self.__cos_2gamma*self.__cos_lambda*self.__sin_2delta*self.__sin_alphi \
            +0.75*self.__sin_2gamma*self.__cos_lambda*self.__cos_lambda*self.__cos_delta*self.__cos_delta
        self.__Gcross = self.__cos_2gamma*self.__sin_lambda*self.__sin_delta*self.__cos_2alphi \
            +0.25*self.__sin_2gamma*self.__3MinusCos2Lambda*self.__sin_delta*self.__sin_2alphi \
            +self.__cos_2gamma*self.__cos_lambda*self.__cos_delta*self.__cos_alphi \
            +0.5*self.__sin_2gamma*self.__sin_2lambda*self.__cos_delta*self.__sin_alphi
        dt_earth = 6371e3/c_SI
        self.__OmegaEarth = 2.*np.pi/(24.*60.*60.)
        self.__delta_t = dt_earth * (self.__cos_lambda*self.__sin_delta + self.__cos_delta*self.__sin_lambda*self.__cos_alphi)

    @property
    def Gplus(self):
        return self.__Gplus
    def Gplus_t(self, t_SI):
        alphi_t = self.__alphi - self.__OmegaEarth * t_SI
        cos_2alphi_t = np.cos(2.*alphi_t)
        sin_2alphi_t = np.sin(2.*alphi_t)
        cos_alphi_t = np.cos(alphi_t)
        sin_alphi_t = np.sin(alphi_t)
        return (1./16.)*self.__sin_2gamma*self.__3MinusCos2Lambda*self.__3MinusCos2Delta*cos_2alphi_t \
            -0.25*self.__cos_2gamma*self.__sin_lambda*self.__3MinusCos2Delta*sin_2alphi_t \
            +0.25*self.__sin_2gamma*self.__sin_2lambda*self.__sin_2delta*cos_alphi_t \
            -0.5*self.__cos_2gamma*self.__cos_lambda*self.__sin_2delta*sin_alphi_t \
            +0.75*self.__sin_2gamma*self.__cos_lambda*self.__cos_lambda*self.__cos_delta*self.__cos_delta

    @property
    def Gcross(self):
        return self.__Gcross
    def Gcross_t(self, t_SI):
        alphi_t = self.__alphi - self.__OmegaEarth * t_SI
        cos_2alphi_t = np.cos(2.*alphi_t)
        sin_2alphi_t = np.sin(2.*alphi_t)
        cos_alphi_t = np.cos(alphi_t)
        sin_alphi_t = np.sin(alphi_t)
        return self.__cos_2gamma*self.__sin_lambda*self.__sin_delta*cos_2alphi_t \
                    +0.25*self.__sin_2gamma*self.__3MinusCos2Lambda*self.__sin_delta*sin_2alphi_t \
                    +self.__cos_2gamma*self.__cos_lambda*self.__cos_delta*cos_alphi_t \
                    +0.5*self.__sin_2gamma*self.__sin_2lambda*self.__cos_delta*sin_alphi_t

    @property
    def Fplus(self):
        return self.__sin_zeta*self.Gplus
    def Fplus_t(self, t_SI):
        return self.__sin_zeta*self.Gplus_t(t_SI)
    
    @property
    def Fcross(self):
        return self.__sin_zeta*self.Gcross
    def Fcross_t(self, t_SI):
        return self.__sin_zeta*self.Gcross_t(t_SI)

    @property
    def barFplus(self):
        return self.Fplus*self.__cos_2psi + self.Fcross*self.__sin_2psi
    def barFplus_t(self, t_SI):
        return self.Fplus_t(t_SI)*self.__cos_2psi + self.Fcross_t(t_SI)*self.__sin_2psi
    
    @property
    def barFcross(self):
        return self.Fcross*self.__cos_2psi - self.Fplus*self.__sin_2psi
    def barFcross_t(self, t_SI):
        return self.Fcross_t(t_SI)*self.__cos_2psi - self.Fplus_t(t_SI)*self.__sin_2psi

    def detector_strain(self, hpc:cTimeSeries, tc = 0.0):
        Fp = self.barFplus
        Fc = self.barFcross
        hp = hpc.real
        hc = -hpc.imag
        tini = tc + self.__delta_t
        h = hp * Fp + hc * Fc
        return rTimeSeries(hpc.time + tini, h)

    def detector_strain_t(self, hpc:cTimeSeries, tc = 0.0):
        hp = hpc.real
        hc = -hpc.imag
        tini = tc + self.__delta_t
        time = hpc.time + tini
        Fp = self.barFplus_t(time)
        Fc = self.barFcross_t(time)
        h = hp * Fp + hc * Fc
        return rTimeSeries(time, h)


class GWDetector(object):
    def __init__(self, name):
        self.__name = name
        self.__choose_psd_and_detector()

    def __choose_psd_and_detector(self):
        for case in switch(self.__name):
            if case('advLIGO-H1'):
                self.__det = Detector('H1')
                self.__xdet = DetectorDict['H1']
                self.__psd = DetectorPSD('advLIGO_zerodethp')
                break
            if case('advLIGO-L1'):
                self.__det = Detector('L1')
                self.__xdet = DetectorDict['L1']
                self.__psd = DetectorPSD('advLIGO_zerodethp')
                break
            if case('H1') or case('LIGO-H1') or case('LIGO-O4-H1'):
                self.__det = Detector('H1')
                self.__xdet = DetectorDict['H1']
                self.__psd = DetectorPSD('advLIGO_O4_high')
                break
            if case('L1') or case('LIGO-L1') or case('LIGO-O4-L1'):
                self.__det = Detector('L1')
                self.__xdet = DetectorDict['L1']
                self.__psd = DetectorPSD('advLIGO_O4_high')
                break
            if case('Virgo') or case('V1') or case('O4-V1'):
                self.__det = Detector('V1')
                self.__xdet = DetectorDict['Virgo']
                self.__psd = DetectorPSD('advVirgo_O4_high')
                break
            if case('KAGRA') or case('K1') or case('KAGRA-80Mpc') or case('K1-80Mpc'):
                self.__det = Detector('K1')
                self.__xdet = DetectorDict['KAGRA']
                self.__psd = DetectorPSD('KAGRA_80Mpc')
                break
            if case('KAGRA-25Mpc') or case('K1-25Mpc'):
                self.__det = Detector('K1')
                self.__xdet = DetectorDict['K1']
                self.__psd = DetectorPSD('KAGRA_25Mpc')
                break
            if case('KAGRA-128Mpc') or case('K1-128Mpc'):
                self.__det = Detector('K1')
                self.__xdet = DetectorDict['K1']
                self.__psd = DetectorPSD('KAGRA_128Mpc')
                break
            if case('ET1') or case('E1'):
                self.__det = Detector('E1')
                self.__xdet = DetectorDict['E1']
                self.__psd = DetectorPSD('ET_D_3')
                break
            if case('ET2') or case('E2'):
                self.__det = Detector('E2')
                self.__xdet = DetectorDict['E2']
                self.__psd = DetectorPSD('ET_D_3')
                break
            if case('ET3') or case('E3'):
                self.__det = Detector('E3')
                self.__xdet = DetectorDict['E3']
                self.__psd = DetectorPSD('ET_D_3')
                break
            else:
                raise Exception(f'cannot find such detector configuration {self.__name}')
    
    @property
    def name(self):
        return self.__name

    def __call__(self, *args, **kwargs):
        return self.__psd(*args, **kwargs)
    
    @property
    def psd(self):
        return self.__psd

    def antenna_pattern_calculator(self, psi:float, ra:float, dec:float):
        return AntennaPatternFCalculator(psi, ra, dec, self.__xdet.latitude, self.__xdet.longitude, self.__xdet.gamma, self.__xdet.zeta)

    def antenna_pattern_gps(self, psi:float, ra:float, dec:float, gps:float):
        return self.__det.antenna_pattern(ra, dec, psi, gps)
    
    def time_delay(self, ra:float, dec:float, gps:float):
        return self.__det.time_delay_from_earth_center(ra, dec, gps)

    def antenna_pattern_gmst(self, psi:float, ra:float, dec:float, gmst:float):
        return self.__det.antenna_pattern_gmst(ra, dec, psi, gmst)
