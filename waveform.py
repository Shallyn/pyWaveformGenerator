#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed, 29 Mar 2023 13:50:54 +0000

@author: Shallyn
"""
import sys, os
import numpy as np
from . import SEOBNRWaveformCaller
from .psd import GWDetector
from .pyUtils import rTimeSeries


def calculate_waveform(params, f_min, Mf_ref = 0.002, srate = 16384, code_version = 1, **kwargs):
    '''
        Calculate EOB waveform
        INPUT
        params: a list of BBH parameters, 
            (m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, e0, dL, zeta_rad, iota_rad, beta_rad, Phi_rad)
            > m1, m2: component mass of BH in solar mass,
            > s1x, s1y, s1z, s2x, s2y, s2z: dimentionless spin vector chi of BH, the norm of these vectors should less than 1
            > e0 initial eccentricity at given reference frequency Mf_ref
            > dL luminosity distance in Mpc 
            > zeta_rad relativistic anomaly zeta in r = p / (1 + e cos(zeta)),
            > iota_rad inclination angle in rad
            > beta_rad represent the initial direction of major axis of the elliptical orbit, 
                which is equivalent to the phiRef parameter in lalsim-inspiral
            > Phi_rad the complex angle of h+ - ihx at the merge stage (default unused other than you set use_coaphase=True)
        NOTE: If you want to generate waveform that both containing eccentricity (e0 != 0) and spin precession (sx !=0 or sy != 0),
              you need to set prec_flag=3 (see arxiv:2310.04552). 
              However, in this case the initial conditions are limited to the periapsis of the equatorial plane,
              which means that zeta_rad is not supported (and also egw_flag=1 doesn't work in this case). 
        f_min: initial (angular) frequency
        Mf_ref: reference (angular) frequency
        srate: the output sample rate
            NOTE: Lowering this value does not significantly improve computational speed, 
                as the time consuming part of the code is the numerical solution of the dynamic system, 
                where adaptive stepping is used. 
                If you want to sacrifice accuracy to improve calculation speed, 
                you can set EPS_REL and EPS_ABS (in kwargs).
                For spin aligned cases, the default setting is EPS_REL=1e-9, EPS_ABS=1e-10
        code_version: if you set code_version=0, this code would equivalent to SEOBNRv4PHM 
                        (sometimes you may need to add an extra argument is_constp=True)
                      if you set code_version=1, this code will use another eccentricity definition
        
        OUTPUT:
        waveform: a collection of polarizations and GW spherical modes, see more details in __init__.py
            > waveform.timeM: the time in M (ndarray)
            > waveform.time: the time in s (ndarray)
            > waveform.hpc: the complex GW polarization waveform h+ - ihx (cTimeSeries)
            > waveform.h22
            > waveform.h21
            > waveform.h33
            > waveform.h44
            > waveform.h55: the GW spherical modes (cTimeSeries)
                Note.1 for the type cTimeSeries, you can use .real .imag to get real part and imaginary part (ndarray)
                    or .amp and .phase to get amplitude and phase (ndarray)
                Note.2 if waveform = None, the generation failed
        
        dynamics: a collection of dynamical variables, see more details in __init__.py
            > dynamics.timeM: the time in M (ndarray)
            > dynamics.m1 .m2 .eta: component masses and symmetric mass ratio
            > dynamics.rVec: Cartesian coordinates of the orbit (ndarray)
            > dynamics.pTVec: (tortoise) momentum vector of the orbit (ndarray)
            > dynamics.vVec: velocity vector of the orbit (ndarray)
            > dynamics.s1Vec .s2Vec: the spin vectors of the two component BH (ndarray)
            > dynamics.chi1Vec .chi2Vec: the dimentionless spin vectors of the two component BH (ndarray)
            > dynamics.omega: the angular velocity of the orbit (ndarray)
            > dynamics.Hreal: the EOB Hamiltonian of the orbit (ndarray)
    '''
    m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, e0, dL, zeta_rad, iota_rad, beta_rad, Phi_rad = params
    ge = SEOBNRWaveformCaller()
    ge.set_params(m1 = m1, m2 = m2, 
                  s1x = s1x, s1y = s1y, s1z = s1z,
                  s2x = s2x, s2y = s2y, s2z = s2z, srate=srate, Mf_ref = Mf_ref, egw_flag=1,
                  ecc = e0, distance = dL, code_version=code_version, zeta_rad = zeta_rad,
                  inc_rad = iota_rad, beta_rad = beta_rad, phiRef_rad = Phi_rad, f_min = f_min, ret_dyn = True, **kwargs)
    waveform, dynamics = ge.run()
    if waveform is None:
        return None, None
    return waveform, dynamics
    
# unused
def calculate_strain(params, f_min, Mf_ref = 0.002, gps0 = 1356566418, detectors = None, is_only22 = False, log_level = 1, code_version = 2, ret_dyn = False, **kwargs) -> rTimeSeries:
    m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, e0, dL, theta, phi, iota, psi, t_c, beta_c, Phi_c = params
    ge = SEOBNRWaveformCaller()
    is_ret_dyn = 0
    if ret_dyn:
        is_ret_dyn = 1
    ge.set_params(m1 = m1, m2 = m2, 
                  s1x = s1x, s1y = s1y, s1z = s1z,
                  s2x = s2x, s2y = s2y, s2z = s2z, Mf_ref = Mf_ref,
                  ecc = e0, distance = dL, code_version=code_version, egw_flag=1,
                  inc_rad = iota, beta_rad = beta_c, phiRef_rad = Phi_c, f_min = f_min,
                  log_level=log_level, is_only22=is_only22, ret_dyn = is_ret_dyn, **kwargs)
    dynamics = None
    waveform, dynamics = ge.run()
    if waveform is None:
        return None, None, None
    if detectors is None:
        # default use L1
        detectors = [GWDetector('L1')]
    # apCalculator = detector.antenna_pattern(psi, phi, theta)
    # strain = apCalculator.detector_strain_t(waveform.hpc, t_c+gps0)
    hp = waveform.hpc.real
    hc = -waveform.hpc.imag
    tpeak = waveform.hpc.time[waveform.h22.argpeak]
    ret = []
    for detector in detectors:
        dt = detector.time_delay(phi, np.pi-theta, t_c + gps0)
        Fplus, Fcross = detector.antenna_pattern_gps(psi, phi, np.pi-theta, t_c + gps0 + dt)
        h = Fplus * hp + Fcross * hc
        t = waveform.hpc.time - tpeak
        ret.append(rTimeSeries(t, h, t0 = t_c + dt, gps0=gps0))
    return ret, waveform, dynamics


def calculate_waveform_ep(params, f_min, Mf_ref = 0.002, srate = 16384, is_coframe = False, code_version=1, **kwargs):
    '''
        Calculate EOB Eccentric-Precession waveform
        INPUT
        params: a list of BBH parameters, 
            (m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, e0, dL, zeta_rad, iota_rad, beta_rad, Phi_rad)
            > m1, m2: component mass of BH in solar mass,
            > s1x, s1y, s1z, s2x, s2y, s2z: dimentionless spin vector chi of BH, the norm of these vectors should less than 1
            > e0 initial eccentricity at given reference frequency Mf_ref
            > dL luminosity distance in Mpc
            > zeta_rad relativistic anomaly zeta in r = p / (1 + e cos(zeta)),
                you need to set egw_flag=1 to turn on it.
            > iota_rad inclination angle in rad
            > beta_rad represent the initial direction of major axis of the elliptical orbit, 
                which is equivalent to the phiRef parameter in lalsim-inspiral
            > Phi_rad the complex angle of h+ - ihx at the merge stage (default unused other than you set use_coaphase=True)
        f_min: initial (angular) frequency
        Mf_ref: reference (angular) frequency
        srate: the output sample rate
            NOTE: Lowering this value does not significantly improve computational speed, 
                as the time consuming part of the code is the numerical solution of the dynamic system, 
                where adaptive stepping is used. 
                If you want to sacrifice accuracy to improve calculation speed, 
                you can set EPS_REL and EPS_ABS (in kwargs).
                For spin aligned cases, the default setting is EPS_REL=1e-9, EPS_ABS=1e-10
        code_version: if you set code_version=0, this code would equivalent to SEOBNRv4PHM 
                        (sometimes you may need to add an extra argument is_constp=True)
                      if you set code_version=1, this code will use another eccentricity definition
        OUTPUT:
        waveform: a collection of polarizations and GW spherical modes, see more details in __init__.py
        NOTE: if you set is_coframe = True, this waveform is in co-precession frame
            > waveform.timeM: the time in M (ndarray)
            > waveform.time: the time in s (ndarray)
            > waveform.hpc: the complex GW polarization waveform h+ - ihx (cTimeSeries)
            > waveform.h22
            > waveform.h21
            > waveform.h33
            > waveform.h44
            > waveform.h55: the GW spherical modes (cTimeSeries)
                Note.1 for the type cTimeSeries, you can use .real .imag to get real part and imaginary part (ndarray)
                    or .amp and .phase to get amplitude and phase (ndarray)
                Note.2 if waveform = None, the generation failed
        
        dynamics: a collection of dynamical variables, see more details in __init__.py
            > dynamics.timeM: the time in M (ndarray)
            > dynamics.m1 .m2 .eta: component masses and symmetric mass ratio
            > dynamics.rVec: Cartesian coordinates of the orbit (ndarray)
            > dynamics.pTVec: (tortoise) momentum vector of the orbit (ndarray)
            > dynamics.vVec: velocity vector of the orbit (ndarray)
            > dynamics.s1Vec .s2Vec: the spin vectors of the two component BH (ndarray)
            > dynamics.chi1Vec .chi2Vec: the dimentionless spin vectors of the two component BH (ndarray)
            > dynamics.omega: the angular velocity of the orbit (ndarray)
            > dynamics.Hreal: the EOB Hamiltonian of the orbit (ndarray)
    '''
    return calculate_waveform(params, f_min, Mf_ref, srate, code_version=code_version, is_coframe=is_coframe, prec_flag=3, **kwargs)