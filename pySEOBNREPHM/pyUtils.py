#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri, 24 Mar 2023 04:11:29 +0000

@author: Shallyn
"""
import os
import sys
from optparse import OptionParser
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy.interpolate import InterpolatedUnivariateSpline, splev, splrep

# pwd = Path(__file__).absolute().parent
# pwd = Path(sys.path[0])

# -----Constants-----#
k_B_SI = 1.3806505e-23  # J/K
h_SI = 6.62606896e-34  # J s
ly_SI = 9.4605284e15  # 1 ly to m
AU_SI = 149597870700  # 1 AU to m
pc_SI = 3.08567758e16  # 1 Parsec to m
c_SI = 299792458  # Speed Of Light [m s^-1]
M_Sun_SI = 1.98892e30  # Mass of Sun [kg]
R_Sun_SI = 6.96342e8  # Radius of Sun [m]
alp_GP_SI = 192.85948 * np.pi / 180  # Direction Of Galactic Polar
det_GP_SI = 27.12825 * np.pi / 180  # Direction Of Galactic Polar
l_CP_SI = 122.932
G_SI = 6.672e-11  # Gravitational Constant [m^3 s^-2 kg^-1]
h_0_SI = 0.678  # Hubble Constant
MRSUN_SI = 1.47662504e3
MTSUN_SI = 4.92549095e-6
# ---------Comm--------#


def dim_t(M):
    return c_SI**3 / (M * M_Sun_SI * G_SI)


def get_fmin_from_Mfmin(Mfmin, Mtotal):
    return Mfmin * pow(c_SI, 3) / (G_SI * M_Sun_SI * Mtotal)


def get_Mfmin_from_fmin(fmin, Mtotal):
    return fmin * G_SI * Mtotal * M_Sun_SI / pow(c_SI, 3)


class bspline1d(object):
    def __init__(self, x, y, **kwargs):
        self.__spl = splrep(x, y, **kwargs)

    def __call__(self, x, **kwargs):
        return splev(x, self.__spl, **kwargs)


class bspline1d_complex(object):
    def __init__(self, x, y, **kwargs):
        yreal = y.real
        yimag = y.imag
        self.__x = x
        self.__y = y
        self.__spl_real = splrep(x, yreal, **kwargs)
        self.__spl_imag = splrep(x, yimag, **kwargs)

    def __call__(self, x, **kwargs):
        yr = splev(x, self.__spl_real, **kwargs)
        yi = splev(x, self.__spl_imag, **kwargs)
        return yr + 1.0j * yi

    @property
    def x(self):
        return self.__x

    @property
    def y(self):
        return self.__y


class gpstime(object):
    def __init__(self, gps0: str):
        self.__gps_str = gps0


class rTimeSeries(object):
    def __init__(self, time, value, srate=None, t0=None, gps0=None):
        if gps0 is None:
            gps0 = 0
        self.__gps0 = gps0
        if time is None:
            if srate is None:
                raise Exception('you should set time or srate')
            time = np.arange(0, len(value) / srate, 1.0 / srate)
        else:
            time = np.asarray(time)
        self.__t0 = time[0] if t0 is None else time[0] + t0
        self.__time = time - time[0]
        self.__mode = np.asarray(value)

    @property
    def gps0(self):
        return self.__gps0

    def reset_gps0(self, gps0):
        gps0_old = self.__gps0
        dgps0 = gps0_old - gps0
        self.__t0 = self.__t0 + dgps0
        self.__gps0 = gps0

    @property
    def t0(self):
        return self.__t0

    @property
    def fs(self):
        return 1.0 / (self.__time[1] - self.__time[0])

    @property
    def length(self):
        return len(self.__time)

    def set_t0(self, t0):
        self.__t0 = t0

    @property
    def tduration(self):
        return self.__time[-1]

    @property
    def time(self):
        return self.__time

    @property
    def time0(self):
        return self.__time + self.t0

    def resample(self, srate: float):
        dt = 1.0 / srate
        t_new = np.arange(0, self.time[-1], dt)
        val_new = self.interpolateFrom0(t_new)
        return rTimeSeries(t_new + self.t0, val_new, gps0=self.gps0)

    @property
    def dot(self):
        # return np.gradient(self.__mode) / np.gradient(self.__time)
        return self.interpolateFrom0(self.__time, der=1)

    @property
    def value(self):
        return self.__mode

    @property
    def duration(self):
        return self.time[-1] - self.time[0]

    @property
    def interpolate(self):
        return bspline1d(self.__time + self.t0, self.__mode)

    @property
    def interpolateFrom0(self):
        return bspline1d(self.__time, self.__mode)

    def copy(self):
        return rTimeSeries(self.__time.copy() + self.t0, self.__mode.copy())

    def __str__(self):
        return '{}'.format(self.__mode)

    def __len__(self):
        return len(self.__mode)

    def __repr__(self):
        return self.__str__()

    def __iter__(self):
        for x in self.__mode:
            yield x

    def __getitem__(self, key):
        if isinstance(key, int) or isinstance(key, np.integer):
            return self.__mode[key]
        return self.__getslice(key)

    def __setitem__(self, key, value):
        self.__mode[key] = value

    def __getslice(self, index):
        if index.start is not None and index.start < 0:
            raise ValueError(('Negative start index ({}) is not supported').format(index.start))
        return rTimeSeries(self.t0 + self.__time[index], self[index])

    def __del__(self):
        del self.__time
        del self.__mode

    def __neg__(self):
        return rTimeSeries(self.time + self.t0, -self.value, gps0=self.gps0)

    def dump(self, fname, **kwargs):
        data = np.stack([self.time, self.value], axis=1)
        np.savetxt(fname, data, **kwargs)
        return

    def pad(self, pad_width, mode, deltaT, **kwargs):
        self.__mode = np.pad(self.__mode, pad_width, mode, **kwargs)
        self.__time = np.arange(self.__time[0], self.__mode.size * deltaT, deltaT)
        self.__t0 = self.__t0 - pad_width[0] * deltaT


class cTimeSeries(object):
    def __init__(self, time, hreal, himag):
        time = np.asarray(time)
        self.__t0 = time[0]
        self.__time = time - time[0]
        self.__mode = np.asarray(hreal) + 1.0j * np.asarray(himag)

    def apply_phic(self, phic):
        self.__mode * np.exp(1.0j * phic)

    @property
    def t0(self):
        return self.__t0

    def set_t0(self, t0):
        self.__t0 = t0

    @property
    def length(self):
        return len(self.__time)

    @property
    def tduration(self):
        return self.__time[-1]

    @property
    def time(self):
        return self.__time

    @property
    def time0(self):
        return self.__time + self.t0

    def resample(self, srate: float):
        dt = 1.0 / srate
        t_new = np.arange(0, self.time[-1], dt)
        val_new = self.interpolateFrom0(t_new)
        return cTimeSeries(t_new + self.t0, val_new.real, val_new.imag)

    @property
    def dot(self):
        # return np.gradient(self.__mode) / np.gradient(self.__time)
        return self.interpolateFrom0(self.__time, der=1)

    @property
    def value(self):
        return self.__mode

    @property
    def real(self):
        return self.__mode.real

    @property
    def imag(self):
        return self.__mode.imag

    @property
    def real_rts(self):
        return rTimeSeries(self.__time + self.t0, self.__mode.real)

    @property
    def imag_rts(self):
        return rTimeSeries(self.__time + self.t0, self.__mode.imag)

    @property
    def amp(self):
        return np.abs(self.__mode)

    @property
    def amp_rts(self):
        return rTimeSeries(self.__time + self.t0, self.amp)

    @property
    def phase(self):
        return np.unwrap(np.angle(self.__mode))

    @property
    def phase_rts(self):
        return rTimeSeries(self.__time + self.t0, self.phase)

    @property
    def phaseFrom0(self):
        phase = self.phase
        return np.abs(phase - phase[0])

    @property
    def phaseFrom0_rts(self):
        return rTimeSeries(self.__time + self.t0, self.phaseFrom0)

    @property
    def frequency(self):
        amp = self.amp
        amp[amp == 0] = np.inf
        h = self.value
        hDot = self.dot
        freq = (h.real * hDot.imag - h.imag * hDot.real) / (amp * amp)
        return np.abs(freq)

    @property
    def frequency_rts(self):
        return rTimeSeries(self.__time + self.t0, self.frequency)

    @property
    def argpeak(self):
        return np.argmax(self.amp)

    @property
    def tpeak(self):
        return self.time[self.argpeak]

    @property
    def duration(self):
        return self.time[-1] - self.time[0]

    @property
    def interpolate(self):
        return bspline1d_complex(self.__time + self.t0, self.__mode)

    @property
    def interpolateFrom0(self):
        return bspline1d_complex(self.__time, self.__mode)

    @property
    def conjvalue(self):
        return self.__mode.conjugate()

    def conjugate(self):
        return cTimeSeries(self.time + self.t0, self.real, -self.imag)

    def copy(self):
        return cTimeSeries(self.__time.copy() + self.t0, self.__mode.copy().real, self.__mode.copy().imag)

    def __str__(self):
        return '{}'.format(self.__mode)

    def __len__(self):
        return len(self.__mode)

    def __repr__(self):
        return self.__str__()

    def __iter__(self):
        for x in self.__mode:
            yield x

    def __getitem__(self, key):
        if isinstance(key, int) or isinstance(key, np.integer):
            return self.__mode[key]
        return self.__getslice(key)

    def __setitem__(self, key, value):
        self.__mode[key] = value

    def __getslice(self, index):
        if index.start is not None and index.start < 0:
            raise ValueError(('Negative start index ({}) is not supported').format(index.start))
        return cTimeSeries(self.__time[index] + self.t0, self.real[index], self.imag[index])

    def __del__(self):
        del self.__time
        del self.__mode

    def dump(self, fname, **kwargs):
        data = np.stack([self.time, self.real, self.imag], axis=1)
        np.savetxt(fname, data, **kwargs)
        return

    def pad(self, pad_width, mode, deltaT, **kwargs):
        self.__mode = np.pad(self.__mode, pad_width, mode, **kwargs)
        self.__time = np.arange(self.__time[0], self.__mode.size * deltaT, deltaT)
        self.__t0 = self.__t0 - pad_width[0] * deltaT


def cTimeSeries_alignment(modeA: cTimeSeries, modeB: cTimeSeries, deltaT=None):
    tA = modeA.time
    tB = modeB.time
    dtA = tA[1] - tA[0]
    dtB = tB[1] - tB[0]
    if dtA != dtB:
        if deltaT is not None:
            dt_final = deltaT
        elif dtA < dtB:
            dt_final = dtB
        else:
            dt_final = dtA
        tA_new = np.arange(tA[0], tA[-1], dt_final)
        tB_new = np.arange(tB[0], tB[-1], dt_final)
        valA = modeA.interpolateFrom0(tA_new)
        valB = modeB.interpolateFrom0(tB_new)
    else:
        tA_new = tA
        tB_new = tB
        valA = modeA.value
        valB = modeB.value
        dt_final = dtA
    wfA = cTimeSeries(tA_new, valA.real, valA.imag)
    wfB = cTimeSeries(tB_new, valB.real, valB.imag)
    if len(valA) == len(valB):
        return wfA, wfB
    ipeak_A = wfA.argpeak
    ipeak_B = wfB.argpeak
    if ipeak_A > ipeak_B:
        idx_A = ipeak_A - ipeak_B
        idx_B = 0
    else:
        idx_A = 0
        idx_B = ipeak_B - ipeak_A
    # tmove = (ipeak_A - ipeak_B) * dt_final
    wfA = wfA[idx_A:]
    wfB = wfB[idx_B:]
    lenA = len(wfA)
    lenB = len(wfB)
    ipeak_A = ipeak_A - idx_A
    ipeak_B = ipeak_B - idx_B
    tail_A = lenA - ipeak_A
    tail_B = lenB - ipeak_B
    if tail_A > tail_B:
        lpad = tail_A - tail_B
        wfB.pad((0, lpad), 'constant', dt_final)
    else:
        lpad = tail_B - tail_A
        wfA.pad((0, lpad), 'constant', dt_final)
    return wfA, wfB


def rTimeSeries_Aligned(h1: rTimeSeries, h2: rTimeSeries, deltaT=None, window_alpha1=None, window_alpha2=None, method: int = 0):
    tA = h1.time
    tB = h2.time
    gps0 = h1.gps0
    if h1.gps0 != h2.gps0:
        h2.reset_gps0(gps0)
    tA0 = h1.t0
    tB0 = h2.t0
    tAe = tA[-1] + tA0
    tBe = tB[-1] + tB0
    dtA = np.min(np.gradient(tA))
    dtB = np.min(np.gradient(tB))
    if deltaT is not None:
        dt_final = deltaT
    elif dtA < dtB:
        dt_final = dtB
    else:
        dt_final = dtA
    if method == 0:
        t0 = min(tA0, tB0)
        te = max(tAe, tBe)
    elif method == 1:
        t0 = tA0
        te = tAe
    else:
        t0 = tB0
        te = tBe
    t_final = np.arange(t0, te, dt_final)
    valA_new = np.zeros(len(t_final))
    valB_new = np.zeros(len(t_final))
    isliceA = np.where((t_final > tA0) & (t_final < tAe))
    isliceB = np.where((t_final > tB0) & (t_final < tBe))
    if window_alpha1 is not None:
        dwA = signal.tukey(len(h1), alpha=window_alpha1)
    else:
        dwA = 1.0

    if window_alpha2 is not None:
        dwB = signal.tukey(len(h2), alpha=window_alpha2)
    else:
        dwB = 1.0
    h1_windowed = rTimeSeries(h1.time + h1.t0, h1.value * dwA)
    h2_windowed = rTimeSeries(h2.time + h2.t0, h2.value * dwB)
    valA_new[isliceA] = h1_windowed.interpolate(t_final[isliceA])
    valB_new[isliceB] = h2_windowed.interpolate(t_final[isliceB])
    wfA = rTimeSeries(t_final, valA_new, gps0=gps0)
    wfB = rTimeSeries(t_final, valB_new, gps0=gps0)
    wfA.set_t0(t0)
    wfB.set_t0(t0)
    return wfA, wfB


def rTimeSeries_overlap(h1: rTimeSeries, h2: rTimeSeries, psd, deltaT=None, fmin=None, fmax=None, **kwargs):
    h1, h2 = rTimeSeries_Aligned(h1, h2, deltaT, **kwargs)
    NFFT = len(h1)
    dt = h1.time[1] - h1.time[0]
    df = 1.0 / NFFT / dt
    h1tilde = np.fft.rfft(h1) * dt
    h2tilde = np.fft.rfft(h2) * dt
    freqs = np.fft.rfftfreq(NFFT, dt)
    if fmin is None:
        fmin = 0.0
    if fmax is None:
        fmax = np.inf
    islice = np.where((freqs > fmin) & (freqs < fmax))
    power_vec = np.abs(psd(freqs[islice]))
    O11 = np.sum(h1tilde[islice] * h1tilde[islice].conjugate() / power_vec).real * df
    O22 = np.sum(h2tilde[islice] * h2tilde[islice].conjugate() / power_vec).real * df
    O12 = h1tilde[islice] * h2tilde[islice].conjugate() / power_vec
    # FF = np.fft.irfft(O12) * len(O12) * df
    return np.sum(O12).real * df / np.sqrt(O11 * O22)


def rTimeSeries_innerproduct(h: rTimeSeries, psd, fmin=None, fmax=None):
    NFFT = len(h)
    dt = h.time[1] - h.time[0]
    df = 1.0 / NFFT / dt
    htilde = np.fft.rfft(h) * dt
    freqs = np.fft.rfftfreq(NFFT, dt)
    if fmin is None:
        fmin = 0.0
    if fmax is None:
        fmax = np.inf
    islice = np.where((freqs > fmin) & (freqs < fmax))
    power_vec = np.abs(psd(freqs[islice]))
    O12 = htilde[islice] * htilde[islice].conjugate() / power_vec
    return 4.0 * np.sum(O12).real * df


def rTimeSeries_filtering(s: rTimeSeries, h: rTimeSeries, psd, deltaT=None, fmin=None, fmax=None, **kwargs):
    h1, h2 = rTimeSeries_Aligned(s, h, deltaT, **kwargs)
    NFFT = len(h1)
    dt = h1.time[1] - h1.time[0]
    df = 1.0 / NFFT / dt
    h1tilde = np.fft.rfft(h1) * dt
    h2tilde = np.fft.rfft(h2) * dt
    freqs = np.fft.rfftfreq(NFFT, dt)
    if fmin is None:
        fmin = 0.0
    if fmax is None:
        fmax = np.inf
    islice = np.where((freqs > fmin) & (freqs < fmax))
    power_vec = np.abs(psd(freqs[islice]))
    O22 = 4.0 * np.sum(h2tilde[islice] * h2tilde[islice].conjugate() / power_vec).real * df
    O12 = h1tilde[islice] * h2tilde[islice].conjugate() / power_vec
    # FF = np.fft.irfft(O12) * len(O12) * df
    return 4.0 * np.sum(O12).real * df / np.sqrt(O22)


def rTimeSeries_add(h1: rTimeSeries, h2: rTimeSeries, srate=None, window_alpha1=None, window_alpha2=None, method: int = 0):
    if srate is None:
        deltaT = None
    else:
        deltaT = 1.0 / srate
    h1_aligned, h2_aligned = rTimeSeries_Aligned(h1, h2, deltaT=deltaT, window_alpha1=window_alpha1, window_alpha2=window_alpha2, method=method)
    return rTimeSeries(h1_aligned.time0, h1_aligned.value + h2_aligned.value, gps0=h1_aligned.gps0)
