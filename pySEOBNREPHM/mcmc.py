#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun, 26 Mar 2023 14:11:28 +0000

@author: Shallyn
"""
import os
import sys
import time as pytime
from optparse import OptionParser
from pathlib import Path

import corner
import emcee
import h5py
import matplotlib.pyplot as plt
import numpy as np
import tqdm

from . import SEOBNRWaveformCaller
from .detector import GWDetector
from .pyUtils import rTimeSeries, rTimeSeries_add, rTimeSeries_Aligned, rTimeSeries_filtering, rTimeSeries_innerproduct

try:
    from pycbc.waveform import get_td_waveform

    is_pycbc_loaded = True
except ImportWarning:
    is_pycbc_loaded = False
# pwd = Path(__file__).absolute().parent
# pwd = Path(sys.path[0])
# pwd = Path(os.getcwd())


class MCMCRunner(object):
    def __init__(self, ndim: 'int > 0', nwalkers: 'int > 0', func_lnprob: 'function', threads: 'int > 0', conv_thresh: 'float > 0'):
        self.__ndim = ndim
        self.__nwalkers = nwalkers
        self.__lnprob = func_lnprob
        self.__threads = threads
        self.__conv_thresh = conv_thresh
        self.__sampler = emcee.EnsembleSampler(nwalkers, ndim, func_lnprob, threads=threads)

    def run(self, p0, max_step: 'int > 100', fsave: str = 'chain.txt', fconv: str = 'converge_status.txt', sigma=1e-4, is_delete=True, min_step: 'int>100' = 200, debug=False, **kwargs):
        if len(p0) != self.__ndim:
            raise Exception('incompatible dimension')
        fchain = Path(fsave)
        if is_delete:
            f = open(fchain, 'w')
            f.close()
            f = open(fconv, 'w')
            f.close()
            p0 = np.asarray(p0)
            sigp0 = sigma * np.abs(p0)
            sigp0[sigp0 < 1e-16] = sigma * np.random.randn(len(sigp0[sigp0 < 1e-16]))
            pini = [(p0 + sigp0 * np.random.randn(self.__ndim)).tolist() for i in range(self.__nwalkers)]
            pos, prob, state = self.__sampler.run_mcmc(pini, 5, skip_initial_state_check=True)
            self.__sampler.reset()
            i_count = 0
        else:
            # try load fchain
            try:
                prev_chains = np.loadtxt(fchain)
                pini, length = get_ele_from_save_chain(prev_chains, self.__ndim, self.__nwalkers)
                print(f'length of previous chain = {length}')
                # print(pini)
                # length, ndim = prev_chains.shape
                # if ndim != self.__ndim:
                #     raise Exception(f'incorrect shape of the chain (should be {self.__ndim}, however {ndim})\n')
                # if length % self.__nwalkers != 0:
                #     raise Exception(f'incompatible walker number (should be {self.__nwalkers})\n')
                # pini = prev_chains[-self.__nwalkers:, :].tolist()
            except FileNotFoundError:
                if debug:
                    sys.stderr.write('cannot read previous chains, return')
                    prev_chains = np.loadtxt(fchain)
                    pini = get_ele_from_save_chain(prev_chains, self.__ndim, self.__nwalkers)
                    return
                sys.stderr.write('cannot read previous chains, so we replace it')
                f = open(fchain, 'w')
                f.close()
                f = open(fconv, 'w')
                f.close()
                p0 = np.asarray(p0)
                length = 0
                sigp0 = sigma * np.abs(p0)
                sigp0[sigp0 < 1e-16] = sigma * np.random.randn(len(sigp0[sigp0 < 1e-16]))
                pini = [(p0 + sigp0 * np.random.randn(self.__ndim)).tolist() for i in range(self.__nwalkers)]
            pos, prob, state = self.__sampler.run_mcmc(pini, 5, skip_initial_state_check=True)
            self.__sampler.reset()
            i_count = length / self.__nwalkers
        start_time = pytime.time()
        R_c = np.linspace(0, 1, self.__ndim)
        converg_state = False
        for position, prob, state in self.__sampler.sample(pos, iterations=max_step, store=False, skip_initial_state_check=True, **kwargs):
            # if np.sum(np.isinf(prob)) > 0:
            #     raise Exception(f'Nan appears')
            f = open(fchain, 'a')
            for k in range(position.shape[0]):
                f.write('1 %.8f ' % (prob[k]))
                for j in range(self.__ndim):
                    f.write('%.8f ' % position[k, j])
                f.write('\n')
            f.close()
            i_count += 1
            if i_count > min_step and i_count % 100 == 0:
                converg_state = self.__CheckConverge(fchain, i_count, fconv)
            if converg_state:
                i_count = max_step
                break
        end_time = pytime.time()
        sys.stderr.write(f'MCMC finished, time cost {end_time - start_time}\n')

    def __CheckConverge(self, fchain, i_count, fconv):
        chains = np.loadtxt(fchain)[:, 2:]
        length, ndim = chains.shape
        if ndim != self.__ndim:
            raise Exception(f'incorrect shape of the chain (should be {self._ndim}, however {ndim})')
        i_c = 0
        thresh_list = []
        n_count = int(length / self.__nwalkers)
        for i in range(ndim):
            para = chains[:, i]
            para_reshape = para.reshape((n_count, self.__nwalkers))
            para_reshape = para_reshape[int(n_count / 2) : n_count, :]
            lenN = para_reshape.shape[0]
            walker_mean = np.mean(para_reshape, axis=0, keepdims=True)
            Bn = np.var(walker_mean, ddof=1)
            walker_var = np.var(para_reshape, axis=0, keepdims=True, ddof=1)
            Wn = np.mean(walker_var)
            R_c = (Wn * (1.0 - 1.0 / lenN) + Bn * (1.0 + 1.0 / self.__nwalkers)) / Wn
            # print(f'{R_c}, {Bn}, {Wn}, {walker_mean}')
            thresh = abs(R_c - 1)
            thresh_list.append(thresh)
            if abs(R_c - 1) < self.__conv_thresh:
                i_c = i_c + 1

        # print('\n')
        f = open(fconv, 'a')
        for k in range(len(thresh_list)):
            f.write('%.8f\t' % thresh_list[k])
        f.write('%.8f\n' % np.max(thresh_list))
        f.close()
        if i_c == self.__ndim:
            return True
        return False


def get_ele_from_save_chain(chain: np.ndarray, ndim, nwalker):
    correct_ndim = ndim + 2
    length, chain_ndim = chain.shape
    if chain_ndim != correct_ndim:
        raise Exception(f'incorrect shape of the chain (should be {correct_ndim}, however {chain_ndim})\n')
    if length % nwalker != 0:
        raise Exception(f'incompatible walker number (should be {nwalker})\n')
    pini = chain[-nwalker:, 2:].tolist()
    return pini, length


def write_paramnames(fparamnames, ndim):
    # q, Mchirp, th1, ph1, chi1, th2, ph2, chi2, e0, dL, theta, phi, iota, psi, t_c, Phi_c
    f = open(fparamnames, 'w')
    f.write('q\tq\n')
    f.write('Mchirp\t\\mathcal{M}\n')
    if ndim == 15 or ndim == 16:
        f.write('th1\t\\theta_{1}\n')
        f.write('ph1\t\\phi_{1}\n')
    f.write('chi1\t\\chi_{1}\n')
    if ndim == 15 or ndim == 16:
        f.write('th2\t\\theta_{2}\n')
        f.write('ph2\t\\phi_{2}\n')
    f.write('chi2\t\\chi_{2}\n')
    if ndim == 12 or ndim == 16:
        f.write('ecc\te_0\n')
    f.write('dL\td_{L}\n')
    f.write('theta\t\\theta\n')
    f.write('phi\t\\phi\n')
    f.write('iota\t\\iota\n')
    f.write('psi\t\\psi\n')
    f.write('tc\tt_c\n')
    f.write('phic\t\\phi_c\n')
    f.close()
    return


def write_ranges(franges, ndim, t_extra, emax):
    # q, Mchirp, th1, ph1, chi1, th2, ph2, chi2, e0, dL, theta, phi, iota, psi, t_c, Phi_c
    f = open(franges, 'w')
    f.write('q\t\t%.8e\t%.8e\n' % (1, 20))
    f.write('Mchirp\t\t%.8e\t%.8e\n' % (1, 200))
    if ndim == 15 or ndim == 16:
        f.write('th1\t\t%.8e\t%.8e\n' % (0, np.pi))
        f.write('ph1\t\t%.8e\t%.8e\n' % (0, 2 * np.pi))
    f.write('chi1\t\t%.8e\t%.8e\n' % (0, 0.98))
    if ndim == 15 or ndim == 16:
        f.write('th2\t\t%.8e\t%.8e\n' % (0, np.pi))
        f.write('ph2\t\t%.8e\t%.8e\n' % (0, 2 * np.pi))
    f.write('chi2\t\t%.8e\t%.8e\n' % (0, 0.98))
    if ndim == 12 or ndim == 16:
        f.write('ecc\t\t%.8e\t%.8e\n' % (0.0, emax))
    f.write('dL\t\t%.8e\tN\n' % 0)
    f.write('theta\t\t%.8e\t%.8e\n' % (0, np.pi))
    f.write('phi\t\t%.8e\t%.8e\n' % (0, 2.0 * np.pi))
    f.write('iota\t\t%.8e\t%.8e\n' % (0, np.pi))
    f.write('psi\t\t%.8e\t%.8e\n' % (0, np.pi))
    f.write('tc\t\t%.8e\t%.8e\n' % (-2.0 * t_extra, 2.0 * t_extra))
    f.write('phic\t\t%.8e\t%.8e\n' % (0, 2.0 * np.pi))
    f.close()
    return


def calculate_m1m2_from_qMchirp(q, Mchirp):
    eta = q / (1.0 + q) / (1.0 + q)
    Mtot = Mchirp * np.power(eta, -3.0 / 5.0)
    m1 = Mtot * q / (1.0 + q)
    m2 = Mtot / (1.0 + q)
    return m1, m2


def calculate_qMchirp_from_m1m2(m1, m2):
    q = m1 / m2
    eta = q / (1.0 + q) / (1.0 + q)
    Mtot = m1 + m2
    Mchirp = Mtot * np.power(eta, 3.0 / 5.0)
    return q, Mchirp


def calculate_chiVec_from_thphchi(th, ph, chi):
    chix = chi * np.sin(th) * np.cos(ph)
    chiy = chi * np.sin(th) * np.sin(ph)
    chiz = chi * np.cos(th)
    return chix, chiy, chiz


def calculate_thphchi_from_chiVec(chix, chiy, chiz):
    chi = np.sqrt(chix * chix + chiy * chiy + chiz * chiz)
    rho = np.sqrt(chix * chix + chiy * chiy)
    ph = np.arctan2(chiy, chix)
    if ph < 0:
        ph = 2.0 * np.pi + ph
    th = np.arctan2(rho, chiz)
    return th, ph, chi


def simulate_waveform_from_detector_list(params, f_min=40.0, gps0=1356566418, detectors=None, is_only22=False, log_level=1, retall=False, **kwargs) -> rTimeSeries:
    if 'debug' in kwargs:
        if kwargs['debug']:
            sys.stderr.write(f'DEBUG: params = {params}, f_min = {f_min}, gps0 = {gps0}\n')
    q, Mchirp, th1, ph1, chi1, th2, ph2, chi2, e0, dL, theta, phi, iota, psi, t_c, Phi_c = params
    s1x, s1y, s1z = calculate_chiVec_from_thphchi(th1, ph1, chi1)
    s2x, s2y, s2z = calculate_chiVec_from_thphchi(th2, ph2, chi2)
    m1, m2 = calculate_m1m2_from_qMchirp(q, Mchirp)
    ge = SEOBNRWaveformCaller()
    ge.set_params(m1=m1, m2=m2, s1x=s1x, s1y=s1y, s1z=s1z, s2x=s2x, s2y=s2y, s2z=s2z, ecc=e0, distance=dL, inc_rad=iota, beta_rad=Phi_c, f_min=f_min, log_level=log_level, is_only22=is_only22, **kwargs)
    waveform, _ = ge.run()
    if waveform is None:
        sys.stderr.write(f'WARNING: generation failed\n')
        return None
    if detectors is None:
        # default use L1
        detectors = [GWDetector('L1')]
    # apCalculator = detector.antenna_pattern(psi, phi, theta)
    # strain = apCalculator.detector_strain_t(waveform.hpc, t_c+gps0)
    hp = waveform.hpc.real
    hc = -waveform.hpc.imag
    tpeak = waveform.hpc.time[waveform.hpc.argpeak]
    ret = []
    for detector in detectors:
        dt = detector.time_delay(phi, np.pi - theta, t_c + gps0)
        Fplus, Fcross = detector.antenna_pattern_gps(psi, phi, np.pi - theta, t_c + gps0 + dt)
        h = Fplus * hp + Fcross * hc
        t = waveform.hpc.time - tpeak
        ret.append(rTimeSeries(t, h, t0=t_c + dt, gps0=gps0))
    if retall:
        return ret, waveform
    return ret


def simulate_waveform_from_detector_list_pycbc(params, f_min=40.0, gps0=1356566418, detectors=None, approx=None, srate=16384, **kwargs):
    if approx is None:
        approx = 'SEOBNRv4PHM'
    q, Mchirp, th1, ph1, chi1, th2, ph2, chi2, e0, dL, theta, phi, iota, psi, t_c, Phi_c = params
    s1x, s1y, s1z = calculate_chiVec_from_thphchi(th1, ph1, chi1)
    s2x, s2y, s2z = calculate_chiVec_from_thphchi(th2, ph2, chi2)
    m1, m2 = calculate_m1m2_from_qMchirp(q, Mchirp)
    hp, hc = get_td_waveform(
        approximant=approx, mass1=m1, mass2=m2, spin1x=s1x, spin1y=s1y, spin1z=s1z, spin2x=s2x, spin2y=s2y, spin2z=s2z, distance=dL, inclination=iota * 180 / np.pi, coa_phase=Phi_c * 180 / np.pi, delta_t=1.0 / srate, f_lower=f_min
    )
    time0 = hp.sample_times
    hp = hp.data
    hc = hc.data
    plt.plot(time0, hp)
    plt.plot(time0, hc)
    plt.show()
    if detectors is None:
        # default use L1
        detectors = [GWDetector('L1')]
    ret = []
    for detector in detectors:
        dt = detector.time_delay(phi, np.pi - theta, t_c + gps0)
        Fplus, Fcross = detector.antenna_pattern_gps(psi, phi, np.pi - theta, t_c + gps0 + dt)
        h = Fplus * hp + Fcross * hc
        t = time0
        ret.append(rTimeSeries(t, h, t0=t_c + dt, gps0=gps0))
    return ret


def parse_save_prefix(prefix):
    lst = prefix.replace(' ', '').split('/')
    lst_re = []
    for val in lst:
        if val != '':
            lst_re.append(val)
    if len(lst_re) == 1:
        datadir = lst[0]
        dprefix = 'default'
        return datadir, dprefix
    datadir = '/'.join(lst_re[:-1])
    dprefix = lst_re[-1]
    return datadir, dprefix


def dump_waveform_to_hdf5(waveform, fh5_fake_signal):
    length = len(waveform)
    fh5 = h5py.File(str(fh5_fake_signal), "a")
    if 'waveform' not in fh5:
        h5g = fh5.create_group('waveform')
        d = h5g.create_dataset('timeM', (length,))
        d[...] = waveform.timeM
        d = h5g.create_dataset('h22_real', (length,))
        d[...] = waveform.h22.real
        d = h5g.create_dataset('h22_imag', (length,))
        d[...] = waveform.h22.imag
        d = h5g.create_dataset('h21_real', (length,))
        d[...] = waveform.h21.real
        d = h5g.create_dataset('h21_imag', (length,))
        d[...] = waveform.h21.imag
        d = h5g.create_dataset('h33_real', (length,))
        d[...] = waveform.h33.real
        d = h5g.create_dataset('h33_imag', (length,))
        d[...] = waveform.h33.imag
        d = h5g.create_dataset('h44_real', (length,))
        d[...] = waveform.h44.real
        d = h5g.create_dataset('h44_imag', (length,))
        d[...] = waveform.h44.imag
        d = h5g.create_dataset('h55_real', (length,))
        d[...] = waveform.h55.real
        d = h5g.create_dataset('h55_imag', (length,))
        d[...] = waveform.h55.imag
    else:
        sys.stderr.write('WARNING: \'waveform\' already exists')
    fh5.close()


def dump_rts_list_to_hdf5(group_name, h_rts_list, detector_name_list, fh5_fake_signal):
    fh5 = h5py.File(str(fh5_fake_signal), "a")
    if group_name not in fh5:
        h5g = fh5.create_group(group_name)
        for i in range(len(h_rts_list)):
            detname = detector_name_list[i]
            h5g_det = h5g.create_group(detname)
            h_rts = h_rts_list[i]
            length = len(h_rts)
            h5g_det['t0'] = h_rts.t0
            h5g_det['gps0'] = h_rts.gps0
            h5g_det['fs'] = h_rts.fs
            d = h5g_det.create_dataset('time', (length,))
            d[...] = h_rts.time
            d = h5g_det.create_dataset('strain', (length,))
            d[...] = h_rts.value
    else:
        sys.stderr.write(f'WARNING: \'{group_name}\' already exists')
    fh5.close()


def load_rts_list_from_hdf5(fh5_fake_signal, group_name, detector_name_list):
    fh5 = h5py.File(str(fh5_fake_signal), "r")
    if group_name not in fh5:
        raise Exception(f'cannot load signal from {group_name}')
    h5g = fh5[group_name]
    ret = []
    for i in range(len(detector_name_list)):
        detname = detector_name_list[i]
        if detname not in h5g:
            raise Exception(f'cannot load signal from {group_name}/{detname}')
        h5g_det = h5g[detname]
        t0 = h5g_det['t0'][()]
        gps0 = h5g_det['gps0'][()]
        time = h5g_det['time'][()]
        strain = h5g_det['strain'][()]
        ret.append(rTimeSeries(time, strain, t0=t0, gps0=gps0))
    fh5.close()
    return ret


def simulate_mcmc(argv):
    parser = OptionParser(description='FIXME')
    parser.add_option('--log-level', type='int', default=1, help='debug level')
    parser.add_option('--debug', action='store_true', help='debug test')
    parser.add_option('--code-debug', action='store_true', help='code debug test')
    parser.add_option('--is-only22', action='store_true', help='only 22 mode will be included')
    parser.add_option('--jobid', type='int', default=1, help='job id')

    # source parameters
    parser.add_option('--m1', type='float', default=20, help='input m1')
    parser.add_option('--m2', type='float', default=10, help='input m2')
    parser.add_option('--s1x', type='float', default=0, help='input s1x')
    parser.add_option('--s1y', type='float', default=0, help='input s1y')
    parser.add_option('--s1z', type='float', default=0, help='input s1z')
    parser.add_option('--s2x', type='float', default=0, help='input s2x')
    parser.add_option('--s2y', type='float', default=0, help='input s2y')
    parser.add_option('--s2z', type='float', default=0, help='input s2z')
    parser.add_option('--code-version', type='int', default=1, help='version of code')
    parser.add_option('--eccentricity', type='float', default=0, help='input ecc')
    parser.add_option('--distance', type='float', default=100.0, help='input distance')
    parser.add_option('--inclination', type='float', default=0, help='input inclination in degree')
    parser.add_option('--phi-ref', type='float', default=0, help='input reference phic in degree')
    parser.add_option('--f-min', type='float', default=40, help='input initial frequency [Hz]')
    parser.add_option('--gps0', type='float', default=1356566418, help='input gps0')
    parser.add_option('--polarization', type='float', default=0, help='input polarization angel in degree')
    parser.add_option('--ra', type='float', default=0, help='input right ascension in degree')
    parser.add_option('--dec', type='float', default=0, help='input declination in degree')
    parser.add_option('--sample-rate', type='float', default=16384, help='sample rate')
    parser.add_option('--t-extra', type='float', default=5.0, help='extra time for fake data')
    parser.add_option('--emax', type='float', default=0.6, help='maximum eccentricity')

    # detector settings
    parser.add_option('--detector', type='str', action='append', help='detector prefix')
    parser.add_option('--prefix', type='str', default='results/test', help='prefix')
    parser.add_option('--fsignal', type='str', help='signal h5 file')
    parser.add_option('--fnoise', type='str', help='noise h5 file')
    parser.add_option('--fstrain', type='str', help='strain h5 file')
    parser.add_option('--model', type='str', default='preccirc', help='model we use, [PrecCirc or SpinAligned]')
    parser.add_option('--adjust', action='store_true', help='adjust eccentricity')

    # mcmc settings
    parser.add_option('--nwalkers', type='int', default=32, help='number of mcmc walkers')
    parser.add_option('--max-steps', type='int', default=10000, help='number of mcmc walkers')
    parser.add_option('--thresh', type='float', default=0.05, help='MCMC thresh')
    parser.add_option('--sigma0', type='float', default=0.01, help='random vairation for initial step')
    parser.add_option('--nthreads', type='int', default=8, help='number threads for mcmc')
    parser.add_option('--progress', action='store_true', help='number threads for mcmc')
    parser.add_option('--delete', action='store_true', help='whether delete previous data')

    args, _ = parser.parse_args(argv)
    pwd = Path(sys.path[0])
    detector_name_list = []
    if args.detector is None or len(args.detector) == 0:
        raise Exception(f'detectors are not specified')
    jobid = args.jobid
    if jobid != 1:
        pytime.sleep(100)
    detector_name_list = args.detector
    detector_list = [GWDetector(detector_name) for detector_name in detector_name_list]
    sample_rate = args.sample_rate
    gps0 = args.gps0
    t_extra = np.abs(args.t_extra)

    m1 = args.m1
    m2 = args.m2
    q, Mchirp = calculate_qMchirp_from_m1m2(m1, m2)
    spinChi1 = [args.s1x, args.s1y, args.s1z]  # initial dimensionless spin chi, which satisfies |chi| < 1
    spinChi2 = [args.s2x, args.s2y, args.s2z]
    th1, ph1, chi1 = calculate_thphchi_from_chiVec(spinChi1[0], spinChi1[1], spinChi1[2])
    th2, ph2, chi2 = calculate_thphchi_from_chiVec(spinChi2[0], spinChi2[1], spinChi2[2])
    ecc = args.eccentricity  # eccentricity, <0.6 is recommended
    f_min = args.f_min  # Hz
    # extrinsic parameters
    dL = args.distance  # Mpc
    theta = np.pi - args.dec * np.pi / 180.0  # pi/2 - declination at the time gps0+t_c, in rad
    phi = args.ra * np.pi / 180.0  # right ascension at the time gps0+t_c, in rad
    psi = args.polarization * np.pi / 180.0  # polarization angle, in rad
    iota = args.inclination * np.pi / 180.0  # inclination angle, in rad
    Phi_c = args.phi_ref * np.pi / 180.0  # reference orbital phase, in rad

    datadirname, prefix = parse_save_prefix(args.prefix)
    datadir = pwd / datadirname
    if not datadir.exists():
        datadir.mkdir(parents=True)

    # 1. generate fake strain and dump to datadir/prefix_fake_signal.h5
    fh5_fake_signal = datadir / f'{prefix}_fake_signal.h5'
    if args.fsignal is not None:
        fh5_fake_signal = Path(args.fsignal)
    is_regen = False
    is_delete = True
    if jobid == 1 and args.delete:
        is_regen = True
        fh5 = h5py.File(fh5_fake_signal, 'w')
        fh5.close()
    elif jobid == 1 and not args.delete and fh5_fake_signal.exists():
        try:
            is_delete = False
            sys.stderr.write(f'try loading signals from {fh5_fake_signal}\n')
            signal_rts_list = load_rts_list_from_hdf5(fh5_fake_signal, 'signal', detector_name_list)
            sys.stderr.write('loading signal success\n')
        except FileNotFoundError:
            sys.stderr.write(f'failed to load signals from {fh5_fake_signal}\n')
            is_regen = True
            is_delete = True
    elif jobid == 1 and not fh5_fake_signal.exists():
        is_regen = True
        fh5 = h5py.File(fh5_fake_signal, 'w')
        fh5.close()
    params0 = (q, Mchirp, th1, ph1, chi1, th2, ph2, chi2, ecc, dL, theta, phi, iota, psi, 0.0, Phi_c)
    sys.stderr.write(f'params0 = {params0}\n')
    if is_regen:
        sys.stderr.write(f'Generate fake strain, dump to {fh5_fake_signal}\n')
        is_gen_strain = True
        if args.fstrain is not None:
            if Path(args.fstrain).exists():
                is_gen_strain = False
        if is_gen_strain:
            f_minE = f_min
            if args.adjust:
                f_minE = f_min * (1.0 - ecc)
            if args.model.lower() == 'spinaligned' or args.model.lower() == 'sa':
                hstrain_list, waveform = simulate_waveform_from_detector_list(params0, gps0=gps0, f_min=f_minE, detectors=detector_list, is_only22=False, retall=True, egw_flag=1, log_level=args.log_level, code_version=args.code_version)
            else:
                hstrain_list, waveform = simulate_waveform_from_detector_list(params0, gps0=gps0, f_min=f_minE, detectors=detector_list, is_only22=False, retall=True, prec_flag=3, log_level=args.log_level, code_version=args.code_version)
            hstrain_rts_list = [hstrain.resample(sample_rate) for hstrain in hstrain_list]
            dump_waveform_to_hdf5(waveform, fh5_fake_signal)
        else:
            sys.stderr.write(f'load strain from {args.fstrain}\n')
            hstrain_rts_list = load_rts_list_from_hdf5(Path(args.fstrain), 'strain', detector_name_list)
        dump_rts_list_to_hdf5('strain', hstrain_rts_list, detector_name_list, fh5_fake_signal)
        # 2. simulate noise and dump noise and signal to datadir/prefix_fake_signal.h5
        sys.stderr.write(f'Simulate noise and signal, dump to {fh5_fake_signal}\n')
        is_gen_noise = True
        if args.fnoise is not None:
            if Path(args.fnoise).exists():
                is_gen_noise = False
        if is_gen_noise:
            sys.stderr.write('simulate noise\n')
            noise_list = [detector_list[i].psd.generate_noise(hstrain_rts_list[i].tduration + 2.0 * t_extra, sample_rate) for i in range(len(detector_name_list))]
            noise_rts_list = [rTimeSeries(None, noise_list[i], srate=sample_rate, t0=hstrain_rts_list[i].t0 - t_extra, gps0=gps0) for i in range(len(detector_name_list))]
            signal_rts_list = [rTimeSeries_add(hstrain_rts_list[i], noise_rts_list[i], srate=sample_rate) for i in range(len(detector_name_list))]
        else:
            sys.stderr.write(f'load noise from {args.fnoise}\n')
            noise_rts_list = load_rts_list_from_hdf5(Path(args.fnoise), 'noise', detector_name_list)
            signal_rts_list = [rTimeSeries_add(hstrain_rts_list[i], noise_rts_list[i], srate=sample_rate, method=2) for i in range(len(detector_name_list))]
        dump_rts_list_to_hdf5('noise', noise_rts_list, detector_name_list, fh5_fake_signal)
        dump_rts_list_to_hdf5('signal', signal_rts_list, detector_name_list, fh5_fake_signal)
        for i in range(len(detector_name_list)):
            hAligned_rts, sAligned_rts = rTimeSeries_Aligned(hstrain_rts_list[i], signal_rts_list[i], method=2)
            fig = plt.figure(figsize=(10, 6))
            ax1 = fig.add_subplot(2, 1, 1)
            ax2 = fig.add_subplot(2, 1, 2)
            ax1.plot(hAligned_rts.time0, hAligned_rts.value, color='black')
            ax2.plot(sAligned_rts.time0, sAligned_rts.value, color='black')
            ax2.set_xlabel(f't[s]')
            fig.savefig(datadir / f'{prefix}_fig_strain_signal_{detector_name_list[i]}.jpg', dpi=200)
            fig.savefig(datadir / f'{prefix}_fig_strain_signal_{detector_name_list[i]}.eps', dpi=200)
            plt.close(fig)
    elif jobid != 1:
        i_wait = 0
        is_found = False
        while i_wait < 100:
            i_wait += 1
            if fh5_fake_signal.exists():
                is_found = True
                pytime.sleep(10 * jobid)
                break
            pytime.sleep(10)
        if is_found:
            sys.stderr.write(f'try loading signals from {fh5_fake_signal}\n')
            signal_rts_list = load_rts_list_from_hdf5(fh5_fake_signal, 'signal', detector_name_list)
            is_delete = args.delete
            sys.stderr.write('loading signal success\n')
        else:
            raise Exception(f'cannot find signal file {fh5_fake_signal}')

    # 3. MCMC lnprob setting
    # Define the matched filtering product
    def residuals_inner_product(residuals: rTimeSeries, detector: GWDetector, fmin=None, fmax=None):
        # model = waveform_model(params, f, psd)
        # residuals = strain - model
        # X.L.: here the input type of rTimeSeries_innerproduct must be rTimeSeries
        #       this function calculate inner product <a|a> = 4\int a^2/psd(f) df,
        #       here I haven't specified the integration range
        SNRF = rTimeSeries_innerproduct(residuals, detector.psd, fmin=None, fmax=None)
        if SNRF < 0 or np.isnan(SNRF):
            return -np.inf
        return SNRF

    def log_likelihood(params):
        q, Mchirp, th1, ph1, chi1, th2, ph2, chi2, ecc, dL, theta, phi, iota, psi, t_c, Phi_c = params
        hmodel_rts_list = simulate_waveform_from_detector_list(
            (q, Mchirp, th1, ph1, chi1, th2, ph2, chi2, ecc, dL, theta, phi, iota, psi, t_c, Phi_c), gps0=gps0, f_min=f_min, detectors=detector_list, is_only22=False, prec_flag=3, code_version=args.code_version, debug=args.code_debug
        )
        if hmodel_rts_list is None:
            return -np.inf
        yres_rts_list = [rTimeSeries_add(signal_rts_list[i], -hmodel_rts_list[i], srate=sample_rate, method=1) for i in range(len(detector_list))]
        lnlike = np.sum([-0.5 * residuals_inner_product(yres_rts_list[i], detector_list[i], fmin=None, fmax=None) for i in range(len(detector_list))])
        return lnlike

    # Define the log likelihood
    def log_likelihood_full(params):
        q, Mchirp, th1, ph1, chi1, th2, ph2, chi2, e0, dL, theta, phi, iota, psi, t_c, Phi_c = params
        hmodel_rts_list = simulate_waveform_from_detector_list(
            (q, Mchirp, th1, ph1, chi1, th2, ph2, chi2, e0, dL, theta, phi, iota, psi, t_c, Phi_c), gps0=gps0, f_min=f_min, detectors=detector_list, is_only22=False, prec_flag=3, code_version=args.code_version, debug=args.code_debug
        )
        if hmodel_rts_list is None:
            return -np.inf
        yres_rts_list = [rTimeSeries_add(signal_rts_list[i], -hmodel_rts_list[i], srate=sample_rate, method=1) for i in range(len(detector_list))]
        lnlike = np.sum([-0.5 * residuals_inner_product(yres_rts_list[i], detector_list[i], fmin=None, fmax=None) for i in range(len(detector_list))])
        return lnlike

    def log_prior_full(params):
        q, Mchirp, th1, ph1, chi1, th2, ph2, chi2, e0, dL, theta, phi, iota, psi, t_c, Phi_c = params
        # we may include the spin parameters here if it is more convenient
        # X.L. It is best not to set the eccentricity over 0.6, otherwise it will be easy to raise problems.
        # chi1 = np.array([s1x, s1y, s1z])
        # chi2 = np.array([s2x, s2y, s2z])
        if (
            1 <= q <= 20
            and (5.0 <= Mchirp <= 200.0)
            and 0 <= chi1 <= 0.95
            and 0 <= chi2 <= 0.95
            and 0.0 <= th1 <= np.pi
            and 0.0 <= ph1 <= 2.0 * np.pi
            and 0.0 <= th2 <= np.pi
            and 0.0 <= ph2 <= 2.0 * np.pi
            and (0.0 <= e0 < args.emax)
            and dL > 0.0
            and 0.0 <= theta <= np.pi
            and 0.0 <= phi <= 2.0 * np.pi
            and 0.0 <= psi <= np.pi
            and 0.0 <= iota <= np.pi
            and -2.0 * t_extra < t_c < 2.0 * t_extra
            and 0.0 <= Phi_c <= 2.0 * np.pi
        ):
            return 0.0
        return -np.inf

    def log_probability_full(params):
        lp = log_prior_full(params)
        if not np.isfinite(lp):
            return -np.inf
        return lp + log_likelihood_full(params)

    def log_likelihood_PrecCirc(params):
        q, Mchirp, th1, ph1, chi1, th2, ph2, chi2, dL, theta, phi, iota, psi, t_c, Phi_c = params
        hmodel_rts_list = simulate_waveform_from_detector_list(
            (q, Mchirp, th1, ph1, chi1, th2, ph2, chi2, 0.0, dL, theta, phi, iota, psi, t_c, Phi_c), gps0=gps0, f_min=f_min, detectors=detector_list, is_only22=args.is_only22, code_version=0, debug=args.code_debug
        )
        if hmodel_rts_list is None:
            return -np.inf
        yres_rts_list = [rTimeSeries_add(signal_rts_list[i], -hmodel_rts_list[i], srate=sample_rate, method=1) for i in range(len(detector_list))]
        lnlike = np.sum([-0.5 * residuals_inner_product(yres_rts_list[i], detector_list[i], fmin=None, fmax=None) for i in range(len(detector_list))])
        return lnlike

    def log_prior_PrecCirc(params):
        q, Mchirp, th1, ph1, chi1, th2, ph2, chi2, dL, theta, phi, iota, psi, t_c, Phi_c = params
        # we may include the spin parameters here if it is more convenient
        # X.L. It is best not to set the eccentricity over 0.6, otherwise it will be easy to raise problems.
        # chi1 = np.array([s1x, s1y, s1z])
        # chi2 = np.array([s2x, s2y, s2z])
        if (
            1 <= q <= 20
            and (5.0 <= Mchirp <= 200.0)
            and 0 <= chi1 <= 0.95
            and 0 <= chi2 <= 0.95
            and 0.0 <= th1 <= np.pi
            and 0.0 <= ph1 <= 2.0 * np.pi
            and 0.0 <= th2 <= np.pi
            and 0.0 <= ph2 <= 2.0 * np.pi
            and dL > 0.0
            and 0.0 <= theta <= np.pi
            and 0.0 <= phi <= 2.0 * np.pi
            and 0.0 <= psi <= np.pi
            and 0.0 <= iota <= np.pi
            and -2.0 * t_extra < t_c < 2.0 * t_extra
            and 0.0 <= Phi_c <= 2.0 * np.pi
        ):
            return 0.0
        return -np.inf

    def log_probability_PrecCirc(params):
        lp = log_prior_PrecCirc(params)
        if not np.isfinite(lp):
            return -np.inf
        return lp + log_likelihood_PrecCirc(params)

    def log_likelihood_SA(params):
        q, Mchirp, chi1, chi2, e0, dL, theta, phi, iota, psi, t_c, Phi_c = params
        hmodel_rts_list = simulate_waveform_from_detector_list(
            (q, Mchirp, 0, 0, chi1, 0, 0, chi2, e0, dL, theta, phi, iota, psi, t_c, Phi_c), gps0=gps0, f_min=f_min, egw_flag=1, detectors=detector_list, is_only22=args.is_only22, debug=args.code_debug
        )
        if hmodel_rts_list is None:
            return -np.inf
        yres_rts_list = [rTimeSeries_add(signal_rts_list[i], -hmodel_rts_list[i], srate=sample_rate, method=1) for i in range(len(detector_list))]
        lnlike = np.sum([-0.5 * residuals_inner_product(yres_rts_list[i], detector_list[i], fmin=None, fmax=None) for i in range(len(detector_list))])
        return lnlike

    def log_prior_SA(params):
        q, Mchirp, chi1, chi2, e0, dL, theta, phi, iota, psi, t_c, Phi_c = params
        # we may include the spin parameters here if it is more convenient
        # X.L. It is best not to set the eccentricity over 0.6, otherwise it will be easy to raise problems.
        if (
            1 <= q <= 20
            and (5.0 <= Mchirp <= 200.0)
            and 0 <= chi1 <= 0.95
            and 0 <= chi2 <= 0.95
            and dL > 0.0
            and 0.0 <= e0 < args.emax
            and 0.0 <= theta <= np.pi
            and 0.0 <= phi <= 2.0 * np.pi
            and 0.0 <= psi <= np.pi
            and 0.0 <= iota <= np.pi
            and -2.0 * t_extra < t_c < 2.0 * t_extra
            and 0.0 <= Phi_c <= 2.0 * np.pi
        ):
            return 0.0
        return -np.inf

    def log_probability_SA(params):
        lp = log_prior_SA(params)
        if not np.isfinite(lp):
            return -np.inf
        return lp + log_likelihood_SA(params)

    if is_pycbc_loaded:
        # Define the log likelihood
        def log_likelihood_pycbc(params):
            q, Mchirp, th1, ph1, chi1, th2, ph2, chi2, dL, theta, phi, iota, psi, t_c, Phi_c = params
            # hmodel_rts_list = simulate_waveform_from_detector_list((m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, 0.0, dL, theta, phi, iota, psi, t_c, Phi_c),
            #             gps0 = gps0, f_min = f_min,
            #             detectors = detector_list, is_only22 = False, code_version = 0)
            hmodel_rts_list = simulate_waveform_from_detector_list_pycbc(
                (q, Mchirp, th1, ph1, chi1, th2, ph2, chi2, 0.0, dL, theta, phi, iota, psi, t_c, Phi_c), gps0=gps0, f_min=f_min, approx=args.model, detectors=detector_list, sample_rate=sample_rate
            )
            if hmodel_rts_list is None:
                return -np.inf
            yres_rts_list = [rTimeSeries_add(signal_rts_list[i], -hmodel_rts_list[i], srate=sample_rate, method=1) for i in range(len(detector_list))]
            lnlike = np.sum([-0.5 * residuals_inner_product(yres_rts_list[i], detector_list[i], fmin=None, fmax=None) for i in range(len(detector_list))])
            return lnlike

        def log_prior_pycbc(params):
            q, Mchirp, th1, ph1, chi1, th2, ph2, chi2, dL, theta, phi, iota, psi, t_c, Phi_c = params
            # we may include the spin parameters here if it is more convenient
            # X.L. It is best not to set the eccentricity over 0.6, otherwise it will be easy to raise problems.
            # chi1 = np.array([s1x, s1y, s1z])
            # chi2 = np.array([s2x, s2y, s2z])
            if (
                1 <= q <= 20
                and (5.0 <= Mchirp <= 200.0)
                and 0 <= chi1 <= 0.95
                and 0 <= chi2 <= 0.95
                and 0.0 <= th1 <= np.pi
                and 0.0 <= ph1 <= 2.0 * np.pi
                and 0.0 <= th2 <= np.pi
                and 0.0 <= ph2 <= 2.0 * np.pi
                and dL > 0.0
                and 0.0 <= theta <= np.pi
                and 0.0 <= phi <= 2.0 * np.pi
                and 0.0 <= psi <= np.pi
                and 0.0 <= iota <= np.pi
                and -2.0 * t_extra < t_c < 2.0 * t_extra
                and 0.0 <= Phi_c <= 2.0 * np.pi
            ):
                return 0.0
            return -np.inf

        def log_probability_pycbc(params):
            lp = log_prior_pycbc(params)
            if not np.isfinite(lp):
                return -np.inf
            return lp + log_likelihood_pycbc(params)

    model_use = args.model.lower()
    max_steps = args.max_steps
    nwalkers = args.nwalkers
    nthreads = args.nthreads
    thresh = args.thresh
    sigma = args.sigma0
    if model_use == 'full':
        ndim = 16
        p0 = np.asarray(params0)
        sys.stderr.write(f'mcmc p0 = {p0}\n')
        if args.debug:
            testprob = log_likelihood(params0)
            sys.stderr.write(f'test run: prob0 = {testprob}\n')
            tmpp0 = p0 + sigma * np.abs(p0) * np.random.randn(ndim)
            testprob = log_probability_full(tmpp0)
            sys.stderr.write(f'test run 1: p = {tmpp0}\n')
            sys.stderr.write(f'test run 1: prob = {testprob}\n')
            tmpp0 = p0 + sigma * np.abs(p0) * np.random.randn(ndim)
            testprob = log_probability_full(tmpp0)
            sys.stderr.write(f'test run 2: p = {tmpp0}\n')
            sys.stderr.write(f'test run 2: prob = {testprob}\n')
            tmpp0 = p0 + sigma * np.abs(p0) * np.random.randn(ndim)
            testprob = log_probability_full(tmpp0)
            sys.stderr.write(f'test run 2: p = {tmpp0}\n')
            sys.stderr.write(f'test run 2: prob = {testprob}\n')
        sampler = MCMCRunner(ndim, nwalkers, log_probability_full, threads=nthreads, conv_thresh=thresh)
    elif model_use == 'preccirc' or model_use == 'pc':
        ndim = 15
        # try call logprob function
        # m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, dL, theta, phi, iota, psi, t_c, Phi_c
        # p0 = np.array([m1, m2, spinChi1[0], spinChi1[1], spinChi1[2], spinChi2[0], spinChi2[1], spinChi2[2], \
        #         dL, theta, phi, iota, psi, 0.0, Phi_c])
        p0 = np.array(
            [
                params0[0],
                params0[1],  # q, Mchirp
                params0[2],
                params0[3],
                params0[4],  # th1, ph1, chi1
                params0[5],
                params0[6],
                params0[7],  # th2, ph2, chi2
                params0[9],
                params0[10],
                params0[11],  # dL, theta, phi
                params0[12],
                params0[13],
                params0[14],
                params0[15],
            ]
        )  # iota, psi, t_c, Phi_c
        sys.stderr.write(f'mcmc p0 = {p0}\n')
        if args.debug:
            testprob = log_likelihood(params0)
            sys.stderr.write(f'test run: prob0 = {testprob}\n')
            tmpp0 = p0 + sigma * np.abs(p0) * np.random.randn(ndim)
            testprob = log_probability_PrecCirc(tmpp0)
            sys.stderr.write(f'test run 1: p = {tmpp0}\n')
            sys.stderr.write(f'test run 1: prob = {testprob}\n')
            tmpp0 = p0 + sigma * np.abs(p0) * np.random.randn(ndim)
            testprob = log_probability_PrecCirc(tmpp0)
            sys.stderr.write(f'test run 2: p = {tmpp0}\n')
            sys.stderr.write(f'test run 2: prob = {testprob}\n')
            tmpp0 = p0 + sigma * np.abs(p0) * np.random.randn(ndim)
            testprob = log_probability_PrecCirc(tmpp0)
            sys.stderr.write(f'test run 2: p = {tmpp0}\n')
            sys.stderr.write(f'test run 2: prob = {testprob}\n')
        sampler = MCMCRunner(ndim, nwalkers, log_probability_PrecCirc, threads=nthreads, conv_thresh=thresh)
        # sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability_PrecCirc)
        # sampler.run_mcmc(pos, 5000, progress=True)
    elif model_use == 'spinaligned' or model_use == 'sa':
        ndim = 12
        # m1, m2, chi1, chi2, e0, dL, theta, phi, iota, psi, t_c, Phi_c
        p0 = np.array([params0[0], params0[1], params0[4], params0[7], params0[8], params0[9], params0[10], params0[11], params0[12], params0[13], params0[14], params0[15]])  # q, Mchirp, chi1, chi2, e0  # dL, theta, phi  # iota, psi, t_c, Phi_c
        sys.stderr.write(f'mcmc p0 = {p0}\n')
        if args.debug:
            testprob = log_likelihood(params0)
            sys.stderr.write(f'test run: prob0 = {testprob}\n')
            tmpp0 = p0 + sigma * np.abs(p0) * np.random.randn(ndim)
            testprob = log_probability_SA(tmpp0)
            sys.stderr.write(f'test run 1: p = {tmpp0}\n')
            sys.stderr.write(f'test run 1: prob = {testprob}\n')
            tmpp0 = p0 + sigma * np.abs(p0) * np.random.randn(ndim)
            testprob = log_probability_SA(tmpp0)
            sys.stderr.write(f'test run 2: p = {tmpp0}\n')
            sys.stderr.write(f'test run 2: prob = {testprob}\n')
            tmpp0 = p0 + sigma * np.abs(p0) * np.random.randn(ndim)
            testprob = log_probability_SA(tmpp0)
            sys.stderr.write(f'test run 2: p = {tmpp0}\n')
            sys.stderr.write(f'test run 2: prob = {testprob}\n')
        sampler = MCMCRunner(ndim, nwalkers, log_probability_SA, nthreads, thresh)
        # sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability_SA)
        # sampler.run_mcmc(pos, 5000, progress=True)
    elif is_pycbc_loaded:
        ndim = 15
        # try call logprob function
        # m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, dL, theta, phi, iota, psi, t_c, Phi_c
        # p0 = np.array([m1, m2, spinChi1[0], spinChi1[1], spinChi1[2], spinChi2[0], spinChi2[1], spinChi2[2], \
        #         dL, theta, phi, iota, psi, 0.0, Phi_c])
        p0 = np.array(
            [
                params0[0],
                params0[1],  # q, Mchirp
                params0[2],
                params0[3],
                params0[4],  # th1, ph1, chi1
                params0[5],
                params0[6],
                params0[7],  # th2, ph2, chi2
                params0[9],
                params0[10],
                params0[11],  # dL, theta, phi
                params0[12],
                params0[13],
                params0[14],
                params0[15],
            ]
        )  # iota, psi, t_c, Phi_c
        sys.stderr.write(f'mcmc p0 = {p0}\n')
        if args.debug:
            testprob = log_likelihood(params0)
            sys.stderr.write(f'test run: prob0 = {testprob}\n')
            tmpp0 = p0 + sigma * np.abs(p0) * np.random.randn(ndim)
            testprob = log_probability_pycbc(tmpp0)
            sys.stderr.write(f'test run 1: p = {tmpp0}\n')
            sys.stderr.write(f'test run 1: prob = {testprob}\n')
            tmpp0 = p0 + sigma * np.abs(p0) * np.random.randn(ndim)
            testprob = log_probability_pycbc(tmpp0)
            sys.stderr.write(f'test run 2: p = {tmpp0}\n')
            sys.stderr.write(f'test run 2: prob = {testprob}\n')
            tmpp0 = p0 + sigma * np.abs(p0) * np.random.randn(ndim)
            testprob = log_probability_pycbc(tmpp0)
            sys.stderr.write(f'test run 2: p = {tmpp0}\n')
            sys.stderr.write(f'test run 2: prob = {testprob}\n')
        sampler = MCMCRunner(ndim, nwalkers, log_probability_pycbc, nthreads, thresh)
        # sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability_SA)
        # sampler.run_mcmc(pos, 5000, progress=True)
    else:
        raise Exception('model is not supported')

    if args.debug:
        return 0
    # 4. dump parameter infos
    fparamnames = datadir / f'{prefix}.paramnames'
    franges = datadir / f'{prefix}.ranges'
    if jobid == 1:
        write_paramnames(fparamnames, ndim)
        write_ranges(franges, ndim, t_extra, args.emax)
    # 5. run MCMC
    fchain = datadir / f'{prefix}_{jobid}.txt'
    fconv = datadir / f'{prefix}_{jobid}.converge_status'
    sys.stderr.write(f'Run MCMC...\n')
    sys.stderr.write(f'dump chain to {fchain}\n')
    sampler.run(p0, max_steps, fchain, fconv, sigma=sigma, progress=args.progress, is_delete=is_delete)
    sys.stderr.write(f'Finished\n')
    return 0
