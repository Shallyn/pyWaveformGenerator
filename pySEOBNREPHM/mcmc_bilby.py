#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue, 18 Apr 2023 09:47:34 +0000

@author: Shallyn
"""
import os
import sys
from optparse import OptionParser
from pathlib import Path

import bilby
import h5py
import matplotlib.pyplot as plt
import numpy as np
from gwpy.timeseries import TimeSeries

from .mcmc import load_rts_list_from_hdf5

# pwd = Path(__file__).absolute().parent
# pwd = Path(sys.path[0])
pwd = Path(os.getcwd())


def run_bilby_mcmc(argv):
    parser = OptionParser(description='FIXME')
    parser.add_option('--detector', type='str', action='append', help='detector prefix')
    parser.add_option('--fsignal', type='str', help='signal file')
    parser.add_option('--sample-rate', type='int', help='sample rate')
    parser.add_option('--outdir', type='str', help='output dir')
    parser.add_option('--label', type='str', help='job label')
    args, _ = parser.parse_args(argv)

    detector_name_list = []
    if args.detector is None or len(args.detector) == 0:
        raise Exception(f'detectors are not specified')
    detector_name_list = args.detector

    fsignal = pwd / args.fsignal
    if not fsignal.exists():
        raise Exception(f'{args.fsignal} do not exists')
    outdir = args.outdir
    label = args.label
    # load parameters from fsignal
    signal_rts_list = load_rts_list_from_hdf5(fsignal, 'signal', detector_name_list)
    duration = signal_rts_list[0].duration
    sampling_frequency = signal_rts_list[0].fs
    gps0 = 1356566418
    # set waveform generator
    # Frequency band
    flow = 10
    fhigh = 512
    prior = bilby.core.prior.PriorDict()
    # prior['chirp_mass'] = bilby.core.prior.Uniform(name='chirp_mass', minimum=10.0,maximum=100.0)
    # prior['mass_ratio'] = bilby.core.prior.Uniform(name='mass_ratio', minimum=0.1, maximum=1)
    prior['mass_1'] = bilby.core.prior.Uniform(name='mass_1', minimum=10, maximum=75)
    prior['mass_2'] = bilby.core.prior.Uniform(name='mass_2', minimum=10, maximum=75)
    # prior['a_1'] =  0.0
    # prior['a_2'] =  0.0
    # prior['tilt_1'] =  0.0
    # prior['tilt_2'] =  0.0
    # prior['phi_12'] =  0.0
    # prior['phi_jl'] =  0.0

    # prior['luminosity_distance'] = bilby.core.prior.Uniform(name="luminosity_distance", minimum=10., maximum=2000.)

    # prior['theta_jn'] = 0.0
    # prior['psi'] = 0.0
    # prior['phase'] = 0.0
    # prior['geocent_time'] = gps0
    # prior['ra'] = 0.0
    # prior['dec'] = 0.0
    # prior['mass_1'] = 30.
    # prior['mass_2'] = 30.
    prior['a_1'] = bilby.core.prior.Uniform(name='a_1', minimum=0.0, maximum=0.95)
    prior['a_2'] = bilby.core.prior.Uniform(name='a_2', minimum=0.0, maximum=0.95)
    prior['tilt_1'] = bilby.core.prior.Uniform(name='tilt_1', minimum=0.0, maximum=np.pi)
    prior['tilt_2'] = bilby.core.prior.Uniform(name='tilt_2', minimum=0.0, maximum=np.pi)
    prior['phi_12'] = bilby.core.prior.Uniform(name='phi_12', minimum=0.0, maximum=2.0 * np.pi, boundary="periodic")
    prior['phi_jl'] = bilby.core.prior.Uniform(name='phi_jl', minimum=0.0, maximum=2.0 * np.pi, boundary="periodic")

    prior['luminosity_distance'] = bilby.core.prior.Uniform(name="luminosity_distance", minimum=10.0, maximum=2000.0)

    prior['theta_jn'] = bilby.core.prior.Uniform(name="theta_jn", minimum=0.0, maximum=np.pi)
    prior['psi'] = bilby.core.prior.Uniform(name="psi", minimum=0, maximum=np.pi)
    prior['phase'] = bilby.core.prior.Uniform(name="phase", minimum=0, maximum=2 * np.pi, boundary="periodic")
    prior['geocent_time'] = bilby.core.prior.Uniform(name="geocent_time", minimum=gps0 - 1.0, maximum=gps0 + 1.0)
    prior['ra'] = bilby.core.prior.Uniform(name="phase", minimum=0, maximum=2 * np.pi, boundary="periodic")
    prior['dec'] = bilby.core.prior.Uniform(name="dec", minimum=-np.pi / 2.0, maximum=np.pi / 2.0)
    # ifos = bilby.gw.detector.InterferometerList(detector_name_list)
    ifos = bilby.gw.detector.InterferometerList([])
    ifo_list = [bilby.gw.detector.get_empty_interferometer(detname) for detname in detector_name_list]
    signal_ts = [TimeSeries(signal.value, t0=signal.t0 + gps0, sample_rate=sampling_frequency) for signal in signal_rts_list]
    for i, ifo in enumerate(ifo_list):
        ifo.set_strain_data_from_gwpy_timeseries(time_series=signal_ts[i])
        ifos.append(ifo)
    # h1_gwRTS = TimeSeries()
    waveform_arguments = dict(
        waveform_approximant="SEOBNRv4PHM",
        reference_frequency=10.0,
        minimum_frequency=10.0,
    )
    waveform_generator = bilby.gw.WaveformGenerator(
        duration=duration,
        sampling_frequency=sampling_frequency,
        frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
        parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
        waveform_arguments=waveform_arguments,
    )
    likelihood = bilby.gw.GravitationalWaveTransient(interferometers=ifos, waveform_generator=waveform_generator, priors=prior)

    result = bilby.run_sampler(
        likelihood=likelihood,
        priors=prior,
        sampler="dynesty",
        nlive=1000,
        walks=32,
        nact=50,
        maxmcmc=2000,
        npool=4,
        outdir=outdir,
        label=label,
        conversion_function=bilby.gw.conversion.generate_all_bbh_parameters,
        result_class=bilby.gw.result.CBCResult,
    )
