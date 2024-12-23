import warnings

warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")

import numpy as np
from lal import MTSUN_SI, LIGOTimeGPS
from pycbc.types import FrequencySeries, TimeSeries

from pySEOBNREPHM.waveform import calculate_waveform_ep


def normalize_spin_parameters(p):
    """
    Normalizes the spin parameters for two objects in the input dictionary.

    Parameters:
        p (dict): A dictionary containing spin components with keys:
                  - 'spin1x', 'spin1y', 'spin1z'
                  - 'spin2x', 'spin2y', 'spin2z'

    Returns:
        dict: A new dictionary with normalized spin components.

    Raises:
        ValueError: If the norm of any spin vector is zero, preventing normalization.
    """
    # Extract spin components for spin1 and spin2
    spin1 = np.array([p.get("spin1x", 0.0), p.get("spin1y", 0.0), p.get("spin1z", 0.0)])
    spin2 = np.array([p.get("spin2x", 0.0), p.get("spin2y", 0.0), p.get("spin2z", 0.0)])

    # Calculate norms
    norm1 = np.linalg.norm(spin1)
    norm2 = np.linalg.norm(spin2)

    # Check for zero norms to avoid division by zero
    if norm1 == 0:
        raise ValueError("The norm of spin1 is zero. Cannot normalize a zero vector.")
    if norm2 == 0:
        raise ValueError("The norm of spin2 is zero. Cannot normalize a zero vector.")

    # Normalize spin vectors
    normalized_spin1 = spin1 / norm1
    normalized_spin2 = spin2 / norm2

    # Create a new dictionary with normalized spins
    normalized_p = p.copy()
    normalized_p["spin1x"], normalized_p["spin1y"], normalized_p["spin1z"] = normalized_spin1
    normalized_p["spin2x"], normalized_p["spin2y"], normalized_p["spin2z"] = normalized_spin2

    return normalized_p


def gen_seobnrephm_td(**p):
    # p.update(
    #     {
    #         "approximant": "SEOBNRv5EHM",  # I call it "SEOBNRv5E" in PyCBC
    #         "ModeArray": [(2, 2)],  # only consider (2,2) mode
    #         "rel_anomaly": p[
    #             "rel_anomaly"
    #         ],  # relativity anomaly,  needed for eccentric waveform
    #         "phi_ref": p["coa_phase"],  # reference phase needed by SEOBNRv5
    #         "f22_start": p["f_lower"],  # starting frequency
    #         "f_ref": p["f_lower"],  # reference frequency
    #         "deltaT": p["delta_t"],
    #         # "postadiabatic": False ,   # turn off postadiabatic correction,
    #         # default is False in SEOBNRv5EHM
    #         # "h_0": 1.0,                # initial time step in the integration of the ODEs.
    #         # Default Value is 1.0
    #         "lmax_nyquist": 1,  # maximum L to be checked against Nyquist frequency
    #     }
    # )
    normalized_p = normalize_spin_parameters(p)
    waveform, dynamics = calculate_waveform_ep(
        (
            p["mass1"],
            p["mass2"],
            normalized_p["spin1x"],
            normalized_p["spin1y"],
            normalized_p["spin1z"],
            normalized_p["spin2x"],
            normalized_p["spin2y"],
            normalized_p["spin2z"],
            p["eccentricity"],
            p["distance"],
            p["rel_anomaly"] if "rel_anomaly" in p else 0,
            p["inclination"],
            0,  # beta_rad,
            p["coa_phase"],
        ),
        p["f_lower"],
        Mf_ref=p["f_lower"] * (p["mass1"] + p["mass2"]) * MTSUN_SI,
        srate=1.0 / p["delta_t"],
        is_coframe=False,
    )

    _hp = waveform.hpc.real
    _hc = -waveform.hpc.imag
    t_peak = waveform.hpc.time[waveform.hpc.argpeak]
    t_vec = waveform.time - t_peak
    epoch = LIGOTimeGPS(t_vec[0])

    # Build the PyCBC TimeSeries format
    hp = TimeSeries(_hp, delta_t=p["delta_t"], epoch=epoch)
    hc = TimeSeries(_hc, delta_t=p["delta_t"], epoch=epoch)

    return hp, hc


# def gen_seobnrv5ehm_td(**p):
#     p.update(
#         {
#             "approximant": "SEOBNRv5EHM",
#             "rel_anomaly": (
#                 p["rel_anomaly"] if "rel_anomaly" in p else 0
#             ),  # relativity anomaly,  needed for eccentric waveform
#             "phi_ref": p["coa_phase"],  # reference phase needed by SEOBNRv5
#             "f22_start": p["f_lower"],  # starting frequency
#             "f_ref": p["f_lower"],  # reference frequency
#             "deltaT": p["delta_t"],
#             "lmax_nyquist": 1,  # maximum L to be checked against Nyquist frequency
#         }
#     )

#     waveform = GenerateWaveform(p)
#     hp, hc = waveform.generate_td_polarizations()

#     # Build the PyCBC TimeSeries format
#     hp_pycbc = TimeSeries(hp.data.data[:], delta_t=hp.deltaT, epoch=hp.epoch)
#     hc_pycbc = TimeSeries(hc.data.data[:], delta_t=hp.deltaT, epoch=hp.epoch)

#     return hp_pycbc, hc_pycbc


# def gen_seobnrv5e_fd(**p):
#     p.update(
#         {
#             "approximant": "SEOBNRv5EHM",
#             "ModeArray": [(2, 2)],  # only consider (2,2) mode
#             "rel_anomaly": (
#                 p["rel_anomaly"] if "rel_anomaly" in p else 0
#             ),  # relativity anomaly,  needed for eccentric waveform
#             "phi_ref": p["coa_phase"],  # reference phase needed by SEOBNRv5
#             "f22_start": p["f_lower"],  # starting frequency
#             "f_ref": p["f_lower"],  # reference frequency
#             "deltaF": p["delta_f"],
#             # "postadiabatic": False ,   # turn off postadiabatic correction,
#             # default is False in SEOBNRv5EHM
#             # "h_0": 1.0,                # initial time step in the integration of the ODEs.
#             # Default Value is 1.0
#             "lmax_nyquist": 1,  # maximum L to be checked against Nyquist frequency
#         }
#     )

#     waveform = GenerateWaveform(p)
#     hp, hc, template_duration = waveform.generate_fd_polarizations()

#     # Build the PyCBC TimeSeries format
#     hp_pycbc = FrequencySeries(hp.data.data[:], delta_f=hp.deltaF, epoch=hp.epoch)
#     hc_pycbc = FrequencySeries(hc.data.data[:], delta_f=hc.deltaF, epoch=hp.epoch)

#     hp_pycbc.eob_template_duration = template_duration
#     hc_pycbc.eob_template_duration = template_duration

#     return hp_pycbc, hc_pycbc


# def gen_seobnrv5ehm_fd(**p):
#     p.update(
#         {
#             "approximant": "SEOBNRv5EHM",
#             "rel_anomaly": (
#                 p["rel_anomaly"] if "rel_anomaly" in p else 0
#             ),  # relativity anomaly,  needed for eccentric waveform
#             "phi_ref": p["coa_phase"],  # reference phase needed by SEOBNRv5
#             "f22_start": p["f_lower"],  # starting frequency
#             "f_ref": p["f_lower"],  # reference frequency
#             "deltaF": p["delta_f"],
#         }
#     )

#     waveform = GenerateWaveform(p)
#     hp, hc, template_duration = waveform.generate_fd_polarizations()

#     # Build the PyCBC TimeSeries format
#     hp_pycbc = FrequencySeries(hp.data.data[:], delta_f=hp.deltaF, epoch=hp.epoch)
#     hc_pycbc = FrequencySeries(hc.data.data[:], delta_f=hc.deltaF, epoch=hp.epoch)

#     hp_pycbc.eob_template_duration = template_duration
#     hc_pycbc.eob_template_duration = template_duration

#     return hp_pycbc, hc_pycbc


# def seobnrv5phm_length_in_time(**kwds):
#     from pycbc.waveform.waveform import get_hm_length_in_time

#     return get_hm_length_in_time("SEOBNRv5", 5, **kwds)


# def seobnrv5e_length_in_time(**kwds):
#    from pycbc.waveform.waveform import get_waveform_filter_length_in_time
#    if "approximant" in kwds:
#        kwds.pop("approximant")
#    return get_waveform_filter_length_in_time(approximant='SEOBNRv5_ROM', **kwds)
