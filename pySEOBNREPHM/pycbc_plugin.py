import numpy as np
from lal import LIGOTimeGPS
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
            p["coa_phase"],  # beta_rad,
            0,  # Phi_rad, this is polarization angle so not needed
        ),
        p["f_lower"],
        Mf_ref=0.003,
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


def gen_seobnrep_td(**p):
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
            p["coa_phase"],  # beta_rad,
            0,  # Phi_rad, this is polarization angle so not needed
        ),
        p["f_lower"],
        Mf_ref=0.003,
        srate=1.0 / p["delta_t"],
        is_coframe=False,
        is_only22=True,
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
