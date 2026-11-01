"""Physics that scales the simulation to the 6.4800's real hardware: field, magnetization, coil, signal, noise.

The simulation kernel works in relative units (proton density 1 = water). These turn that into volts at the coil:

    B0       = larmor / (gamma / 2 pi)                        107.3 mT at 4.567 MHz
    M0       = rho_H gamma^2 hbar^2 B0 / (4 k T)              equilibrium magnetization of water, A/m
    B1/I     = mu0 N / sqrt(L^2 + 4 a^2)                      field per amp at the centre of an N-turn solenoid
    EMF      = omega0 (B1/I) M0 V                             peak voltage a precessing volume V induces

(the last is reciprocity: the coil's field per amp tells you how well it hears each voxel). At low field the signal
is small because both M0 and omega0 scale with B0, so the signal goes as B0^2.
"""

import math

GAMMA = 2.675_221_874e8          # proton gyromagnetic ratio, rad/s/T
GAMMA_BAR = GAMMA / (2 * math.pi)  # Hz/T (42.577 MHz/T)
HBAR = 1.054_571_817e-34         # J s
K_B = 1.380_649e-23              # J/K
MU0 = 4e-7 * math.pi             # T m/A
RHO_H_WATER = 6.69e28            # hydrogen nuclei per m^3 in water


def b0_tesla(larmor_hz):
    return larmor_hz / GAMMA_BAR


def magnetization(b0_t, temperature_c, proton_density=1.0):
    """Equilibrium magnetization M0 (A/m) of water (scaled by proton_density) at field b0_t and temperature."""
    t_k = temperature_c + 273.15
    return proton_density * RHO_H_WATER * GAMMA**2 * HBAR**2 * b0_t / (4 * K_B * t_k)


def solenoid_field_per_amp(turns, radius_m, length_m):
    """B1/I (T/A) at the centre of a finite solenoid."""
    return MU0 * turns / math.sqrt(length_m**2 + 4 * radius_m**2)


def emf_per_volume(larmor_hz, temperature_c, coil_t_per_a):
    """Peak coil voltage (V) per m^3 of water fully tipped into the transverse plane."""
    b0 = b0_tesla(larmor_hz)
    return 2 * math.pi * larmor_hz * coil_t_per_a * magnetization(b0, temperature_c)


def thermal_noise_density(resistance_ohm, temperature_c):
    """Johnson noise of the coil's resistance, V/sqrt(Hz)."""
    return math.sqrt(4 * K_B * (temperature_c + 273.15) * resistance_ohm)


def db_to_gain(db):
    return 10 ** (db / 20)


def summary(design, coil_t_per_a=None):
    """The derived numbers for a design, as a dict (shown by `python -m sim6480 numbers`).

    coil_t_per_a: the coil's field per amp at the centre, from its map; default: the generated solenoid's formula.
    """
    m, c, r = design["magnet"], design["fields"]["coil"], design["receive"]
    b0 = b0_tesla(m["larmor_hz"])
    if coil_t_per_a is None:
        sol = c["solenoid"]
        coil_t_per_a = solenoid_field_per_amp(sol["turns"], sol["radius_mm"] * 1e-3, sol["length_mm"] * 1e-3)
    coil = coil_t_per_a
    gain_db = r["preamp_db"] + r["lna_db"] + r["pga_db"]
    return {
        "B0 (mT)": b0 * 1e3,
        "M0 of water (A/m)": magnetization(b0, m["temperature_c"]),
        "coil field per amp (mT/A)": coil * 1e3,
        "EMF per mm^3 of water (nV)": emf_per_volume(m["larmor_hz"], m["temperature_c"], coil) * 1e-9 * 1e9,
        "coil thermal noise (nV/rtHz)": thermal_noise_density(c["resistance_ohm"], m["temperature_c"]) * 1e9,
        "noise used (nV/rtHz)": r["noise_nv_per_rthz"],
        "analog gain (dB)": gain_db,
        "volts per ADC count (mV)": r["adc_vpp"] / 2 ** r["adc_bits"] * 1e3,
    }
