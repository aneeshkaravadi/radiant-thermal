"""Radiative cooling surface physics.

A surface outdoors gains sunlight and downwelling sky infrared, and loses its
own thermal emission plus convection to the air:

    q_net_out = eps*sigma*Ts^4 - eps*[F_sky*L_sky + (1-F_sky)*sigma*Tg^4]
                - alpha_sol*G_sun + h*(Ts - Ta)

The atmosphere is partly transparent between about 8 and 13 um, so the clear
sky radiates like a body at Ta with emissivity below 1. A surface that reflects
sunlight and emits strongly in the infrared can sit below air temperature in
full sun. That is passive daytime radiative cooling.

``net_heat_out`` treats the infrared as one gray band. ``net_heat_out_spectral``
splits it at the 8-13 um window: the sky is taken as black outside the window
and partly transparent inside it, with the window's emissivity set so the total
still matches Berdahl & Martin. A gray surface gets exactly the same answer
either way; a selective emitter (one that emits mainly inside the window) only
gets the right answer from the two-band version.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

SIGMA = 5.670374419e-8  # W/m^2/K^4


def sky_emissivity(dew_point_c, cloud_fraction=0.0):
    """Clear-sky emissivity from dew point (Berdahl & Martin 1984), with a cloud correction.

    Clear sky:  eps = 0.711 + 0.56 (Tdp/100) + 0.73 (Tdp/100)^2,  Tdp in C.
    Clouds are treated as near-black emitters at air temperature covering a
    fraction c of the sky: eps = c + (1 - c) eps_clear (Crawford & Duchon 1999 form).
    """
    t = np.asarray(dew_point_c, float) / 100.0
    clear = 0.711 + 0.56 * t + 0.73 * t**2
    c = np.clip(np.asarray(cloud_fraction, float), 0.0, 1.0)
    return c + (1.0 - c) * clear


def convection_coefficient(wind_speed):
    """External forced/natural convection, h = 5.7 + 3.8 v  (McAdams), W/m^2/K, v in m/s."""
    return 5.7 + 3.8 * np.asarray(wind_speed, float)


C2 = 14387.77  # um K, second radiation constant hc/k
WINDOW_UM = (8.0, 13.0)


def blackbody_fraction(lam_t):
    """Fraction of blackbody emission below wavelength x temperature ``lam_t`` (um K).

    The exact series  F = (15/pi^4) sum_n e^(-n z)/n (z^3 + 3z^2/n + 6z/n^2 + 6/n^3),  z = C2/(lam T).
    """
    z = C2 / np.asarray(lam_t, float)
    total = 0.0
    for n in range(1, 40):
        total = total + np.exp(-n * z) / n * (z**3 + 3 * z**2 / n + 6 * z / n**2 + 6 / n**3)
    return 15.0 / np.pi**4 * total


def window_fraction(T):
    """Share of a blackbody's emission at temperature T (K) that falls inside the 8-13 um window."""
    T = np.asarray(T, float)
    return blackbody_fraction(WINDOW_UM[1] * T) - blackbody_fraction(WINDOW_UM[0] * T)


_T_TABLE = np.arange(150.0, 450.0001, 0.05)
_FW_TABLE = window_fraction(_T_TABLE)


def window_fraction_fast(T: float) -> float:
    """``window_fraction`` from a 0.05 K table, for the inner loop of the enclosure model."""
    return float(np.interp(T, _T_TABLE, _FW_TABLE))


def sky_band_emissivities(eps_sky_total, Ta):
    """Sky emissivity (inside, outside) the 8-13 um window, keeping the total at ``eps_sky_total``.

    Outside the window water vapor and CO2 make the sky nearly black, so it's taken as 1 there,
    and the window gets whatever emissivity makes the band-weighted total match. In very dry air
    that would go below zero; then the window is fully clear and the outside band takes the rest.
    """
    fw = window_fraction(Ta)
    eps = np.asarray(eps_sky_total, float)
    inside = np.clip(1.0 - (1.0 - eps) / fw, 0.0, 1.0)
    outside = (eps - fw * inside) / (1.0 - fw)
    return inside, outside


@dataclass(frozen=True)
class Coating:
    name: str
    alpha_solar: float  # solar absorptance (0.3-2.5 um)
    eps_ir: float  # hemispherical thermal emissivity (~5-25 um), as a total-emissivity measurement at ~300 K reports
    source: str = ""
    bands: tuple[float, float] | None = None  # (emissivity inside the 8-13 um window, outside it); None = gray

    @property
    def eps_bands(self) -> tuple[float, float]:
        return (self.eps_ir, self.eps_ir) if self.bands is None else self.bands


def selective_coating(name: str, alpha_solar: float, eps_window: float, eps_outside: float, source: str = "",
                      T_ref: float = 300.0) -> Coating:
    """A coating that emits ``eps_window`` inside 8-13 um and ``eps_outside`` elsewhere.

    Its gray ``eps_ir`` is the total emissivity at ``T_ref``, what a total-emissivity
    measurement would report, which is all the gray model gets to see.
    """
    fw = float(window_fraction(T_ref))
    return Coating(name, alpha_solar, fw * eps_window + (1 - fw) * eps_outside, source, (eps_window, eps_outside))


@dataclass(frozen=True)
class Film:
    """A semi-transparent laminate, wrap or film: solar R/T/A and IR emissivity."""

    name: str
    R_solar: float
    T_solar: float
    eps_ir: float
    source: str = ""

    @property
    def A_solar(self) -> float:
        return 1.0 - self.R_solar - self.T_solar

    def on(self, substrate: Coating) -> Coating:
        """Effective coating of this film laminated on a substrate.

        Solar light transmitted through the film bounces between substrate and film.
        Summing the geometric series of bounces (film treated as symmetric):

            alpha_eff = A_f + T_f * [alpha_s + (1-alpha_s) * A_f] / (1 - R_f (1-alpha_s))

        In the IR the film is assumed optically thick (true for most polymer
        films over ~50 um), so the emissivity is the film's.
        """
        a_s = substrate.alpha_solar
        Rf, Tf, Af = self.R_solar, self.T_solar, self.A_solar
        alpha = Af + Tf * (a_s + (1 - a_s) * Af) / (1 - Rf * (1 - a_s))
        return Coating(f"{self.name} on {substrate.name}", float(alpha), self.eps_ir, self.source)


def net_heat_out(Ts, coating: Coating, Ta, G, eps_sky, h, F_sky=1.0, Tg=None, sun_factor=1.0):
    """Net heat leaving 1 m^2 of surface (W/m^2). Temperatures in K. Positive = cooling the surface."""
    Tg = Ta if Tg is None else Tg
    emit = coating.eps_ir * SIGMA * Ts**4
    absorb_ir = coating.eps_ir * (F_sky * eps_sky * SIGMA * Ta**4 + (1 - F_sky) * SIGMA * Tg**4)
    absorb_sun = coating.alpha_solar * G * sun_factor
    return emit - absorb_ir - absorb_sun + h * (Ts - Ta)


def net_heat_out_spectral(Ts, coating: Coating, Ta, G, eps_sky, h, F_sky=1.0, Tg=None, sun_factor=1.0):
    """``net_heat_out`` with the infrared split at the 8-13 um window (two bands)."""
    Tg = Ta if Tg is None else Tg
    e_in, e_out = coating.eps_bands
    s_in, s_out = sky_band_emissivities(eps_sky, Ta)
    fs, fa, fg = window_fraction(Ts), window_fraction(Ta), window_fraction(Tg)
    emit = SIGMA * Ts**4 * (e_in * fs + e_out * (1 - fs))
    sky = F_sky * SIGMA * Ta**4 * (e_in * s_in * fa + e_out * s_out * (1 - fa))
    ground = (1 - F_sky) * SIGMA * Tg**4 * (e_in * fg + e_out * (1 - fg))
    return emit - sky - ground - coating.alpha_solar * G * sun_factor + h * (Ts - Ta)


def steady_surface_temperature(coating: Coating, Ta, G, eps_sky, h, spectral: bool = False, **kw) -> float:
    """Temperature (K) of an insulated surface (no heat from behind) in steady state."""
    f = net_heat_out_spectral if spectral else net_heat_out
    return brentq(lambda T: f(T, coating, Ta, G, eps_sky, h, **kw), Ta - 60, Ta + 120)


def cooling_power_at_ambient(coating: Coating, Ta, G, eps_sky, F_sky=1.0, spectral: bool = False) -> float:
    """Radiative cooling power (W/m^2) with the surface held at air temperature (convection drops out)."""
    f = net_heat_out_spectral if spectral else net_heat_out
    return float(f(Ta, coating, Ta, G, eps_sky, 0.0, F_sky=F_sky))


# Representative coatings for comparison. Replace with datasheet values for real designs.
DARK_PAINT = Coating("dark gray paint", 0.75, 0.90, "representative")
WHITE_PAINT = Coating("white paint", 0.25, 0.90, "representative")
IDEAL_COOLER = Coating("ideal radiative cooler", 0.04, 0.95, "upper-bound reference")
