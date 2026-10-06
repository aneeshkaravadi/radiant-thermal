"""Radiative cooling surface physics.

A surface outdoors gains sunlight and downwelling sky infrared, and loses its
own thermal emission plus convection to the air:

    q_net_out = eps*sigma*Ts^4 - eps*[F_sky*L_sky + (1-F_sky)*sigma*Tg^4]
                - alpha_sol*G_sun + h*(Ts - Ta)

The atmosphere is partly transparent between about 8 and 13 um, so the clear
sky radiates like a body at Ta with emissivity below 1. A surface that reflects
sunlight and emits strongly in the infrared can sit below air temperature in
full sun. That is passive daytime radiative cooling.
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


@dataclass(frozen=True)
class Coating:
    name: str
    alpha_solar: float  # solar absorptance (0.3-2.5 um)
    eps_ir: float  # hemispherical thermal emissivity (~5-25 um)
    source: str = ""


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


def steady_surface_temperature(coating: Coating, Ta, G, eps_sky, h, **kw) -> float:
    """Temperature (K) of an insulated surface (no heat from behind) in steady state."""
    return brentq(lambda T: net_heat_out(T, coating, Ta, G, eps_sky, h, **kw), Ta - 60, Ta + 120)


def cooling_power_at_ambient(coating: Coating, Ta, G, eps_sky, F_sky=1.0) -> float:
    """Radiative cooling power (W/m^2) with the surface held at air temperature (convection drops out)."""
    return float(net_heat_out(Ta, coating, Ta, G, eps_sky, 0.0, F_sky=F_sky))


# Representative coatings for comparison. Replace with datasheet values for real designs.
DARK_PAINT = Coating("dark gray paint", 0.75, 0.90, "representative")
WHITE_PAINT = Coating("white paint", 0.25, 0.90, "representative")
IDEAL_COOLER = Coating("ideal radiative cooler", 0.04, 0.95, "upper-bound reference")
