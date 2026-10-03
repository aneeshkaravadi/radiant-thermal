"""Spacecraft radiator sizing.

A radiator at temperature T rejects eps*sigma*T^4 per m^2 to deep space but
absorbs sunlight, Earth infrared and albedo:

    q_abs = alpha * S * cos(theta_sun) + F_earth * (eps * q_earth_ir + alpha * a * S)
    A     = Q / (eps * sigma * T^4 - q_abs)

Here alpha/eps is the figure of merit for a hot environment, and eps alone
for a cold one. The same surface physics as a rooftop cooler, minus the atmosphere.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .radiative import SIGMA, Coating

SOLAR_CONSTANT = 1361.0  # W/m^2
EARTH_IR = 237.0  # W/m^2, orbit-average outgoing longwave
ALBEDO = 0.30


@dataclass(frozen=True)
class Environment:
    name: str
    sun_cos: float  # cosine of sun incidence on the radiator (0 = edge-on / shaded)
    earth_view: float  # view factor to Earth
    distance_au: float = 1.0


DEEP_SPACE = Environment("sun-shaded, no planet", 0.0, 0.0)
LEO_SHADED = Environment("LEO, sun-shaded, Earth view 0.3", 0.0, 0.3)
LEO_HOT = Environment("LEO hot case, sun 60 deg off normal, Earth view 0.3", 0.5, 0.3)


def absorbed_flux(c: Coating, env: Environment) -> float:
    S = SOLAR_CONSTANT / env.distance_au**2
    return c.alpha_solar * S * env.sun_cos + env.earth_view * (c.eps_ir * EARTH_IR + c.alpha_solar * ALBEDO * S)


def area_per_kw(c: Coating, env: Environment, T_kelvin) -> np.ndarray:
    """Radiator area (m^2) to reject 1 kW at radiator temperature T. inf where it cannot reject heat."""
    net = c.eps_ir * SIGMA * np.asarray(T_kelvin, float) ** 4 - absorbed_flux(c, env)
    return np.where(net > 0, 1000.0 / np.where(net > 0, net, 1.0), np.inf)


def equilibrium_temperature(c: Coating, env: Environment, q_internal_w_m2: float = 0.0) -> float:
    """Temperature a radiator panel settles at with a given internal heat load per m^2."""
    return float(((absorbed_flux(c, env) + q_internal_w_m2) / (c.eps_ir * SIGMA)) ** 0.25)


# Representative spacecraft surfaces (typical beginning-of-life values from thermal-control handbooks).
OSR = Coating("optical solar reflector", 0.08, 0.80, "representative BOL")
WHITE_SPACE_PAINT = Coating("white thermal-control paint", 0.20, 0.90, "representative BOL")
BLACK_PAINT = Coating("black paint", 0.95, 0.88, "representative")
