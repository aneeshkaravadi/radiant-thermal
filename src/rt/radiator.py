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


# ---------------------------------------------------------------- through an orbit, in and out of eclipse
#
# A radiator sized for the hot, sunlit case has nothing coming in but Earth's
# infrared once the spacecraft passes into Earth's shadow, so it cools. How
# far depends on its heat capacity per square meter and how long the eclipse is.

MU_EARTH = 3.986004418e14  # m^3/s^2
R_EARTH = 6.371e6  # m


def orbit_period(altitude_m: float) -> float:
    """Circular orbit period (s), Kepler's third law."""
    return float(2 * np.pi * np.sqrt((R_EARTH + altitude_m) ** 3 / MU_EARTH))


def eclipse_fraction(altitude_m: float, beta_deg: float = 0.0) -> float:
    """Share of a circular orbit spent in Earth's (cylindrical) shadow.

    beta is the angle between the orbit plane and the sun direction. The orbit is in
    shadow for  cos(theta) < -sqrt(1 - (R/r)^2) / cos(beta),  which gives
    f = acos(sqrt(h^2 + 2 R h) / (r cos beta)) / pi, and no eclipse at all past the beta
    where that argument reaches 1.
    """
    r = R_EARTH + altitude_m
    arg = np.sqrt(altitude_m**2 + 2 * R_EARTH * altitude_m) / (r * np.cos(np.radians(beta_deg)))
    return float(np.arccos(arg) / np.pi) if arg < 1 else 0.0


def orbit_transient(c: Coating, sunlit_env: Environment, areal_heat_capacity: float, q_internal_w_m2: float,
                    altitude_m: float = 400e3, beta_deg: float = 0.0, orbits: int = 6, steps_per_orbit: int = 4000,
                    T0: float | None = None):
    """Temperature of a radiator panel (per square meter) over several orbits.

    In sunlight it absorbs ``absorbed_flux(c, sunlit_env)``; in eclipse only Earth's infrared,
    F_earth eps q_IR. The panel is one thermal mass, c dT/dt = q_int + q_abs - eps sigma T^4,
    stepped with RK4, starting at the sunlit equilibrium unless T0 is given.
    Returns time (s), temperature (K) and a sunlit flag for each step.
    """
    period = orbit_period(altitude_m)
    f_ecl = eclipse_fraction(altitude_m, beta_deg)
    dt = period / steps_per_orbit
    q_sun = absorbed_flux(c, sunlit_env) + q_internal_w_m2
    q_dark = sunlit_env.earth_view * c.eps_ir * EARTH_IR + q_internal_w_m2
    n = orbits * steps_per_orbit
    t = np.arange(n + 1) * dt
    phase = (t % period) / period
    sunlit = phase < 1 - f_ecl  # the eclipse is the last part of each orbit
    T = np.empty(n + 1)
    T[0] = equilibrium_temperature(c, sunlit_env, q_internal_w_m2) if T0 is None else T0

    def rate(temp, q_in):
        return (q_in - c.eps_ir * SIGMA * temp**4) / areal_heat_capacity

    for i in range(n):
        q = q_sun if sunlit[i] else q_dark
        k1 = rate(T[i], q)
        k2 = rate(T[i] + 0.5 * dt * k1, q)
        k3 = rate(T[i] + 0.5 * dt * k2, q)
        k4 = rate(T[i] + dt * k3, q)
        T[i + 1] = T[i] + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return t, T, sunlit
