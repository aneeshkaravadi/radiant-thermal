"""Outdoor battery enclosure: a two-node thermal model driven by real weather.

    skin:     C_s dTs/dt = A [alpha f_sun G + eps (F L_sky + (1-F) sigma Ta^4) - eps sigma Ts^4 - h (Ts - Ta)]
                           + UA (Tb - Ts)
    battery:  C_b dTb/dt = Q_gen(t) - UA (Tb - Ts)

The coating only enters through alpha and eps on the skin. Everything else is
the same for each design, so differences in battery temperature come from the
coating alone. Parameters describe a generic ~20 kWh wall/ground unit and are
illustrative, not any company's product.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .radiative import SIGMA, Coating, convection_coefficient, sky_emissivity


@dataclass
class Enclosure:
    area: float = 3.0  # m^2 of exterior skin
    skin_heat_capacity: float = 1.5e4  # J/K  (~30 kg steel)
    battery_heat_capacity: float = 2.0e5  # J/K  (~200 kg of cells and modules)
    ua_internal: float = 12.0  # W/K, battery <-> skin conductance
    sky_view: float = 0.75  # box: top sees all sky, walls half
    sun_factor: float = 0.55  # mean projected-area fraction of the skin facing the sun, x GHI


@dataclass
class Weather:
    time: pd.DatetimeIndex
    Ta: np.ndarray  # K
    dew_c: np.ndarray
    ghi: np.ndarray  # W/m^2
    wind: np.ndarray  # m/s
    cloud: np.ndarray  # 0..1

    @classmethod
    def from_csv(cls, path) -> "Weather":
        d = pd.read_csv(path, parse_dates=["time"])
        return cls(pd.DatetimeIndex(d["time"]), d["temperature_2m"].to_numpy() + 273.15, d["dew_point_2m"].to_numpy(),
                   d["shortwave_radiation"].to_numpy(), d["wind_speed_10m"].to_numpy(), d["cloud_cover"].to_numpy() / 100.0)


@dataclass
class ThermalResult:
    time: pd.DatetimeIndex
    T_skin: np.ndarray  # K, hourly
    T_batt: np.ndarray  # K, hourly
    T_air: np.ndarray

    def celsius(self, which="batt"):
        return (self.T_batt if which == "batt" else self.T_skin) - 273.15


def simulate(coating: Coating, weather: Weather, q_gen_w: np.ndarray, enc: Enclosure = Enclosure(),
             substeps: int = 60) -> ThermalResult:
    """Explicit integration with 1-minute substeps (stable: skin time constant is ~20 min)."""
    n = len(weather.Ta)
    dt = 3600.0 / substeps
    eps_sky = sky_emissivity(weather.dew_c, weather.cloud)
    h = convection_coefficient(weather.wind)
    Ts = Tb = weather.Ta[0]
    out_s, out_b = np.empty(n), np.empty(n)
    a, e, F = coating.alpha_solar, coating.eps_ir, enc.sky_view
    for i in range(n):
        Ta, G, es, hi, q = weather.Ta[i], weather.ghi[i], eps_sky[i], h[i], q_gen_w[i]
        incoming = a * enc.sun_factor * G + e * (F * es + (1 - F)) * SIGMA * Ta**4
        for _ in range(substeps):
            q_int = enc.ua_internal * (Tb - Ts)
            dTs = (enc.area * (incoming - e * SIGMA * Ts**4 - hi * (Ts - Ta)) + q_int) / enc.skin_heat_capacity
            dTb = (q - q_int) / enc.battery_heat_capacity
            Ts += dTs * dt
            Tb += dTb * dt
        out_s[i], out_b[i] = Ts, Tb
    return ThermalResult(weather.time, out_s, out_b, weather.Ta)


def arrhenius_factor(T_kelvin, Ea_j_mol: float = 50e3, T_ref: float = 298.15):
    """Calendar-aging rate relative to 25 C. Literature Ea for Li-ion calendar aging is roughly 30-60 kJ/mol."""
    return np.exp(Ea_j_mol / 8.314462 * (1.0 / T_ref - 1.0 / np.asarray(T_kelvin)))


def derate_fraction(T_kelvin, start_c: float = 45.0, stop_c: float = 55.0):
    """Fraction of rated power available: 1 below start_c, linear to 0 at stop_c (a typical BMS policy)."""
    Tc = np.asarray(T_kelvin) - 273.15
    return np.clip((stop_c - Tc) / (stop_c - start_c), 0.0, 1.0)
