"""Outdoor battery enclosure: a two-node thermal model driven by real weather.

    skin:     C_s dTs/dt = A [alpha f_sun G + eps (F L_sky + (1-F) sigma Ta^4) - eps sigma Ts^4 - h (Ts - Ta)]
                           + UA (Tb - Ts)
    battery:  C_b dTb/dt = Q_gen(t) - UA (Tb - Ts)

The coating only enters through alpha and eps on the skin. Everything else is
the same for each design, so differences in battery temperature come from the
coating alone. Parameters describe a generic ~20 kWh wall/ground unit and are
illustrative, not any company's product.

An optional ventilation fan pulls outside air past the cells when they get warm:

    battery:  ... - eff * rho cp Vdot * (Tb - Ta)     while the fan runs
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .radiative import (SIGMA, Coating, convection_coefficient, sky_band_emissivities, sky_emissivity, window_fraction,
                        window_fraction_fast)


@dataclass
class Enclosure:
    area: float = 3.0  # m^2 of exterior skin
    skin_heat_capacity: float = 1.5e4  # J/K  (~30 kg steel)
    battery_heat_capacity: float = 2.0e5  # J/K  (~200 kg of cells and modules)
    ua_internal: float = 12.0  # W/K, battery <-> skin conductance
    sky_view: float = 0.75  # box: top sees all sky, walls half
    sun_factor: float = 0.55  # mean projected-area fraction of the skin facing the sun, x GHI


@dataclass
class Fan:
    """Ventilation fan pulling outside air past the cells, switched by a thermostat on the cells.

    With the fan running, the cells lose  eff * rho cp Vdot * (Tb - Ta)  to the air stream,
    where eff is how close the exhaust air gets to the cell temperature. It turns on above
    ``on_c`` and off below ``off_c`` (the gap stops it chattering), and only while the
    outside air is cooler than the cells, since blowing in hotter air would heat them.
    It draws ``power_w`` while running.
    """

    flow_m3_s: float = 0.05  # about 100 cfm
    effectiveness: float = 0.5
    power_w: float = 30.0
    on_c: float = 35.0
    off_c: float = 32.0

    @property
    def ua(self) -> float:
        """Conductance from the cells to outside air while running, W/K (air at about 30 C)."""
        return self.effectiveness * 1.16 * 1007.0 * self.flow_m3_s


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
    fan_on: np.ndarray | None = None  # fraction of each hour the fan ran

    def fan_energy_kwh(self, fan: Fan) -> np.ndarray:
        """Fan electricity in each hour, kWh."""
        return np.zeros_like(self.T_air) if self.fan_on is None else self.fan_on * fan.power_w / 1000.0

    def celsius(self, which="batt"):
        return (self.T_batt if which == "batt" else self.T_skin) - 273.15


def simulate(coating: Coating, weather: Weather, q_gen_w: np.ndarray, enc: Enclosure = Enclosure(),
             substeps: int = 60, fan: Fan | None = None, spectral: bool = False) -> ThermalResult:
    """Explicit integration with 1-minute substeps (stable: skin time constant is ~20 min).

    ``spectral=True`` splits the skin's infrared exchange at the 8-13 um window
    (see radiative.net_heat_out_spectral). It changes nothing for a gray coating.
    """
    n = len(weather.Ta)
    dt = 3600.0 / substeps
    eps_sky = sky_emissivity(weather.dew_c, weather.cloud)
    h = convection_coefficient(weather.wind)
    Ts = Tb = weather.Ta[0]
    out_s, out_b, out_fan = np.empty(n), np.empty(n), np.zeros(n)
    a, e, F = coating.alpha_solar, coating.eps_ir, enc.sky_view
    if spectral:
        e_in, e_out = coating.eps_bands
        s_in, s_out = sky_band_emissivities(eps_sky, weather.Ta)
        fa = window_fraction(weather.Ta)
        ir_in = (e_in * (F * s_in + 1 - F) * fa + e_out * (F * s_out + 1 - F) * (1 - fa)) * SIGMA * weather.Ta**4
    running = False
    for i in range(n):
        Ta, G, es, hi, q = weather.Ta[i], weather.ghi[i], eps_sky[i], h[i], q_gen_w[i]
        if spectral:
            incoming = a * enc.sun_factor * G + ir_in[i]
        else:
            incoming = a * enc.sun_factor * G + e * (F * es + (1 - F)) * SIGMA * Ta**4
        on_steps = 0
        for _ in range(substeps):
            if fan is not None:
                Tc = Tb - 273.15
                running = (Tc > fan.off_c) if running else (Tc > fan.on_c)
                running = running and Ta < Tb
                on_steps += running
            q_int = enc.ua_internal * (Tb - Ts)
            q_fan = fan.ua * (Tb - Ta) if running else 0.0
            if spectral:
                fw = window_fraction_fast(Ts)
                emit = SIGMA * Ts**4 * (e_in * fw + e_out * (1 - fw))
            else:
                emit = e * SIGMA * Ts**4
            dTs = (enc.area * (incoming - emit - hi * (Ts - Ta)) + q_int) / enc.skin_heat_capacity
            dTb = (q - q_int - q_fan) / enc.battery_heat_capacity
            Ts += dTs * dt
            Tb += dTb * dt
        out_s[i], out_b[i], out_fan[i] = Ts, Tb, on_steps / substeps
    return ThermalResult(weather.time, out_s, out_b, weather.Ta, None if fan is None else out_fan)


def arrhenius_factor(T_kelvin, Ea_j_mol: float = 50e3, T_ref: float = 298.15):
    """Calendar-aging rate relative to 25 C. Literature Ea for Li-ion calendar aging is roughly 30-60 kJ/mol."""
    return np.exp(Ea_j_mol / 8.314462 * (1.0 / T_ref - 1.0 / np.asarray(T_kelvin)))


def derate_fraction(T_kelvin, start_c: float = 45.0, stop_c: float = 55.0):
    """Fraction of rated power available: 1 below start_c, linear to 0 at stop_c (a typical BMS policy)."""
    Tc = np.asarray(T_kelvin) - 273.15
    return np.clip((stop_c - Tc) / (stop_c - start_c), 0.0, 1.0)
