"""Outdoor battery enclosure: a two-node thermal model driven by real weather.

    skin:     C_s dTs/dt = A [alpha f_sun G + eps (F L_sky + (1-F) sigma Ta^4) - eps sigma Ts^4 - h (Ts - Ta)]
                           + UA (Tb - Ts) + Q_inv(t)
    battery:  C_b dTb/dt = Q_cell(t) - UA (Tb - Ts)

Q_cell is the heat released in the cells and Q_inv the inverter's, which I take as heat-sunk to the skin.

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

from .radiative import (
    CLOUD_EFFECT,
    SIGMA,
    Coating,
    combined_coefficient,
    convection_coefficient,
    local_wind,
    sky_band_emissivities,
    sky_emissivity,
    window_fraction,
    window_fraction_fast,
)


@dataclass
class Enclosure:
    area: float = 3.0  # m^2 of exterior skin
    skin_heat_capacity: float = 1.5e4  # J/K  (~30 kg steel)
    battery_heat_capacity: float = 2.0e5  # J/K  (~200 kg of cells and modules)
    ua_internal: float = 12.0  # W/K, battery <-> skin conductance
    sky_view: float = 0.75  # box: top sees all sky, walls half
    sun_factor: float = 0.55  # mean projected-area fraction of the skin facing the sun, x GHI
    wind_height_m: float = 1.0  # middle of a ground-standing box
    terrain: str = "suburbs"
    convection: str = "convection_only"  # or "combined_10m", my first setup (a combined coefficient on the 10 m wind)

    def h(self, wind_10m) -> np.ndarray:
        """Outside convection coefficient (W/m^2/K) from the 10 m wind."""
        if self.convection == "combined_10m":
            return combined_coefficient(wind_10m)
        if self.convection != "convection_only":
            raise ValueError(self.convection)
        return convection_coefficient(local_wind(wind_10m, self.wind_height_m, self.terrain))


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


def instant_to_hour(x) -> np.ndarray:
    """Average of instantaneous values at each hour's start and end (the last hour keeps its last value)."""
    x = np.asarray(x, float)
    return np.r_[0.5 * (x[:-1] + x[1:]), x[-1]]


@dataclass
class Weather:
    """Weather for each simulated hour. Row i holds the conditions over [time[i], time[i] + 1 h).

    ``ghi`` is the mean sunlight over that hour and the other fields are hour averages too.
    """

    time: pd.DatetimeIndex  # hour starts
    Ta: np.ndarray  # K
    dew_c: np.ndarray
    ghi: np.ndarray  # W/m^2
    wind: np.ndarray  # m/s at 10 m
    cloud: np.ndarray  # (n, 3): low, mid and high cloud cover, 0..1
    cloud_effect: tuple[float, float, float] = CLOUD_EFFECT  # see radiative.sky_emissivity

    @property
    def eps_sky(self) -> np.ndarray:
        return sky_emissivity(self.dew_c, *np.asarray(self.cloud, float).T, effect=self.cloud_effect)

    @classmethod
    def from_csv(cls, path) -> Weather:
        """Read an Open-Meteo hourly file and put every variable on the model's hour-start convention.

        Open-Meteo reports shortwave_radiation at time t as the mean over the hour before t, and
        temperature, dew point, wind and cloud cover as instantaneous values at t. The model steps
        the hour starting at t_i, so it needs the sunlight stamped t_(i+1) and the midpoint of the
        instantaneous values at t_i and t_(i+1). The files end at local midnight, so the last hour
        gets no sunlight and keeps its last instantaneous values.
        """
        d = pd.read_csv(path, parse_dates=["time"])
        midpoint = instant_to_hour
        ghi = np.r_[d["shortwave_radiation"].to_numpy(float)[1:], 0.0]
        cloud = np.column_stack([midpoint(d[f"cloud_cover_{k}"]) / 100.0 for k in ("low", "mid", "high")])
        return cls(pd.DatetimeIndex(d["time"]), midpoint(d["temperature_2m"]) + 273.15, midpoint(d["dew_point_2m"]),
                   ghi, midpoint(d["wind_speed_10m"]), cloud)


@dataclass
class ThermalResult:
    time: pd.DatetimeIndex
    T_skin: np.ndarray  # K, mean over each hour
    T_batt: np.ndarray  # K, mean over each hour
    T_air: np.ndarray
    fan_on: np.ndarray | None = None  # fraction of each hour the fan ran

    def fan_energy_kwh(self, fan: Fan) -> np.ndarray:
        """Fan electricity in each hour, kWh."""
        return np.zeros_like(self.T_air) if self.fan_on is None else self.fan_on * fan.power_w / 1000.0

    def celsius(self, which="batt"):
        return (self.T_batt if which == "batt" else self.T_skin) - 273.15


def simulate(coating: Coating, weather: Weather, q_gen_w: np.ndarray, enc: Enclosure | None = None,
             substeps: int = 60, fan: Fan | None = None, spectral: bool = False,
             q_skin_w: np.ndarray | None = None) -> ThermalResult:
    """Explicit integration with 1-minute substeps.

    The skin's time constant C_s / (A (h + 4 eps sigma T^3) + UA) is 3 to 7 minutes over 2025's
    winds, so Euler steps of 1 minute are well inside the stability limit of twice that.

    ``q_gen_w`` is the heat released in the cells each hour and ``q_skin_w`` (optional) heat put
    straight into the skin, such as an inverter's losses.
    ``spectral=True`` splits the skin's infrared exchange at the 8-13 um window
    (see radiative.net_heat_out_spectral). It changes nothing for a gray coating.
    """
    enc = Enclosure() if enc is None else enc
    n = len(weather.Ta)
    dt = 3600.0 / substeps
    eps_sky = weather.eps_sky
    h = enc.h(weather.wind)
    Ts = Tb = weather.Ta[0]
    out_s, out_b, out_fan = np.empty(n), np.empty(n), np.zeros(n)
    a, e, F = coating.alpha_solar, coating.eps_ir, enc.sky_view
    if spectral:
        e_in, e_out = coating.eps_bands
        s_in, s_out = sky_band_emissivities(eps_sky, weather.Ta)
        fa = window_fraction(weather.Ta)
        ir_in = (e_in * (F * s_in + 1 - F) * fa + e_out * (F * s_out + 1 - F) * (1 - fa)) * SIGMA * weather.Ta**4
    q_skin = np.zeros(n) if q_skin_w is None else np.asarray(q_skin_w, float)
    running = False
    for i in range(n):
        Ta, G, es, hi, q, qs = weather.Ta[i], weather.ghi[i], eps_sky[i], h[i], q_gen_w[i], q_skin[i]
        if spectral:
            incoming = a * enc.sun_factor * G + ir_in[i]
        else:
            incoming = a * enc.sun_factor * G + e * (F * es + (1 - F)) * SIGMA * Ta**4
        on_steps, sum_s, sum_b = 0, 0.0, 0.0
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
            dTs = (enc.area * (incoming - emit - hi * (Ts - Ta)) + q_int + qs) / enc.skin_heat_capacity
            dTb = (q - q_int - q_fan) / enc.battery_heat_capacity
            Ts += dTs * dt
            Tb += dTb * dt
            sum_s += Ts
            sum_b += Tb
        out_s[i], out_b[i], out_fan[i] = sum_s / substeps, sum_b / substeps, on_steps / substeps
    return ThermalResult(weather.time, out_s, out_b, weather.Ta, None if fan is None else out_fan)


# Activation energy for capacity fade of large-format stationary-storage LFP cells, 37.3 kJ/mol over all of
# their calendar and cycling protocols at 35 and 50 C (Yagci, Richter, Behmann and Bessler, J. Energy
# Storage 124, 116774, 2025).
EA_LFP_J_MOL = 37.3e3


def arrhenius_factor(T_kelvin, Ea_j_mol: float = EA_LFP_J_MOL, T_ref: float = 298.15):
    """Aging rate relative to 25 C, exp[Ea/R (1/T_ref - 1/T)]."""
    return np.exp(Ea_j_mol / 8.314462 * (1.0 / T_ref - 1.0 / np.asarray(T_kelvin)))


def derate_fraction(T_kelvin, start_c: float = 45.0, stop_c: float = 55.0):
    """Fraction of rated power available: 1 below start_c, linear to 0 at stop_c.

    This cell-temperature policy is my assumption, not a datasheet's. For scale, the Powerwall 3
    datasheet says performance "may be de-rated at operating temperatures above 40 C", meaning air temperature.
    """
    Tc = np.asarray(T_kelvin) - 273.15
    return np.clip((stop_c - Tc) / (stop_c - start_c), 0.0, 1.0)
