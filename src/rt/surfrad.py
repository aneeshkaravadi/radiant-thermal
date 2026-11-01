"""Check the sky's infrared against measurement: NOAA SURFRAD pyrgeometers at two sites.

``data/surfrad/surfrad_2025_hourly.csv`` holds hourly downwelling longwave, air temperature and
relative humidity measured at Goodwin Creek, Mississippi (gwn, humid) and Desert Rock, Nevada
(dra, dry) in January and July 2025, next to Open-Meteo's weather and cloud layers for the same
hours (examples/build_surfrad_fixture.py). The model's sky radiates eps_sky sigma Ta^4.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from .radiative import CLOUD_EFFECT, SIGMA, sky_emissivity

FIXTURE = Path(__file__).resolve().parents[2] / "data/surfrad/surfrad_2025_hourly.csv"


def dew_point_c(temp_c, rh_pct):
    """Magnus formula with b = 17.625 and c = 243.04 C (Lawrence 2005, BAMS 86:225)."""
    b, c = 17.625, 243.04
    g = np.log(np.asarray(rh_pct, float) / 100.0) + b * np.asarray(temp_c, float) / (c + np.asarray(temp_c, float))
    return c * g / (b - g)


def load(path=FIXTURE) -> pd.DataFrame:
    d = pd.read_csv(path, parse_dates=["time_utc"])
    d["station_dew_c"] = dew_point_c(d["station_temp_c"], d["station_rh_pct"])
    d["month"] = d["time_utc"].dt.month
    for k in ("cloud", "cloud_low", "cloud_mid", "cloud_high"):
        d[k] = d[f"om_{k}_pct"] / 100.0
    return d


def modeled_longwave(d: pd.DataFrame, effect=CLOUD_EFFECT, inputs: str = "station", clouds: str = "layers") -> np.ndarray:
    """Downwelling longwave (W/m^2) the model would give for each hour.

    inputs="station" uses the station's own temperature and humidity (a check of the sky model alone);
    inputs="open-meteo" uses the reanalysis values, the way the enclosure model is driven.
    clouds="layers" is the current model; "black_total" is my first one (total cover as black cloud at
    air temperature); "none" ignores clouds.
    """
    if inputs == "station":
        T, dew = d["station_temp_c"].to_numpy() + 273.15, d["station_dew_c"].to_numpy()
    else:
        T, dew = d["om_temp_c"].to_numpy() + 273.15, d["om_dew_c"].to_numpy()
    if clouds == "layers":
        eps = sky_emissivity(dew, d["cloud_low"], d["cloud_mid"], d["cloud_high"], effect=effect)
    elif clouds == "black_total":
        eps = sky_emissivity(dew, d["cloud"], effect=(1.0, 0.0, 0.0))
    elif clouds == "none":
        eps = sky_emissivity(dew)
    else:
        raise ValueError(clouds)
    return np.asarray(eps) * SIGMA * T**4


def fit_cloud_effects(d: pd.DataFrame) -> np.ndarray:
    """Least-squares (k_low, k_mid, k_high) in [0, 1] against the measured longwave, station inputs."""
    def resid(k):
        return modeled_longwave(d, tuple(k)) - d["lw_down_wm2"].to_numpy()
    return least_squares(resid, [0.9, 0.5, 0.3], bounds=([0, 0, 0], [1, 1, 1])).x


def error_stats(d: pd.DataFrame, model: np.ndarray) -> dict:
    """Bias and RMS error (model minus measured, W/m^2) for each station and month, overall and by sky."""
    e = model - d["lw_down_wm2"].to_numpy()
    out = {}
    for (st, mo), g in d.groupby(["station", "month"]):
        cover = g["cloud"].to_numpy()
        ee = e[g.index]
        row = {}
        for name, sel in (("all", cover >= 0), ("clear", cover < 0.1), ("overcast", cover > 0.9)):
            row[name] = {"n": int(sel.sum()), "bias": round(float(ee[sel].mean()), 1) if sel.any() else None,
                         "rms": round(float(np.sqrt(np.mean(ee[sel] ** 2))), 1) if sel.any() else None}
        out[f"{st}_{mo:02d}"] = row
    return out
