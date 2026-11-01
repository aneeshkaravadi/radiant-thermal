"""Dallas vs Houston: same battery, same skins, each city's own weather and ERCOT load-zone prices,
passive and with the ventilation fan, and which inputs explain the difference.

    python examples/compare_cities.py                     # both parts
    python examples/compare_cities.py --part cities       # or one at a time: cities, attribution

Each part rewrites only its own keys in docs/city_results.json.
"""
from __future__ import annotations

import argparse
import dataclasses
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from rt import coupled
from rt import dispatch as dp
from rt import enclosure as enc
from rt import radiative as rad

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "docs/city_results.json"
CITIES = {
    "Dallas (LZ_NORTH)": ("data/weather/dallas_2025_hourly.csv", "data/ercot/dam_2025_lz_north.csv"),
    "Houston (LZ_HOUSTON)": ("data/weather/houston_2025_hourly.csv", "data/ercot/dam_2025_lz_houston.csv"),
}
SKINS = [rad.DARK_PAINT, rad.WHITE_PAINT, rad.IDEAL_COOLER]


def load(city):
    weather_csv, price_csv = CITIES[city]
    w = enc.Weather.from_csv(ROOT / weather_csv)
    prices = dp.load_prices(ROOT / price_csv)
    assert (w.time == prices.index).all()
    return w, prices


def summer_mask(w):
    local = w.time.tz_convert("America/Chicago")
    return (local.month >= 6) & (local.month <= 8)


def box_numbers(r: coupled.CoupledResult, summer) -> dict:
    return {"summer_mean_batt_C": round(float(r.thermal.celsius()[summer].mean()), 2),
            "hours_derated": r.hours_derated,
            "aging_rate_vs_25C": round(float(enc.arrhenius_factor(r.thermal.T_batt).mean()), 3),
            "revenue_usd": round(r.revenue_usd, 1),
            "lost_usd": round(r.uncapped.revenue_usd - r.revenue_usd, 1)}


def run_city(city, battery=None):
    battery = dp.Battery() if battery is None else battery
    w, prices = load(city)
    summer = summer_mask(w)
    base = dp.dispatch(prices, battery)
    out = {"upper_bound_usd": round(base.revenue_usd, 1),
           "summer_mean_air_C": round(float(w.Ta[summer].mean() - 273.15), 2),
           "annual_mean_air_C": round(float(w.Ta.mean() - 273.15), 2),
           "summer_mean_dew_point_C": round(float(w.dew_c[summer].mean()), 2),
           "summer_mean_wind_10m_m_s": round(float(w.wind[summer].mean()), 2),
           "summer_mean_ghi_W_m2": round(float(w.ghi[summer].mean()), 1),
           "summer_mean_sky_emissivity": round(float(w.eps_sky[summer].mean()), 3),
           "skins": {}}
    fan = enc.Fan()
    for c in SKINS:
        out["skins"][c.name] = box_numbers(coupled.run(c, w, prices, battery), summer)
        rf = coupled.run(c, w, prices, battery, fan=fan)
        kwh = rf.thermal.fan_energy_kwh(fan)
        out["skins"][c.name]["with_fan"] = box_numbers(rf, summer) | {
            "fan_hours": round(float(rf.thermal.fan_on.sum())), "fan_kwh": round(float(kwh.sum()), 1),
            "fan_cost_usd": round(float(np.sum(prices.to_numpy() * kwh) / 1000), 2)}
    return out


def swap(w_to, w_from, what):
    """``w_to`` with one group of inputs taken from ``w_from`` (same hours of the same year)."""
    if what == "dew point":
        return dataclasses.replace(w_to, dew_c=w_from.dew_c)
    if what == "air temperature":  # keep the dew point at or below the new air temperature
        return dataclasses.replace(w_to, Ta=w_from.Ta, dew_c=np.minimum(w_to.dew_c, w_from.Ta - 273.15))
    if what == "wind":
        return dataclasses.replace(w_to, wind=w_from.wind)
    if what == "sunlight":
        return dataclasses.replace(w_to, ghi=w_from.ghi)
    if what == "clouds":
        return dataclasses.replace(w_to, cloud=w_from.cloud)
    raise ValueError(what)


def part_attribution(results):
    """Dark box. Start from one city and swap in one group of the other city's inputs at a time."""
    names = list(CITIES)
    (wd, pd_), (wh, ph) = load(names[0]), load(names[1])
    summer = summer_mask(wd)
    out = {"skin": rad.DARK_PAINT.name, "note": "each row: the starting city with one input group from the other city"}
    for label, (w0, p0, w1, p1) in {"Houston with Dallas's": (wh, ph, wd, pd_), "Dallas with Houston's": (wd, pd_, wh, ph)}.items():
        rows = {"nothing swapped": box_numbers(coupled.run(rad.DARK_PAINT, w0, p0), summer)}
        for what in ("dew point", "wind", "air temperature", "sunlight", "clouds"):
            rows[what] = box_numbers(coupled.run(rad.DARK_PAINT, swap(w0, w1, what), p0), summer)
        rows["prices"] = box_numbers(coupled.run(rad.DARK_PAINT, w0, p1), summer)
        out[label] = rows
    results["attribution_dark_box"] = out


def part_cities(results):
    for city in CITIES:
        results[city] = run_city(city)
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.8))
    x = np.arange(len(SKINS))
    for k, city in enumerate(CITIES):
        r = results[city]
        lost = [r["skins"][c.name]["lost_usd"] for c in SKINS]
        hrs = [r["skins"][c.name]["hours_derated"] for c in SKINS]
        fan_hrs = [r["skins"][c.name]["with_fan"]["fan_hours"] for c in SKINS]
        axes[0].bar(x + (k - 0.5) * 0.38, lost, 0.38, label=city)
        axes[1].bar(x + (k - 0.5) * 0.38, hrs, 0.38, label=city)
        axes[2].bar(x + (k - 0.5) * 0.38, fan_hrs, 0.38, label=city)
    for ax in axes:
        ax.set_xticks(x, [c.name.replace(" radiative", "\nradiative") for c in SKINS], fontsize=8)
        ax.grid(axis="y", alpha=0.3)
    axes[0].set_ylabel("revenue lost to derating, 2025 ($)")
    axes[1].set_ylabel("hours derated, 2025")
    axes[2].set_ylabel("hours the fan ran, 2025 (with a fan)")
    axes[0].legend(fontsize=8)
    fig.suptitle("Same battery, two cities, each with its own weather and ERCOT load-zone prices", fontsize=10)
    fig.tight_layout()
    fig.savefig(ROOT / "docs/figures/city_comparison.png", dpi=140)
    plt.close(fig)


PARTS = {"cities": part_cities, "attribution": part_attribution}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", choices=[*PARTS, "all"], default="all")
    args = ap.parse_args()
    results = json.loads(RESULTS.read_text()) if RESULTS.exists() else {}
    for name, fn in PARTS.items():
        if args.part in (name, "all"):
            fn(results)
            RESULTS.write_text(json.dumps(results, indent=2))
    print(json.dumps(results, indent=2))
