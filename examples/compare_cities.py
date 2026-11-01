"""Dallas vs Houston: same battery, same skins, each city's own weather and ERCOT load-zone prices,
passive and with the ventilation fan.

    python examples/compare_cities.py
"""
from __future__ import annotations

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
CITIES = {
    "Dallas (LZ_NORTH)": ("data/weather/dallas_2025_hourly.csv", "data/ercot/dam_2025_lz_north.csv"),
    "Houston (LZ_HOUSTON)": ("data/weather/houston_2025_hourly.csv", "data/ercot/dam_2025_lz_houston.csv"),
}
SKINS = [rad.DARK_PAINT, rad.WHITE_PAINT, rad.IDEAL_COOLER]


def run_city(weather_csv, price_csv, battery=None):
    battery = dp.Battery() if battery is None else battery
    w = enc.Weather.from_csv(ROOT / weather_csv)
    prices = dp.load_prices(ROOT / price_csv)
    assert (w.time == prices.index).all()
    local = w.time.tz_convert("America/Chicago")
    summer = (local.month >= 6) & (local.month <= 8)
    base = dp.dispatch(prices, battery)
    out = {"upper_bound_usd": round(base.revenue_usd, 1), "summer_mean_air_C": round(float(w.Ta[summer].mean() - 273.15), 2),
           "summer_mean_dew_point_C": round(float(w.dew_c[summer].mean()), 2), "skins": {}}
    fan = enc.Fan()
    for c in SKINS:
        r = coupled.run(c, w, prices, battery)
        out["skins"][c.name] = {"summer_mean_batt_C": round(float(r.thermal.celsius()[summer].mean()), 2),
                                "hours_derated": r.hours_derated,
                                "aging_rate_vs_25C": round(float(enc.arrhenius_factor(r.thermal.T_batt).mean()), 3),
                                "revenue_usd": round(r.revenue_usd, 1)}
        rf = coupled.run(c, w, prices, battery, fan=fan)
        kwh = rf.thermal.fan_energy_kwh(fan)
        out["skins"][c.name]["with_fan"] = {
            "summer_mean_batt_C": round(float(rf.thermal.celsius()[summer].mean()), 2),
            "hours_derated": rf.hours_derated,
            "aging_rate_vs_25C": round(float(enc.arrhenius_factor(rf.thermal.T_batt).mean()), 3),
            "revenue_usd": round(rf.revenue_usd, 1), "fan_hours": round(float(rf.thermal.fan_on.sum())),
            "fan_kwh": round(float(kwh.sum()), 1), "fan_cost_usd": round(float(np.sum(prices.to_numpy() * kwh) / 1000), 2)}
    return out


if __name__ == "__main__":
    results = {city: run_city(*files) for city, files in CITIES.items()}
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.8))
    x = np.arange(len(SKINS))
    for k, (city, r) in enumerate(results.items()):
        lost = [r["upper_bound_usd"] - r["skins"][c.name]["revenue_usd"] for c in SKINS]
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
    fig.suptitle("Same battery, two cities: how much heat costs depends on where you are", fontsize=10)
    fig.tight_layout()
    fig.savefig(ROOT / "docs/figures/city_comparison.png", dpi=140)
    (ROOT / "docs/city_results.json").write_text(json.dumps(results, indent=2))
    print(json.dumps(results, indent=2))
