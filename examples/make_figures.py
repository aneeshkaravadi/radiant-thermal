"""Regenerate every figure and number in the README.   python examples/make_figures.py"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rt import dispatch as dp
from rt import enclosure as enc
from rt import radiative as rad, radiator

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "docs" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"figure.dpi": 140, "axes.grid": True, "grid.alpha": 0.3, "axes.spines.top": False,
                     "axes.spines.right": False, "font.size": 10})
results: dict = {}


def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIG / name)
    plt.close(fig)


w = enc.Weather.from_csv(ROOT / "data/weather/dallas_2025_hourly.csv")
prices = dp.load_prices(ROOT / "data/ercot/dam_2025_lz_north.csv")
assert (w.time == prices.index).all()
tag = ""
summer = (w.time.tz_convert("America/Chicago").month >= 6) & (w.time.tz_convert("America/Chicago").month <= 8)

coatings = [rad.DARK_PAINT, rad.WHITE_PAINT, rad.IDEAL_COOLER]

# ------------------------------------------------------------------ 1. a hot day, surface by surface
local = w.time.tz_convert("America/Chicago")
day = local.normalize() == local[np.argmax(w.Ta)].normalize()
eps_sky = rad.sky_emissivity(w.dew_c, w.cloud)
h = rad.convection_coefficient(w.wind)
fig, ax = plt.subplots(figsize=(7, 4))
for c in coatings:
    dT = [rad.steady_surface_temperature(c, w.Ta[i], w.ghi[i], eps_sky[i], h[i]) - w.Ta[i] for i in np.where(day)[0]]
    ax.plot(local[day].hour, dT, marker="o", ms=3, label=f"{c.name} (a={c.alpha_solar:.2f}, e={c.eps_ir:.2f})")
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("hour (Central time)")
ax.set_ylabel("surface minus air temperature (K)")
ax.set_title(f"Hottest day of 2025 in Dallas ({local[day][0]:%b %d}), insulated surface{tag}", fontsize=10)
ax.legend(fontsize=7)
save(fig, "hot_day_surface.png")

# ------------------------------------------------------------------ 2. a year in a battery enclosure
battery = dp.Battery()
base_dispatch = dp.dispatch(prices, battery)
q = base_dispatch.heat_w(battery.eta_one_way)
results["dispatch_unconstrained_usd"] = round(base_dispatch.revenue_usd, 1)
results["coatings"] = {}
fig, ax = plt.subplots(figsize=(7, 4))
bins = np.arange(25, 62, 1.0)
for c in coatings:
    th = enc.simulate(c, w, q)
    # Feed thermal derating back into dispatch once: power caps from the first thermal pass.
    cap = enc.derate_fraction(th.T_batt) * battery.power_kw
    d2 = dp.dispatch(prices, battery, power_cap_kw=cap)
    th2 = enc.simulate(c, w, d2.heat_w(battery.eta_one_way))
    Tb = th2.celsius()
    results["coatings"][c.name] = {
        "alpha": round(c.alpha_solar, 3), "eps": c.eps_ir,
        "summer_mean_batt_C": round(float(Tb[summer].mean()), 2),
        "max_batt_C": round(float(Tb.max()), 1),
        "aging_rate_vs_25C": round(float(enc.arrhenius_factor(th2.T_batt).mean()), 3),
        "hours_derated": int((enc.derate_fraction(th2.T_batt) < 1).sum()),
        "revenue_usd": round(d2.revenue_usd, 1),
    }
    ax.hist(Tb[summer], bins=bins, histtype="step", lw=1.6, label=f"{c.name}: mean {Tb[summer].mean():.1f} C")
ax.axvline(45, color="k", ls=":", lw=1)
ax.text(45.3, ax.get_ylim()[1] * 0.9, "derating\nstarts", fontsize=8)
ax.set_xlabel("battery temperature, June-August hours (C)")
ax.set_ylabel("hours")
ax.set_title(f"Passive 20 kWh enclosure, Dallas summer 2025, ERCOT-dispatched{tag}", fontsize=10)
ax.legend(fontsize=7)
save(fig, "battery_summer_hist.png")

dark, white = results["coatings"][rad.DARK_PAINT.name], results["coatings"][rad.WHITE_PAINT.name]
results["white_vs_dark"] = {
    "summer_mean_drop_C": round(dark["summer_mean_batt_C"] - white["summer_mean_batt_C"], 2),
    "aging_rate_reduction_pct": round(100 * (1 - white["aging_rate_vs_25C"] / dark["aging_rate_vs_25C"]), 1),
    "derated_hours": [dark["hours_derated"], white["hours_derated"]],
}

# ------------------------------------------------------------------ 3. film design map
A_vals = np.linspace(0.0, 0.40, 9)
R_vals = np.linspace(0.0, 0.40, 9)
ref = results["coatings"][rad.WHITE_PAINT.name]["summer_mean_batt_C"]
grid = np.zeros((len(R_vals), len(A_vals)))
for i, Rf in enumerate(R_vals):
    for j, Af in enumerate(A_vals):
        if Rf + Af > 0.95:
            grid[i, j] = np.nan
            continue
        c = rad.Film("f", Rf, 1 - Rf - Af, 0.93).on(rad.WHITE_PAINT)
        grid[i, j] = enc.simulate(c, w, q, substeps=20).celsius()[summer].mean() - ref
fig, ax = plt.subplots(figsize=(6.5, 4.8))
cs = ax.contourf(A_vals, R_vals, grid, levels=np.arange(-3, 6.5, 0.5), cmap="RdBu_r")
ax.contour(A_vals, R_vals, grid, levels=[0], colors="k", linewidths=1.5)
fig.colorbar(cs, label="summer mean battery temp vs plain white paint (K)")
ax.set_xlabel("wrap / laminate solar absorptance A")
ax.set_ylabel("wrap / laminate solar reflectance R")
ax.set_title("Laminating a film or wrap over white paint: only ones that\nbarely absorb sunlight help (black line = break-even)", fontsize=10)
save(fig, "laminate_design_map.png")
results["laminate_map_note"] = "black contour = no change vs plain white paint; laminate emissivity 0.93"

# ------------------------------------------------------------------ 4. ERCOT week with derating
wk = slice(int(np.argmax(local >= pd.Timestamp("2025-08-04", tz="America/Chicago"))), int(np.argmax(local >= pd.Timestamp("2025-08-11", tz="America/Chicago"))))
th_dark = enc.simulate(rad.DARK_PAINT, w, q)
fig, ax = plt.subplots(2, 1, figsize=(8, 5), sharex=True)
ax[0].plot(local[wk], prices.to_numpy()[wk], color="C7")
ax[0].set_ylabel("ERCOT North DAM ($/MWh)")
ax[1].plot(local[wk], base_dispatch.discharge_kw[wk] - base_dispatch.charge_kw[wk], color="C0", label="dispatch (kW, + = export)")
ax2 = ax[1].twinx()
ax2.plot(local[wk], th_dark.celsius()[wk], color="C3", label="cell temp, dark paint (C)")
ax2.axhline(45, color="C3", ls=":", lw=1)
ax[1].set_ylabel("battery power (kW)")
ax2.set_ylabel("cell temperature (C)", color="C3")
ax[0].set_title("First week of August 2025: the battery is hottest exactly when prices peak")
save(fig, "ercot_week.png")

# ------------------------------------------------------------------ 5. spacecraft radiators
T = np.linspace(250, 400, 151)
space_coats = [radiator.OSR, radiator.WHITE_SPACE_PAINT, radiator.BLACK_PAINT]
fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
for ax, env in zip(axes, (radiator.LEO_SHADED, radiator.LEO_HOT)):
    for c in space_coats:
        ax.plot(T - 273.15, radiator.area_per_kw(c, env, T), label=f"{c.name} (a={c.alpha_solar}, e={c.eps_ir})")
    ax.set_title(env.name, fontsize=9)
    ax.set_xlabel("radiator temperature (C)")
    ax.set_ylim(0, 12)
axes[0].set_ylabel("radiator area per kW rejected (m^2)")
axes[0].legend(fontsize=7)
fig.suptitle("Radiator area depends on the same two numbers as a rooftop cooler: alpha and epsilon")
save(fig, "radiator_area.png")
results["radiator_m2_per_kW_at_300K"] = {env.name: {c.name: (lambda v: round(v, 2) if np.isfinite(v) else "cannot reject heat")(float(radiator.area_per_kw(c, env, 300.0)))
                                                     for c in space_coats} for env in (radiator.LEO_SHADED, radiator.LEO_HOT)}

(ROOT / "docs" / "results.json").write_text(json.dumps(results, indent=2))
print(json.dumps(results, indent=2))
