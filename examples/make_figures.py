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
from rt import radiative as rad
from rt import radiator

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

# ------------------------------------------------------------------ 2b. the same skins with a ventilation fan
fan = enc.Fan()
results["fan"] = {"flow_m3_s": fan.flow_m3_s, "power_w": fan.power_w, "on_off_C": [fan.on_c, fan.off_c], "skins": {}}
for c in coatings:
    th = enc.simulate(c, w, q, fan=fan)
    d2 = dp.dispatch(prices, battery, power_cap_kw=enc.derate_fraction(th.T_batt) * battery.power_kw)
    th2 = enc.simulate(c, w, d2.heat_w(battery.eta_one_way), fan=fan)
    fan_kwh = th2.fan_energy_kwh(fan)
    results["fan"]["skins"][c.name] = {
        "summer_mean_batt_C": round(float(th2.celsius()[summer].mean()), 2),
        "hours_derated": int((enc.derate_fraction(th2.T_batt) < 1).sum()),
        "aging_rate_vs_25C": round(float(enc.arrhenius_factor(th2.T_batt).mean()), 3),
        "revenue_usd": round(d2.revenue_usd, 1),
        "fan_hours": round(float(th2.fan_on.sum())),
        "fan_kwh": round(float(fan_kwh.sum()), 1),
        "fan_cost_usd_at_dam": round(float(np.sum(prices.to_numpy() * fan_kwh) / 1000.0), 2),
    }
fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
x = np.arange(len(coatings))
names = [c.name for c in coatings]
passive = [results["coatings"][n] for n in names]
fanned = [results["fan"]["skins"][n] for n in names]
upper = results["dispatch_unconstrained_usd"]
axes[0].bar(x - 0.2, [upper - r["revenue_usd"] for r in passive], 0.4, label="passive")
axes[0].bar(x + 0.2, [upper - f["revenue_usd"] + f["fan_cost_usd_at_dam"] for f in fanned], 0.4,
            label="with fan (lost revenue + fan electricity)")
axes[0].set_ylabel("$ lost in 2025")
axes[1].bar(x - 0.2, [r["aging_rate_vs_25C"] for r in passive], 0.4, label="passive")
axes[1].bar(x + 0.2, [f["aging_rate_vs_25C"] for f in fanned], 0.4, label="with fan")
axes[1].set_ylabel("calendar aging rate vs 25 C")
axes[2].bar(x, [f["fan_hours"] for f in fanned], 0.5, color="C1")
for k, f in enumerate(fanned):
    axes[2].text(k, f["fan_hours"] + 40, f"{f['fan_kwh']:.0f} kWh", ha="center", fontsize=8)
axes[2].set_ylabel("hours the fan ran in 2025")
for ax in axes:
    ax.set_xticks(x, [n.replace(" radiative", "\nradiative") for n in names], fontsize=8)
axes[0].legend(fontsize=7)
axes[1].legend(fontsize=7)
fig.suptitle(f"A {fan.power_w:.0f} W, {fan.flow_m3_s * 2119:.0f} cfm fan (on at {fan.on_c:.0f} C, off at {fan.off_c:.0f} C): "
             "the paint matters much less once it's running", fontsize=10)
save(fig, "fan_cooling.png")

# ------------------------------------------------------------------ 2c. selective emitters, gray vs two-band
sel_white = rad.selective_coating("selective white", rad.WHITE_PAINT.alpha_solar, 0.95, 0.10, "hypothetical")
sel_ideal = rad.selective_coating("selective ideal cooler", rad.IDEAL_COOLER.alpha_solar, 0.95, 0.05, "hypothetical")


def run_box(c, spectral):
    th = enc.simulate(c, w, q, spectral=spectral)
    d2 = dp.dispatch(prices, battery, power_cap_kw=enc.derate_fraction(th.T_batt) * battery.power_kw)
    th2 = enc.simulate(c, w, d2.heat_w(battery.eta_one_way), spectral=spectral)
    return {"summer_mean_batt_C": round(float(th2.celsius()[summer].mean()), 2),
            "hours_derated": int((enc.derate_fraction(th2.T_batt) < 1).sum()), "revenue_usd": round(d2.revenue_usd, 1)}


results["spectral"] = {}
pairs = [(rad.WHITE_PAINT, sel_white), (rad.IDEAL_COOLER, sel_ideal)]
for broad, sel in pairs:
    results["spectral"][sel.name] = {"bands": list(sel.bands), "total_eps_at_300K": round(sel.eps_ir, 3),
                                     "broadband": results["coatings"][broad.name] | {},
                                     "gray_model": run_box(sel, False), "two_band_model": run_box(sel, True)}
day_idx = np.where(day)[0]
surf = {c.name: [rad.steady_surface_temperature(c, w.Ta[i], w.ghi[i], eps_sky[i], h[i], spectral=True) - w.Ta[i]
                 for i in day_idx] for c in (rad.WHITE_PAINT, sel_white)}
fig, axes = plt.subplots(1, 2, figsize=(11, 3.9))
for name, dT in surf.items():
    axes[0].plot(local[day].hour, dT, marker="o", ms=3, label=name)
axes[0].axhline(0, color="k", lw=0.8)
axes[0].set_xlabel("hour (Central time)")
axes[0].set_ylabel("insulated surface minus air (K)")
axes[0].set_title(f"{local[day][0]:%b %d}, both with solar absorptance 0.25", fontsize=9)
axes[0].legend(fontsize=8)
labels, vals, colors = [], [], []
short = {rad.WHITE_PAINT.name: "white\npaint", rad.IDEAL_COOLER.name: "ideal\ncooler"}
for broad, sel in pairs:
    r = results["spectral"][sel.name]
    labels += [short[broad.name], "selective\n(gray\nmodel)", "selective\n(two-band)"]
    vals += [results["coatings"][broad.name]["hours_derated"], r["gray_model"]["hours_derated"],
             r["two_band_model"]["hours_derated"]]
    colors += ["C0", "C7", "C1"]
axes[1].bar(range(len(vals)), vals, color=colors)
axes[1].set_xticks(range(len(vals)), labels, fontsize=7.5)
axes[1].set_ylabel("battery hours derated, 2025")
axes[1].set_title("a box that runs warm wants to emit in every band", fontsize=9)
fig.suptitle("Selective emitters (0.95 inside 8-13 um, 0.05-0.10 outside): colder at night, worse for a warm box,\n"
             "and the one-band (gray) model badly underrates them", fontsize=10)
save(fig, "selective_emitters.png")
results["spectral"]["hottest_day_surface_minus_air_K"] = {k: [round(v, 1) for v in dT] for k, dT in surf.items()}

# ------------------------------------------------------------------ 2d. what a tenth of solar absorptance is worth
alphas = np.round(np.arange(0.05, 0.96, 0.10), 2)
sweep = []
for a in alphas:
    c = rad.Coating(f"alpha {a:.2f}", float(a), 0.90)
    th = enc.simulate(c, w, q)
    d2 = dp.dispatch(prices, battery, power_cap_kw=enc.derate_fraction(th.T_batt) * battery.power_kw)
    th2 = enc.simulate(c, w, d2.heat_w(battery.eta_one_way))
    sweep.append((d2.revenue_usd, int((enc.derate_fraction(th2.T_batt) < 1).sum())))
lost = results["dispatch_unconstrained_usd"] - np.array([s_[0] for s_ in sweep])
hours = np.array([s_[1] for s_ in sweep])
per_tenth = np.diff(lost)  # alphas are 0.1 apart
results["absorptance_value"] = {"emissivity": 0.90, "alpha": alphas.tolist(), "revenue_lost_usd": np.round(lost, 1).tolist(),
                                "hours_derated": hours.tolist(), "usd_per_tenth": np.round(per_tenth, 1).tolist()}
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(alphas, lost, "o-", color="C3", label="revenue lost to derating")
ax.set_xlabel("skin solar absorptance (infrared emissivity 0.90)")
ax.set_ylabel("\\$ lost in 2025", color="C3")
ax2 = ax.twinx()
ax2.plot(alphas, hours, "s--", color="C0", label="hours derated")
ax2.set_ylabel("hours derated in 2025", color="C0")
for c in (rad.WHITE_PAINT, rad.DARK_PAINT):
    ax.axvline(c.alpha_solar, color="k", ls=":", lw=1)
    ax.text(c.alpha_solar + 0.01, lost.max() * 0.95, c.name, fontsize=8)
ax.set_title(f"Each tenth of absorptance costs \\${per_tenth[0]:.0f}/yr near white and \\${per_tenth[-1]:.0f}/yr near black",
             fontsize=10)  # \\$ so matplotlib doesn't read a pair of dollar signs as math
save(fig, "absorptance_value.png")

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


def area_or_note(v):
    return round(v, 2) if np.isfinite(v) else "cannot reject heat"


results["radiator_m2_per_kW_at_300K"] = {env.name: {c.name: area_or_note(float(radiator.area_per_kw(c, env, 300.0)))
                                                     for c in space_coats} for env in (radiator.LEO_SHADED, radiator.LEO_HOT)}

(ROOT / "docs" / "results.json").write_text(json.dumps(results, indent=2))
print(json.dumps(results, indent=2))
