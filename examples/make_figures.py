"""Regenerate every figure and number in the README.

    python examples/make_figures.py                  # everything (a few minutes)
    python examples/make_figures.py --part box       # or one part at a time: box, sweep,
                                                     # sens-enclosure, sens-physics, radiator

Each part rewrites only its own keys in docs/results.json, so the parts can run separately.
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
import pandas as pd

from rt import coupled
from rt import dispatch as dp
from rt import enclosure as enc
from rt import radiative as rad
from rt import radiator

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "docs" / "figures"
RESULTS = ROOT / "docs" / "results.json"
plt.rcParams.update({"figure.dpi": 140, "axes.grid": True, "grid.alpha": 0.3, "axes.spines.top": False,
                     "axes.spines.right": False, "font.size": 10})

w = enc.Weather.from_csv(ROOT / "data/weather/dallas_2025_hourly.csv")
prices = dp.load_prices(ROOT / "data/ercot/dam_2025_lz_north.csv")
assert (w.time == prices.index).all()
local = w.time.tz_convert("America/Chicago")
summer = (local.month >= 6) & (local.month <= 8)
battery = dp.Battery()
coatings = [rad.DARK_PAINT, rad.WHITE_PAINT, rad.IDEAL_COOLER]
DARK, WHITE = rad.DARK_PAINT.name, rad.WHITE_PAINT.name


def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIG / name)
    plt.close(fig)


def summary(r: coupled.CoupledResult) -> dict:
    Tb = r.thermal.celsius()
    return {"summer_mean_batt_C": round(float(Tb[summer].mean()), 2), "max_batt_C": round(float(Tb.max()), 1),
            "aging_rate_vs_25C": round(float(enc.arrhenius_factor(r.thermal.T_batt).mean()), 3),
            "hours_derated": r.hours_derated, "revenue_usd": round(r.revenue_usd, 1), "coupling_passes": r.passes}


def white_advantage(weather=w, enclosure=None, b=battery, inverter_heat="skin"):
    rev = {c.name: coupled.run(c, weather, prices, b, enclosure, inverter_heat=inverter_heat).revenue_usd
           for c in (rad.DARK_PAINT, rad.WHITE_PAINT)}
    return round(rev[WHITE] - rev[DARK], 1)


# ====================================================================== box: the passive and fan-cooled enclosure
def part_box(results):
    # ---------------------------------------------------------------- a hot day, surface by surface
    day = local.normalize() == local[np.argmax(w.Ta)].normalize()
    eps_sky = w.eps_sky
    h = enc.Enclosure().h(w.wind)  # the box's own convection, for the insulated-surface comparison
    hours = local[day].hour + 0.5  # each value is the mean over the hour
    fig, ax = plt.subplots(figsize=(7, 4))
    for c in coatings:
        dT = [rad.steady_surface_temperature(c, w.Ta[i], w.ghi[i], eps_sky[i], h[i]) - w.Ta[i] for i in np.where(day)[0]]
        ax.plot(hours, dT, marker="o", ms=3, label=f"{c.name} (a={c.alpha_solar:.2f}, e={c.eps_ir:.2f})")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("hour (Central time)")
    ax.set_ylabel("surface minus air temperature (K)")
    ax.set_title(f"Hottest day of 2025 in Dallas ({local[day][0]:%b %d}), insulated surface", fontsize=10)
    ax.legend(fontsize=7)
    save(fig, "hot_day_surface.png")

    # ---------------------------------------------------------------- a year in a passive enclosure
    upper = dp.dispatch(prices, battery)
    results["dispatch_unconstrained_usd"] = round(upper.revenue_usd, 1)
    results["battery"] = {"eta_one_way": battery.eta_one_way, "cell_round_trip": battery.cell_round_trip,
                          "eta_cell": round(battery.eta_cell, 4), "eta_inverter": round(battery.eta_inverter, 4),
                          "cell_share_of_loss": round(float(upper.cell_heat_w(battery).sum()
                                                            / (upper.cell_heat_w(battery) + upper.inverter_heat_w(battery)).sum()), 3)}
    runs = {c.name: coupled.run(c, w, prices, battery) for c in coatings}
    results["coatings"] = {c.name: {"alpha": round(c.alpha_solar, 3), "eps": c.eps_ir} | summary(runs[c.name]) for c in coatings}
    fig, ax = plt.subplots(figsize=(7, 4))
    bins = np.arange(25, 62, 1.0)
    for c in coatings:
        Tb = runs[c.name].thermal.celsius()
        ax.hist(Tb[summer], bins=bins, histtype="step", lw=1.6, label=f"{c.name}: mean {Tb[summer].mean():.1f} C")
    ax.axvline(45, color="k", ls=":", lw=1)
    ax.text(45.3, ax.get_ylim()[1] * 0.9, "derating\nstarts", fontsize=8)
    ax.set_xlabel("battery temperature, June-August hours (C)")
    ax.set_ylabel("hours")
    ax.set_title("Passive 20 kWh enclosure, Dallas summer 2025, ERCOT-dispatched", fontsize=10)
    ax.legend(fontsize=7)
    save(fig, "battery_summer_hist.png")

    dark, white = results["coatings"][DARK], results["coatings"][WHITE]
    results["white_vs_dark"] = {
        "revenue_usd": round(white["revenue_usd"] - dark["revenue_usd"], 1),
        "lost_usd": [round(results["dispatch_unconstrained_usd"] - dark["revenue_usd"], 1),
                     round(results["dispatch_unconstrained_usd"] - white["revenue_usd"], 1)],
        "summer_mean_drop_C": round(dark["summer_mean_batt_C"] - white["summer_mean_batt_C"], 2),
        "aging_rate_reduction_pct": round(100 * (1 - white["aging_rate_vs_25C"] / dark["aging_rate_vs_25C"]), 1),
        "derated_hours": [dark["hours_derated"], white["hours_derated"]],
    }
    # the aging rate with the activation energy I used first, unsourced (no effect on revenue)
    results["aging_rate_vs_25C_at_50kJ_mol"] = {
        c.name: round(float(enc.arrhenius_factor(runs[c.name].thermal.T_batt, 50e3).mean()), 3) for c in coatings}

    # ---------------------------------------------------------------- when the box is hot, and what power is worth then
    th_d, d_d = runs[DARK].thermal, runs[DARK].dispatch
    p = prices.to_numpy()
    sdays = pd.Series(np.arange(len(p)))[summer].groupby(local[summer].date)
    peak_T = [local[i[np.argmax(th_d.T_batt[i])]].hour for i in (g.to_numpy() for _, g in sdays)]
    peak_p = [local[i[np.argmax(p[i])]].hour for i in (g.to_numpy() for _, g in sdays)]
    derated = enc.derate_fraction(th_d.T_batt) < 1
    at_price_peak = [i[np.argmax(p[i])] for i in (g.to_numpy() for _, g in sdays)]
    results["timing_dark_box_summer"] = {
        "share_of_days_derated_in_the_peak_price_hour": round(float(np.mean(derated[at_price_peak])), 3),
        "median_hour_of_daily_peak_cell_temp": float(np.median(peak_T)),
        "median_hour_of_daily_peak_price": float(np.median(peak_p)),
        "corr_cell_temp_vs_price": round(float(np.corrcoef(th_d.T_batt[summer], p[summer])[0, 1]), 3),
        "mean_price_all_summer_hours": round(float(p[summer].mean()), 1),
        "mean_price_in_derated_hours": round(float(p[derated & summer].mean()), 1),
        "derated_hours_in_summer": int((derated & summer).sum()),
    }

    # ---------------------------------------------------------------- the same skins with a ventilation fan
    fan = enc.Fan()
    results["fan"] = {"flow_m3_s": fan.flow_m3_s, "power_w": fan.power_w, "on_off_C": [fan.on_c, fan.off_c], "skins": {}}
    fan_runs = {}
    for c in coatings:
        r = fan_runs[c.name] = coupled.run(c, w, prices, battery, fan=fan)
        kwh = r.thermal.fan_energy_kwh(fan)
        night = (local.hour >= 21) | (local.hour < 9)
        results["fan"]["skins"][c.name] = summary(r) | {
            "fan_hours": round(float(r.thermal.fan_on.sum())), "fan_kwh": round(float(kwh.sum()), 1),
            "fan_cost_usd_at_dam": round(float(np.sum(p * kwh) / 1000.0), 2),
            "summer_fan_share_9pm_to_9am": round(float(r.thermal.fan_on[summer & night].sum() / r.thermal.fan_on[summer].sum()), 3)}
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
    x = np.arange(len(coatings))
    names = [c.name for c in coatings]
    passive = [results["coatings"][n] for n in names]
    fanned = [results["fan"]["skins"][n] for n in names]
    top = results["dispatch_unconstrained_usd"]
    axes[0].bar(x - 0.2, [top - r["revenue_usd"] for r in passive], 0.4, label="passive")
    axes[0].bar(x + 0.2, [top - f["revenue_usd"] + f["fan_cost_usd_at_dam"] for f in fanned], 0.4,
                label="with fan (lost revenue + fan electricity)")
    axes[0].set_ylabel("$ lost in 2025")
    axes[1].bar(x - 0.2, [r["aging_rate_vs_25C"] for r in passive], 0.4, label="passive")
    axes[1].bar(x + 0.2, [f["aging_rate_vs_25C"] for f in fanned], 0.4, label="with fan")
    axes[1].set_ylabel("aging rate vs 25 C")
    axes[2].bar(x, [f["fan_hours"] for f in fanned], 0.5, color="C1")
    for k, f in enumerate(fanned):
        axes[2].text(k, f["fan_hours"] + 40, f"{f['fan_kwh']:.0f} kWh", ha="center", fontsize=8)
    axes[2].set_ylabel("hours the fan ran in 2025")
    for ax in axes:
        ax.set_xticks(x, [n.replace(" radiative", "\nradiative") for n in names], fontsize=8)
    axes[0].legend(fontsize=7)
    axes[1].legend(fontsize=7)
    fig.suptitle(f"A {fan.power_w:.0f} W, {fan.flow_m3_s * 2119:.0f} cfm fan (on at {fan.on_c:.0f} C, off at {fan.off_c:.0f} C)",
                 fontsize=10)
    save(fig, "fan_cooling.png")

    # ---------------------------------------------------------------- selective emitters, gray vs two-band
    sel_white = rad.selective_coating("selective white", rad.WHITE_PAINT.alpha_solar, 0.95, 0.10, "hypothetical")
    sel_ideal = rad.selective_coating("selective ideal cooler", rad.IDEAL_COOLER.alpha_solar, 0.95, 0.05, "hypothetical")

    def run_box(c, spectral):
        r = coupled.run(c, w, prices, battery, spectral=spectral)
        return {"summer_mean_batt_C": round(float(r.thermal.celsius()[summer].mean()), 2),
                "hours_derated": r.hours_derated, "revenue_usd": round(r.revenue_usd, 1)}

    results["spectral"] = {}
    pairs = [(rad.WHITE_PAINT, sel_white), (rad.IDEAL_COOLER, sel_ideal)]
    for broad, sel in pairs:
        results["spectral"][sel.name] = {"bands": list(sel.bands), "total_eps_at_300K": round(sel.eps_ir, 3),
                                         "broadband": dict(results["coatings"][broad.name]),
                                         "gray_model": run_box(sel, False), "two_band_model": run_box(sel, True)}
    surf = {c.name: [rad.steady_surface_temperature(c, w.Ta[i], w.ghi[i], eps_sky[i], h[i], spectral=True) - w.Ta[i]
                     for i in np.where(day)[0]] for c in (rad.WHITE_PAINT, sel_white)}
    results["spectral"]["hottest_day_surface_minus_air_K"] = {k: [round(v, 1) for v in dT] for k, dT in surf.items()}
    diff = np.array(surf[sel_white.name]) - np.array(surf[rad.WHITE_PAINT.name])
    results["spectral"]["hottest_day_selective_minus_white_K"] = {"coldest": round(float(diff.min()), 1),
                                                                  "coldest_at_hour_start": int(local[day][np.argmin(diff)].hour),
                                                                  "warmest": round(float(diff.max()), 1),
                                                                  "warmest_at_hour_start": int(local[day][np.argmax(diff)].hour)}
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.9))
    for name, dT in surf.items():
        axes[0].plot(hours, dT, marker="o", ms=3, label=name)
    axes[0].axhline(0, color="k", lw=0.8)
    axes[0].set_xlabel("hour (Central time)")
    axes[0].set_ylabel("insulated surface minus air (K)")
    axes[0].set_title(f"{local[day][0]:%b %d}, both with solar absorptance 0.25", fontsize=9)
    axes[0].legend(fontsize=8)
    labels, vals, colors = [], [], []
    short = {WHITE: "white\npaint", rad.IDEAL_COOLER.name: "ideal\ncooler"}
    for broad, sel in pairs:
        r = results["spectral"][sel.name]
        labels += [short[broad.name], "selective\n(gray\nmodel)", "selective\n(two-band)"]
        vals += [results["coatings"][broad.name]["hours_derated"], r["gray_model"]["hours_derated"],
                 r["two_band_model"]["hours_derated"]]
        colors += ["C0", "C7", "C1"]
    axes[1].bar(range(len(vals)), vals, color=colors)
    axes[1].set_xticks(range(len(vals)), labels, fontsize=7.5)
    axes[1].set_ylabel("battery hours derated, 2025")
    axes[1].set_title("battery box, Dallas 2025", fontsize=9)
    fig.suptitle("Selective emitters (0.95 inside 8-13 um, 0.05-0.10 outside) vs broadband skins, and the one-band (gray) "
                 "model of the same coatings", fontsize=10)
    save(fig, "selective_emitters.png")

    # ---------------------------------------------------------------- film design map
    # Both the films and the plain-white reference use the uncapped schedule's heat, so the map
    # compares skins under the same load.
    q_cell, q_skin = coupled.heat_inputs(upper, battery)
    ref = enc.simulate(rad.WHITE_PAINT, w, q_cell, substeps=20, q_skin_w=q_skin).celsius()[summer].mean()
    A_vals = np.linspace(0.0, 0.40, 9)
    R_vals = np.linspace(0.0, 0.40, 9)
    grid = np.zeros((len(R_vals), len(A_vals)))
    for i, Rf in enumerate(R_vals):
        for j, Af in enumerate(A_vals):
            if Rf + Af > 0.95:
                grid[i, j] = np.nan
                continue
            c = rad.Film("f", Rf, 1 - Rf - Af, 0.93).on(rad.WHITE_PAINT)
            grid[i, j] = enc.simulate(c, w, q_cell, substeps=20, q_skin_w=q_skin).celsius()[summer].mean() - ref
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    cs = ax.contourf(A_vals, R_vals, grid, levels=np.arange(-3, 6.5, 0.5), cmap="RdBu_r")
    ax.contour(A_vals, R_vals, grid, levels=[0], colors="k", linewidths=1.5)
    fig.colorbar(cs, label="summer mean battery temp vs plain white paint (K)")
    ax.set_xlabel("wrap / laminate solar absorptance A")
    ax.set_ylabel("wrap / laminate solar reflectance R")
    ax.set_title("Laminating a film or wrap over white paint: only ones that\nbarely absorb sunlight help (black line = break-even)",
                 fontsize=10)
    save(fig, "laminate_design_map.png")
    results["laminate_map"] = {"note": "black contour = no change vs plain white paint; laminate emissivity 0.93",
                               "clear_film_change_K": round(float(grid[0, 0]), 3),
                               "best_change_K": round(float(np.nanmin(grid)), 2), "worst_change_K": round(float(np.nanmax(grid)), 2)}

    # ---------------------------------------------------------------- one week of August
    wk = slice(int(np.argmax(local >= pd.Timestamp("2025-08-04", tz="America/Chicago"))),
               int(np.argmax(local >= pd.Timestamp("2025-08-11", tz="America/Chicago"))))
    fig, ax = plt.subplots(2, 1, figsize=(8, 5), sharex=True)
    ax[0].plot(local[wk], p[wk], color="C7")
    ax[0].set_ylabel("ERCOT North DAM ($/MWh)")
    ax[1].plot(local[wk], d_d.discharge_kw[wk] - d_d.charge_kw[wk], color="C0", label="dispatch (kW, + = export)")
    ax2 = ax[1].twinx()
    ax2.plot(local[wk], th_d.celsius()[wk], color="C3", label="cell temp, dark paint (C)")
    ax2.axhline(45, color="C3", ls=":", lw=1)
    ax[1].set_ylabel("battery power (kW)")
    ax2.set_ylabel("cell temperature (C)", color="C3")
    ax[0].set_title("First week of August 2025: prices, the dark box's dispatch and its cell temperature", fontsize=10)
    save(fig, "ercot_week.png")


# ====================================================================== sweep: what a tenth of solar absorptance is worth
def part_sweep(results):
    alphas = np.round(np.arange(0.05, 0.96, 0.10), 2)
    rows = [coupled.run(rad.Coating(f"alpha {a:.2f}", float(a), 0.90), w, prices, battery) for a in alphas]
    top = dp.dispatch(prices, battery).revenue_usd
    lost = top - np.array([r.revenue_usd for r in rows])
    hours = np.array([r.hours_derated for r in rows])
    per_tenth = np.diff(lost)  # per_tenth[k] is the step from alphas[k] to alphas[k + 1]
    k_white = int(np.argmin(abs(alphas - rad.WHITE_PAINT.alpha_solar)))
    results["absorptance_value"] = {
        "emissivity": 0.90, "alpha": alphas.tolist(), "revenue_lost_usd": np.round(lost, 1).tolist(),
        "hours_derated": hours.tolist(), "usd_per_tenth": np.round(per_tenth, 1).tolist(),
        "white_0.25_to_0.35": {"usd": round(float(per_tenth[k_white]), 1), "hours": int(hours[k_white + 1] - hours[k_white])}}
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
    ax.set_title(f"Each tenth of absorptance costs about \\${per_tenth[k_white]:.0f}/yr going darker from white paint "
                 f"and \\${per_tenth[-1]:.0f}/yr at the black end", fontsize=9)  # \\$ so it isn't read as math
    save(fig, "absorptance_value.png")


# ====================================================================== how much the answer leans on my assumptions
ENCLOSURE_ROWS = {"sun_factor": (0.4, 0.7, "share of the skin facing the sun"),
                  "ua_internal": (8.0, 20.0, "cells-to-skin conductance (W/K)"),
                  "battery_heat_capacity": (1.0e5, 4.0e5, "battery heat capacity (J/K)"),
                  "sky_view": (0.5, 1.0, "view of the sky"),
                  "area": (2.0, 4.0, "skin area (m^2)")}
PHYSICS_ROWS = {"terrain": ("country", "city", "terrain for the wind at the box"),
                "cell_round_trip": (0.978, 0.95**2, "cells' round-trip efficiency"),
                "inverter_heat": ("outside", None, "where the inverter's heat goes"),
                "convection": ("combined_10m", None, "convection coefficient"),
                "clouds": ("none", "black_total", "cloud model")}
BASE_LABEL = {"terrain": "suburbs", "cell_round_trip": battery.cell_round_trip, "inverter_heat": "skin",
              "convection": "2.8+3.0V at the box", "clouds": "layers"}
VALUE_LABEL = {"country": "open country", "city": "city", 0.978: "0.978 (best HTW)", 0.95**2: "0.9025 (all loss\nin cells, my first)",
               "outside": "outside the box", "combined_10m": "5.7+3.8V on the\n10 m wind (my first)",
               "none": "no clouds", "black_total": "all cloud black\n(my first)"}


def physics_case(field, value):
    if field == "terrain":
        return white_advantage(enclosure=enc.Enclosure(terrain=value))
    if field == "convection":
        return white_advantage(enclosure=enc.Enclosure(convection=value))
    if field == "cell_round_trip":
        return white_advantage(b=dp.Battery(cell_round_trip=value))
    if field == "inverter_heat":
        return white_advantage(inverter_heat=value)
    if field == "clouds":
        if value == "none":
            return white_advantage(weather=dataclasses.replace(w, cloud=np.zeros_like(w.cloud)))
        total = enc.instant_to_hour(pd.read_csv(ROOT / "data/weather/dallas_2025_hourly.csv")["cloud_cover"]) / 100.0
        black = np.column_stack([total, np.zeros_like(total), np.zeros_like(total)])
        return white_advantage(weather=dataclasses.replace(w, cloud=black, cloud_effect=(1.0, 0.0, 0.0)))
    raise ValueError(field)


def part_sensitivity(results, which):
    sens = results.setdefault("assumption_sensitivity", {})
    sens["white_minus_dark_usd_base"] = white_advantage()
    base_enc = enc.Enclosure()
    if which == "enclosure":
        for field, (lo, hi, label) in ENCLOSURE_ROWS.items():
            sens[field] = {"label": label, "base": getattr(base_enc, field),
                           "low": [lo, white_advantage(enclosure=dataclasses.replace(base_enc, **{field: lo}))],
                           "high": [hi, white_advantage(enclosure=dataclasses.replace(base_enc, **{field: hi}))]}
    else:
        for field, (lo, hi, label) in PHYSICS_ROWS.items():
            sens[field] = {"label": label, "base": BASE_LABEL[field], "low": [lo, physics_case(field, lo)]}
            if hi is not None:
                sens[field]["high"] = [hi, physics_case(field, hi)]
    rows = [k for k in list(ENCLOSURE_ROWS) + list(PHYSICS_ROWS) if k in sens]
    vals = [v for k in rows for v in (sens[k]["low"][1], sens[k].get("high", [None, sens[k]["low"][1]])[1])]
    sens["range_usd"] = [min(vals + [sens["white_minus_dark_usd_base"]]), max(vals + [sens["white_minus_dark_usd_base"]])]
    draw_sensitivity(sens, rows)


def draw_sensitivity(sens, rows):
    base = sens["white_minus_dark_usd_base"]

    def spread(k):
        r = sens[k]
        return abs(r.get("high", r["low"])[1] - r["low"][1]) + abs(r["low"][1] - base)

    order = sorted(rows, key=spread)
    fig, ax = plt.subplots(figsize=(9, 0.42 * len(order) + 1.6))
    for k, field in enumerate(order):
        r = sens[field]
        for end, color in (("low", "C0"), ("high", "C1")):
            if end not in r:
                continue
            value, adv = r[end]
            ax.barh(k, adv - base, left=base, color=color if "high" in r else "C7", height=0.5)
            name = VALUE_LABEL.get(value, f"{value:g}" if isinstance(value, float) else str(value))
            ax.text(adv + (1 if adv >= base else -1), k, name, ha="left" if adv >= base else "right", va="center", fontsize=7)
    ax.axvline(base, color="k", lw=1)
    def base_text(f):
        b = sens[f]["base"]
        return f"{b:g}" if isinstance(b, float) else str(b)
    ax.set_yticks(range(len(order)), [f"{sens[f]['label']} (base {base_text(f)})" for f in order], fontsize=7.5)
    ax.set_xlabel("white paint's advantage over dark gray, 2025 (\\$ per battery)")
    lo, hi = sens["range_usd"]
    ax.set_xlim(min(0, lo) - 30, hi + 30)
    ax.set_xticks(np.arange(0, hi + 30, 20))
    ax.set_title(f"One assumption at a time: white paint's advantage runs from \\${lo:.0f} to \\${hi:.0f}", fontsize=10)
    save(fig, "assumption_sensitivity.png")


# ====================================================================== spacecraft radiators
def part_radiator(results):
    T = np.linspace(250, 400, 151)
    space_coats = [radiator.OSR, radiator.WHITE_SPACE_PAINT, radiator.BLACK_PAINT]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, env in zip(axes, (radiator.LEO_SHADED, radiator.LEO_HOT), strict=True):
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

    results["radiator_environment"] = {"solar_flux_W_m2": radiator.HOT_SOLAR_FLUX, "albedo": radiator.HOT_ALBEDO,
                                       "earth_ir_W_m2": radiator.HOT_EARTH_IR}
    results["radiator_m2_per_kW_at_300K"] = {env.name: {c.name: area_or_note(float(radiator.area_per_kw(c, env, 300.0)))
                                                         for c in space_coats} for env in (radiator.LEO_SHADED, radiator.LEO_HOT)}
    area = float(radiator.area_per_kw(radiator.OSR, radiator.LEO_HOT, 300.0))  # m^2 to reject 1 kW at 300 K, sunlit hot case
    q_int = 1000.0 / area
    period = radiator.orbit_period(400e3)
    fig, ax = plt.subplots(figsize=(8, 3.8))
    results["radiator_orbit"] = {"area_m2_per_kW": round(area, 2), "period_min": round(period / 60, 1),
                                 "eclipse_min": round(radiator.eclipse_fraction(400e3) * period / 60, 1), "swing_K": {}}
    for cap, col in ((2000.0, "C3"), (5000.0, "C1"), (20000.0, "C0")):
        t, T, sun = radiator.orbit_transient(radiator.OSR, radiator.LEO_HOT, cap, q_int, orbits=6)
        last = t >= t[-1] - 3 * period
        ax.plot((t[last] - t[last][0]) / 60, T[last] - 273.15, color=col, label=f"{cap / 1000:g} kJ/m^2K")
        one = t >= t[-1] - period
        results["radiator_orbit"]["swing_K"][f"{cap / 1000:g}_kJ_m2K"] = [round(float(T[one].min()), 1), round(float(T[one].max()), 1)]
    f_ecl = radiator.eclipse_fraction(400e3)
    for k in range(3):  # the eclipse is the last part of each orbit
        ax.axvspan((k + 1 - f_ecl) * period / 60, (k + 1) * period / 60, color="k", alpha=0.08)
    ax.set_xlabel("time (min), three orbits at 400 km (gray: eclipse)")
    ax.set_ylabel("radiator temperature (C)")
    ax.set_title(f"Optical solar reflector radiator sized to reject 1 kW at 300 K in the hot case ({area:.1f} m^2):\n"
                 "how far it cools in Earth's shadow depends on its heat capacity", fontsize=9)
    ax.legend(fontsize=8, title="heat capacity per m^2", title_fontsize=8)
    save(fig, "radiator_orbit.png")


def redraw_sensitivity(results):
    sens = results["assumption_sensitivity"]
    draw_sensitivity(sens, [k for k in list(ENCLOSURE_ROWS) + list(PHYSICS_ROWS) if k in sens])


PARTS = {"box": part_box, "sweep": part_sweep, "sens-enclosure": lambda r: part_sensitivity(r, "enclosure"),
         "sens-physics": lambda r: part_sensitivity(r, "physics"), "radiator": part_radiator}
REDRAW = {"sens-figure": redraw_sensitivity}  # redraw from docs/results.json without rerunning

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", choices=[*PARTS, *REDRAW, "all"], default="all")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    results = json.loads(RESULTS.read_text()) if RESULTS.exists() else {}
    if args.part in REDRAW:
        REDRAW[args.part](results)
    for name, fn in PARTS.items():
        if args.part in (name, "all"):
            fn(results)
            RESULTS.write_text(json.dumps(results, indent=2))
    print(json.dumps(results, indent=2))
