"""Compare the sky's infrared with NOAA SURFRAD measurements.   python examples/validate_sky.py

Writes docs/sky_validation.json and docs/figures/sky_validation.png.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from rt import radiative as rad
from rt import surfrad

ROOT = Path(__file__).resolve().parents[1]
d = surfrad.load()
jan = d[d["month"] == 1].reset_index(drop=True)
models = {"clear sky only": "none", "first model: total cloud, black": "black_total", "cloud layers (fit on January)": "layers"}
out = {"cloud_effect_low_mid_high": list(rad.CLOUD_EFFECT),
       "refit_on_january": [round(float(k), 3) for k in surfrad.fit_cloud_effects(jan)],
       "refit_on_all_hours": [round(float(k), 3) for k in surfrad.fit_cloud_effects(d)],
       "hours": {f"{st}_{mo:02d}": int(len(g)) for (st, mo), g in d.groupby(["station", "month"])},
       "station_inputs": {}, "open_meteo_inputs": {}}
for label, clouds in models.items():
    out["station_inputs"][label] = surfrad.error_stats(d, surfrad.modeled_longwave(d, clouds=clouds))
    out["open_meteo_inputs"][label] = surfrad.error_stats(d, surfrad.modeled_longwave(d, inputs="open-meteo", clouds=clouds))
(ROOT / "docs/sky_validation.json").write_text(json.dumps(out, indent=2))

fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
keys = ["gwn_01", "gwn_07", "dra_01", "dra_07"]
names = {"gwn": "Goodwin Creek", "dra": "Desert Rock"}
x = np.arange(len(keys))
for ax, stat in zip(axes, ("bias", "rms"), strict=True):
    for k, label in enumerate(models):
        vals = [out["open_meteo_inputs"][label][key]["all"][stat] for key in keys]
        ax.bar(x + (k - 1) * 0.27, vals, 0.27, label=label)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(x, [f"{names[key[:3]]}\n{'Jan' if key.endswith('01') else 'Jul'} 2025" for key in keys], fontsize=8)
    ax.set_title(f"modeled minus measured downwelling infrared, {stat} (W/m$^2$)", fontsize=9)
axes[0].legend(fontsize=7)
fig.suptitle("The sky model against SURFRAD pyrgeometers, driven by Open-Meteo weather (July was not used in the fit)", fontsize=10)
fig.tight_layout()
fig.savefig(ROOT / "docs/figures/sky_validation.png", dpi=140)
print(json.dumps({k: v for k, v in out.items() if k not in ("station_inputs", "open_meteo_inputs")}, indent=2))
for label in models:
    print(label, {key: (out["open_meteo_inputs"][label][key]["all"]["bias"], out["open_meteo_inputs"][label][key]["all"]["rms"],
                        out["open_meteo_inputs"][label][key]["overcast"]["bias"]) for key in keys})
