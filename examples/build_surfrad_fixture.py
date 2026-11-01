"""Build data/surfrad/surfrad_2025_hourly.csv from raw downloads (not needed to run the tests).

    python examples/build_surfrad_fixture.py <download_dir>

<download_dir> holds
  surfrad/{gwn,dra}25{001..031,182..212}.dat   NOAA SURFRAD daily files for January and July 2025,
      from https://gml.noaa.gov/aftp/data/radiation/surfrad/<station>/2025/
  om_gwn_2025.json, om_dra_2025.json           Open-Meteo archive API, 2025-01-01..2025-07-31, timezone GMT,
      hourly=temperature_2m,dew_point_2m,cloud_cover,cloud_cover_low,cloud_cover_mid,cloud_cover_high,
      shortwave_radiation,wind_speed_10m and wind_speed_unit=ms, at Goodwin Creek (34.25, -89.87)
      and Desert Rock (36.624, -116.019)

SURFRAD reports 1-minute means with a quality flag after each value (0 = passed every check). I keep
minutes where downwelling infrared, air temperature and relative humidity are all flagged 0, average
them over the hour centered on each Open-Meteo timestamp (Open-Meteo's temperature, dew point and
cloud cover are instantaneous at the timestamp), and keep hours with at least 50 good minutes.
"""
from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
NAMES = ["dw_solar", "uw_solar", "direct_n", "diffuse", "dw_ir", "dw_casetemp", "dw_dometemp", "uw_ir", "uw_casetemp",
         "uw_dometemp", "uvb", "par", "netsolar", "netir", "totalnet", "temp", "rh", "windspd", "winddir", "pressure"]


def read_station(src: Path, station: str) -> pd.DataFrame:
    rows = []
    for f in sorted(glob.glob(str(src / "surfrad" / f"{station}25*.dat"))):
        for line in Path(f).read_text().splitlines()[2:]:
            v = line.split()
            if len(v) < 48:
                continue
            y, _, mo, dy, hh, mm = (int(x) for x in v[:6])
            r = dict(zip(NAMES, (float(x) for x in v[8::2]), strict=False))  # after 8 date fields, value and flag alternate
            q = dict(zip(NAMES, (int(x) for x in v[9::2]), strict=False))
            good = q["dw_ir"] == 0 and q["temp"] == 0 and q["rh"] == 0 and r["dw_ir"] > 0 and r["rh"] > 0
            rows.append((pd.Timestamp(year=y, month=mo, day=dy, hour=hh, minute=mm, tz="UTC"),
                         *((r["dw_ir"], r["temp"], r["rh"]) if good else (np.nan,) * 3)))
    d = pd.DataFrame(rows, columns=["t", "dw_ir", "temp", "rh"]).set_index("t")
    d.index = d.index + pd.Timedelta(30, unit="min")  # so each hourly bin is centered on its label
    g = d.resample("1h", label="left")
    h = g.mean()
    h["n_minutes"] = g["dw_ir"].count()
    return h[h["n_minutes"] >= 50]


def main(src: Path):
    frames = []
    for station in ("gwn", "dra"):
        h = read_station(src, station)
        om = json.loads((src / f"om_{station}_2025.json").read_text())["hourly"]
        o = pd.DataFrame(om)
        o.index = pd.to_datetime(o.pop("time"), utc=True)
        h = h.join(o, how="inner")
        h.insert(0, "station", station)
        frames.append(h)
    df = pd.concat(frames)
    df.index.name = "time_utc"
    df = df.rename(columns={"dw_ir": "lw_down_wm2", "temp": "station_temp_c", "rh": "station_rh_pct",
                            "temperature_2m": "om_temp_c", "dew_point_2m": "om_dew_c", "cloud_cover": "om_cloud_pct",
                            "cloud_cover_low": "om_cloud_low_pct", "cloud_cover_mid": "om_cloud_mid_pct",
                            "cloud_cover_high": "om_cloud_high_pct"})
    df = df.drop(columns=["shortwave_radiation", "wind_speed_10m"])
    df = df.round({"lw_down_wm2": 1, "station_temp_c": 2, "station_rh_pct": 1})
    out = ROOT / "data/surfrad/surfrad_2025_hourly.csv"
    df.to_csv(out, date_format="%Y-%m-%dT%H:%MZ")
    print(df.groupby("station").size())


if __name__ == "__main__":
    main(Path(sys.argv[1]))
