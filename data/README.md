# Data sources

| File | Source | Notes |
|---|---|---|
| `weather/dallas_2025_hourly.csv` | [Open-Meteo historical weather API](https://open-meteo.com/) (ERA5-based reanalysis), 32.78 N, 96.80 W | Hourly air temperature, dew point, global horizontal irradiance, 10 m wind, cloud cover. UTC timestamps, 2025 Central-time calendar year. |
| `ercot/dam_2025_lz_north.csv`, `ercot/dam_2025_lz_houston.csv` | ERCOT day-ahead market settlement point prices, load zones North and Houston, fetched with [gridstatus](https://github.com/gridstatus/gridstatus) | $/MWh, hour-start in UTC |

Weather and prices are both 8,760 hourly rows for 2025, aligned in UTC.
