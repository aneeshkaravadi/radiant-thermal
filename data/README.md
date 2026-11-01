# Data sources

| File | Source | Notes |
|---|---|---|
| `weather/dallas_2025_hourly.csv`, `weather/houston_2025_hourly.csv` | [Open-Meteo historical weather API](https://open-meteo.com/) (ERA5-based reanalysis), Dallas 32.78 N 96.80 W, Houston 29.76 N 95.37 W | Hourly air temperature, dew point, global horizontal irradiance, 10 m wind, total cloud cover and low, mid and high cloud cover. UTC timestamps, 2025 Central-time calendar year. |
| `surfrad/surfrad_2025_hourly.csv` | [NOAA GML SURFRAD](https://gml.noaa.gov/aftp/data/radiation/surfrad/) 1-minute data at Goodwin Creek MS (gwn) and Desert Rock NV (dra), January and July 2025 (public domain, CC0), next to the Open-Meteo archive at each station | Hourly means of measured downwelling infrared, air temperature and relative humidity, keeping minutes with quality flag 0, centered on each Open-Meteo timestamp. Built by `examples/build_surfrad_fixture.py`. SURFRAD: Augustine, DeLuisi and Long (2000), BAMS 81, 2341-2357. |
| `ercot/dam_2025_lz_north.csv`, `ercot/dam_2025_lz_houston.csv` | ERCOT day-ahead market settlement point prices, load zones North and Houston, fetched with [gridstatus](https://github.com/gridstatus/gridstatus) | $/MWh, hour-start in UTC |

Weather and prices are both 8,760 hourly rows for 2025, aligned in UTC.

**Timestamps.** The [Open-Meteo docs](https://open-meteo.com/en/docs/historical-weather-api) list `shortwave_radiation` as the mean over the hour *before* its timestamp, while temperature, dew point, wind and cloud cover are instantaneous values at the timestamp. ERCOT prices are stamped with the start of the hour they apply to. The model steps hour by hour from each hour's start, so `Weather.from_csv` gives each hour the sunlight stamped one hour after its start, and the average of the instantaneous values at its start and end. `test_sunlight_lands_in_the_hour_it_belongs_to` checks this against the sun's position on cloudless days.
