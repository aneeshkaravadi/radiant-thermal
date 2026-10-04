# radiant-thermal

[![tests](https://github.com/aneeshkaravadi/radiant-thermal/actions/workflows/ci.yml/badge.svg)](https://github.com/aneeshkaravadi/radiant-thermal/actions/workflows/ci.yml)

Does the paint colour of an outdoor home battery matter, in dollars? I live in North Texas, where batteries sit outside in summer sun and the grid's prices spike on the hottest afternoons, so I wanted to see whether those two things interact. They do.

![Battery temperature over the summer](docs/figures/battery_summer_hist.png)

## What I simulated

A generic 20 kWh, 5 kW battery in a passive outdoor enclosure, run through every hour of 2025 in Dallas. All of the inputs are real data:

- **weather:** hourly temperature, dew point, sunlight, wind and cloud cover from Open-Meteo
- **prices:** ERCOT's 2025 day-ahead prices for the North load zone, pulled with gridstatus

Each day a small linear program decides when to charge and discharge to make the most money. The enclosure is a two-node thermal model (the steel skin and the battery inside), and the skin trades heat with the sun, the sky and the air. When the cells pass 45 °C the battery management system starts cutting power, which means less to sell, so I feed that derating back into the dispatch and run it again.

## What came out

| Skin | Summer avg cell temp | Hours derated | Aging rate vs 25 °C | 2025 revenue |
|---|---|---|---|---|
| dark grey paint | 42.3 °C | 1,195 | 2.29× | $264 |
| white paint | 40.7 °C | 339 | 2.01× | $304 |
| ideal radiative cooler | 39.7 °C | 208 | 1.88× | $312 |

If the battery never had to derate it would make $319, so the dark box loses about $56 of that and the white box about $16. White paint is worth roughly $40 per battery per year here, plus 12% slower calendar aging, which isn't nothing across a fleet.

The reason is easiest to see in one week of August: the cells are hottest right when prices peak, which is exactly when you most want the battery at full power.

![ERCOT week](docs/figures/ercot_week.png)

## Dallas vs Houston

I reran everything for Houston with its own 2025 weather and its own ERCOT load-zone prices (`examples/compare_cities.py`). I expected Houston to be easier on the batteries, since its summer air was actually 0.7 °C cooler than Dallas's. It wasn't: the dark box derated for 1,308 hours against 1,195 in Dallas, and aged faster too. The difference is humidity. Houston's summer dew point averages 23.6 °C against 21.0 °C in Dallas, and water vapour closes the sky's infrared window, so the enclosure can't radiate heat away as well at night or during the day. Air temperature alone would have told the wrong story.

<img src="docs/figures/city_comparison.png" width="80%">

## Wraps and laminates

I also expected that laminating a film or wrap over a white enclosure could help, and mostly it can't. Sunlight that passes through the film bounces between it and the paint, so I summed those bounces as a geometric series. The result is that on an opaque box, only a film that absorbs almost no sunlight beats plain white paint, and anything that absorbs UV or near-infrared makes it worse.

<img src="docs/figures/laminate_design_map.png" width="55%">

## Same equation, in orbit

A spacecraft radiator is the same energy balance with no atmosphere and no air, just sunlight, Earth's infrared and reflected sunlight coming in against $\varepsilon\sigma T^4$ going out. At 300 K in a hot low-Earth-orbit case, an optical solar reflector needs about 4.1 m² per kW, white paint 5.3, and black paint can't get rid of heat at all, which my first version reported as a literal infinity in the results file.

<img src="docs/figures/radiator_area.png" width="85%">

## Mistakes and fixes

- I first pulled the weather in local time while the ERCOT prices came timezone-aware, and around daylight saving those two don't line up hour for hour. Everything is in UTC now, and a test checks the two files match row for row.
- My first pass didn't feed derating back into dispatch, so it counted revenue the battery couldn't actually earn while it was hot.
- My first "energy balance" test only checked that the temperatures came out finite, which doesn't prove anything. I replaced it with one that holds the weather constant and checks that the battery settles exactly where the heat it conducts out equals the heat it generates.

## Running it

```bash
pip install -e ".[dev]"
pytest -q
python examples/make_figures.py
```

The enclosure numbers describe a generic passive unit, not anyone's product. Real units with fans or active cooling run cooler, so the differences between skins are the result, not the absolute temperatures. Assumptions and derivations are in [DERIVATIONS.md](DERIVATIONS.md), and data sources are in [data/README.md](data/README.md).

---

Aneesh Karavadi, engineering at UNT (TAMS). Separately from this repo, I'm an undergraduate researcher in Dr. Zihao Richard Zhang's lab at UNT, testing transparent UV/IR-blocking radiative-cooling films. A window solar-heat-gain model would be a natural place to apply that kind of film here someday. I used Claude Code to write a lot of the implementation, but the questions and conclusions are mine.
