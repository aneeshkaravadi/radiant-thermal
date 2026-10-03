# radiant-thermal

**One surface energy balance, three uses: coatings and wraps, outdoor batteries dispatched on real ERCOT prices, and spacecraft radiators.**

[![tests](https://github.com/aneeshkaravadi/radiant-thermal/actions/workflows/ci.yml/badge.svg)](https://github.com/aneeshkaravadi/radiant-thermal/actions/workflows/ci.yml)

![Battery temperature distribution by coating](docs/figures/battery_summer_hist.png)

## What this shows

I simulated a passive 20 kWh / 5 kW outdoor battery enclosure through **every hour of 2025 in Dallas**. It was dispatched against **real ERCOT North-zone day-ahead prices** by a daily linear program. Its cell temperature fed back into a battery-management derating rule (full power below 45 °C, none at 55 °C).

| Enclosure skin | Summer mean cell temp | Hours derated | Calendar-aging rate vs 25 °C | 2025 arbitrage revenue |
|---|---|---|---|---|
| dark grey paint (α 0.75) | 42.3 °C | 1,195 | 2.29× | $264 |
| white paint (α 0.25) | 40.7 °C | 339 | 2.01× | $304 |
| ideal radiative cooler (α 0.04, ε 0.95) | 39.7 °C | 208 | 1.88× | $312 |

The unconstrained arbitrage upper bound is **$319**. Painting the box white instead of dark grey recovers $40 of that per unit per year. It also slows Arrhenius calendar aging by 12% and cuts derated hours by 72%.

The reason is in one plot: the battery runs hottest in exactly the hours when ERCOT prices spike.

![ERCOT week](docs/figures/ercot_week.png)

### Does a wrap or laminate help?

For an **opaque** enclosure, a film or wrap laminated over white paint only helps if it reflects more sunlight than it absorbs. The stack model sums the light bouncing between film and substrate, and the design map shows the break-even line against plain white paint:

<img src="docs/figures/laminate_design_map.png" width="60%">

Only laminates that barely absorb sunlight beat plain white. Anything that absorbs UV or near-IR adds heat to an opaque box.

### Same physics in orbit

A spacecraft radiator balances $\varepsilon\sigma T^4$ against absorbed sunlight, Earth infrared and albedo, the same two numbers (α, ε) with no atmosphere. At 300 K in the LEO hot case, an optical solar reflector needs 4.1 m² per kW and white paint 5.3 m². Black paint cannot reject heat at all.

<img src="docs/figures/radiator_area.png" width="90%">

## Checks

9 tests in CI, including:
- the Berdahl–Martin sky emissivity formula
- a black body in overcast darkness sitting exactly at air temperature
- the laminate stack's limits (invisible film, black film and mirror film behave correctly)
- the enclosure settling to the analytic steady state
- dispatch respecting state-of-charge and power limits, with a power cap reducing revenue
- the radiator closed form
- the real weather and price files aligning hour-for-hour in UTC

## Quick start

```bash
git clone https://github.com/aneeshkaravadi/radiant-thermal && cd radiant-thermal
pip install -e ".[dev]"
pytest -q
python examples/make_figures.py      # ~30 s, writes docs/figures and docs/results.json
```

## How it works

| Module | What it does |
|---|---|
| `radiative.py` | Sky emissivity from dew point and cloud cover, surface energy balance, laminate-on-substrate optics |
| `enclosure.py` | Two-node (skin, battery) transient model on hourly weather; Arrhenius aging; BMS derating |
| `dispatch.py` | Daily perfect-foresight arbitrage LP on ERCOT prices, with optional per-hour power caps |
| `radiator.py` | Spacecraft radiator area and equilibrium temperature in LEO or deep space |

Derivations and assumptions are in [DERIVATIONS.md](DERIVATIONS.md); data sources are in [data/README.md](data/README.md).

> [!NOTE]
> The enclosure parameters describe a generic passive unit, not any company's product. Real units with active cooling run cooler. The *differences* between coatings are the result, not the absolute temperatures.

## About

Built by **Aneesh Karavadi**, an engineering student at the University of North Texas (TAMS). I used **Claude Code** as a pair programmer. The modelling choices and conclusions are mine to defend.

Separately, I'm an undergraduate researcher in Dr. Zihao Richard Zhang's lab at UNT, testing transparent UV/IR-blocking radiative-cooling films. That work isn't part of this repo. A potential application would be adding a window solar-heat-gain model here, where a transparent film's UV and near-IR blocking is the whole point.
