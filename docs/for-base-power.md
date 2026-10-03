# For the Base Power team

You put batteries outside Texas homes and run them against ERCOT prices. This repo couples exactly those two things:

- **Real 2025 inputs:** hourly Dallas weather (Open-Meteo/ERA5) and ERCOT North-zone day-ahead prices.
- **Thermal ↔ dispatch coupling:**
  1. a daily arbitrage LP picks the charge and discharge hours
  2. a two-node enclosure model turns those hours plus the weather into cell temperature
  3. a derating rule (45 °C to 55 °C) caps power
  4. the battery is re-dispatched with those caps
- **Result (generic passive 20 kWh unit):** dark vs white skin is worth **$40 per unit-year** of arbitrage revenue (of a $319 upper bound), **12%** slower calendar aging and **72%** fewer derated hours. The battery is hottest exactly when prices peak. ([figure](figures/ercot_week.png))

**What it doesn't know:** your enclosure, cooling system, cell chemistry or BMS policy. The parameters are generic and the LP has perfect foresight (an upper bound).

**A 3-week project I could do for you, remote:** calibrate the two-node model to logged cell and ambient temperatures from a few units. Then quantify, fleet-wide, the revenue and aging value of skin coating, orientation (sun-facing vs shaded wall) and fan setpoints across your service area. The deliverable is the model, a validation plot and a ranked list of cheap thermal fixes.

— Aneesh Karavadi · aneesh.karavadi@gmail.com
