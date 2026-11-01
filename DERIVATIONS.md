# Derivations

## 1. Surface energy balance (`radiative.net_heat_out`)

Per square meter of surface at temperature $T_s$, these are the energy flows:

| term | expression |
|---|---|
| own thermal emission (out) | $\varepsilon\sigma T_s^4$ |
| sky infrared absorbed (in) | $\varepsilon\,F_\text{sky}\,\varepsilon_\text{sky}\sigma T_a^4$ |
| ground infrared absorbed (in) | $\varepsilon\,(1-F_\text{sky})\,\sigma T_g^4$ |
| sunlight absorbed (in) | $\alpha_\text{sol}\, f\, G$ |
| convection (out) | $h\,(T_s - T_a)$ |

**Kirchhoff's law.** In the thermal infrared, a surface's absorptivity equals its emissivity, so the same ε appears in emission and in absorption of sky radiation. Sunlight sits in a different wavelength band (0.3–2.5 µm), so its absorptance $\alpha_\text{sol}$ is independent of ε. That independence is what makes radiative cooling possible: reflect sunlight (low α) while emitting strongly in the IR (high ε).

**Cooling power at ambient.** Setting $T_s = T_a$ removes convection and gives the standard figure of merit:

$$P_\text{cool} = \varepsilon\sigma T_a^4\,(1 - F\varepsilon_\text{sky}) - \alpha G$$

At night under a clear sky this is about 100 W/m², checked in `test_reflective_emitter_cools_below_air_under_clear_sky`.

## 2. Sky emissivity

The atmosphere is mostly transparent between 8 and 13 µm (the "window"), so the clear sky radiates less than a black body at air temperature. Berdahl & Martin (1984) fitted US measurements against dew point:

$$\varepsilon_\text{clear} = 0.711 + 0.56\left(\frac{T_{dp}}{100}\right) + 0.73\left(\frac{T_{dp}}{100}\right)^2 \qquad (T_{dp}\text{ in } ^\circ\text{C})$$

Humid air closes the window (more water vapor emission), so the formula rises with dew point.

**Clouds.** A cloud fills part of the gap between the clear sky and a black sky at air temperature. Low water clouds are warm and nearly black, while thin, cold high cloud adds much less, so I treat Open-Meteo's low, mid and high cloud cover ($c_l, c_m, c_h$) as randomly overlapping layers that each close a share $k$ of the gap:

$$\varepsilon_\text{sky} = \varepsilon_\text{clear} + (1 - \varepsilon_\text{clear})\left[1 - (1 - k_l c_l)(1 - k_m c_m)(1 - k_h c_h)\right]$$

My first version used the total cloud cover as black cloud, $\varepsilon_\text{sky} = c + (1-c)\,\varepsilon_\text{clear}$, and never checked it against a measurement. So I did: NOAA's SURFRAD network measures downwelling infrared with pyrgeometers, and I took hourly means for January and July 2025 at Goodwin Creek, Mississippi (humid) and Desert Rock, Nevada (dry), next to Open-Meteo's weather for the same hours (`data/surfrad/`). Fitting the three $k$ to January only gives $(k_l, k_m, k_h) = (0.80, 0.57, 0.29)$. July, which the fit never saw, then checks it (`examples/validate_sky.py`, `docs/sky_validation.json`):

| driven by Open-Meteo inputs, modeled minus measured (W/m²) | Goodwin Creek, July | Desert Rock, July |
|---|---|---|
| clear sky only | bias −2.7, RMS 11.7 | bias −7.2, RMS 15.7 |
| my first model (total cloud, black) | bias +16.9, RMS 25.9 | bias +2.0, RMS 22.1 |
| cloud layers | bias +5.8, RMS 12.6 | bias −2.8, RMS 13.7 |

Under overcast skies the first model read 43 W/m² high at Goodwin Creek in July and 72 W/m² high at Desert Rock, so it understated how well the box could cool to the sky. The layered model still reads 12 and 23 W/m² high there (only 21 overcast hours at Desert Rock). On clear hours, Berdahl–Martin with each station's own temperature and humidity reads 4 to 15 W/m² low, most in winter. These are two sites and two months, so the error bars are indicative, not a full validation. The tests bound the July error and check that refitting January reproduces the constants.

**Splitting out the window (`radiative.net_heat_out_spectral`).** The share of a blackbody's emission at temperature $T$ that falls between 8 and 13 µm is $f_w(T) = F(13T) - F(8T)$, where $F(\lambda T)$ is the exact series

$$F(\lambda T) = \frac{15}{\pi^4}\sum_{n=1}^{\infty} \frac{e^{-nz}}{n}\left(z^3 + \frac{3z^2}{n} + \frac{6z}{n^2} + \frac{6}{n^3}\right), \qquad z = \frac{C_2}{\lambda T}$$

(checked against direct integration of Planck's law). At 300 K the window holds 32% of the emission. Outside the window, water vapor and CO₂ make the sky nearly black, so I take its emissivity there as 1. Inside, it gets the emissivity that keeps the total at the Berdahl–Martin value:

$$\varepsilon_\text{in} = 1 - \frac{1 - \varepsilon_\text{sky}}{f_w(T_a)}$$

In very dry air that would go negative. Then the window is fully clear ($\varepsilon_\text{in} = 0$) and the outside band takes the rest, $\varepsilon_\text{out} = \varepsilon_\text{sky}/(1 - f_w)$. A surface with emissivities $e_\text{in}$ and $e_\text{out}$ then emits $\sigma T_s^4[e_\text{in} f_w(T_s) + e_\text{out}(1 - f_w(T_s))]$ and absorbs from the sky $\sigma T_a^4[e_\text{in}\varepsilon_\text{in} f_w(T_a) + e_\text{out}\varepsilon_\text{out}(1 - f_w(T_a))]$.

**Checks:** with $e_\text{in} = e_\text{out}$ both sums collapse to the gray formula, because the bands add back up to the total. The tests confirm that for dew points from −15 to 24 °C and in a full enclosure run. A selective surface stagnates colder than a broadband one at night but rejects less heat once it's warmer than the air (`tests/test_models.py`).

## 3. A laminate or wrap over a substrate (`Film.on`)

A film with solar reflectance $R_f$, transmittance $T_f$ and absorptance $A_f = 1 - R_f - T_f$ sits on a substrate with absorptance $\alpha_s$. Follow one unit of sunlight:

1. The film absorbs $A_f$ on the way in. $T_f$ reaches the substrate.
2. The substrate absorbs $\alpha_s T_f$ and reflects $(1-\alpha_s)T_f$ back up into the film.
3. On the way up, the film absorbs $A_f$ of that, transmits $T_f$ out, and reflects $R_f$ back down. Then the cycle repeats.

Each round trip multiplies by $r = R_f(1-\alpha_s)$, a geometric series $\sum r^k = 1/(1-r)$:

$$\alpha_\text{eff} = A_f + \frac{T_f\,[\alpha_s + (1-\alpha_s)A_f]}{1 - R_f(1-\alpha_s)}$$

**Checks:** a perfectly clear film ($T_f=1$) leaves $\alpha_s$ unchanged; a black film gives 1; a mirror gives 0 (`test_film_stack_limits`).

In the IR, most polymer films thicker than about 50 µm are close to opaque, so the stack's emissivity is the film's.

## 4. Enclosure transient (`enclosure.simulate`)

Two lumped nodes with heat capacities $C_s$ (skin) and $C_b$ (cells), coupled by conductance $UA$:

$$C_s\dot T_s = A\,[\text{surface balance}] + UA\,(T_b - T_s) + Q_\text{inv}, \qquad C_b \dot T_b = Q_\text{cell} - UA\,(T_b - T_s)$$

**Battery heat.** The system's one-way efficiency $\eta = 0.95$ (about 90% round trip) is split between the cells and the inverter, $\eta = \eta_\text{cell}\,\eta_\text{inv}$. Following the energy through a charge $c$ (AC in) and a discharge $d$ (AC out):

| | inverter loss | cell loss |
|---|---|---|
| charging | $c\,(1-\eta_\text{inv})$ | $c\,\eta_\text{inv}(1-\eta_\text{cell})$ |
| discharging | $d\,(1/\eta_\text{inv}-1)$ | $(d/\eta_\text{inv})(1/\eta_\text{cell}-1)$ |

The four add up to $c(1-\eta) + d(1/\eta - 1)$, which is everything bought minus everything sold over a day that ends at the charge it started with (`test_battery_heat_is_everything_in_minus_everything_out`). For the cells I use the mean battery round-trip efficiency HTW Berlin measured on 22 home batteries, 95.7% ([Stromspeicher-Inspektion 2025](https://solar.htw-berlin.de/wp-content/uploads/HTW-Stromspeicher-Inspektion-2025.pdf), p. 25), so $\eta_\text{cell} = \sqrt{0.957} = 0.978$ and the inverter gets the rest, $\eta_\text{inv} = 0.971$. That sits between the same study's mean measured discharge conversion efficiency (95.7%, p. 28) and the 97.5% CEC-weighted inverter figure on the [Powerwall 3 datasheet](https://digitalassets.tesla.com/tesla-contents/image/upload/powerwall-3-datasheet-na-en.pdf). About 42% of the loss is then released in the cells. The inverter's share goes into the skin node, as if it were heat-sunk to the enclosure wall.

**Convection.** The skin loses $h\,(T_s - T_a)$ per m² to the air with $h = 2.8 + 3.0\,V$ W/m²K, the convection-only correlation that the [EnergyPlus Engineering Reference](https://bigladdersoftware.com/epx/docs/23-1/engineering-reference/solar-collectors.html) uses for "the convection heat transfer coefficient from the outer cover to the surrounding air". My first version used $h = 5.7 + 3.8\,V$ on the 10 m wind, but the same reference calls that one "the combined conductance from the outer cover to the surrounding" (citing Duffie & Beckman), meaning it already includes radiation, which the model also adds explicitly. So it counted radiation twice.

$V$ is the wind at the box, not the 10 m reanalysis wind. I bring it down with the ASHRAE power law that EnergyPlus uses for every exterior surface ([Local Wind Speed Calculation](https://bigladdersoftware.com/epx/docs/23-1/engineering-reference/outside-surface-heat-balance.html)):

$$V_z = V_\text{met}\left(\frac{\delta_\text{met}}{z_\text{met}}\right)^{\alpha_\text{met}}\left(\frac{z}{\delta}\right)^{\alpha}$$

with the weather wind treated as a 10 m open-field reading ($\alpha_\text{met} = 0.14$, $\delta_\text{met} = 270$ m) and the box at $z = 1$ m in EnergyPlus's "suburbs" terrain ($\alpha = 0.22$, $\delta = 370$ m). That makes the wind at the box 0.43 times the 10 m wind. Open country would make it 0.72 times and a city 0.21 times; both are in the sensitivity study. Treating a reanalysis grid-cell wind as an open-field station reading is my assumption.

**Time step.** The skin time constant $C_s / [A(h + 4\varepsilon\sigma T^3) + UA]$ is 3 to 7 minutes over 2025's winds. Explicit Euler with 1-minute substeps is stable (it needs $\Delta t < 2\tau$).

**Check:** with constant inputs, the steady state must satisfy $UA(T_b - T_s) = Q_\text{gen}$ (`test_enclosure_settles_to_steady_state`).

**Fan.** A ventilation fan moving $\dot V$ of outside air past the cells removes

$$Q_\text{fan} = \epsilon\,\rho c_p \dot V\,(T_b - T_a)$$

from the battery node. Here $\epsilon$ is how close the exhaust air gets to the cell temperature (0.5), and $\rho c_p \approx 1168$ J/m³K for air near 30 °C. At 0.05 m³/s that's 29 W/K, against 12 W/K through the skin. A thermostat on the cells switches it on above 35 °C and off below 32 °C. The gap between them keeps it from switching every substep. It also only runs while $T_a < T_b$, because otherwise ventilation would heat the cells. It draws 30 W while running, priced at the hour's day-ahead price.

**Checks:** a fan that never starts leaves every temperature identical. An always-on fan settles where $UA(T_b - T_s) + Q_\text{fan} = Q_\text{gen}$. With cool air and a load that would take a passive box to about 45 °C, the cells cycle inside the 32–35 °C band. With the air hotter than the cells, the fan stays off even above its set point (`tests/test_models.py`).

## 5. Aging and derating

**Arrhenius aging.** Calendar aging accelerates with temperature: $k(T)/k(25^\circ\text{C}) = \exp\!\left[\frac{E_a}{R}\left(\frac{1}{298.15} - \frac{1}{T}\right)\right]$. I use $E_a = 50$ kJ/mol; the literature for Li-ion calendar fade spans about 30–60 kJ/mol. I report the time-averaged *relative* rate, not absolute capacity fade, because that needs cell-specific data.

**Derating.** Power is available in proportion to $(55 - T)/(55 - 45)$, clipped to [0, 1]. That is a generic BMS policy, not a specific product's.

## 6. Dispatch LP (`dispatch.py`)

For each Central-time day, with hourly prices $p_h$:

$$\max \sum_h p_h (d_h - c_h) \quad \text{s.t.}\quad s_{k} = \tfrac{E}{2} + \sum_{j\le k}\left(\eta c_j - \tfrac{d_j}{\eta}\right) \in [0, E],\; s_{24} = \tfrac{E}{2},\; 0 \le c_h, d_h \le P_h$$

It is linear, so the HiGHS solver returns the global optimum in milliseconds per day. Because the model has perfect foresight, this is an **upper bound** on energy-arbitrage revenue.

**Coupling to temperature (`coupled.run`).** The battery management system caps power at $P\,\text{derate}(T_b)$, and the power it allows changes how much heat the cells make. I want a schedule that never asks for more than the cap implied by its own temperatures:

1. dispatch with no caps and simulate the enclosure
2. tighten the caps, $\text{cap}_k = \min(\text{cap}_{k-1},\ P\,\text{derate}(T_k))$
3. re-dispatch under $\text{cap}_k$ and re-simulate, and stop once every hour satisfies $\max(c_h, d_h) \le P\,\text{derate}(T_h)$ to within 1 W

$T_h$ here is the cells' mean temperature over hour $h$. My first version did a single re-dispatch with caps from the uncapped schedule's temperatures and claimed that was enough. It wasn't: repeating it flips between two schedules, and the one I reported asked for more power than its own temperatures allowed in 182 hours. Because these caps only tighten, the loop can't oscillate, and it settles in 3 to 8 passes. The result is a feasible schedule rather than necessarily the best one, so its revenue is a lower bound on what a perfect-foresight battery could earn under derating, with the uncapped dispatch as the upper bound. `test_coupled_dispatch_never_asks_for_more_than_its_own_temperatures_allow` checks the stopping condition on two hot weeks.

## 7. Spacecraft radiator (`radiator.py`)

There is no atmosphere or convection. The radiator emits $\varepsilon\sigma T^4$ to 3 K space and absorbs:

$$q_\text{abs} = \alpha S\cos\theta + F_E\,(\varepsilon\, q_{IR} + \alpha\, a\, S)$$

with hot-case values from NASA's TFAWS 2014 course [Introduction to On-Orbit Thermal Environments](https://tfaws.nasa.gov/wp-content/uploads/On-Orbit_Thermal_Environments_TFAWS_2014.pdf): $S = 1423$ W/m², "a typical hot case solar flux design parameter" (slide 34), and the 30-minute "combined maximum" pair of albedo $a = 0.27$ and Earth IR $q_{IR} = 260$ W/m² (slide 78), since a radiator's time constant is minutes. My first version used the mean solar constant (1361 W/m²) and called it a hot case. The area to reject heat $Q$ is $A = Q/(\varepsilon\sigma T^4 - q_\text{abs})$.

**Coatings.** The optical solar reflector ($\alpha = 0.07$, $\varepsilon = 0.79$, a silver-backed second-surface mirror), white paint (S-13G, 0.20 and 0.90) and black paint (Chemglaze Z306, 0.96 and 0.91) are fresh-sample values from NASA RP-1121 (Henninger 1984). They're beginning-of-life numbers: paints darken under UV over a mission, which makes the radiators bigger.

**Black paint.** In the hot case its absorbed sunlight exceeds its emission at 300 K, so no area is large enough. That is why radiators use low-α, high-ε surfaces.

**Through an orbit (`radiator.orbit_transient`).** A circular orbit at radius $r = R_E + h$ takes $T = 2\pi\sqrt{r^3/\mu}$, which is 92.4 min at 400 km. With the sun at angle $\beta$ to the orbit plane, the spacecraft is in Earth's cylindrical shadow while $\cos\theta < -\sqrt{1 - (R_E/r)^2}/\cos\beta$, so the eclipse fraction is

$$f = \frac{1}{\pi}\arccos\frac{\sqrt{h^2 + 2R_E h}}{r\cos\beta}$$

or zero once that argument passes 1. At 400 km and $\beta = 0$ that's 39%, 36 minutes. A panel with heat capacity $c$ per m² then follows

$$c\,\dot T = q_\text{int} + q_\text{abs}(t) - \varepsilon\sigma T^4$$

with $q_\text{abs}$ switching between the sunlit value and Earth's infrared alone, stepped with RK4.

In sunlight the panel absorbs the whole hot-case input, albedo included, for the entire sunlit arc. Real albedo fades toward the terminator, so this overstates the sunlit input near sunrise and sunset and makes the hot end of the swing slightly too hot.

**Checks:** the eclipse fraction against counting the points of a tilted orbit that fall inside the shadow cylinder, to $10^{-4}$, and the β above which a 408 km orbit never enters shadow (70.0°) against the "about 71°" in the TFAWS course (slide 132). With no input and no load, the transient against the closed-form radiative cooldown $T = (T_0^{-3} + 3\varepsilon\sigma t/c)^{-1/3}$, to $10^{-9}$. Energy in over an orbit against energy radiated, once the orbits repeat, within 0.2% (`tests/test_models.py`).

## Limitations

- Two infrared bands, not a full spectrum. Taking the sky as black outside the 8–13 µm window ignores the partly clear 16–25 µm region in dry air, so it slightly understates the benefit of emitting there.
- The enclosure is two lumped nodes with no internal gradients. The fan is a simple on/off ventilation fan, and there's no air conditioning or liquid cooling.
- The orbit transient holds the sunlit input, albedo included, constant over the whole sunlit arc.
- No soiling or UV degradation of coatings. Organic films degrade under space UV and atomic oxygen, so the radiator section uses standard spacecraft coatings.
- Two cities (Dallas and Houston), one year of weather and prices.
