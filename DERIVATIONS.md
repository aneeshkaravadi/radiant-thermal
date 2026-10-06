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

Humid air closes the window (more water vapor emission), so the formula rises with dew point. Clouds are treated as near-black emitters at air temperature covering a fraction $c$ of the sky: $\varepsilon_\text{sky} = c + (1-c)\,\varepsilon_\text{clear}$.

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

$$C_s\dot T_s = A\,[\text{surface balance}] + UA\,(T_b - T_s), \qquad C_b \dot T_b = Q_\text{gen} - UA\,(T_b - T_s)$$

**Battery heat.** $Q_\text{gen} = (1-\eta)(P_\text{charge} + P_\text{discharge})$, the one-way conversion loss.

**Time step.** The skin time constant is $C_s / [A(h + 4\varepsilon\sigma T^3)] \approx 5$ min. Explicit Euler with 1-minute substeps is stable (it needs $\Delta t < 2\tau$).

**Check:** with constant inputs, the steady state must satisfy $UA(T_b - T_s) = Q_\text{gen}$ (`test_enclosure_settles_to_steady_state`).

## 5. Aging and derating

**Arrhenius aging.** Calendar aging accelerates with temperature: $k(T)/k(25^\circ\text{C}) = \exp\!\left[\frac{E_a}{R}\left(\frac{1}{298.15} - \frac{1}{T}\right)\right]$. I use $E_a = 50$ kJ/mol; the literature for Li-ion calendar fade spans about 30–60 kJ/mol. I report the time-averaged *relative* rate, not absolute capacity fade, because that needs cell-specific data.

**Derating.** Power is available in proportion to $(55 - T)/(55 - 45)$, clipped to [0, 1]. That is a generic BMS policy, not a specific product's.

## 6. Dispatch LP (`dispatch.py`)

For each Central-time day, with hourly prices $p_h$:

$$\max \sum_h p_h (d_h - c_h) \quad \text{s.t.}\quad s_{k} = \tfrac{E}{2} + \sum_{j\le k}\left(\eta c_j - \tfrac{d_j}{\eta}\right) \in [0, E],\; s_{24} = \tfrac{E}{2},\; 0 \le c_h, d_h \le P_h$$

It is linear, so the HiGHS solver returns the global optimum in milliseconds per day. Because the model has perfect foresight, this is an **upper bound** on energy-arbitrage revenue. The thermal coupling runs in three steps:
1. dispatch unconstrained
2. simulate temperature, then turn temperature into per-hour power caps
3. re-dispatch and re-simulate

One fixed-point pass is enough here because the caps change heat generation only a little.

## 7. Spacecraft radiator (`radiator.py`)

There is no atmosphere or convection. The radiator emits $\varepsilon\sigma T^4$ to 3 K space and absorbs:

$$q_\text{abs} = \alpha S\cos\theta + F_E\,(\varepsilon\, q_{IR} + \alpha\, a\, S)$$

with $S = 1361$ W/m², Earth IR $q_{IR} \approx 237$ W/m², albedo $a \approx 0.30$. The area to reject heat $Q$ is $A = Q/(\varepsilon\sigma T^4 - q_\text{abs})$.

**Black paint.** In the hot case its absorbed sunlight exceeds its emission at 300 K, so no area is large enough. That is why radiators use low-α, high-ε surfaces.

## Limitations

- Band-averaged (gray) optics, not spectral. A spectral model would resolve the 8–13 µm window explicitly.
- The enclosure is two lumped nodes: no internal gradients, no fans or HVAC.
- No soiling or UV degradation of coatings. Organic films degrade under space UV and atomic oxygen, so the radiator section uses standard spacecraft coatings.
- One weather location, one year.
