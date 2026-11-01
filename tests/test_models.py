from itertools import pairwise
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from rt import coupled
from rt import dispatch as dp
from rt import enclosure as enc
from rt import radiative as rad
from rt import radiator
from rt import surfrad

ROOT = Path(__file__).resolve().parents[1]


def test_cloud_layers_limits():
    """No cloud gives the clear sky; layers that close the whole gap give a black sky; one black layer
    of cover c gives the black-cloud form c + (1 - c) eps_clear; more cloud never lowers the sky."""
    clear = rad.clear_sky_emissivity(15.0)
    assert rad.sky_emissivity(15.0) == pytest.approx(clear)
    assert rad.sky_emissivity(15.0, 1.0, 1.0, 1.0, effect=(1.0, 1.0, 1.0)) == pytest.approx(1.0)
    assert rad.sky_emissivity(15.0, 0.4, effect=(1.0, 0.0, 0.0)) == pytest.approx(0.4 + 0.6 * clear)
    c = np.linspace(0, 1, 6)
    for layer in range(3):
        cover = [np.zeros(6)] * 3
        cover[layer] = c
        assert np.all(np.diff(rad.sky_emissivity(15.0, *cover)) > 0)


def test_blackbody_surface_in_dark_overcast_sits_at_air_temperature():
    c = rad.Coating("black", 1.0, 1.0)
    T = rad.steady_surface_temperature(c, 300.0, 0.0, 1.0, 10.0)
    assert T == pytest.approx(300.0, abs=1e-6)


def test_reflective_emitter_cools_below_air_under_clear_sky():
    T = rad.steady_surface_temperature(rad.IDEAL_COOLER, 303.0, 900.0, rad.sky_emissivity(10.0), 8.0)
    assert T < 303.0
    assert rad.cooling_power_at_ambient(rad.IDEAL_COOLER, 303.0, 0.0, 0.75) > 80  # night, clear: ~100 W/m^2


def test_blackbody_fraction_matches_planck_integration():
    from scipy.integrate import quad
    c1 = 3.741771852e8  # W um^4 / m^2

    def planck(x):
        return c1 / (x**5 * np.expm1(rad.C2 / x))

    for lam_t in (1500.0, 2898.0, 5000.0, 10000.0):
        edges = [100.0] + [e for e in (1000.0, 2000.0, 4000.0) if e < lam_t] + [lam_t]
        direct = sum(quad(planck, lo, hi, epsabs=0, epsrel=1e-12, limit=500)[0] for lo, hi in pairwise(edges))
        assert rad.blackbody_fraction(lam_t) == pytest.approx(direct / rad.SIGMA, abs=1e-6)
    assert rad.blackbody_fraction(2898.0) == pytest.approx(0.25, abs=1e-3)  # a quarter lies below the peak


@pytest.mark.parametrize("dew", [-15.0, 0.0, 12.0, 24.0])
def test_two_band_model_gives_a_gray_surface_the_gray_answer(dew):
    """The window split keeps the sky's total, so a gray surface sees no difference, even in dry air."""
    eps_sky = rad.sky_emissivity(dew, 0.2)
    for Ts in (285.0, 300.0, 330.0):
        gray = rad.net_heat_out(Ts, rad.WHITE_PAINT, 300.0, 500.0, eps_sky, 6.0, F_sky=0.75, Tg=305.0)
        two = rad.net_heat_out_spectral(Ts, rad.WHITE_PAINT, 300.0, 500.0, eps_sky, 6.0, F_sky=0.75, Tg=305.0)
        assert two == pytest.approx(gray, rel=1e-9)


def test_selective_emitter_runs_colder_below_ambient_but_cools_less_above_it():
    """The classic trade: emitting only through the window keeps out the sky's other bands (colder stagnation),
    but rejects less heat once the surface is warmer than the air."""
    broad = rad.Coating("broadband", 0.05, 0.95)
    sel = rad.selective_coating("selective", 0.05, 0.95, 0.05)
    eps_sky = rad.sky_emissivity(10.0)
    T_broad = rad.steady_surface_temperature(broad, 300.0, 0.0, eps_sky, 2.0, spectral=True)
    T_sel = rad.steady_surface_temperature(sel, 300.0, 0.0, eps_sky, 2.0, spectral=True)
    assert T_sel < T_broad < 300.0
    hot = 320.0
    assert (rad.net_heat_out_spectral(hot, sel, 300.0, 0.0, eps_sky, 0.0)
            < rad.net_heat_out_spectral(hot, broad, 300.0, 0.0, eps_sky, 0.0))


def test_film_stack_limits():
    sub = rad.Coating("s", 0.25, 0.9)
    clear = rad.Film("clear", 0.0, 1.0, 0.9).on(sub)
    assert clear.alpha_solar == pytest.approx(0.25)  # invisible film changes nothing
    opaque = rad.Film("black", 0.0, 0.0, 0.9).on(sub)
    assert opaque.alpha_solar == pytest.approx(1.0)
    mirror = rad.Film("mirror", 1.0, 0.0, 0.9).on(sub)
    assert mirror.alpha_solar == pytest.approx(0.0)


def test_enclosure_settles_to_steady_state():
    """With constant weather and heat load, the battery must settle where UA*(Tb - Ts) equals the load."""
    n = 48
    t = pd.date_range("2025-07-01", periods=n, freq="h", tz="UTC")
    w = enc.Weather(t, np.full(n, 305.0), np.full(n, 15.0), np.full(n, 400.0),
                    np.full(n, 2.0), np.zeros((n, 3)))
    q = np.full(n, 200.0)
    e = enc.Enclosure()
    r = enc.simulate(rad.WHITE_PAINT, w, q, e)
    assert np.all(np.isfinite(r.T_batt))
    # battery time constant C_b/UA ~ 4.6 h, so after 48 h it is at steady state
    assert e.ua_internal * (r.T_batt[-1] - r.T_skin[-1]) == pytest.approx(200.0, rel=0.01)


def test_heat_into_the_skin_acts_like_extra_absorbed_sunlight():
    """Q watts put into the skin enter the same balance as Q / (A alpha f_sun) more W/m^2 of sunlight."""
    n, Q = 48, 150.0
    e, c = enc.Enclosure(), rad.WHITE_PAINT
    base = constant_weather(n, 305.0, ghi=400.0)
    brighter = constant_weather(n, 305.0, ghi=400.0 + Q / (e.area * c.alpha_solar * e.sun_factor))
    q = np.full(n, 100.0)
    a = enc.simulate(c, base, q, e, q_skin_w=np.full(n, Q))
    b = enc.simulate(c, brighter, q, e)
    assert a.T_batt == pytest.approx(b.T_batt, abs=1e-9) and a.T_skin == pytest.approx(b.T_skin, abs=1e-9)


def constant_weather(n, Ta, ghi=0.0, dew=15.0):
    t = pd.date_range("2025-07-01", periods=n, freq="h", tz="UTC")
    return enc.Weather(t, np.full(n, Ta), np.full(n, dew), np.full(n, ghi), np.full(n, 2.0), np.zeros((n, 3)))


def test_fan_that_never_starts_changes_nothing():
    w = constant_weather(48, 305.0, ghi=400.0)
    q = np.full(48, 200.0)
    plain = enc.simulate(rad.WHITE_PAINT, w, q)
    idle = enc.simulate(rad.WHITE_PAINT, w, q, fan=enc.Fan(on_c=200.0, off_c=199.0))
    assert np.array_equal(plain.T_batt, idle.T_batt) and not idle.fan_on.any()


def test_running_fan_balances_conduction_and_ventilation():
    """A fan that's always on: at steady state the load leaves through the skin and the air stream."""
    w = constant_weather(72, 300.0)
    fan = enc.Fan(on_c=-100.0, off_c=-101.0)
    e = enc.Enclosure()
    r = enc.simulate(rad.WHITE_PAINT, w, np.full(72, 400.0), e, fan=fan)
    out = e.ua_internal * (r.T_batt[-1] - r.T_skin[-1]) + fan.ua * (r.T_batt[-1] - 300.0)
    assert out == pytest.approx(400.0, rel=0.01)
    assert r.fan_on[-1] == 1.0 and r.fan_energy_kwh(fan)[-1] == pytest.approx(fan.power_w / 1000)


def test_fan_thermostat_holds_the_cells_in_its_band():
    """Cool air and a load that would take the cells to ~45 C passive: the fan cycles between its set points."""
    w = constant_weather(96, 293.15)
    fan = enc.Fan()
    r = enc.simulate(rad.WHITE_PAINT, w, np.full(96, 300.0), fan=fan)
    late = slice(48, None)
    assert r.celsius()[late].min() > fan.off_c - 0.1 and r.celsius()[late].max() < fan.on_c + 0.1
    assert 0.1 < r.fan_on[late].mean() < 0.9


def test_fan_does_not_blow_in_hotter_air():
    """Air at 47 C and no load: the cells stay below the air, so the fan never runs even above its set point."""
    w = constant_weather(24, 320.0)
    r = enc.simulate(rad.WHITE_PAINT, w, np.zeros(24), fan=enc.Fan())
    assert r.celsius().max() > enc.Fan().on_c and not r.fan_on.any()


def test_skin_temperature_responds_to_absorptance_like_the_linearized_balance():
    """With no heat from the battery, the steady skin warms by f G / (h + 4 eps sigma T^3) per unit of absorptance."""
    w = constant_weather(72, 305.0, ghi=800.0)
    e, a, da = enc.Enclosure(), 0.40, 0.01
    T = [enc.simulate(rad.Coating("c", aa, 0.9), w, np.zeros(72), e).T_skin[-1] for aa in (a, a + da)]
    h = e.h(2.0)
    expected = e.sun_factor * 800.0 / (h + 4 * 0.9 * rad.SIGMA * T[0] ** 3)
    assert (T[1] - T[0]) / da == pytest.approx(expected, rel=0.02)


def test_two_band_enclosure_matches_the_gray_one_for_a_gray_skin():
    w = constant_weather(48, 305.0, ghi=400.0, dew=18.0)
    q = np.full(48, 200.0)
    gray = enc.simulate(rad.WHITE_PAINT, w, q)
    two = enc.simulate(rad.WHITE_PAINT, w, q, spectral=True)
    assert two.T_batt == pytest.approx(gray.T_batt, abs=1e-9)


def test_local_wind_profile_and_convection():
    """The ASHRAE power law returns the 10 m reading in open country at 10 m, and gets slower lower down
    and over rougher ground. The convection-only coefficient is below the combined one at every speed."""
    assert rad.local_wind(3.0, 10.0, "country") == pytest.approx(3.0)
    assert rad.local_wind(3.0, 1.0, "city") < rad.local_wind(3.0, 1.0, "suburbs") < rad.local_wind(3.0, 1.0, "country") < 3.0
    v = np.linspace(0, 10, 11)
    assert np.all(rad.convection_coefficient(v) < rad.combined_coefficient(v))
    # EnergyPlus's suburban profile at 1 m: (270/10)^0.14 * (1/370)^0.22
    assert rad.local_wind(1.0, 1.0, "suburbs") == pytest.approx(27**0.14 * (1 / 370) ** 0.22)


def test_arrhenius_reference():
    assert enc.arrhenius_factor(298.15) == pytest.approx(1.0)
    # 37.3 kJ/mol: exp(37300 / 8.314462 * (1/298.15 - 1/308.15)) = 1.63 for 10 K
    assert enc.arrhenius_factor(308.15) == pytest.approx(1.630, abs=0.001)


def test_dispatch_matches_the_hand_worked_optimum():
    """$20 for 00-11 UTC and $100 for 12-23 UTC, two UTC days. In Central daylight time that's three local days:
    Jul 31 19:00-23:00 (cheap only), all of Aug 1 (cheap 0-6, dear 7-18, cheap 19-23) and Aug 2 00:00-18:00.
    Each local day starts and ends half full (10 of 20 kWh).
      Jul 31: nothing to gain.
      Aug 1:  fill 10 -> 20 kWh at $20, empty 20 -> 0 at $100, refill 0 -> 10 at $20.
      Aug 2:  fill 10 -> 20 at $20, sell 10 kWh of storage at $100.
    Buying x kWh of storage costs x / eta at the meter and selling it returns x * eta."""
    t = pd.date_range("2025-08-01", periods=48, freq="h", tz="UTC")
    price = pd.Series(np.tile(np.r_[np.full(12, 20.0), np.full(12, 100.0)], 2), index=t)
    b = dp.Battery(20.0, 5.0, 0.95)
    eta = b.eta_one_way
    r = dp.dispatch(price, b)
    aug1 = (100 * 20 * eta - 20 * 20 / eta) / 1000
    aug2 = (100 * 10 * eta - 20 * 10 / eta) / 1000
    assert r.revenue_usd == pytest.approx(aug1 + aug2, rel=1e-6)
    assert r.soc_kwh.min() >= -1e-6 and r.soc_kwh.max() <= b.energy_kwh + 1e-6
    assert r.charge_kw.max() <= 5 + 1e-9 and r.discharge_kw.max() <= 5 + 1e-9
    # Capped at 1 kW the cheap hours run out first: Aug 1 has 12 cheap hours (12 kWh in, 12 eta^2 out),
    # Aug 2 has 7 before the dear block (7 kWh in, 7 eta^2 out).
    capped = dp.dispatch(price, b, power_cap_kw=np.full(48, 1.0))
    expected = sum((100 * n * eta**2 - 20 * n) / 1000 for n in (12, 7))
    assert capped.revenue_usd == pytest.approx(expected, rel=1e-6)


@pytest.mark.parametrize("cell_round_trip", [0.95**2, 0.957, 0.99])
def test_battery_heat_is_everything_in_minus_everything_out(cell_round_trip):
    """Over days that start and end at the same charge, cell heat plus inverter heat must equal
    the AC energy bought minus the AC energy sold."""
    t = pd.date_range("2025-08-01 05:00", periods=72, freq="h", tz="UTC")  # three whole Central-time days
    rng = np.random.default_rng(1)
    price = pd.Series(rng.uniform(10, 200, 72), index=t)  # synthetic prices, just to make the battery cycle
    b = dp.Battery(cell_round_trip=cell_round_trip)
    r = dp.dispatch(price, b)
    assert r.charge_kw.sum() > 10
    heat_kwh = (r.cell_heat_w(b).sum() + r.inverter_heat_w(b).sum()) / 1000
    assert heat_kwh == pytest.approx(r.charge_kw.sum() - r.discharge_kw.sum(), rel=1e-9)
    if cell_round_trip == pytest.approx(0.95**2):
        assert r.inverter_heat_w(b).max() == pytest.approx(0.0, abs=1e-6)  # all the loss in the cells


@pytest.mark.parametrize("city, zone", [("dallas", "lz_north"), ("houston", "lz_houston")])
def test_real_data_files_align(city, zone):
    w = enc.Weather.from_csv(ROOT / f"data/weather/{city}_2025_hourly.csv")
    p = dp.load_prices(ROOT / f"data/ercot/dam_2025_{zone}.csv")
    assert len(w.Ta) == len(p) == 8760 and (w.time == p.index).all()


def cos_zenith(times: pd.DatetimeIndex, lat: float, lon: float) -> np.ndarray:
    """Cosine of the solar zenith angle (0 below the horizon), from NOAA's General Solar Position Calculations."""
    doy = times.dayofyear.to_numpy()
    hr = times.hour.to_numpy() + times.minute.to_numpy() / 60 + times.second.to_numpy() / 3600
    g = 2 * np.pi / 365 * (doy - 1 + (hr - 12) / 24)
    eqtime = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g) - 0.014615 * np.cos(2 * g)
                       - 0.040849 * np.sin(2 * g))
    decl = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g) - 0.006758 * np.cos(2 * g)
            + 0.000907 * np.sin(2 * g) - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    ha = np.radians((hr * 60 + eqtime + 4 * lon) / 4 - 180)  # times are UTC, so no timezone offset
    la = np.radians(lat)
    return np.clip(np.sin(la) * np.sin(decl) + np.cos(la) * np.cos(decl) * np.cos(ha), 0.0, None)


@pytest.mark.parametrize("city, lat, lon", [("dallas", 32.78, -96.80), ("houston", 29.76, -95.37)])
def test_sunlight_lands_in_the_hour_it_belongs_to(city, lat, lon):
    """On cloudless days, the sunlight the model applies to the hour starting at t should follow the sun's
    height averaged over [t, t + 1 h), not over the hour before or the instant t."""
    path = ROOT / f"data/weather/{city}_2025_hourly.csv"
    w = enc.Weather.from_csv(path)
    raw = pd.read_csv(path)
    local_day = pd.Series(w.time.tz_convert("America/Chicago").date)
    clear = raw["cloud_cover"].groupby(local_day).transform("max").to_numpy() < 10
    clear &= np.r_[clear[1:], False]  # the sunlight for hour i comes from row i + 1, so both must be clear
    minutes = [pd.Timedelta(m + 0.5, unit="min") for m in range(60)]
    same_hour = np.mean([cos_zenith(w.time + m, lat, lon) for m in minutes], axis=0)
    hour_before = np.mean([cos_zenith(w.time - m, lat, lon) for m in minutes], axis=0)
    instant = cos_zenith(w.time, lat, lon)
    sel = clear & ((w.ghi > 0) | (same_hour > 0))
    r = {name: np.corrcoef(w.ghi[sel], geo[sel])[0, 1]
         for name, geo in (("same", same_hour), ("before", hour_before), ("instant", instant))}
    assert sel.sum() > 250
    assert r["same"] > 0.99
    assert r["same"] > r["instant"] + 0.02 and r["same"] > r["before"] + 0.05


def test_radiator_closed_form():
    c = rad.Coating("c", 0.2, 0.9)
    A = radiator.area_per_kw(c, radiator.DEEP_SPACE, 300.0)
    assert A == pytest.approx(1000.0 / (0.9 * rad.SIGMA * 300.0**4))
    T_eq = radiator.equilibrium_temperature(c, radiator.Environment("sun normal", 1.0, 0.0, solar_flux=1361.0))
    assert T_eq == pytest.approx((0.2 * 1361 / (0.9 * rad.SIGMA)) ** 0.25)


def test_no_eclipse_above_the_critical_beta_nasa_quotes():
    """NASA's TFAWS 2014 environments course (slide 132) says its 408 km circular orbit sees no eclipse above
    |beta| of about 71 deg. My cylindrical shadow gives asin(R/r) = 70.0 deg."""
    critical = np.degrees(np.arcsin(radiator.R_EARTH / (radiator.R_EARTH + 408e3)))
    assert critical == pytest.approx(71.0, abs=1.5)
    assert radiator.eclipse_fraction(408e3, critical - 0.5) > 0.0
    assert radiator.eclipse_fraction(408e3, critical + 0.5) == 0.0



# ---------------------------------------------------------------- radiator through an orbit

def test_eclipse_fraction_matches_a_numerical_shadow_check():
    """Walk around the orbit and count the points inside Earth's cylindrical shadow (sun along +x)."""
    for alt, beta in ((400e3, 0.0), (400e3, 50.0), (1000e3, 20.0)):
        r = radiator.R_EARTH + alt
        th = np.linspace(0, 2 * np.pi, 200001)[:-1]
        b = np.radians(beta)
        # orbit plane tilted by beta out of the sun direction
        x, y, z = r * np.cos(th) * np.cos(b), r * np.sin(th), r * np.cos(th) * np.sin(b)
        shadow = (x < 0) & (np.hypot(y, z) < radiator.R_EARTH)
        assert radiator.eclipse_fraction(alt, beta) == pytest.approx(shadow.mean(), abs=1e-4)
    assert radiator.eclipse_fraction(400e3, 75.0) == 0.0  # past the critical beta the orbit never enters shadow


def test_radiator_cools_down_like_the_closed_form_in_the_dark():
    """No load and nothing coming in: c dT/dt = -eps sigma T^4 gives T = (T0^-3 + 3 eps sigma t / c)^(-1/3)."""
    dark = radiator.Environment("nothing", 0.0, 0.0)
    c = rad.Coating("panel", 0.1, 0.8)
    t, T, _ = radiator.orbit_transient(c, dark, 5000.0, 0.0, orbits=1, T0=300.0)
    exact = (300.0**-3 + 3 * 0.8 * rad.SIGMA * t / 5000.0) ** (-1 / 3)
    assert T == pytest.approx(exact, rel=1e-9)


def test_radiator_balances_energy_over_an_orbit():
    """Once the orbits repeat, what comes in over one orbit equals what is radiated."""
    c, env, cap, q = radiator.OSR, radiator.LEO_HOT, 3000.0, 150.0
    t, T, sun = radiator.orbit_transient(c, env, cap, q, orbits=8)
    last = t >= t[-1] - radiator.orbit_period(400e3)
    q_in = np.where(sun, radiator.absorbed_flux(c, env), env.earth_view * c.eps_ir * env.earth_ir) + q
    assert np.mean(q_in[last]) == pytest.approx(np.mean(c.eps_ir * rad.SIGMA * T[last] ** 4), rel=2e-3)


# ---------------------------------------------------------------- dispatch and temperature together

def august_slice(days=14):
    w = enc.Weather.from_csv(ROOT / "data/weather/dallas_2025_hourly.csv")
    p = dp.load_prices(ROOT / "data/ercot/dam_2025_lz_north.csv")
    local = w.time.tz_convert("America/Chicago")
    i0 = int(np.argmax(local >= pd.Timestamp("2025-08-01", tz="America/Chicago")))
    sl = slice(i0, i0 + 24 * days)
    ws = enc.Weather(w.time[sl], w.Ta[sl], w.dew_c[sl], w.ghi[sl], w.wind[sl], w.cloud[sl])
    return ws, p.iloc[sl]


def test_coupled_dispatch_never_asks_for_more_than_its_own_temperatures_allow():
    """Two hot weeks with every watt of loss in the cells, so derating bites. The converged schedule must stay
    under derate(T) * P in every hour, using the temperatures that schedule itself produces."""
    w, p = august_slice()
    b = dp.Battery()
    r = coupled.run(rad.DARK_PAINT, w, p, b, inverter_heat="cells")
    assert r.hours_derated > 20  # the case actually exercises the caps
    allowed = enc.derate_fraction(r.thermal.T_batt) * b.power_kw
    assert np.all(np.maximum(r.dispatch.charge_kw, r.dispatch.discharge_kw) <= allowed + 1e-3)
    assert r.revenue_usd < r.uncapped.revenue_usd
    # the temperatures reported are the ones this dispatch produces
    q_cell, q_skin = coupled.heat_inputs(r.dispatch, b, "cells")
    again = enc.simulate(rad.DARK_PAINT, w, q_cell, q_skin_w=q_skin)
    assert again.T_batt == pytest.approx(r.thermal.T_batt, abs=1e-12)


def test_reusing_unchanged_days_gives_the_same_dispatch():
    w, p = august_slice(days=6)
    b = dp.Battery()
    free = dp.dispatch(p, b)
    cap = np.full(len(p), b.power_kw)
    cap[30:40] = 2.0  # only the second local day changes
    fresh = dp.dispatch(p, b, power_cap_kw=cap)
    reused = dp.dispatch(p, b, power_cap_kw=cap, reuse=(None, free))
    assert reused.revenue_usd == pytest.approx(fresh.revenue_usd, abs=1e-9)


# ---------------------------------------------------------------- the sky against measured infrared (NOAA SURFRAD)

@pytest.fixture(scope="module")
def surfrad_hours():
    return surfrad.load()


def test_clear_sky_emissivity_matches_measured_longwave(surfrad_hours):
    """Berdahl & Martin with each station's own temperature and humidity, on hours Open-Meteo calls clear.
    Measured: biases of -4 to -15 W/m^2 (the fit reads a little low, most in winter), RMS 9 to 18."""
    d = surfrad_hours
    stats = surfrad.error_stats(d, surfrad.modeled_longwave(d, clouds="none"))
    for key, row in stats.items():
        assert row["clear"]["n"] > 200, key
        assert -16.0 < row["clear"]["bias"] < 5.0, key
        assert row["clear"]["rms"] < 20.0, key


def test_cloud_layers_fitted_in_january_hold_up_in_july(surfrad_hours):
    """Driven by Open-Meteo inputs the way the enclosure is, July (never used in the fit) comes out within
    8 W/m^2 of the pyrgeometers on average at both sites, with an RMS under 15 W/m^2 and well below what my
    first model (total cloud as black cloud at air temperature) gave."""
    d = surfrad_hours
    new = surfrad.error_stats(d, surfrad.modeled_longwave(d, inputs="open-meteo"))
    old = surfrad.error_stats(d, surfrad.modeled_longwave(d, inputs="open-meteo", clouds="black_total"))
    for key in ("gwn_07", "dra_07"):
        assert abs(new[key]["all"]["bias"]) < 8.0, key
        assert new[key]["all"]["rms"] < 15.0, key
        assert new[key]["all"]["rms"] < 0.7 * old[key]["all"]["rms"], key


def test_cloud_constants_are_the_january_fit(surfrad_hours):
    d = surfrad_hours
    k = surfrad.fit_cloud_effects(d[d["month"] == 1].reset_index(drop=True))
    assert k == pytest.approx(rad.CLOUD_EFFECT, abs=0.01)
