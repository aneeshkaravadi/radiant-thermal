from itertools import pairwise
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from rt import dispatch as dp
from rt import enclosure as enc
from rt import radiative as rad
from rt import radiator

ROOT = Path(__file__).resolve().parents[1]


def test_berdahl_martin_values():
    # Tdp = 0 C -> 0.711 ; Tdp = 20 C -> 0.711 + 0.112 + 0.0292
    assert rad.sky_emissivity(0.0) == pytest.approx(0.711)
    assert rad.sky_emissivity(20.0) == pytest.approx(0.711 + 0.56 * 0.2 + 0.73 * 0.04)
    assert rad.sky_emissivity(10.0, cloud_fraction=1.0) == pytest.approx(1.0)


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
                    np.full(n, 2.0), np.zeros(n))
    q = np.full(n, 200.0)
    e = enc.Enclosure()
    r = enc.simulate(rad.WHITE_PAINT, w, q, e)
    assert np.all(np.isfinite(r.T_batt))
    # battery time constant C_b/UA ~ 4.6 h, so after 48 h it is at steady state
    assert e.ua_internal * (r.T_batt[-1] - r.T_skin[-1]) == pytest.approx(200.0, rel=0.01)


def constant_weather(n, Ta, ghi=0.0, dew=15.0):
    t = pd.date_range("2025-07-01", periods=n, freq="h", tz="UTC")
    return enc.Weather(t, np.full(n, Ta), np.full(n, dew), np.full(n, ghi), np.full(n, 2.0), np.zeros(n))


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


def test_two_band_enclosure_matches_the_gray_one_for_a_gray_skin():
    w = constant_weather(48, 305.0, ghi=400.0, dew=18.0)
    q = np.full(48, 200.0)
    gray = enc.simulate(rad.WHITE_PAINT, w, q)
    two = enc.simulate(rad.WHITE_PAINT, w, q, spectral=True)
    assert two.T_batt == pytest.approx(gray.T_batt, abs=1e-9)


def test_arrhenius_reference():
    assert enc.arrhenius_factor(298.15) == pytest.approx(1.0)
    assert enc.arrhenius_factor(308.15) > 1.5  # ~1.9x per 10 K at 50 kJ/mol


def test_dispatch_respects_limits_and_is_profitable():
    t = pd.date_range("2025-08-01", periods=48, freq="h", tz="UTC")
    price = pd.Series(np.tile(np.r_[np.full(12, 20.0), np.full(12, 100.0)], 2), index=t)
    b = dp.Battery(20.0, 5.0, 0.95)
    r = dp.dispatch(price, b)
    assert r.soc_kwh.min() >= -1e-6 and r.soc_kwh.max() <= b.energy_kwh + 1e-6
    assert r.charge_kw.max() <= 5 + 1e-9 and r.discharge_kw.max() <= 5 + 1e-9
    assert r.revenue_usd > 0
    capped = dp.dispatch(price, b, power_cap_kw=np.full(48, 1.0))
    assert capped.revenue_usd < r.revenue_usd


@pytest.mark.parametrize("city, zone", [("dallas", "lz_north"), ("houston", "lz_houston")])
def test_real_data_files_align(city, zone):
    w = enc.Weather.from_csv(ROOT / f"data/weather/{city}_2025_hourly.csv")
    p = dp.load_prices(ROOT / f"data/ercot/dam_2025_{zone}.csv")
    assert len(w.Ta) == len(p) == 8760 and (w.time == p.index).all()


def test_radiator_closed_form():
    c = rad.Coating("c", 0.2, 0.9)
    A = radiator.area_per_kw(c, radiator.DEEP_SPACE, 300.0)
    assert A == pytest.approx(1000.0 / (0.9 * rad.SIGMA * 300.0**4))
    T_eq = radiator.equilibrium_temperature(c, radiator.Environment("sun normal", 1.0, 0.0))
    assert T_eq == pytest.approx((0.2 * 1361 / (0.9 * rad.SIGMA)) ** 0.25)
