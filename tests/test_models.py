from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from rt import dispatch as dp
from rt import enclosure as enc
from rt import radiative as rad, radiator

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
