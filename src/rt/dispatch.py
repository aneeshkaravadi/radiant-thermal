"""Price-arbitrage dispatch of a battery against real ERCOT day-ahead prices.

For each day, a linear program picks hourly charge c_h and discharge d_h in [0, P]:

    maximize   sum_h price_h * (d_h - c_h)                      ($ per MWh x MWh)
    subject to soc_{h+1} = soc_h + eta * c_h - d_h / eta,   0 <= soc_h <= E
               soc_0 = soc_24 = E / 2

eta is the one-way efficiency (round trip eta^2). This is perfect-foresight
arbitrage: an upper bound on energy-only revenue, not a trading strategy.

The one-way loss is split between the cells and the inverter, eta = eta_cell * eta_inv,
because only the cells' share heats the cells directly:

    charging c (AC in):     inverter loses c (1 - eta_inv),   cells lose c eta_inv (1 - eta_cell)
    discharging d (AC out): inverter loses d (1/eta_inv - 1), cells lose (d / eta_inv) (1/eta_cell - 1)

which adds up to c (1 - eta) + d (1/eta - 1), everything that goes in minus everything that comes out.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.optimize import linprog


# Mean battery round-trip efficiency (energy out over energy in at the battery's own terminals) of the
# 22 home batteries HTW Berlin measured in its Stromspeicher-Inspektion 2025, p. 25: 95.7 %.
CELL_ROUND_TRIP_HTW_2025 = 0.957


@dataclass
class Battery:
    energy_kwh: float = 20.0
    power_kw: float = 5.0
    eta_one_way: float = 0.95  # whole system, AC to stored energy; round trip about 90 %
    cell_round_trip: float = CELL_ROUND_TRIP_HTW_2025  # the cells alone; the rest of the loss is the inverter's

    def __post_init__(self):
        if not self.eta_one_way**2 - 1e-12 <= self.cell_round_trip <= 1.0:
            raise ValueError("cell round trip must lie between the system's round trip and 1")

    @property
    def eta_cell(self) -> float:
        return float(np.sqrt(self.cell_round_trip))

    @property
    def eta_inverter(self) -> float:
        return min(1.0, self.eta_one_way / self.eta_cell)


@dataclass
class DispatchResult:
    time: pd.DatetimeIndex
    charge_kw: np.ndarray
    discharge_kw: np.ndarray
    soc_kwh: np.ndarray
    price: np.ndarray  # $/MWh

    @property
    def revenue_usd(self) -> float:
        return float(np.sum(self.price * (self.discharge_kw - self.charge_kw)) / 1000.0)

    def cell_heat_w(self, b: Battery) -> np.ndarray:
        """Heat released inside the cells in each hour, W."""
        ei, ec = b.eta_inverter, b.eta_cell
        return (self.charge_kw * ei * (1 - ec) + self.discharge_kw / ei * (1 / ec - 1)) * 1000.0

    def inverter_heat_w(self, b: Battery) -> np.ndarray:
        """Heat released in the inverter and power electronics in each hour, W."""
        ei = b.eta_inverter
        return (self.charge_kw * (1 - ei) + self.discharge_kw * (1 / ei - 1)) * 1000.0


def _solve_day(price: np.ndarray, b: Battery, power_cap: np.ndarray | None = None):
    n = len(price)
    cap = np.full(n, b.power_kw) if power_cap is None else np.minimum(b.power_kw, power_cap)
    # variables: c (n), d (n);  soc_k = E/2 + sum_{j<k} (eta c_j - d_j/eta)
    cost = np.concatenate([price, -price])  # minimize price*c - price*d
    L = np.tril(np.ones((n, n)))
    G = np.hstack([b.eta_one_way * L, -L / b.eta_one_way])  # soc_k - E/2 for k = 1..n
    A_ub = np.vstack([G, -G])
    b_ub = np.concatenate([np.full(n, b.energy_kwh / 2), np.full(n, b.energy_kwh / 2)])
    A_eq = G[-1:, :]
    res = linprog(cost, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[0.0],
                  bounds=[(0, c) for c in cap] * 2, method="highs")
    c, d = res.x[:n], res.x[n:]
    soc = b.energy_kwh / 2 + G @ res.x
    return c, d, soc


def dispatch(prices: pd.Series, b: Battery | None = None, power_cap_kw: np.ndarray | None = None,
             reuse: tuple[np.ndarray | None, DispatchResult] | None = None) -> DispatchResult:
    """``prices`` indexed by hour start (tz-aware), $/MWh. ``power_cap_kw`` optionally limits each hour.

    ``reuse=(previous_cap, previous_result)`` copies the previous solution for every day whose caps
    haven't changed, since the LP for that day is then identical. ``previous_cap`` None means uncapped.
    """
    b = Battery() if b is None else b
    t = pd.DatetimeIndex(prices.index)
    p = prices.to_numpy(float)
    n = len(p)
    C, D, S = np.zeros_like(p), np.zeros_like(p), np.zeros_like(p)

    def caps(c):
        return np.full(n, b.power_kw) if c is None else np.minimum(b.power_kw, np.asarray(c, float))

    cap = caps(power_cap_kw)
    prev_cap = None if reuse is None else caps(reuse[0])
    local_day = t.tz_convert("America/Chicago").date if t.tz is not None else t.date
    days = pd.Series(np.arange(n)).groupby(local_day)
    for _, idx in days:
        i = idx.to_numpy()
        if reuse is not None and np.array_equal(cap[i], prev_cap[i]):
            prev = reuse[1]
            C[i], D[i], S[i] = prev.charge_kw[i], prev.discharge_kw[i], prev.soc_kwh[i]
        else:
            C[i], D[i], S[i] = _solve_day(p[i], b, cap[i])
    return DispatchResult(t, C, D, S, p)


def load_prices(path) -> pd.Series:
    d = pd.read_csv(path)
    return pd.Series(d["price_usd_mwh"].to_numpy(float), index=pd.DatetimeIndex(pd.to_datetime(d["hour_start_utc"], utc=True)))
