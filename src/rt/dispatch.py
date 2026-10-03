"""Price-arbitrage dispatch of a battery against real ERCOT day-ahead prices.

For each day, a linear program picks hourly charge c_h and discharge d_h in [0, P]:

    maximize   sum_h price_h * (d_h - c_h)                      ($ per MWh x MWh)
    subject to soc_{h+1} = soc_h + eta * c_h - d_h / eta,   0 <= soc_h <= E
               soc_0 = soc_24 = E / 2

eta is the one-way efficiency (round trip eta^2). This is perfect-foresight
arbitrage: an upper bound on energy-only revenue, not a trading strategy.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.optimize import linprog


@dataclass
class Battery:
    energy_kwh: float = 20.0
    power_kw: float = 5.0
    eta_one_way: float = 0.95


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

    def heat_w(self, eta_one_way: float) -> np.ndarray:
        """Heat released in the cells: the one-way loss on whatever flows in or out."""
        return (1 - eta_one_way) * (self.charge_kw + self.discharge_kw) * 1000.0


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


def dispatch(prices: pd.Series, b: Battery = Battery(), power_cap_kw: np.ndarray | None = None) -> DispatchResult:
    """``prices`` indexed by hour start (tz-aware), $/MWh. ``power_cap_kw`` optionally limits each hour."""
    t = pd.DatetimeIndex(prices.index)
    p = prices.to_numpy(float)
    C, D, S = np.zeros_like(p), np.zeros_like(p), np.zeros_like(p)
    local_day = t.tz_convert("America/Chicago").date if t.tz is not None else t.date
    days = pd.Series(np.arange(len(p))).groupby(local_day)
    for _, idx in days:
        i = idx.to_numpy()
        cap = None if power_cap_kw is None else power_cap_kw[i]
        C[i], D[i], S[i] = _solve_day(p[i], b, cap)
    return DispatchResult(t, C, D, S, p)


def load_prices(path) -> pd.Series:
    d = pd.read_csv(path)
    return pd.Series(d["price_usd_mwh"].to_numpy(float), index=pd.DatetimeIndex(pd.to_datetime(d["hour_start_utc"], utc=True)))
