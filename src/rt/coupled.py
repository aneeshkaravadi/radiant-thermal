"""Battery dispatch and enclosure temperature, solved together.

The battery management system cuts power when the cells are hot, and the power it allows changes
how much heat the cells make. So I look for a dispatch that never asks for more power than the cap
implied by its own temperatures:

    1. dispatch with no caps and simulate the enclosure
    2. cap_k = min(cap_(k-1), derate(T_k) * P)        (caps only ever tighten)
    3. re-dispatch under cap_k, re-simulate, and stop once no hour asks for more than derate(T) * P

Plain re-dispatching (cap_k = derate(T_k) * P) can flip back and forth between two schedules
instead of settling. With caps that only tighten it can't, and it stops at a schedule that respects
its own temperatures. That schedule is feasible but not necessarily the best one, so its revenue is a
lower bound on what a perfect-foresight battery could earn, and the uncapped dispatch is the upper bound.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from . import dispatch as dp
from . import enclosure as enc
from .radiative import Coating

INVERTER_HEAT = ("skin", "cells", "outside")


@dataclass
class CoupledResult:
    dispatch: dp.DispatchResult
    thermal: enc.ThermalResult
    uncapped: dp.DispatchResult
    passes: int

    @property
    def revenue_usd(self) -> float:
        return self.dispatch.revenue_usd

    @property
    def hours_derated(self) -> int:
        return int((enc.derate_fraction(self.thermal.T_batt) < 1).sum())


def heat_inputs(d: dp.DispatchResult, b: dp.Battery, inverter_heat: str = "skin"):
    """(heat into the cells, heat into the skin), W, with the inverter's loss sent where ``inverter_heat`` says."""
    if inverter_heat not in INVERTER_HEAT:
        raise ValueError(f"inverter_heat must be one of {INVERTER_HEAT}")
    cells, inv = d.cell_heat_w(b), d.inverter_heat_w(b)
    if inverter_heat == "cells":
        return cells + inv, None
    return cells, (inv if inverter_heat == "skin" else None)


def run(coating: Coating, weather: enc.Weather, prices: pd.Series, battery: dp.Battery | None = None,
        enclosure: enc.Enclosure | None = None, fan: enc.Fan | None = None, spectral: bool = False,
        inverter_heat: str = "skin", substeps: int = 60, max_passes: int = 30, tol_kw: float = 1e-3) -> CoupledResult:
    b = dp.Battery() if battery is None else battery
    uncapped = dp.dispatch(prices, b)
    d, cap = uncapped, None
    for k in range(1, max_passes + 1):
        q_cell, q_skin = heat_inputs(d, b, inverter_heat)
        th = enc.simulate(coating, weather, q_cell, enclosure, substeps=substeps, fan=fan, spectral=spectral,
                          q_skin_w=q_skin)
        own = enc.derate_fraction(th.T_batt) * b.power_kw
        if np.all(np.maximum(d.charge_kw, d.discharge_kw) <= own + tol_kw):
            return CoupledResult(d, th, uncapped, k)
        new_cap = own if cap is None else np.minimum(cap, own)
        d, cap = dp.dispatch(prices, b, power_cap_kw=new_cap, reuse=(cap, d)), new_cap
    raise RuntimeError(f"thermal-dispatch coupling did not settle in {max_passes} passes")
