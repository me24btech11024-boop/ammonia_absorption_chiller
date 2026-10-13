"""
Task 4 Fresh-Air Sensitivity Handler.

Compares 30% and 40% outdoor-air fractions for all six ambient conditions.
Reports changes in cooling load, cycle duties, refrigerant/flue-gas flow,
and condensate generation.

Condenser temperature uses the first-pass assignment estimate:
T_cond = outdoor dry-bulb temperature + 12 K.
"""

import math
import pandas as pd

from models.zone_psychrometrics import ZonePsychrometrics
from models.absorption_cycle import AbsorptionCycle
from models.decanter import DecanterModel


FRESH_AIR_BASELINE = 0.30
FRESH_AIR_INCREASED = 0.40
ROOM_RH_LIMIT_PERCENT = 75.0
CONDENSER_APPROACH_K = 12.0

AMBIENT_CONDITIONS = [
    {"name": "Summer Day", "T_db": 40.0, "rh": 50.0, "solar": True},
    {"name": "Summer Night", "T_db": 30.0, "rh": 50.0, "solar": False},
    {"name": "Monsoon Day", "T_db": 35.0, "rh": 85.0, "solar": True},
    {"name": "Monsoon Night", "T_db": 30.0, "rh": 85.0, "solar": False},
    {"name": "Winter Day", "T_db": 28.0, "rh": 50.0, "solar": True},
    {"name": "Winter Night", "T_db": 16.0, "rh": 50.0, "solar": False},
]


def _percent_change(baseline, increased):
    """Return percentage change, or NaN if the baseline is effectively zero."""
    if not math.isfinite(baseline) or not math.isfinite(increased):
        return float("nan")
    if abs(baseline) < 1e-12:
        return float("nan")
    return 100.0 * (increased - baseline) / abs(baseline)


def _solve_case(condition, fresh_air_frac):
    """Calculate the zone load, cycle performance, and condensate rate."""
    psych = ZonePsychrometrics.solve_condition(
        condition["T_db"],
        condition["rh"],
        condition["solar"],
        fresh_air_frac=fresh_air_frac,
    )

    T_cond_C = condition["T_db"] + CONDENSER_APPROACH_K

    cycle = AbsorptionCycle.solve_cycle(
        psych["Q_evap_chiller_kW"],
        T_cond_C=T_cond_C,
        T_evap_C=2.0,
    )

    condensate = DecanterModel.calculate_condensate_rates(
        psych["m_condensate_district_kg_s"]
    )

    return {
        "psych": psych,
        "cycle": cycle,
        "condensate": condensate,
        "T_cond_C": T_cond_C,
    }


def run_task4():
    """Run the 30%-versus-40% fresh-air sensitivity analysis."""
    results = []

    for condition in AMBIENT_CONDITIONS:
        baseline = _solve_case(condition, FRESH_AIR_BASELINE)
        increased = _solve_case(condition, FRESH_AIR_INCREASED)

        p30 = baseline["psych"]
        p40 = increased["psych"]
        c30 = baseline["cycle"]
        c40 = increased["cycle"]
        d30 = baseline["condensate"]
        d40 = increased["condensate"]

        results.append({
            "Condition": condition["name"],
            "Outdoor DB (°C)": condition["T_db"],
            "Outdoor RH (%)": condition["rh"],

            "Q_evap 30% (kW)": p30["Q_evap_chiller_kW"],
            "Q_evap 40% (kW)": p40["Q_evap_chiller_kW"],
            "ΔQ_evap (%)": _percent_change(
                p30["Q_evap_chiller_kW"],
                p40["Q_evap_chiller_kW"],
            ),

            "Q_gen 30% (kW)": c30["Q_gen_kW"],
            "Q_gen 40% (kW)": c40["Q_gen_kW"],
            "ΔQ_gen (%)": _percent_change(
                c30["Q_gen_kW"], c40["Q_gen_kW"]
            ),

            "Refrigerant flow 30% (kg/s)": c30["m_r_kg_s"],
            "Refrigerant flow 40% (kg/s)": c40["m_r_kg_s"],
            "Δm_r (%)": _percent_change(
                c30["m_r_kg_s"], c40["m_r_kg_s"]
            ),

            "Flue gas 30% (kg/s)": c30["m_flue_gas_kg_s"],
            "Flue gas 40% (kg/s)": c40["m_flue_gas_kg_s"],
            "ΔFlue gas (%)": _percent_change(
                c30["m_flue_gas_kg_s"], c40["m_flue_gas_kg_s"]
            ),

            "Condensate 30% (kg/h)": d30["mass_rate_kg_h"],
            "Condensate 40% (kg/h)": d40["mass_rate_kg_h"],
            "ΔCondensate (%)": _percent_change(
                d30["mass_rate_kg_h"], d40["mass_rate_kg_h"]
            ),

            "Room RH 30% (%)": p30.get("room_rh_percent", float("nan")),
            "Room RH 40% (%)": p40.get("room_rh_percent", float("nan")),
            "RH violation 30%": p30.get("rh_violated", None),
            "RH violation 40%": p40.get("rh_violated", None),
            "Room RH limit (%)": ROOM_RH_LIMIT_PERCENT,
            "T_cond assumed (°C)": baseline["T_cond_C"],
        })

    return pd.DataFrame(results)
