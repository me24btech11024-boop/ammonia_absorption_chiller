
"""
Task 2: Condenser design sizing and off-design performance.

This version reports infeasible cycle trial points instead of allowing
a pressure-limit exception to terminate the whole simulation.
"""

import math
import pandas as pd

from models.zone_psychrometrics import ZonePsychrometrics
from models.absorption_cycle import AbsorptionCycle
from models.condenser import condenser


CONDITIONS = [
    {"name": "Summer Day", "T_db": 40.0, "rh": 50.0, "solar": True},
    {"name": "Summer Night", "T_db": 30.0, "rh": 50.0, "solar": False},
    {"name": "Monsoon Day", "T_db": 35.0, "rh": 85.0, "solar": True},
    {"name": "Monsoon Night", "T_db": 30.0, "rh": 85.0, "solar": False},
    {"name": "Winter Day", "T_db": 28.0, "rh": 50.0, "solar": True},
    {"name": "Winter Night", "T_db": 16.0, "rh": 50.0, "solar": False},
]


def run_task2():
    """Return condenser design information and seasonal operating results."""
    seasonal = []

    for condition in CONDITIONS:
        psych = ZonePsychrometrics.solve_condition(
            T_db_out=condition["T_db"],
            rh_out=condition["rh"],
            solar_on=condition["solar"],
            fresh_air_frac=0.30,
        )
        seasonal.append({
            **condition,
            "Q_evap_kW": float(psych["Q_evap_chiller_kW"]),
        })

    # Select a provisional design load from the highest calculated load.
    design_condition = max(seasonal, key=lambda x: x["Q_evap_kW"])
    design_T_cond = design_condition["T_db"] + 8.0

    # The condenser module's exact public API may differ between versions.
    # Do not guess a sizing method: verify the module before invoking it.
    try:
        design_sizing = Condenser.size_condenser(
            Q_cond_kW=design_condition["Q_evap_kW"],
            T_ambient_C=design_condition["T_db"],
            T_cond_C=design_T_cond,
            condenser_type="evaporative",
        )
    except (AttributeError, TypeError, ValueError) as exc:
        raise RuntimeError(
            "Task 2 needs its condenser API verified in models/condenser.py. "
            f"The sizing call could not be completed: {exc}"
        ) from exc

    rows = []

    for condition in seasonal:
        T_cond = condition["T_db"] + 8.0
        row = {
            "Condition": condition["name"],
            "Q_evap_kW": condition["Q_evap_kW"],
            "T_cond_C": T_cond,
            "P_high_bar": float("nan"),
            "T_gen_C": float("nan"),
            "COP": float("nan"),
            "m_flue_gas_kg_s": float("nan"),
            "Feasible": False,
            "Status": "",
        }

        try:
            cycle = AbsorptionCycle.solve_cycle(
                Q_evap_kW=condition["Q_evap_kW"],
                T_cond_C=T_cond,
                T_evap_C=2.0,
            )
            row.update({
                "P_high_bar": cycle["P_high_bar"],
                "T_gen_C": cycle["T_gen_C"],
                "COP": cycle["COP"],
                "m_flue_gas_kg_s": cycle["m_flue_gas_kg_s"],
                "Feasible": cycle["check_p_high_pass"],
                "Status": cycle.get("model_status", "Calculated"),
            })
        except (ValueError, ArithmeticError, KeyError) as exc:
            row["Status"] = f"INFEASIBLE: {exc}"

        rows.append(row)

    return design_sizing, pd.DataFrame(rows)
