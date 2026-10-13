
"""Task 1: Seasonal load analysis and ammonia circuit sizing."""

import pandas as pd
from models.zone_psychrometrics import ZonePsychrometrics
from models.absorption_cycle import AbsorptionCycle


def run_task1():
    conditions = [
        {"name": "Summer Day", "T_db": 40.0, "rh": 50.0, "solar": True},
        {"name": "Summer Night", "T_db": 30.0, "rh": 50.0, "solar": False},
        {"name": "Monsoon Day", "T_db": 35.0, "rh": 85.0, "solar": True},
        {"name": "Monsoon Night", "T_db": 30.0, "rh": 85.0, "solar": False},
        {"name": "Winter Day", "T_db": 28.0, "rh": 50.0, "solar": True},
        {"name": "Winter Night", "T_db": 16.0, "rh": 50.0, "solar": False},
    ]

    results = []
    psych_results = []

    for cond in conditions:
        psych = ZonePsychrometrics.solve_condition(
            T_db_out=cond["T_db"],
            rh_out=cond["rh"],
            solar_on=cond["solar"],
            fresh_air_frac=0.30,
        )
        psych_results.append((cond, psych))

    max_load = max(
        (float(p["Q_evap_chiller_kW"]) for _, p in psych_results),
        default=0.0,
    )

    for cond, psych in psych_results:
        q_load = float(psych["Q_evap_chiller_kW"])
        t_cond = cond["T_db"] + 12.0
        row = {
            "Condition": cond["name"],
            "T_cond_est_C": t_cond,
            "Q_evap_load_kW": q_load,
            "P_high_bar": float("nan"),
            "Q_gen_kW": float("nan"),
            "COP": float("nan"),
            "m_refrigerant_kg_s": float("nan"),
            "m_flue_gas_kg_s": float("nan"),
            "Cycle_feasible": False,
            "Cycle_status": "",
            "Turndown_percent": (
                round(100.0 * q_load / max_load, 1)
                if max_load > 0 else 0.0
            ),
        }

        try:
            cycle = AbsorptionCycle.solve_cycle(
                Q_evap_kW=q_load,
                T_cond_C=t_cond,
                T_evap_C=2.0,
            )
            row.update({
                "P_high_bar": cycle["P_high_bar"],
                "Q_gen_kW": cycle["Q_gen_kW"],
                "COP": cycle["COP"],
                "m_refrigerant_kg_s": cycle["m_r_kg_s"],
                "m_flue_gas_kg_s": cycle["m_flue_gas_kg_s"],
                "Cycle_feasible": cycle["check_p_high_pass"],
                "Cycle_status": cycle.get("model_status", "Calculated"),
            })
        except (ValueError, ArithmeticError, KeyError) as exc:
            row["Cycle_status"] = f"INFEASIBLE: {exc}"

        results.append(row)

    df = pd.DataFrame(results)
    feasible = df[df["Cycle_feasible"]]

    governing = (
        feasible.loc[feasible["Q_evap_load_kW"].idxmax(), "Condition"]
        if not feasible.empty
        else "No feasible cycle condition"
    )

    return df, governing
