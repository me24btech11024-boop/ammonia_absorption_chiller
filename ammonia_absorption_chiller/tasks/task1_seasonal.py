"""
Task 1 Execution Handler: Seasonal Load Analysis & Ammonia Circuit Sizing.
"""

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
    max_Q_evap = 0.0

    # First pass: solve psychrometrics to find peak load for turndown ratio
    temp_results = []
    for cond in conditions:
        psych = ZonePsychrometrics.solve_condition(
            T_db_out=cond["T_db"],
            rh_out=cond["rh"],
            solar_on=cond["solar"],
            fresh_air_frac=0.30,
        )
        if psych["Q_evap_chiller_kW"] > max_Q_evap:
            max_Q_evap = psych["Q_evap_chiller_kW"]
        temp_results.append((cond, psych))

    # Second pass: evaluate unconstrained absorption cycle balances
    for cond, psych in temp_results:
        Q_evap = psych["Q_evap_chiller_kW"]
        T_cond_est = cond["T_db"] + 12.0
        
        cycle = AbsorptionCycle.solve_cycle(
            Q_evap_kW=Q_evap,
            T_cond_C=T_cond_est,
            T_evap_C=2.0,
        )

        turndown_ratio = (Q_evap / max_Q_evap) * 100.0 if max_Q_evap > 0 else 0.0
        turndown_ok = turndown_ratio >= 25.0 or Q_evap == 0.0

        results.append({
            "Condition": cond["name"],
            "T_mix (°C)": round(psych["T_mix"], 2),
            "ADP (°C)": round(psych["ADP"], 2),
            "w_supply (g/kg)": round(psych["w_supply"] * 1000, 3),
            "Room RH (%)": round(psych["rh_room"], 1),
            "RH Viol": psych["rh_violated"],
            "Q_evap (kW)": round(cycle["Q_evap_actual_kW"], 1),
            "Q_gen (kW)": round(cycle["Q_gen_kW"], 1),
            "COP": round(cycle["COP"], 3),
            "m_r (kg/s)": round(cycle["m_r_kg_s"], 2),
            "m_flue (kg/s)": round(cycle["m_flue_gas_kg_s"], 2),
            "Flue OK": cycle["check_flue_pass"],
            "Turndown (%)": round(turndown_ratio, 1),
            "Turndown OK": turndown_ok,
        })

    df_results = pd.DataFrame(results)
    governing_row = df_results.loc[df_results["Q_evap (kW)"].idxmax()]

    return df_results, governing_row["Condition"]
