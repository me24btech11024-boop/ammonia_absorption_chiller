"""
Task 2 Condenser Selection, Sizing & Off-Design Performance Mapper.
"""

import pandas as pd
from models.zone_psychrometrics import ZonePsychrometrics
from models.absorption_cycle import AbsorptionCycle
from models.condenser import CondenserModel


def run_task2():
    conditions = [
        {"name": "Summer Day", "T_db": 40.0, "rh": 50.0, "solar": True},
        {"name": "Summer Night", "T_db": 30.0, "rh": 50.0, "solar": False},
        {"name": "Monsoon Day", "T_db": 35.0, "rh": 85.0, "solar": True},
        {"name": "Monsoon Night", "T_db": 30.0, "rh": 85.0, "solar": False},
        {"name": "Winter Day", "T_db": 28.0, "rh": 50.0, "solar": True},
        {"name": "Winter Night", "T_db": 16.0, "rh": 50.0, "solar": False},
    ]

    # Find governing condition (Summer Day)
    gov_cond = conditions[0]
    psych_gov = ZonePsychrometrics.solve_condition(gov_cond["T_db"], gov_cond["rh"], gov_cond["solar"])
    cycle_gov = AbsorptionCycle.solve_cycle(Q_evap_kW=psych_gov["Q_evap_chiller_kW"], T_cond_C=gov_cond["T_db"] + 12.0)
    Q_cond_gov_kW = psych_gov["Q_evap_chiller_kW"] + cycle_gov["Q_gen_kW"] + cycle_gov["W_pump_kW"]

    # Size Evaporative Condenser (U = 600 W/m^2*K)
    design_sizing = CondenserModel.size_design_condenser(
        Q_cond_kW=Q_cond_gov_kW,
        T_db_C=gov_cond["T_db"],
        rh_percent=gov_cond["rh"],
        condenser_type="evaporative",
        approach_K=8.0,
    )

    UA_design_kW_K = design_sizing["UA_kW_K"]

    # Off-design mapping
    off_design_results = []
    for cond in conditions:
        psych = ZonePsychrometrics.solve_condition(cond["T_db"], cond["rh"], cond["solar"])
        Q_evap = psych["Q_evap_chiller_kW"]
        
        T_cond_est = cond["T_db"] + 10.0
        for _ in range(10):
            cycle = AbsorptionCycle.solve_cycle(Q_evap_kW=Q_evap, T_cond_C=T_cond_est)
            Q_cond_kW = Q_evap + cycle["Q_gen_kW"] + cycle["W_pump_kW"]
            
            off_design = CondenserModel.solve_off_design(
                Q_cond_kW=Q_cond_kW,
                UA_kW_K=UA_design_kW_K,
                T_db_C=cond["T_db"],
                rh_percent=cond["rh"],
                condenser_type="evaporative",
            )
            T_cond_est = off_design["T_cond_C"]

        off_design_results.append({
            "Condition": cond["name"],
            "T_cond (°C)": round(off_design["T_cond_C"], 2),
            "P_high (bar)": round(off_design["P_high_bar"], 2),
            "T_gen (°C)": round(cycle["T_gen_C"], 2),
            "COP": round(cycle["COP"], 3),
            "P_high <= 20bar": off_design["check_p_high_pass"],
            "T_cond >= 25C": off_design["check_t_cond_min_pass"],
        })

    return design_sizing, pd.DataFrame(off_design_results)
