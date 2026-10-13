"""
Task 3 Execution Handler: Water Decanter (Condensate) System Sizing.

Calculates:
1. Seasonal condensate generation per zone and for the district across all 6 conditions[cite: 10, 11].
2. Identifies the governing condition for condensate collection[cite: 11].
3. Central decanter/holding tank volume sizing for water recovery[cite: 10, 11].
"""

import pandas as pd
from models.zone_psychrometrics import ZonePsychrometrics
from models.decanter import DecanterModel


def run_task3():
    """
    Executes Task 3 condensate generation analysis and decanter tank sizing[cite: 10, 11].
    """
    conditions = [
        {"name": "Summer Day", "T_db": 50.0, "rh": 40.0, "solar": True},
        {"name": "Summer Night", "T_db": 50.0, "rh": 30.0, "solar": False},
        {"name": "Monsoon Day", "T_db": 85.0, "rh": 35.0, "solar": True},
        {"name": "Monsoon Night", "T_db": 85.0, "rh": 30.0, "solar": False},
        {"name": "Winter Day", "T_db": 50.0, "rh": 28.0, "solar": True},
        {"name": "Winter Night", "T_db": 50.0, "rh": 16.0, "solar": False},
    ]

    results = []
    max_daily_vol = 0.0
    gov_cond_decanter = ""

    for cond in conditions:
        psych = ZonePsychrometrics.solve_condition(cond["T_db"], cond["rh"], cond["solar"])
        m_condensate = psych["m_condensate_district_kg_s"]
        
        rates = DecanterModel.calculate_condensate_rates(m_condensate)
        
        if rates["daily_vol_m3_day"] > max_daily_vol:
            max_daily_vol = rates["daily_vol_m3_day"]
            gov_cond_decanter = cond["name"]

        results.append({
            "Condition": cond["name"],
            "m_cond_zone (kg/h)": round(psych["m_condensate_zone_kg_s"] * 3600.0, 2),
            "m_cond_district (kg/h)": round(rates["mass_rate_kg_h"], 1),
            "Vol Rate (m3/h)": round(rates["vol_rate_m3_h"], 2),
            "Daily Vol (m3/day)": round(rates["daily_vol_m3_day"], 1),
        })

    df_results = pd.DataFrame(results)

    # Size central decanter storage tank for 12 hours of retention buffer[cite: 10, 11]
    tank_sizing = DecanterModel.size_decanter_tank(max_daily_vol_m3_day=max_daily_vol, retention_hours=12.0)

    return df_results, gov_cond_decanter, tank_sizing