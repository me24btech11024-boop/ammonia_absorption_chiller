"""
Task 4 Fresh-Air Sensitivity Handler (30% vs 40%).
"""

import pandas as pd
from models.zone_psychrometrics import ZonePsychrometrics
from models.absorption_cycle import AbsorptionCycle
from models.decanter import DecanterModel


def run_task4():
    conditions = [
        {"name": "Summer Day", "T_db": 50.0, "rh": 40.0, "solar": True},
        {"name": "Summer Night", "T_db": 50.0, "rh": 30.0, "solar": False},
        {"name": "Monsoon Day", "T_db": 35.0, "rh": 85.0, "solar": True},
        {"name": "Monsoon Night", "T_db": 30.0, "rh": 85.0, "solar": False},
        {"name": "Winter Day", "T_db": 28.0, "rh": 50.0, "solar": True},
        {"name": "Winter Night", "T_db": 16.0, "rh": 50.0, "solar": False},
    ]

    sensitivity_results = []

    for cond in conditions:
        p30 = ZonePsychrometrics.solve_condition(cond["T_db"], cond["rh"], cond["solar"], fresh_air_frac=0.30)
        c30 = AbsorptionCycle.solve_cycle(p30["Q_evap_chiller_kW"], T_cond_C=cond["T_db"] + 12.0)
        d30 = DecanterModel.calculate_condensate_rates(p30["m_condensate_district_kg_s"])

        p40 = ZonePsychrometrics.solve_condition(cond["T_db"], cond["rh"], cond["solar"], fresh_air_frac=0.40)
        c40 = AbsorptionCycle.solve_cycle(p40["Q_evap_chiller_kW"], T_cond_C=cond["T_db"] + 12.0)
        d40 = DecanterModel.calculate_condensate_rates(p40["m_condensate_district_kg_s"])

        # Evaporator & Generator percentage changes
        if p30["Q_evap_chiller_kW"] > 0:
            pct_Q_evap = ((p40["Q_evap_chiller_kW"] - p30["Q_evap_chiller_kW"]) / p30["Q_evap_chiller_kW"]) * 100.0
            pct_Q_gen = ((c40["Q_gen_kW"] - c30["Q_gen_kW"]) / c30["Q_gen_kW"]) * 100.0
            pct_m_r = ((c40["m_r_kg_s"] - c30["m_r_kg_s"]) / c30["m_r_kg_s"]) * 100.0
            pct_flue = ((c40["m_flue_gas_kg_s"] - c30["m_flue_gas_kg_s"]) / c30["m_flue_gas_kg_s"]) * 100.0
        else:
            pct_Q_evap = 0.0
            pct_Q_gen = 0.0
            pct_m_r = 0.0
            pct_flue = 0.0

        # Condensate percentage change
        if d30["mass_rate_kg_h"] > 1.0:
            pct_cond = ((d40["mass_rate_kg_h"] - d30["mass_rate_kg_h"]) / d30["mass_rate_kg_h"]) * 100.0
        else:
            pct_cond = 0.0

        sensitivity_results.append({
            "Condition": cond["name"],
            "ΔQ_evap (%)": round(pct_Q_evap, 2),
            "ΔQ_gen (%)": round(pct_Q_gen, 2),
            "Δm_r (%)": round(pct_m_r, 2),
            "ΔFlue Gas (%)": round(pct_flue, 2),
            "ΔCondensate (%)": round(pct_cond, 2),
            "RH Viol (40%)": p40["rh_violated"],
        })

    return pd.DataFrame(sensitivity_results)