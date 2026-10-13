"""Task 3 Execution Handler: Water Decanter (Condensate) System Sizing.

Calculates:
1. Seasonal condensate generation per zone and district for all six conditions.
2. The governing condition for condensate collection.
3. Central decanter holding-tank volume using a 12-hour retention assumption.
"""

import pandas as pd

from models.zone_psychrometrics import ZonePsychrometrics
from models.decanter import DecanterModel


N_ZONES = 400
FRESH_AIR_FRAC = 0.30
TANK_RETENTION_HOURS = 12.0

AMBIENT_CONDITIONS = [
    {"name": "Summer Day", "T_db": 40.0, "rh": 50.0, "solar": True},
    {"name": "Summer Night", "T_db": 30.0, "rh": 50.0, "solar": False},
    {"name": "Monsoon Day", "T_db": 35.0, "rh": 85.0, "solar": True},
    {"name": "Monsoon Night", "T_db": 30.0, "rh": 85.0, "solar": False},
    {"name": "Winter Day", "T_db": 28.0, "rh": 50.0, "solar": True},
    {"name": "Winter Night", "T_db": 16.0, "rh": 50.0, "solar": False},
]


def run_task3():
    """Calculate condensate rates and size the district decanter tank."""

    results = []

    for cond in AMBIENT_CONDITIONS:
        psych = ZonePsychrometrics.solve_condition(
            cond["T_db"],
            cond["rh"],
            cond["solar"],
            fresh_air_frac=FRESH_AIR_FRAC,
        )

        # ZonePsychrometrics supplies the district total for the 400-zone case.
        m_zone_kg_s = psych["m_condensate_zone_kg_s"]
        m_district_kg_s = psych["m_condensate_district_kg_s"]

        rates = DecanterModel.calculate_condensate_rates(m_district_kg_s)

        results.append({
            "Condition": cond["name"],
            "Outdoor DB (°C)": cond["T_db"],
            "Outdoor RH (%)": cond["rh"],
            "Condensate per zone (kg/h)": m_zone_kg_s * 3600.0,
            "District condensate (kg/h)": rates["mass_rate_kg_h"],
            "District flow (m³/h)": rates["vol_rate_m3_h"],
            "District volume (m³/day)": rates["daily_vol_m3_day"],
        })

    df_results = pd.DataFrame(results)

    # Select the governing condition using the largest daily condensate volume.
    governing_index = df_results["District volume (m³/day)"].idxmax()
    gov_cond_decanter = df_results.loc[governing_index, "Condition"]
    max_daily_vol = float(
        df_results.loc[governing_index, "District volume (m³/day)"]
    )

    # Explicitly state that 12 hours is the selected design retention period.
    tank_sizing = DecanterModel.size_decanter_tank(
        max_daily_vol_m3_day=max_daily_vol,
        retention_hours=TANK_RETENTION_HOURS,
    )

    # Round only for presentation; retain unrounded values in the returned table.
    display_columns = [
        "Condensate per zone (kg/h)",
        "District condensate (kg/h)",
        "District flow (m³/h)",
        "District volume (m³/day)",
    ]
    df_results[display_columns] = df_results[display_columns].round(3)

    return df_results, gov_cond_decanter, tank_sizing
