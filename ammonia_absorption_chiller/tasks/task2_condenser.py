"""
Task 2: Condenser selection, sizing, and off-design performance.

Compares air-cooled and evaporative condenser design estimates and
evaluates off-design operating conditions across all six seasons.

The condenser model uses a simplified LMTD approximation; results
should be treated as provisional until the cooling-medium energy
balance is independently closed.
"""

import math
import pandas as pd

from models.zone_psychrometrics import ZonePsychrometrics
from models.absorption_cycle import AbsorptionCycle
from models.condenser import CondenserModel


AMBIENT_CONDITIONS = [
    {"name": "Summer Day", "T_db": 40.0, "rh": 50.0, "solar": True},
    {"name": "Summer Night", "T_db": 30.0, "rh": 50.0, "solar": False},
    {"name": "Monsoon Day", "T_db": 35.0, "rh": 85.0, "solar": True},
    {"name": "Monsoon Night", "T_db": 30.0, "rh": 85.0, "solar": False},
    {"name": "Winter Day", "T_db": 28.0, "rh": 50.0, "solar": True},
    {"name": "Winter Night", "T_db": 16.0, "rh": 50.0, "solar": False},
]


def _solve_zone_load(condition):
    """Return the psychrometric result for one ambient condition."""
    return ZonePsychrometrics.solve_condition(
        T_db_out=condition["T_db"],
        rh_out=condition["rh"],
        solar_on=condition["solar"],
        fresh_air_frac=0.30,
    )


def _solve_cycle_at_temperature(Q_evap_kW, T_cond_C):
    """Solve the absorption cycle at a specified condensing temperature."""
    return AbsorptionCycle.solve_cycle(
        Q_evap_kW=Q_evap_kW,
        T_cond_C=T_cond_C,
        T_evap_C=2.0,
    )


def _get_condenser_duty(Q_evap_kW, cycle):
    """Estimate condenser duty using the overall cycle energy balance."""
    return (
        Q_evap_kW
        + cycle["Q_gen_kW"]
        + cycle["W_pump_kW"]
    )


def _solve_coupled_off_design(
    condition,
    Q_evap_kW,
    UA_design_kW_K,
    condenser_type,
    max_iterations=50,
    tolerance_C=0.05,
):
    """
    Iterate between cycle duty and condenser off-design temperature.

    Convergence requires the outer condensing-temperature estimate to
    stabilize and the condenser solver to report convergence.
    """
    T_cond_est = condition["T_db"] + 12.0
    previous_T_cond = None
    cycle = None
    condenser = None

    for iteration in range(1, max_iterations + 1):
        cycle = _solve_cycle_at_temperature(
            Q_evap_kW=Q_evap_kW,
            T_cond_C=T_cond_est,
        )

        Q_cond_kW = _get_condenser_duty(Q_evap_kW, cycle)

        condenser = CondenserModel.solve_off_design(
            Q_cond_kW=Q_cond_kW,
            UA_kW_K=UA_design_kW_K,
            T_db_C=condition["T_db"],
            rh_percent=condition["rh"],
            condenser_type=condenser_type,
        )

        if not condenser["converged"]:
            raise RuntimeError(
                "Condenser solver did not converge for "
                f"{condition['name']} ({condenser_type})."
            )

        T_cond_new = condenser["T_cond_C"]

        if (
            previous_T_cond is not None
            and abs(T_cond_new - T_cond_est) <= tolerance_C
        ):
            # Re-solve the cycle at the final reported temperature.
            T_cond_est = T_cond_new
            cycle = _solve_cycle_at_temperature(
                Q_evap_kW=Q_evap_kW,
                T_cond_C=T_cond_est,
            )

            return cycle, condenser, iteration

        previous_T_cond = T_cond_est
        T_cond_est = T_cond_new

    raise RuntimeError(
        f"Cycle/condenser iteration did not converge for "
        f"{condition['name']} ({condenser_type}) after "
        f"{max_iterations} iterations."
    )


def run_task2():
    """Run condenser design comparisons and six-condition off-design map."""
    zone_results = []

    for condition in AMBIENT_CONDITIONS:
        psych = _solve_zone_load(condition)
        zone_results.append((condition, psych))

    # Use the maximum calculated evaporator duty as the preliminary
    # chiller governing condition for condenser sizing.
    governing_condition, governing_psych = max(
        zone_results,
        key=lambda item: item[1]["Q_evap_chiller_kW"],
    )

    Q_evap_design_kW = governing_psych["Q_evap_chiller_kW"]

    # First-pass condenser temperature per the assignment.
    T_cond_design_est = governing_condition["T_db"] + 12.0

    design_cycle = _solve_cycle_at_temperature(
        Q_evap_kW=Q_evap_design_kW,
        T_cond_C=T_cond_design_est,
    )

    Q_cond_design_kW = _get_condenser_duty(
        Q_evap_design_kW,
        design_cycle,
    )

    # Compare both condenser types at the governing cooling-load condition.
    design_results = {}

    for condenser_type, approach_K in (
        ("air_cooled", 12.0),
        ("evaporative", 8.0),
    ):
        design_results[condenser_type] = (
            CondenserModel.size_design_condenser(
                Q_cond_kW=Q_cond_design_kW,
                T_db_C=governing_condition["T_db"],
                rh_percent=governing_condition["rh"],
                condenser_type=condenser_type,
                approach_K=approach_K,
            )
        )

    # Select the evaporative option for the off-design map.
    # The final condenser-type selection should be justified from
    # area, operating constraints, water use, and site requirements.
    selected_type = "evaporative"
    selected_design = design_results[selected_type]
    UA_design_kW_K = selected_design["UA_kW_K"]

    off_design_results = []

    for condition, psych in zone_results:
        Q_evap_kW = psych["Q_evap_chiller_kW"]

        if Q_evap_kW <= 0:
            off_design_results.append({
                "Condition": condition["name"],
                "Condenser type": selected_type,
                "Q_evap (kW)": Q_evap_kW,
                "T_cond (°C)": None,
                "P_high (bar)": None,
                "T_gen (°C)": None,
                "COP": None,
                "Condenser converged": False,
                "Condenser energy residual (kW)": None,
                "P_high <= 20 bar": None,
                "T_cond >= 25 C": None,
                "Flue gas <= 80 kg/s": None,
                "Status": "No positive cooling load",
            })
            continue

        cycle, condenser, iterations = _solve_coupled_off_design(
            condition=condition,
            Q_evap_kW=Q_evap_kW,
            UA_design_kW_K=UA_design_kW_K,
            condenser_type=selected_type,
        )

        m_flue = cycle["m_flue_gas_kg_s"]

        off_design_results.append({
            "Condition": condition["name"],
            "Condenser type": selected_type,
            "Q_evap (kW)": round(Q_evap_kW, 2),
            "Q_cond (kW)": round(
                _get_condenser_duty(Q_evap_kW, cycle), 2
            ),
            "T_cond (°C)": round(condenser["T_cond_C"], 2),
            "P_high (bar)": round(condenser["P_high_bar"], 3),
            "T_gen (°C)": round(cycle["T_gen_C"], 2),
            "COP": round(cycle["COP"], 3),
            "Condenser converged": condenser["converged"],
            "Condenser iterations": condenser["iterations"],
            "Outer iterations": iterations,
            "Condenser energy residual (kW)": (
                round(condenser["energy_residual_kW"], 4)
                if condenser["energy_residual_kW"] is not None
                else None
            ),
            "P_high <= 20 bar": condenser["check_p_high_pass"],
            "T_cond >= 25 C": condenser["check_t_cond_min_pass"],
            "Flue gas (kg/s)": round(m_flue, 3),
            "Flue gas <= 80 kg/s": bool(m_flue <= 80.0),
            "Status": "Provisional simplified condenser model",
        })

    return (
        design_results,
        pd.DataFrame(off_design_results),
    )
