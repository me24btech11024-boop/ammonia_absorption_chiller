
def run_task1():
    conditions = [
        {"name": "Summer Day", "T_db": 40.0, "rh": 50.0, "solar": True},
        {"name": "Summer Night", "T_db": 30.0, "rh": 50.0, "solar": False},
        {"name": "Monsoon Day", "T_db": 35.0, "rh": 85.0, "solar": True},
        {"name": "Monsoon Night", "T_db": 30.0, "rh": 85.0, "solar": False},
        {"name": "Winter Day", "T_db": 28.0, "rh": 50.0, "solar": True},
        {"name": "Winter Night", "T_db": 16.0, "rh": 50.0, "solar": False},
    ]

    psych_results = []
    for cond in conditions:
        psych = ZonePsychrometrics.solve_condition(
            T_db_out=cond["T_db"],
            rh_out=cond["rh"],
            solar_on=cond["solar"],
            fresh_air_frac=0.30,
        )
        psych_results.append((cond, psych))

    max_Q_evap = max(
        (p["Q_evap_chiller_kW"] for _, p in psych_results),
        default=0.0,
    )

    results = []

    for cond, psych in psych_results:
        Q_evap = float(psych["Q_evap_chiller_kW"])
        T_cond_est = cond["T_db"] + 12.0

        turndown_ratio = (
            100.0 * Q_evap / max_Q_evap
            if max_Q_evap > 0.0 else 0.0
        )
        turndown_ok = turndown_ratio >= 25.0 or Q_evap == 0.0

        row = {
            "Condition": cond["name"],
            "T_cond_est (°C)": T_cond_est,
            "T_mix (°C)": round(psych["T_mix"], 2),
            "ADP (°C)": round(psych["ADP"], 2),
            "w_supply (g/kg)": round(psych["w_supply"] * 1000, 3),
            "Room RH (%)": round(
                psych.get("room_rh_percent", psych.get("rh_room", float("nan"))),
                1,
            ),
            "RH Viol": psych["rh_violated"],
            "Q_evap_load (kW)": round(Q_evap, 2),
            "Turndown (%)": round(turndown_ratio, 1),
            "Turndown OK": turndown_ok,
            "Cycle feasible": False,
            "P_high (bar)": float("nan"),
            "Q_evap (kW)": float("nan"),
            "Q_gen (kW)": float("nan"),
            "COP": float("nan"),
            "m_r (kg/s)": float("nan"),
            "m_flue (kg/s)": float("nan"),
            "Flue OK": False,
            "Cycle status": "",
        }

        try:
            cycle = AbsorptionCycle.solve_cycle(
                Q_evap_kW=Q_evap,
                T_cond_C=T_cond_est,
                T_evap_C=2.0,
            )

            row.update({
                "Cycle feasible": bool(cycle["check_p_high_pass"]),
                "P_high (bar)": round(cycle["P_high_bar"], 3),
                "Q_evap (kW)": round(cycle["Q_evap_actual_kW"], 2),
                "Q_gen (kW)": round(cycle["Q_gen_kW"], 2),
                "COP": round(cycle["COP"], 3),
                "m_r (kg/s)": round(cycle["m_r_kg_s"], 5),
                "m_flue (kg/s)": round(cycle["m_flue_gas_kg_s"], 3),
                "Flue OK": bool(cycle["check_flue_pass"]),
                "Cycle status": cycle.get("model_status", "Calculated"),
            })

        except (ValueError, ArithmeticError) as exc:
            row["Cycle status"] = f"INFEASIBLE: {exc}"

        results.append(row)

    df_results = pd.DataFrame(results)

    feasible = df_results[df_results["Cycle feasible"]]
    if not feasible.empty:
        governing_condition = feasible.loc[
            feasible["Q_evap (kW)"].idxmax(), "Condition"
        ]
    else:
        governing_condition = "No feasible cycle condition"

    return df_results, governing_condition
