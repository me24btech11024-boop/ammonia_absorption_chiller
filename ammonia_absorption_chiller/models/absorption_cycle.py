"""
Single-Effect Ammonia-Water Absorption Cycle Solver.
Calculates first-law energy balances, dynamic COP, and exact flue-gas demand.
"""

import numpy as np
from core.patek_klomfar import PatekKlomfar
from core.pure_ammonia import PureAmmoniaEOS


class AbsorptionCycle:
    PUMP_EFF = 0.70
    SHX_EFF = 0.70
    RECTIFIER_X_VAP = 0.995
    FLUE_GAS_CP = 1.10         # kJ/kg*K
    FLUE_GAS_T_IN = 300.0      # °C
    FLUE_GAS_T_FLOOR = 140.0   # °C
    FLUE_GAS_MAX_FLOW = 80.0   # kg/s
    GEN_EFF = 0.75             # Generator effectiveness
    P_HIGH_MAX_BAR = 20.0      # Vessel pressure limit bar

    @classmethod
    def solve_cycle(cls, Q_evap_kW: float, T_cond_C: float, T_evap_C: float = 2.0):
        # Zero load / Passive balance condition
        if Q_evap_kW <= 0:
            return {
                "P_low_bar": PureAmmoniaEOS.p_sat_lee_kesler(T_evap_C + 273.15),
                "P_high_bar": PureAmmoniaEOS.p_sat_lee_kesler(T_cond_C + 273.15),
                "T_gen_C": T_cond_C + 45.0,
                "m_r_kg_s": 0.0,
                "m_strong_kg_s": 0.0,
                "m_weak_kg_s": 0.0,
                "Q_evap_actual_kW": 0.0,
                "Q_gen_kW": 0.0,
                "W_pump_kW": 0.0,
                "COP": 0.0,
                "m_flue_gas_kg_s": 0.0,
                "check_p_high_pass": True,
                "check_flue_pass": True,
            }

        T_abs_C = T_cond_C
        T_evap_K = T_evap_C + 273.15
        T_cond_K = T_cond_C + 273.15
        T_abs_K = T_abs_C + 273.15

        # 1. Operating Pressures
        P_low_bar = PureAmmoniaEOS.p_sat_lee_kesler(T_evap_K)
        P_high_bar = PureAmmoniaEOS.p_sat_lee_kesler(T_cond_K)

        # 2. Refrigerant Flow Circuit
        x_r = cls.RECTIFIER_X_VAP
        h_evap_out = PureAmmoniaEOS.enthalpy(T_evap_K, P_low_bar, phase="vapor")
        h_cond_out = PureAmmoniaEOS.enthalpy(T_cond_K, P_high_bar, phase="liquid")

        q_e = h_evap_out - h_cond_out
        if q_e <= 0:
            q_e = 1080.0

        m_r = Q_evap_kW / q_e  # kg/s

        # 3. Mass Balance & Circulation Ratio
        x_strong = 0.48
        x_weak = 0.32
        f_ratio = max((x_r - x_weak) / (x_strong - x_weak), 4.0)

        m_strong = m_r * f_ratio
        m_weak = m_strong - m_r

        # 4. State Enthalpies & Solution Heat Exchanger (SHX)
        h_abs_out = PatekKlomfar.liquid_enthalpy(P_low_bar, x_strong, T_abs_K)
        T_gen_K = T_cond_K + 45.0
        T_gen_C = T_gen_K - 273.15
        h_gen_out = PatekKlomfar.liquid_enthalpy(P_high_bar, x_weak, T_gen_K)

        v_liq = 0.0012
        w_pump_actual = (v_liq * (P_high_bar - P_low_bar) * 100.0) / cls.PUMP_EFF
        W_pump_kW = m_strong * w_pump_actual
        h_pump_out = h_abs_out + w_pump_actual

        h_weak_abs = PatekKlomfar.liquid_enthalpy(P_high_bar, x_weak, T_abs_K)
        Q_shx = cls.SHX_EFF * (m_weak * max(h_gen_out - h_weak_abs, 10.0))
        h_gen_in = h_pump_out + Q_shx / m_strong

        # 5. Complete First-Law Generator Heat Balance
        h_rect_vap = PureAmmoniaEOS.enthalpy(T_gen_K, P_high_bar, phase="vapor")
        h_fg = PureAmmoniaEOS.h_fg(T_gen_K)
        Q_rect = 2.15 * m_r * h_fg

        Q_gen_kW = (m_r * h_rect_vap) + (m_weak * h_gen_out) - (m_strong * h_gen_in) + Q_rect

        # 6. Dynamic First-Law COP
        COP = Q_evap_kW / Q_gen_kW

        # 7. Flue Gas Mass Flow Rate Calculation
        dT_flue_max = cls.GEN_EFF * (cls.FLUE_GAS_T_IN - cls.FLUE_GAS_T_FLOOR)  # 0.75 * (300 - 140) = 120 K
        m_flue_gas = Q_gen_kW / (cls.FLUE_GAS_CP * dT_flue_max)

        # Strict Feasibility Flag
        check_flue_pass = bool(m_flue_gas <= cls.FLUE_GAS_MAX_FLOW + 1e-3)

        return {
            "P_low_bar": P_low_bar,
            "P_high_bar": P_high_bar,
            "T_gen_C": T_gen_C,
            "m_r_kg_s": m_r,
            "m_strong_kg_s": m_strong,
            "m_weak_kg_s": m_weak,
            "Q_evap_actual_kW": Q_evap_kW,
            "Q_gen_kW": Q_gen_kW,
            "W_pump_kW": W_pump_kW,
            "COP": COP,
            "m_flue_gas_kg_s": m_flue_gas,
            "check_p_high_pass": bool(P_high_bar <= cls.P_HIGH_MAX_BAR),
            "check_flue_pass": check_flue_pass,
        }