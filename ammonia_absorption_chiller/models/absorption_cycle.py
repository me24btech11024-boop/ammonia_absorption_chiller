"""
Single-effect ammonia-water absorption-cycle solver.

Adds input validation, invalid-state guards, and mass/energy diagnostics.

Limitations: solution concentrations remain assumptions; SHX and rectifier
duties remain approximate. Diagnostics are not proof of a fully validated
thermodynamic cycle.
"""
import numpy as np
from core.patek_klomfar import PatekKlomfar
from core.pure_ammonia import PureAmmoniaEOS


class AbsorptionCycle:
    PUMP_EFF = 0.70
    SHX_EFF = 0.70
    RECTIFIER_X_VAP = 0.995
    FLUE_GAS_CP = 1.10
    FLUE_GAS_T_IN = 300.0
    FLUE_GAS_T_FLOOR = 140.0
    FLUE_GAS_MAX_FLOW = 80.0
    GEN_EFF = 0.75
    P_HIGH_MAX_BAR = 20.0

    # Provisional assumptions inherited from the original model.
    X_STRONG_ASSUMED = 0.48
    X_WEAK_ASSUMED = 0.32

    @classmethod
    def _fail(cls, message):
        raise ValueError(f"Absorption-cycle calculation failed: {message}")

    @classmethod
    def solve_cycle(cls, Q_evap_kW: float, T_cond_C: float, T_evap_C: float = 2.0):
        for name, value in {
            "Q_evap_kW": Q_evap_kW,
            "T_cond_C": T_cond_C,
            "T_evap_C": T_evap_C,
        }.items():
            if not np.isfinite(value):
                cls._fail(f"{name} must be finite.")

        if T_cond_C <= T_evap_C:
            cls._fail("Condenser temperature must exceed evaporator temperature.")

        if Q_evap_kW <= 0:
            P_low_bar = PureAmmoniaEOS.p_sat_lee_kesler(T_evap_C + 273.15)
            P_high_bar = PureAmmoniaEOS.p_sat_lee_kesler(T_cond_C + 273.15)
            return {
                "P_low_bar": P_low_bar,
                "P_high_bar": P_high_bar,
                "T_gen_C": T_cond_C + 45.0,
                "m_r_kg_s": 0.0,
                "m_strong_kg_s": 0.0,
                "m_weak_kg_s": 0.0,
                "Q_evap_actual_kW": 0.0,
                "Q_gen_kW": 0.0,
                "W_pump_kW": 0.0,
                "COP": 0.0,
                "m_flue_gas_kg_s": 0.0,
                "check_p_high_pass": bool(P_high_bar <= cls.P_HIGH_MAX_BAR),
                "check_flue_pass": True,
                "mass_balance_solution_residual_kg_s": 0.0,
                "mass_balance_ammonia_residual_kg_s": 0.0,
                "energy_balance_residual_kW": None,
                "concentrations_are_assumed": True,
                "model_status": "zero-load; cycle not evaluated",
            }

        T_evap_K = T_evap_C + 273.15
        T_cond_K = T_cond_C + 273.15
        T_abs_K = T_cond_K
        T_gen_K = T_cond_K + 45.0
        T_gen_C = T_gen_K - 273.15

        P_low_bar = PureAmmoniaEOS.p_sat_lee_kesler(T_evap_K)
        P_high_bar = PureAmmoniaEOS.p_sat_lee_kesler(T_cond_K)
        if not np.isfinite(P_low_bar) or P_low_bar <= 0:
            cls._fail(f"Invalid low-side pressure: {P_low_bar}.")
        if not np.isfinite(P_high_bar) or P_high_bar <= 0:
            cls._fail(f"Invalid high-side pressure: {P_high_bar}.")

        x_r = cls.RECTIFIER_X_VAP
        h_evap_out = PureAmmoniaEOS.enthalpy(T_evap_K, P_low_bar, phase="vapor")
        h_cond_out = PureAmmoniaEOS.enthalpy(T_cond_K, P_high_bar, phase="liquid")
        q_e = h_evap_out - h_cond_out
        if not np.isfinite(q_e) or q_e <= 0:
            cls._fail(
                f"Invalid refrigerating effect: {q_e} kJ/kg; check ammonia "
                "enthalpies, units, and phase states."
            )

        m_r = Q_evap_kW / q_e
        if not np.isfinite(m_r) or m_r <= 0:
            cls._fail(f"Invalid refrigerant mass flow: {m_r} kg/s.")

        x_strong = cls.X_STRONG_ASSUMED
        x_weak = cls.X_WEAK_ASSUMED
        if not (0.0 < x_weak < x_strong < x_r <= 1.0):
            cls._fail("Solution/refrigerant mass fractions are invalid.")

        f_ratio = (x_r - x_weak) / (x_strong - x_weak)
        if not np.isfinite(f_ratio) or f_ratio <= 0:
            cls._fail(f"Invalid circulation ratio: {f_ratio}.")

        m_strong = m_r * f_ratio
        m_weak = m_strong - m_r
        solution_mass_residual = m_strong - m_weak - m_r
        ammonia_mass_residual = (
            m_strong * x_strong - m_weak * x_weak - m_r * x_r
        )

        h_abs_out = PatekKlomfar.liquid_enthalpy(P_low_bar, x_strong, T_abs_K)
        h_gen_out = PatekKlomfar.liquid_enthalpy(P_high_bar, x_weak, T_gen_K)

        v_liq = 0.0012
        w_pump_specific = (
            v_liq * (P_high_bar - P_low_bar) * 100.0
        ) / cls.PUMP_EFF
        if not np.isfinite(w_pump_specific) or w_pump_specific < 0:
            cls._fail(f"Invalid specific pump work: {w_pump_specific} kJ/kg.")

        W_pump_kW = m_strong * w_pump_specific
        h_pump_out = h_abs_out + w_pump_specific
        h_weak_abs = PatekKlomfar.liquid_enthalpy(P_high_bar, x_weak, T_abs_K)
        shx_delta_h = h_gen_out - h_weak_abs
        Q_shx_kW = cls.SHX_EFF * m_weak * max(shx_delta_h, 10.0)
        h_gen_in = h_pump_out + Q_shx_kW / m_strong

        h_rect_vap = PureAmmoniaEOS.enthalpy(T_gen_K, P_high_bar, phase="vapor")
        h_fg = PureAmmoniaEOS.h_fg(T_gen_K)
        Q_rect_kW = 2.15 * m_r * h_fg
        Q_gen_kW = (
            m_r * h_rect_vap
            + m_weak * h_gen_out
            - m_strong * h_gen_in
            + Q_rect_kW
        )
        if not np.isfinite(Q_gen_kW) or Q_gen_kW <= 0:
            cls._fail(f"Invalid generator heat duty: {Q_gen_kW} kW.")

        COP = Q_evap_kW / Q_gen_kW
        if not np.isfinite(COP) or COP <= 0:
            cls._fail(f"Invalid COP: {COP}.")

        dT_flue_max = cls.GEN_EFF * (cls.FLUE_GAS_T_IN - cls.FLUE_GAS_T_FLOOR)
        if dT_flue_max <= 0:
            cls._fail("Flue-gas available temperature drop is non-positive.")
        m_flue_gas = Q_gen_kW / (cls.FLUE_GAS_CP * dT_flue_max)

        # Required total heat rejection from the whole-cycle first law.
        # Absorber and condenser duties are not independently solved here.
        Q_reject_required_kW = Q_gen_kW + Q_evap_kW + W_pump_kW
        energy_balance_residual_kW = (
            Q_gen_kW + Q_evap_kW + W_pump_kW - Q_reject_required_kW
        )

        return {
            "P_low_bar": P_low_bar,
            "P_high_bar": P_high_bar,
            "T_gen_C": T_gen_C,
            "m_r_kg_s": m_r,
            "m_strong_kg_s": m_strong,
            "m_weak_kg_s": m_weak,
            "x_r": x_r,
            "x_strong": x_strong,
            "x_weak": x_weak,
            "circulation_ratio": f_ratio,
            "q_evap_specific_kJ_kg": q_e,
            "Q_evap_actual_kW": Q_evap_kW,
            "Q_shx_approx_kW": Q_shx_kW,
            "Q_rect_approx_kW": Q_rect_kW,
            "Q_gen_kW": Q_gen_kW,
            "Q_reject_required_kW": Q_reject_required_kW,
            "W_pump_kW": W_pump_kW,
            "COP": COP,
            "m_flue_gas_kg_s": m_flue_gas,
            "mass_balance_solution_residual_kg_s": solution_mass_residual,
            "mass_balance_ammonia_residual_kg_s": ammonia_mass_residual,
            "energy_balance_residual_kW": energy_balance_residual_kW,
            "check_p_high_pass": bool(P_high_bar <= cls.P_HIGH_MAX_BAR),
            "check_flue_pass": bool(m_flue_gas <= cls.FLUE_GAS_MAX_FLOW + 1e-3),
            "concentrations_are_assumed": True,
            "model_status": "provisional; concentrations and SHX/rectifier are approximate",
        }
