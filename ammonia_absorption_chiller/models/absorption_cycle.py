"""
Single-effect ammonia-water absorption-cycle solver.

Includes:
- Input validation and invalid-state guards
- Refrigerant and solution mass-flow calculations
- Approximate generator, SHX, and rectifier duties
- High-side pressure and flue-gas capacity checks
- Mass-balance diagnostics

Limitations:
- Strong and weak solution concentrations are assumed.
- Generator temperature is estimated as T_cond + 45 K.
- SHX and rectifier duties are approximate.
- Absorber and condenser duties are not independently solved.
- Energy-balance closure is not claimed by this model.
"""

import math

from core.patek_klomfar import PatekKlomfar
from core.pure_ammonia import PureAmmoniaEOS


class AbsorptionCycle:
    PUMP_EFF = 0.70
    SHX_EFF = 0.70
    RECTIFIER_X_VAP = 0.995

    FLUE_GAS_CP = 1.10       # kJ/(kg K)
    FLUE_GAS_T_IN = 300.0    # °C
    FLUE_GAS_T_FLOOR = 140.0 # °C
    FLUE_GAS_MAX_FLOW = 80.0 # kg/s
    GEN_EFF = 0.75

    P_HIGH_MAX_BAR = 20.0

    # Provisional mass fractions; not calculated from equilibrium.
    X_STRONG_ASSUMED = 0.48
    X_WEAK_ASSUMED = 0.32

    # Approximate solution specific volume for pump work.
    SOLUTION_SPECIFIC_VOLUME_M3_KG = 0.0012

    @classmethod
    def _fail(cls, message):
        raise ValueError(f"Absorption-cycle calculation failed: {message}")

    @classmethod
    def solve_cycle(
        cls,
        Q_evap_kW: float,
        T_cond_C: float,
        T_evap_C: float = 2.0,
    ):
        """Calculate a provisional single-effect cycle operating point."""

        for name, value in {
            "Q_evap_kW": Q_evap_kW,
            "T_cond_C": T_cond_C,
            "T_evap_C": T_evap_C,
        }.items():
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
            ):
                cls._fail(f"{name} must be a finite number.")

        if Q_evap_kW < 0.0:
            cls._fail("Q_evap_kW cannot be negative.")

        if T_cond_C <= T_evap_C:
            cls._fail("Condenser temperature must exceed evaporator temperature.")

        T_evap_K = T_evap_C + 273.15
        T_cond_K = T_cond_C + 273.15

        if T_cond_K >= PureAmmoniaEOS.Tc:
            cls._fail("Condenser temperature must be below ammonia's critical temperature.")

        if T_evap_K >= PureAmmoniaEOS.Tc:
            cls._fail("Evaporator temperature must be below ammonia's critical temperature.")

        P_low_bar = PureAmmoniaEOS.p_sat_lee_kesler(T_evap_K)
        P_high_bar = PureAmmoniaEOS.p_sat_lee_kesler(T_cond_K)

        if not math.isfinite(P_low_bar) or P_low_bar <= 0.0:
            cls._fail(f"Invalid low-side pressure: {P_low_bar} bar.")

        if not math.isfinite(P_high_bar) or P_high_bar <= 0.0:
            cls._fail(f"Invalid high-side pressure: {P_high_bar} bar.")

        # Approximate generator temperature assumption.
        T_gen_C = T_cond_C + 45.0
        T_gen_K = T_gen_C + 273.15

        if T_gen_K >= PureAmmoniaEOS.Tc:
            cls._fail(
                f"Estimated generator temperature {T_gen_C:.2f} °C "
                "is at or above the ammonia critical temperature."
            )

        # Zero-load result: pressure checks remain meaningful, but the cycle
        # mass flows and duties are zero.
        if Q_evap_kW == 0.0:
            return {
                "P_low_bar": P_low_bar,
                "P_high_bar": P_high_bar,
                "T_gen_C": T_gen_C,
                "m_r_kg_s": 0.0,
                "m_strong_kg_s": 0.0,
                "m_weak_kg_s": 0.0,
                "x_r": cls.RECTIFIER_X_VAP,
                "x_strong": cls.X_STRONG_ASSUMED,
                "x_weak": cls.X_WEAK_ASSUMED,
                "circulation_ratio": 0.0,
                "q_evap_specific_kJ_kg": 0.0,
                "Q_evap_kW": 0.0,
                "Q_evap_actual_kW": 0.0,
                "Q_shx_approx_kW": 0.0,
                "Q_rect_approx_kW": 0.0,
                "Q_gen_kW": 0.0,
                "Q_reject_required_kW": 0.0,
                "W_pump_kW": 0.0,
                "COP": 0.0,
                "m_flue_gas_kg_s": 0.0,
                "mass_balance_solution_residual_kg_s": 0.0,
                "mass_balance_ammonia_residual_kg_s": 0.0,
                "energy_balance_residual_kW": None,
                "check_p_high_pass": P_high_bar <= cls.P_HIGH_MAX_BAR,
                "check_flue_pass": True,
                "check_generator_temperature_pass": (
                    T_gen_C < cls.FLUE_GAS_T_IN
                    and T_gen_C < cls.FLUE_GAS_T_FLOOR
                ),
                "concentrations_are_assumed": True,
                "model_status": "zero load; cycle not evaluated",
            }

        if P_high_bar > cls.P_HIGH_MAX_BAR:
            cls._fail(
                f"High-side pressure {P_high_bar:.3f} bar exceeds "
                f"the {cls.P_HIGH_MAX_BAR:.1f} bar limit."
            )

        # Refrigerating effect from pure-ammonia property estimates.
        h_evap_out = PureAmmoniaEOS.enthalpy(
            T_evap_K, P_low_bar, phase="vapor"
        )
        h_cond_out = PureAmmoniaEOS.enthalpy(
            T_cond_K, P_high_bar, phase="liquid"
        )
        q_evap_specific = h_evap_out - h_cond_out

        if not math.isfinite(q_evap_specific) or q_evap_specific <= 0.0:
            cls._fail(
                f"Invalid specific refrigerating effect: "
                f"{q_evap_specific} kJ/kg."
            )

        m_r = Q_evap_kW / q_evap_specific

        if not math.isfinite(m_r) or m_r <= 0.0:
            cls._fail(f"Invalid refrigerant mass flow: {m_r} kg/s.")

        # Solution mass balances using the assumed mass fractions.
        x_r = cls.RECTIFIER_X_VAP
        x_strong = cls.X_STRONG_ASSUMED
        x_weak = cls.X_WEAK_ASSUMED

        if not (0.0 < x_weak < x_strong < x_r <= 1.0):
            cls._fail("Solution and refrigerant mass fractions are invalid.")

        circulation_ratio = (x_r - x_weak) / (x_strong - x_weak)
        m_strong = m_r * circulation_ratio
        m_weak = m_strong - m_r

        solution_mass_residual = m_strong - m_weak - m_r
        ammonia_mass_residual = (
            m_strong * x_strong
            - m_weak * x_weak
            - m_r * x_r
        )

        # Approximate solution enthalpies from the Patek-Klomfar module.
        h_abs_out = PatekKlomfar.liquid_enthalpy(
            P_low_bar, x_strong, T_cond_K
        )
        h_gen_out = PatekKlomfar.liquid_enthalpy(
            P_high_bar, x_weak, T_gen_K
        )
        h_weak_at_abs_temp = PatekKlomfar.liquid_enthalpy(
            P_high_bar, x_weak, T_cond_K
        )

        # Pump specific work: v * delta-P / efficiency.
        delta_P_kPa = (P_high_bar - P_low_bar) * 100.0
        if delta_P_kPa < 0.0:
            cls._fail("High-side pressure is lower than low-side pressure.")

        w_pump_specific = (
            cls.SOLUTION_SPECIFIC_VOLUME_M3_KG
            * delta_P_kPa
            / cls.PUMP_EFF
        )

        W_pump_kW = m_strong * w_pump_specific
        h_pump_out = h_abs_out + w_pump_specific

        # Approximate solution heat exchanger duty.
        shx_delta_h = h_gen_out - h_weak_at_abs_temp
        Q_shx_kW = (
            cls.SHX_EFF
            * m_weak
            * max(shx_delta_h, 0.0)
        )
        h_gen_in = h_pump_out + Q_shx_kW / m_strong

        # Approximate rectifier and generator duties.
        h_rect_vapor = PureAmmoniaEOS.enthalpy(
            T_gen_K, P_high_bar, phase="vapor"
        )
        h_fg = PureAmmoniaEOS.h_fg(T_gen_K)
        Q_rect_kW = 2.15 * m_r * h_fg

        Q_gen_kW = (
            m_r * h_rect_vapor
            + m_weak * h_gen_out
            - m_strong * h_gen_in
            + Q_rect_kW
        )

        if not math.isfinite(Q_gen_kW) or Q_gen_kW <= 0.0:
            cls._fail(f"Invalid generator heat duty: {Q_gen_kW} kW.")

        COP = Q_evap_kW / Q_gen_kW

        if not math.isfinite(COP) or COP <= 0.0:
            cls._fail(f"Invalid COP: {COP}.")

        # Flue-gas heat availability, using the assignment's effectiveness.
        available_flue_temperature_drop_K = (
            cls.GEN_EFF
            * (cls.FLUE_GAS_T_IN - cls.FLUE_GAS_T_FLOOR)
        )

        if available_flue_temperature_drop_K <= 0.0:
            cls._fail("Available flue-gas temperature drop is non-positive.")

        m_flue_gas = Q_gen_kW / (
            cls.FLUE_GAS_CP * available_flue_temperature_drop_K
        )

        if not math.isfinite(m_flue_gas) or m_flue_gas < 0.0:
            cls._fail(f"Invalid flue-gas flow: {m_flue_gas} kg/s.")

        if m_flue_gas > cls.FLUE_GAS_MAX_FLOW + 1e-3:
            cls._fail(
                f"Required flue-gas flow {m_flue_gas:.3f} kg/s exceeds "
                f"the {cls.FLUE_GAS_MAX_FLOW:.1f} kg/s limit."
            )

        T_flue_out_C = (
            cls.FLUE_GAS_T_IN
            - Q_gen_kW / (m_flue_gas * cls.FLUE_GAS_CP)
            if m_flue_gas > 0.0
            else cls.FLUE_GAS_T_IN
        )

        # Whole-cycle rejection requirement from a first-law estimate.
        # This is not an independent condenser/absorber energy balance.
        Q_reject_required_kW = Q_gen_kW + Q_evap_kW + W_pump_kW

        return {
            "P_low_bar": P_low_bar,
            "Q_evap_kW": Q_evap_kW,
            "P_high_bar": P_high_bar,
            "T_gen_C": T_gen_C,
            "m_r_kg_s": m_r,
            "m_strong_kg_s": m_strong,
            "m_weak_kg_s": m_weak,
            "x_r": x_r,
            "x_strong": x_strong,
            "x_weak": x_weak,
            "circulation_ratio": circulation_ratio,
            "q_evap_specific_kJ_kg": q_evap_specific,
            "Q_evap_actual_kW": Q_evap_kW,
            "Q_shx_approx_kW": Q_shx_kW,
            "Q_rect_approx_kW": Q_rect_kW,
            "Q_gen_kW": Q_gen_kW,
            "Q_reject_required_kW": Q_reject_required_kW,
            "W_pump_kW": W_pump_kW,
            "COP": COP,
            "m_flue_gas_kg_s": m_flue_gas,
            "T_flue_out_C": T_flue_out_C,
            "mass_balance_solution_residual_kg_s": solution_mass_residual,
            "mass_balance_ammonia_residual_kg_s": ammonia_mass_residual,
            "energy_balance_residual_kW": None,
            "check_p_high_pass": P_high_bar <= cls.P_HIGH_MAX_BAR,
            "check_flue_pass": (
                m_flue_gas <= cls.FLUE_GAS_MAX_FLOW + 1e-3
                and T_flue_out_C >= cls.FLUE_GAS_T_FLOOR - 1e-6
            ),
            "check_generator_temperature_pass": (
                T_gen_C < cls.FLUE_GAS_T_IN
            ),
            "concentrations_are_assumed": True,
            "model_status": (
                "provisional; solution concentrations, SHX, and rectifier "
                "are approximate; energy balance not independently closed"
            ),
        }
