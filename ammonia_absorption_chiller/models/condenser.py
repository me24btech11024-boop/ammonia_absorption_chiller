"""
Condenser model for air-cooled and evaporative condensers.

Includes:
- Input validation
- Design-point UA and area estimates
- Off-design condensing-temperature iteration
- High-pressure and minimum-condensing-temperature checks

Limitations:
- Cooling-stream outlet temperature is estimated.
- Off-design results use a simplified LMTD model.
- A full cooling-medium energy balance is not implemented.
"""

import math

from core.moist_air import MoistAir
from core.pure_ammonia import PureAmmoniaEOS


class CondenserModel:
    T_COND_MIN_C = 25.0
    P_HIGH_MAX_BAR = 20.0

    # Overall heat-transfer coefficients converted to kW/(m^2 K).
    U_AIR_KW_M2_K = 0.030
    U_EVAP_KW_M2_K = 0.600

    @classmethod
    def _validate_inputs(
        cls,
        Q_cond_kW,
        T_db_C,
        rh_percent,
        condenser_type,
        UA_kW_K=None,
    ):
        values = [Q_cond_kW, T_db_C, rh_percent]

        if UA_kW_K is not None:
            values.append(UA_kW_K)

        if not all(math.isfinite(float(value)) for value in values):
            raise ValueError("Condenser inputs must be finite.")

        if Q_cond_kW <= 0:
            raise ValueError("Q_cond_kW must be positive.")

        if not 0.0 <= rh_percent <= 100.0:
            raise ValueError(
                "Relative humidity must be between 0 and 100%."
            )

        if condenser_type not in ("air_cooled", "evaporative"):
            raise ValueError(
                "condenser_type must be 'air_cooled' or 'evaporative'."
            )

        if UA_kW_K is not None and UA_kW_K <= 0:
            raise ValueError("UA_kW_K must be positive.")

    @classmethod
    def _ambient_temperature(cls, T_db_C, rh_percent, condenser_type):
        """Return dry-bulb or wet-bulb temperature for the selected type."""
        if condenser_type == "air_cooled":
            return float(T_db_C)

        return float(
            MoistAir.wet_bulb_temperature(T_db_C, rh_percent)
        )

    @staticmethod
    def _lmtd(dt1, dt2):
        """Calculate logarithmic mean temperature difference."""
        if not math.isfinite(dt1) or not math.isfinite(dt2):
            raise ValueError("LMTD temperature differences must be finite.")

        if dt1 <= 0 or dt2 <= 0:
            raise ValueError(
                "LMTD requires both terminal temperature differences "
                "to be positive."
            )

        if abs(dt1 - dt2) < 1e-9:
            return 0.5 * (dt1 + dt2)

        ratio = dt1 / dt2

        if ratio <= 0:
            raise ValueError("Invalid temperature-difference ratio.")

        return (dt1 - dt2) / math.log(ratio)

    @classmethod
    def size_design_condenser(
        cls,
        Q_cond_kW: float,
        T_db_C: float,
        rh_percent: float,
        condenser_type: str = "evaporative",
        approach_K: float = 8.0,
    ):
        """Estimate design condensing temperature, UA, and area."""
        cls._validate_inputs(
            Q_cond_kW,
            T_db_C,
            rh_percent,
            condenser_type,
        )

        if not math.isfinite(float(approach_K)) or approach_K <= 0:
            raise ValueError("approach_K must be finite and positive.")

        T_ambient = cls._ambient_temperature(
            T_db_C,
            rh_percent,
            condenser_type,
        )

        T_cond_design = T_ambient + approach_K

        # Simplified assumed cooling-stream outlet temperature.
        # This is an assumption, not an independently solved energy balance.
        T_out = T_ambient + 0.5 * approach_K

        dt1 = T_cond_design - T_ambient
        dt2 = T_cond_design - T_out
        lmtd = cls._lmtd(dt1, dt2)

        UA_design_kW_K = Q_cond_kW / lmtd

        if condenser_type == "air_cooled":
            U_coeff = cls.U_AIR_KW_M2_K
        else:
            U_coeff = cls.U_EVAP_KW_M2_K

        area_m2 = UA_design_kW_K / U_coeff

        T_cond_K = T_cond_design + 273.15
        P_high_bar = PureAmmoniaEOS.p_sat_lee_kesler(T_cond_K)

        return {
            "condenser_type": condenser_type,
            "T_ambient_C": T_ambient,
            "T_cond_design_C": T_cond_design,
            "T_cooling_out_assumed_C": T_out,
            "approach_K": approach_K,
            "LMTD_K": lmtd,
            "U_kW_m2_K": U_coeff,
            "UA_kW_K": UA_design_kW_K,
            "area_m2": area_m2,
            "Q_design_kW": Q_cond_kW,
            "P_high_design_bar": P_high_bar,
            "check_p_high_pass": bool(
                P_high_bar <= cls.P_HIGH_MAX_BAR
            ),
            "model_status": (
                "provisional: cooling outlet temperature is assumed"
            ),
        }

    @classmethod
    def solve_off_design(
        cls,
        Q_cond_kW: float,
        UA_kW_K: float,
        T_db_C: float,
        rh_percent: float,
        condenser_type: str = "evaporative",
    ):
        """Solve a simplified off-design condensing-temperature estimate."""
        cls._validate_inputs(
            Q_cond_kW,
            T_db_C,
            rh_percent,
            condenser_type,
            UA_kW_K=UA_kW_K,
        )

        T_ambient = cls._ambient_temperature(
            T_db_C,
            rh_percent,
            condenser_type,
        )

        # Initial approach estimate, with the specified minimum floor.
        T_cond = max(T_ambient + 8.0, cls.T_COND_MIN_C)

        converged = False
        iterations = 0
        residual_kW = None

        for iteration in range(1, 101):
            iterations = iteration

            T_out = T_ambient + 0.5 * (T_cond - T_ambient)

            dt1 = T_cond - T_ambient
            dt2 = T_cond - T_out
            lmtd = cls._lmtd(dt1, dt2)

            Q_calc = UA_kW_K * lmtd
            residual_kW = Q_calc - Q_cond_kW

            if abs(residual_kW) < 1e-3:
                converged = True
                break

            # For this simplified outlet-temperature model, LMTD is
            # proportional to the condensing-to-ambient temperature gap.
            # Update the gap directly rather than using an arbitrary step.
            gap = T_cond - T_ambient

            if gap <= 0:
                raise ValueError(
                    "Condensing temperature must exceed ambient temperature."
                )

            new_gap = gap * Q_cond_kW / Q_calc

            if not math.isfinite(new_gap) or new_gap <= 0:
                raise ValueError(
                    "Off-design iteration produced an invalid temperature."
                )

            T_cond = T_ambient + new_gap

        # Apply the minimum head-pressure control floor.
        if T_cond < cls.T_COND_MIN_C:
            T_cond = cls.T_COND_MIN_C

            T_out = T_ambient + 0.5 * (T_cond - T_ambient)
            dt1 = T_cond - T_ambient
            dt2 = T_cond - T_out

            if dt1 > 0 and dt2 > 0:
                lmtd = cls._lmtd(dt1, dt2)
                Q_calc = UA_kW_K * lmtd
                residual_kW = Q_calc - Q_cond_kW

            # The floor may prevent the requested duty from being met.
            converged = (
                residual_kW is not None
                and abs(residual_kW) < 1e-3
            )

        T_cond_K = T_cond + 273.15
        P_high_bar = PureAmmoniaEOS.p_sat_lee_kesler(T_cond_K)

        return {
            "T_cond_C": T_cond,
            "T_ambient_C": T_ambient,
            "P_high_bar": P_high_bar,
            "Q_required_kW": Q_cond_kW,
            "Q_calculated_kW": (
                Q_cond_kW + residual_kW
                if residual_kW is not None
                else None
            ),
            "energy_residual_kW": residual_kW,
            "iterations": iterations,
            "converged": converged,
            "check_p_high_pass": bool(
                P_high_bar <= cls.P_HIGH_MAX_BAR
            ),
            "check_t_cond_min_pass": bool(
                T_cond >= cls.T_COND_MIN_C
            ),
            "model_status": (
                "provisional: simplified LMTD and assumed cooling outlet"
            ),
        }
