"""
Condenser Model with Head-Pressure Control Floor (T_cond >= 25.0°C).
"""

from core.moist_air import MoistAir
from core.pure_ammonia import PureAmmoniaEOS


class CondenserModel:
    T_COND_MIN_C = 25.0      # Minimum head pressure floor [°C]
    P_HIGH_MAX_BAR = 20.0    # Vessel pressure safety limit [bar]

    @classmethod
    def size_design_condenser(cls, Q_cond_kW: float, T_db_C: float, rh_percent: float, 
                              condenser_type: str = "evaporative", approach_K: float = 8.0):
        if condenser_type == "air_cooled":
            T_ambient = T_db_C
        else:
            T_ambient = MoistAir.wet_bulb_temperature(T_db_C, rh_percent)

        T_cond_design = T_ambient + approach_K
        T_out = T_ambient + approach_K * 0.5

        # LMTD calculation
        dt1 = T_cond_design - T_ambient
        dt2 = T_cond_design - T_out
        if abs(dt1 - dt2) < 1e-4:
            lmtd = dt1
        else:
            import numpy as np
            lmtd = (dt1 - dt2) / np.log(dt1 / dt2)

        UA_design_kW_K = Q_cond_kW / lmtd
        if condenser_type == "air_cooled":
            U_coeff = 0.030  # kW/(m^2 K): 30 W/(m^2 K)
        elif condenser_type == "evaporative":
            U_coeff = 0.600  # kW/(m^2 K): 600 W/(m^2 K)
        else:
            raise ValueError(
                "condenser_type must be 'air_cooled' or 'evaporative'"
            )
        area_m2 = UA_design_kW_K / U_coeff

        return {
            "condenser_type": condenser_type,
            "T_cond_design_C": T_cond_design,
            "UA_kW_K": UA_design_kW_K,
            "area_m2": area_m2,
            "Q_design_kW": Q_cond_kW,
        }

    @classmethod
    def solve_off_design(cls, Q_cond_kW: float, UA_kW_K: float, T_db_C: float, 
                         rh_percent: float, condenser_type: str = "evaporative"):
        if condenser_type == "air_cooled":
            T_ambient = T_db_C
        else:
            T_ambient = MoistAir.wet_bulb_temperature(T_db_C, rh_percent)

        T_cond = max(T_ambient + 8.0, cls.T_COND_MIN_C)

        for _ in range(50):
            T_out = T_ambient + (T_cond - T_ambient) * 0.5
            dt1 = T_cond - T_ambient
            dt2 = T_cond - T_out
            if abs(dt1 - dt2) < 1e-4:
                lmtd = dt1
            else:
                import numpy as np
                lmtd = (dt1 - dt2) / np.log(dt1 / dt2)

            Q_calc = UA_kW_K * lmtd
            residual = Q_calc - Q_cond_kW
            if abs(residual) < 1e-3:
                break
            T_cond += (Q_cond_kW - Q_calc) / (UA_kW_K + 1e-5)

        # Enforce minimum head pressure control floor
        if T_cond < cls.T_COND_MIN_C:
            T_cond = cls.T_COND_MIN_C

        T_cond_K = T_cond + 273.15
        P_high_bar = PureAmmoniaEOS.p_sat_lee_kesler(T_cond_K)

        return {
            "T_cond_C": T_cond,
            "P_high_bar": P_high_bar,
            "check_p_high_pass": bool(P_high_bar <= cls.P_HIGH_MAX_BAR),
            "check_t_cond_min_pass": bool(T_cond >= cls.T_COND_MIN_C),
        }
