"""
Ideal-Gas Moist Air Psychrometric Calculations and Wet-Bulb Solver.

Implements:
1. Vapor pressure / saturation pressure correlations for H2O.
2. Ideal-gas moist air humidity ratio (w) and enthalpy (h).
3. Psychrometric mixing and cooling coil process analysis.
4. Iterative adiabatic-saturation energy balance solver for Wet-Bulb Temperature (T_wb).
"""

import numpy as np


class MoistAir:
    """
    Self-contained ideal-gas moist air psychrometric property evaluator.
    """

    P_ATM_KPA = 101.325
    CP_DRY_AIR = 1.006
    CP_VAPOR = 1.86
    H_FG_0 = 2501.0

    @classmethod
    def p_sat_water_antoine(cls, T_C: float) -> float:
        """
        Calculates water saturation pressure P_sat (kPa) for T in Celsius (°C).
        Matches standard steam tables (e.g. 22°C ~ 2.645 kPa, 40°C = 7.384 kPa).
        """
        # High-accuracy vapor pressure correlation for water
        T_K = T_C + 273.15
        tc = 647.096
        pc = 22064.0  # kPa
        c = [-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502]
        tau = 1.0 - (T_K / tc)
        
        num = c[0]*tau + c[1]*(tau**1.5) + c[2]*(tau**3) + c[3]*(tau**3.5) + c[4]*(tau**4) + c[5]*(tau**7.5)
        p_sat = pc * np.exp((tc / T_K) * num)
        return p_sat

    @classmethod
    def humidity_ratio(cls, T_C: float, rh_percent: float, P_kPa: float = P_ATM_KPA) -> float:
        """
        Calculates humidity ratio w (kg H2O / kg dry air) from T (°C) and RH (%).
        w = 0.622 * P_v / (P - P_v)
        """
        p_sat = cls.p_sat_water_antoine(T_C)
        p_v = (rh_percent / 100.0) * p_sat

        if p_v >= P_kPa:
            p_v = P_kPa - 1e-5

        w = 0.622 * p_v / (P_kPa - p_v)
        return w

    @classmethod
    def vapor_pressure(cls, w: float, P_kPa: float = P_ATM_KPA) -> float:
        """Calculates vapor pressure P_v (kPa) from humidity ratio w."""
        p_v = P_kPa * w / (0.622 + w)
        return p_v

    @classmethod
    def enthalpy(cls, T_C: float, w: float) -> float:
        """
        Calculates moist air specific enthalpy h (kJ/kg dry air).
        h = Cp_da * T + w * (h_fg_0 + Cp_v * T)
        """
        h = cls.CP_DRY_AIR * T_C + w * (cls.H_FG_0 + cls.CP_VAPOR * T_C)
        return h

    @classmethod
    def relative_humidity(cls, T_C: float, w: float, P_kPa: float = P_ATM_KPA) -> float:
        """Calculates Relative Humidity RH (%) from T (°C) and humidity ratio w."""
        p_v = cls.vapor_pressure(w, P_kPa)
        p_sat = cls.p_sat_water_antoine(T_C)
        rh = (p_v / p_sat) * 100.0
        return min(rh, 100.0)
        
    @classmethod
    def wet_bulb_temperature(
        cls,
        T_db_C: float,
        rh_percent: float,
        P_kPa: float = P_ATM_KPA,
    ) -> float:
        """Solve the adiabatic-saturation energy balance for wet-bulb temperature."""
        if not np.isfinite(T_db_C) or not np.isfinite(rh_percent):
            raise ValueError("Dry-bulb temperature and RH must be finite.")
        if not 0.0 <= rh_percent <= 100.0:
            raise ValueError("Relative humidity must be between 0 and 100%.")
        if P_kPa <= 0.0:
            raise ValueError("Pressure must be positive.")

        if rh_percent == 100.0:
            return float(T_db_C)

        w = cls.humidity_ratio(T_db_C, rh_percent, P_kPa)
        h_actual = cls.enthalpy(T_db_C, w)

        def residual(T_wb_trial):
            w_sat = cls.humidity_ratio(T_wb_trial, 100.0, P_kPa)
            h_sat = cls.enthalpy(T_wb_trial, w_sat)
            h_liquid = 4.184 * T_wb_trial
            return h_sat - (w_sat - w) * h_liquid - h_actual

        # For unsaturated air, wet-bulb temperature lies below dry-bulb.
        low = min(-50.0, T_db_C - 50.0)
        high = T_db_C
        f_low = residual(low)
        f_high = residual(high)

        # Expand the lower bracket if necessary.
        for _ in range(20):
            if f_low * f_high <= 0.0:
                break
            low -= 50.0
            f_low = residual(low)
        else:
            raise ValueError("Could not bracket the wet-bulb energy-balance root.")

        # Bisection: slower than secant, but bounded and robust.
        for _ in range(200):
            mid = 0.5 * (low + high)
            f_mid = residual(mid)

            if abs(f_mid) < 1e-7 or abs(high - low) < 1e-5:
                return float(mid)

            if f_low * f_mid <= 0.0:
                high = mid
                f_high = f_mid
            else:
                low = mid
                f_low = f_mid

        raise ValueError("Wet-bulb solver did not converge.")
