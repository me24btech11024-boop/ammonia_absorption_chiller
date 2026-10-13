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
    def wet_bulb_temperature(cls, T_db_C: float, rh_percent: float, P_kPa: float = P_ATM_KPA) -> float:
        """Solves wet-bulb temperature T_wb (°C) iteratively using adiabatic saturation balance."""
        w = cls.humidity_ratio(T_db_C, rh_percent, P_kPa)
        h_actual = cls.enthalpy(T_db_C, w)

        def residual(T_wb_trial):
            w_sat_wb = cls.humidity_ratio(T_wb_trial, 100.0, P_kPa)
            h_sat_wb = cls.enthalpy(T_wb_trial, w_sat_wb)
            cp_liquid = 4.184
            h_f = cp_liquid * T_wb_trial
            
            h_calc = h_sat_wb - (w_sat_wb - w) * h_f
            return h_calc - h_actual

        t1 = T_db_C - 5.0
        t2 = T_db_C

        for _ in range(50):
            r1 = residual(t1)
            r2 = residual(t2)
            if abs(r2 - r1) < 1e-7:
                break
            t_next = t2 - r2 * (t2 - t1) / (r2 - r1)
            t1, t2 = t2, t_next
            if abs(t2 - t1) < 1e-4:
                break

        return t2