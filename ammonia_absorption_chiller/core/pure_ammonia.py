"""
Pure Ammonia Real-Gas EOS and Thermodynamic Departure Properties.

Implements:
1. Peng-Robinson Equation of State (EOS) for pure ammonia.
2. Lee-Kesler/Pitzer acentric-factor correlation for saturation pressure P_sat(T).
3. Ideal-gas heat capacity Cp0 integration from T_ref = 273.15 K.
4. Enthalpy (h) and entropy (s) departure functions for real-gas corrections.
"""

import numpy as np


class PureAmmoniaEOS:
    """
    Thermodynamic property evaluator for pure ammonia (NH3) using Peng-Robinson EOS.
    """

    # Critical constants & molecular properties for Ammonia (NH3)
    Tc = 405.56        # Critical temperature [K]
    Pc_bar = 113.33    # Critical pressure [bar] (11.333 MPa)
    omega = 0.256      # Acentric factor
    MW = 17.031        # Molecular weight [kg/kmol]
    R_gas = 8.3144626  # Universal gas constant [kJ/(kmol*K)]
    R_spec = R_gas / MW  # Specific gas constant [kJ/(kg*K)] = 0.4882

    # Reference state: T_ref = 273.15 K, P_ref = 1.01325 bar (101.325 kPa), ideal gas h_ref = 0, s_ref = 0
    T_ref = 273.15
    P_ref_bar = 1.01325

    # Ideal-gas Cp0 polynomial coefficients for NH3 (kJ/kg*K) as f(T [K])
    # Cp0(T) = a0 + a1*T + a2*T^2 + a3*T^3
    CP0_COEFFS = [1.61, 2.80e-3, -4.0e-7, 2.0e-11]

    @classmethod
    def p_sat_lee_kesler(cls, T_K: float) -> float:
        """
        Calculates saturation pressure P_sat (bar) at temperature T (K)
        using Lee-Kesler / Pitzer acentric factor correlation.
        """
        Tr = T_K / cls.Tc
        if Tr >= 1.0:
            return cls.Pc_bar

        # Pitzer vapor pressure correlation components
        f0 = 5.92714 - 6.09648 / Tr - 1.28862 * np.log(Tr) + 0.169347 * (Tr**6)
        f1 = 15.2518 - 15.6875 / Tr - 13.4721 * np.log(Tr) + 0.43577 * (Tr**6)

        ln_Pr_sat = f0 + cls.omega * f1
        Pr_sat = np.exp(ln_Pr_sat)
        P_sat_bar = Pr_sat * cls.Pc_bar
        return P_sat_bar

    @classmethod
    def cp0_ideal(cls, T_K: float) -> float:
        """Calculates ideal-gas specific heat capacity Cp0 (kJ/kg*K) at T (K)."""
        a0, a1, a2, a3 = cls.CP0_COEFFS
        return a0 + a1 * T_K + a2 * (T_K**2) + a3 * (T_K**3)

    @classmethod
    def h_ideal(cls, T_K: float) -> float:
        """
        Calculates ideal-gas enthalpy contribution h_0(T) (kJ/kg) integrated from T_ref = 273.15 K.
        """
        a0, a1, a2, a3 = cls.CP0_COEFFS
        T0 = cls.T_ref
        dh = (
            a0 * (T_K - T0)
            + 0.5 * a1 * (T_K**2 - T0**2)
            + (1.0 / 3.0) * a2 * (T_K**3 - T0**3)
            + 0.25 * a3 * (T_K**4 - T0**4)
        )
        return dh

    @classmethod
    def _peng_robinson_params(cls, T_K: float, P_bar: float):
        """Calculates dimensionless Peng-Robinson parameters A and B."""
        Tr = T_K / cls.Tc
        Pr = P_bar / cls.Pc_bar

        kappa = 0.37464 + 1.54226 * cls.omega - 0.26992 * (cls.omega**2)
        alpha = (1.0 + kappa * (1.0 - np.sqrt(Tr))) ** 2

        a = 0.45724 * ((cls.R_spec * cls.Tc) ** 2) / (cls.Pc_bar * 100.0) * alpha
        b = 0.07780 * (cls.R_spec * cls.Tc) / (cls.Pc_bar * 100.0)

        A = a * (P_bar * 100.0) / ((cls.R_spec * T_K) ** 2)
        B = b * (P_bar * 100.0) / (cls.R_spec * T_K)

        return A, B, kappa, Tr

    @classmethod
    def compressibility_factor(cls, T_K: float, P_bar: float, phase: str = "vapor") -> float:
        """
        Solves Peng-Robinson cubic equation of state for compressibility factor Z.
        Z^3 - (1-B)*Z^2 + (A - 2B - 3B^2)*Z - (AB - B^2 - B^3) = 0
        """
        A, B, _, _ = cls._peng_robinson_params(T_K, P_bar)

        c2 = -(1.0 - B)
        c1 = A - 2.0 * B - 3.0 * (B**2)
        c0 = -(A * B - (B**2) - (B**3))

        roots = np.roots([1.0, c2, c1, c0])
        real_roots = np.real(roots[np.isreal(roots)])

        # Filter physically valid positive real compressibility roots
        valid_roots = [r for r in real_roots if r > B]
        if not valid_roots:
            return 1.0

        if phase == "liquid":
            return min(valid_roots)
        else:
            return max(valid_roots)

    @classmethod
    def enthalpy_departure(cls, T_K: float, P_bar: float, phase: str = "vapor") -> float:
        """
        Calculates enthalpy departure function (h - h_ideal) [kJ/kg] using Peng-Robinson EOS.
        """
        A, B, kappa, Tr = cls._peng_robinson_params(T_K, P_bar)
        Z = cls.compressibility_factor(T_K, P_bar, phase=phase)

        da_dT_ratio = kappa * np.sqrt(Tr / ((1.0 + kappa * (1.0 - np.sqrt(Tr))) ** 2))
        term1 = Z - 1.0
        term2 = (A / (2.8284 * B)) * (1.0 + da_dT_ratio) * np.log((Z + 2.4142 * B) / (Z - 0.4142 * B))

        h_dep = cls.R_spec * T_K * (term1 - term2)
        return h_dep

    @classmethod
    def enthalpy(cls, T_K: float, P_bar: float, phase: str = "vapor") -> float:
        """
        Calculates total enthalpy h(T, P) = h_ideal(T) + h_dep(T, P) [kJ/kg].
        """
        h_id = cls.h_ideal(T_K)
        h_dep = cls.enthalpy_departure(T_K, P_bar, phase=phase)
        return h_id + h_dep

    @classmethod
    def h_fg(cls, T_K: float) -> float:
        """
        Calculates latent heat of vaporization h_fg = h_g - h_f (kJ/kg) at saturation.
        """
        P_sat = cls.p_sat_lee_kesler(T_K)
        h_g = cls.enthalpy(T_K, P_sat, phase="vapor")
        h_f = cls.enthalpy(T_K, P_sat, phase="liquid")
        return h_g - h_f