"""
Pure ammonia thermodynamic properties.

Provides:
- Lee-Kesler/Pitzer saturation-pressure correlation
- Peng-Robinson compressibility factor
- Ideal-gas enthalpy and real-gas enthalpy departure
- Saturated vaporization latent heat

Units:
- Temperature: K
- Pressure: bar
- Enthalpy and latent heat: kJ/kg
- Specific gas constant: kJ/(kg K)
"""

import math
import numpy as np


class PureAmmoniaEOS:
    """Pure-ammonia property model using a Peng-Robinson EOS."""

    Tc = 405.56
    Pc_bar = 113.33
    omega = 0.256
    MW = 17.031

    R_gas = 8.3144626  # kJ/(kmol K)
    R_spec = R_gas / MW  # kJ/(kg K)

    T_ref = 273.15
    P_ref_bar = 1.01325

    CP0_COEFFS = [1.61, 2.80e-3, -4.0e-7, 2.0e-11]

    @staticmethod
    def _validate_temperature(T_K):
        if (
            isinstance(T_K, bool)
            or not isinstance(T_K, (int, float))
            or not math.isfinite(T_K)
            or T_K <= 0.0
        ):
            raise ValueError("Temperature must be a finite positive value in K.")

    @staticmethod
    def _validate_pressure(P_bar):
        if (
            isinstance(P_bar, bool)
            or not isinstance(P_bar, (int, float))
            or not math.isfinite(P_bar)
            or P_bar <= 0.0
        ):
            raise ValueError("Pressure must be a finite positive value in bar.")

    @classmethod
    def p_sat_lee_kesler(cls, T_K: float) -> float:
        """Return saturation pressure in bar using the Pitzer correlation."""
        cls._validate_temperature(T_K)

        if T_K >= cls.Tc:
            raise ValueError(
                "Saturation pressure is undefined above the critical temperature."
            )

        Tr = T_K / cls.Tc
        f0 = (
            5.92714
            - 6.09648 / Tr
            - 1.28862 * math.log(Tr)
            + 0.169347 * Tr**6
        )
        f1 = (
            15.2518
            - 15.6875 / Tr
            - 13.4721 * math.log(Tr)
            + 0.43577 * Tr**6
        )

        ln_Pr_sat = f0 + cls.omega * f1
        P_sat_bar = cls.Pc_bar * math.exp(ln_Pr_sat)

        if not math.isfinite(P_sat_bar) or P_sat_bar <= 0.0:
            raise ArithmeticError("Saturation-pressure calculation failed.")

        return P_sat_bar

    @classmethod
    def cp0_ideal(cls, T_K: float) -> float:
        """Ideal-gas specific heat capacity in kJ/(kg K)."""
        cls._validate_temperature(T_K)
        a0, a1, a2, a3 = cls.CP0_COEFFS
        return a0 + a1 * T_K + a2 * T_K**2 + a3 * T_K**3

    @classmethod
    def h_ideal(cls, T_K: float) -> float:
        """Ideal-gas enthalpy relative to T_ref, in kJ/kg."""
        cls._validate_temperature(T_K)
        a0, a1, a2, a3 = cls.CP0_COEFFS
        T0 = cls.T_ref

        return (
            a0 * (T_K - T0)
            + 0.5 * a1 * (T_K**2 - T0**2)
            + (a2 / 3.0) * (T_K**3 - T0**3)
            + 0.25 * a3 * (T_K**4 - T0**4)
        )

    @classmethod
    def _peng_robinson_params(cls, T_K: float, P_bar: float):
        """Return dimensionless PR parameters A, B and alpha parameters."""
        cls._validate_temperature(T_K)
        cls._validate_pressure(P_bar)

        Tr = T_K / cls.Tc
        Pr = P_bar / cls.Pc_bar

        kappa = (
            0.37464
            + 1.54226 * cls.omega
            - 0.26992 * cls.omega**2
        )
        sqrt_Tr = math.sqrt(Tr)
        alpha = (1.0 + kappa * (1.0 - sqrt_Tr)) ** 2

        # R_spec is kJ/(kg K); 1 bar = 100 kPa.
        a = (
            0.45724
            * (cls.R_spec * cls.Tc) ** 2
            / (cls.Pc_bar * 100.0)
            * alpha
        )
        b = (
            0.07780
            * cls.R_spec
            * cls.Tc
            / (cls.Pc_bar * 100.0)
        )

        A = a * (P_bar * 100.0) / (cls.R_spec * T_K) ** 2
        B = b * (P_bar * 100.0) / (cls.R_spec * T_K)

        return A, B, kappa, Tr

    @classmethod
    def compressibility_factor(
        cls, T_K: float, P_bar: float, phase: str = "vapor"
    ) -> float:
        """Solve the PR cubic and select the requested physical root."""
        if phase not in ("vapor", "liquid"):
            raise ValueError("phase must be 'vapor' or 'liquid'.")

        A, B, _, _ = cls._peng_robinson_params(T_K, P_bar)

        coefficients = [
            1.0,
            -(1.0 - B),
            A - 2.0 * B - 3.0 * B**2,
            -(A * B - B**2 - B**3),
        ]

        roots = np.roots(coefficients)
        real_roots = [
            float(root.real)
            for root in roots
            if abs(root.imag) < 1e-8 and root.real > B
        ]

        if not real_roots:
            raise ArithmeticError(
                f"No physical Peng-Robinson root at T={T_K} K, P={P_bar} bar."
            )

        return min(real_roots) if phase == "liquid" else max(real_roots)

    @classmethod
    def enthalpy_departure(
        cls, T_K: float, P_bar: float, phase: str = "vapor"
    ) -> float:
        """Calculate the PR enthalpy departure in kJ/kg."""
        A, B, kappa, Tr = cls._peng_robinson_params(T_K, P_bar)
        Z = cls.compressibility_factor(T_K, P_bar, phase=phase)

        sqrt_Tr = math.sqrt(Tr)
        denominator = 1.0 + kappa * (1.0 - sqrt_Tr)
        da_dT_ratio = (
            kappa * sqrt_Tr / denominator
            if abs(denominator) > 1e-14
            else 0.0
        )

        log_argument_denominator = Z - 0.41421356237 * B
        log_argument_numerator = Z + 2.41421356237 * B

        if (
            log_argument_denominator <= 0.0
            or log_argument_numerator <= 0.0
            or B <= 0.0
        ):
            raise ArithmeticError("Invalid logarithm in PR departure function.")

        log_term = math.log(
            log_argument_numerator / log_argument_denominator
        )

        h_departure = cls.R_spec * T_K * (
            Z - 1.0
            - (A / (2.0 * math.sqrt(2.0) * B))
            * (1.0 + da_dT_ratio)
            * log_term
        )

        if not math.isfinite(h_departure):
            raise ArithmeticError("PR enthalpy-departure calculation failed.")

        return h_departure

    @classmethod
    def enthalpy(
        cls, T_K: float, P_bar: float, phase: str = "vapor"
    ) -> float:
        """Return ideal-gas enthalpy plus PR departure, in kJ/kg."""
        return (
            cls.h_ideal(T_K)
            + cls.enthalpy_departure(T_K, P_bar, phase=phase)
        )

    @classmethod
    def h_fg(cls, T_K: float) -> float:
        """Return saturated latent heat of vaporization in kJ/kg."""
        P_sat_bar = cls.p_sat_lee_kesler(T_K)

        h_vapor = cls.enthalpy(T_K, P_sat_bar, phase="vapor")
        h_liquid = cls.enthalpy(T_K, P_sat_bar, phase="liquid")
        latent_heat = h_vapor - h_liquid

        if not math.isfinite(latent_heat) or latent_heat <= 0.0:
            raise ArithmeticError(
                f"Invalid latent heat at T={T_K} K: {latent_heat} kJ/kg."
            )

        return latent_heat
