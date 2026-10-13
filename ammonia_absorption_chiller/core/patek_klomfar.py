"""
Ammonia-water mixture property utilities for Project 2.

Important model-status note:
- Bubble/dew temperature functions use the coefficient tables currently
  present in this project. Those coefficients have not been independently
  verified against the published Pátek-Klomfar tables.
- Liquid and vapor enthalpies are simplified engineering approximations,
  not a verified implementation of the complete Pátek-Klomfar enthalpy
  correlations.

Units:
- Pressure: bar
- Temperature: K
- Composition: ammonia mass fraction, [0, 1]
- Enthalpy: kJ/kg
"""

import math
import numpy as np


class PatekKlomfar:
    """Ammonia-water property approximations with explicit limitations."""

    T0 = 100.0
    P0 = 2.0  # MPa

    BUBBLE_COEFFS = [
        (0, 0, 3.223020e+00),
        (0, 1, -3.842060e-01),
        (0, 2, 4.609650e-02),
        (0, 3, -3.789450e-03),
        (0, 4, 1.356100e-04),
        (1, 0, 4.877550e-01),
        (1, 1, -1.201080e-01),
        (1, 2, 1.061540e-02),
        (2, 3, -5.335890e-04),
        (4, 0, 7.850410e+00),
        (5, 0, -1.159410e+01),
        (5, 1, -5.231500e-02),
        (6, 0, 4.895960e-02),
        (13, 1, 4.210590e-02),
    ]

    DEW_COEFFS = [
        (0, 0, 3.240040e+00),
        (0, 1, -3.959200e-01),
        (0, 2, -4.356240e-02),
        (0, 3, -2.189430e-03),
        (0, 1, -1.435260e+00),
        (1, 1, 1.052560e+00),
        (1, 2, -7.192810e-02),
        (0, 2, 1.223620e+01),
        (2, 1, -2.243680e+00),
        (0, 3, -2.107800e+01),
        (3, 1, 1.108340e+00),
        (4, 0, 1.453990e+01),
        (4, 2, 6.443120e-01),
        (0, 5, -2.212640e+00),
        (5, 2, -7.562660e-01),
        (0, 6, -1.355290e+00),
        (7, 2, 1.835410e-01),
    ]

    @staticmethod
    def _validate_pressure(p_bar):
        if (
            isinstance(p_bar, bool)
            or not isinstance(p_bar, (int, float))
            or not math.isfinite(p_bar)
            or p_bar <= 0.0
        ):
            raise ValueError("Pressure must be a finite positive value in bar.")

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
    def _validate_fraction(x, name="ammonia mass fraction"):
        if (
            isinstance(x, bool)
            or not isinstance(x, (int, float))
            or not math.isfinite(x)
            or not 0.0 <= x <= 1.0
        ):
            raise ValueError(f"{name} must be between 0 and 1.")

    @classmethod
    def bubble_temperature(cls, p_bar: float, x: float) -> float:
        """Estimate bubble temperature in K from the project's coefficient table."""
        cls._validate_pressure(p_bar)
        cls._validate_fraction(x, "liquid ammonia fraction")

        p_MPa = p_bar / 10.0
        ln_ratio = math.log(cls.P0 / p_MPa)

        value = sum(
            a_i * (1.0 - x) ** m_i * ln_ratio ** n_i
            for m_i, n_i, a_i in cls.BUBBLE_COEFFS
        )
        T_bubble_K = cls.T0 * value

        if not math.isfinite(T_bubble_K) or T_bubble_K <= 0.0:
            raise ArithmeticError("Bubble-temperature estimate is invalid.")

        return T_bubble_K

    @classmethod
    def dew_temperature(cls, p_bar: float, y: float) -> float:
        """Estimate dew temperature in K from the project's coefficient table."""
        cls._validate_pressure(p_bar)
        cls._validate_fraction(y, "vapor ammonia fraction")

        p_MPa = p_bar / 10.0
        ln_ratio = math.log(cls.P0 / p_MPa)

        value = sum(
            a_i * (1.0 - y) ** (m_i / 4.0) * ln_ratio ** n_i
            for m_i, n_i, a_i in cls.DEW_COEFFS
        )
        T_dew_K = cls.T0 * value

        if not math.isfinite(T_dew_K) or T_dew_K <= 0.0:
            raise ArithmeticError("Dew-temperature estimate is invalid.")

        return T_dew_K

    @classmethod
    def liquid_enthalpy(cls, p_bar: float, x: float, T_K: float) -> float:
        """
        Approximate liquid-solution enthalpy in kJ/kg.

        This is a simplified estimate, not the full Pátek-Klomfar correlation.
        Pressure is validated for API consistency but is not used in this fit.
        """
        cls._validate_pressure(p_bar)
        cls._validate_fraction(x, "liquid ammonia fraction")
        cls._validate_temperature(T_K)

        T_ref = 273.15
        h_liq = (
            (1.0 - x) * 4.18 * (T_K - T_ref)
            + x * 4.6 * (T_K - T_ref)
            - 800.0 * x * (1.0 - x)
        )

        if not math.isfinite(h_liq):
            raise ArithmeticError("Liquid-solution enthalpy estimate is invalid.")

        return h_liq

    @classmethod
    def vapor_enthalpy(cls, p_bar: float, y: float, T_K: float) -> float:
        """
        Approximate vapor-mixture enthalpy in kJ/kg.

        This is a simplified estimate, not the full Pátek-Klomfar correlation.
        Pressure is validated for API consistency but is not used in this fit.
        """
        cls._validate_pressure(p_bar)
        cls._validate_fraction(y, "vapor ammonia fraction")
        cls._validate_temperature(T_K)

        T_ref = 273.15
        h_vap = (
            1260.0 * y
            + 2501.0 * (1.0 - y)
            + (1.9 * y + 1.86 * (1.0 - y)) * (T_K - T_ref)
        )

        if not math.isfinite(h_vap):
            raise ArithmeticError("Vapor-mixture enthalpy estimate is invalid.")

        return h_vap
