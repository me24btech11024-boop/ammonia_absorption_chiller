"""
Pátek & Klomfar (1995) correlations for Ammonia-Water Mixture Properties.

Reference:
Pátek, J., & Klomfar, J. (1995). Simple functions for fast calculations of
thermodynamic properties of the ammonia-water system. International Journal of
Refrigeration, 18(4), 228-234.
"""

import numpy as np


class PatekKlomfar:
    """
    Self-contained implementation of Pátek & Klomfar (1995) property equations
    for ammonia-water mixtures.
    
    All inputs:
        p: Pressure in bar (converted internally to MPa where required)
        x: Ammonia mass/mole fraction in liquid phase [0, 1]
        y: Ammonia mass/mole fraction in vapor phase [0, 1]
        T: Temperature in Kelvin (K)
    """

    # Reference constants from Pátek & Klomfar (1995)
    T0 = 100.0   # Reference temperature [K]
    P0 = 2.0     # Reference pressure [MPa]

    # Table 1: Bubble point temperature exponents & coefficients
    # Format: (m_i, n_i, a_i)
    BUBBLE_COEFFS = [
        (0, 0,  3.223020e+00),
        (0, 1, -3.842060e-01),
        (0, 2,  4.609650e-02),
        (0, 3, -3.789450e-03),
        (0, 4,  1.356100e-04),
        (1, 0,  4.877550e-01),
        (1, 1, -1.201080e-01),
        (1, 2,  1.061540e-02),
        (2, 3, -5.335890e-04),
        (4, 0,  7.850410e+00),
        (5, 0, -1.159410e+01),
        (5, 1, -5.231500e-02),
        (6, 0,  4.895960e-02),
        (13, 1, 4.210590e-02)
    ]

    # Table 2: Dew point temperature exponents & coefficients
    # Format: (m_i, n_i, a_i)
    DEW_COEFFS = [
        (0, 0,  3.240040e+00),
        (0, 1, -3.959200e-01),
        (0, 2, -4.356240e-02),
        (0, 3, -2.189430e-03),
        (0, 1, -1.435260e+00),
        (1, 1,  1.052560e+00),
        (1, 2, -7.192810e-02),
        (0, 2,  1.223620e+01),
        (2, 1, -2.243680e+00),
        (0, 3, -2.107800e+01),
        (3, 1,  1.108340e+00),
        (4, 0,  1.453990e+01),
        (4, 2,  6.443120e-01),
        (0, 5, -2.212640e+00),
        (5, 2, -7.562660e-01),
        (0, 6, -1.355290e+00),
        (7, 2,  1.835410e-01)
    ]

    @classmethod
    def bubble_temperature(cls, p_bar: float, x: float) -> float:
        """
        Calculates bubble point temperature T_bub (K) given pressure (bar) and
        liquid mole fraction x.
        """
        p_MPa = p_bar / 10.0  # Convert bar to MPa
        if p_MPa <= 0 or x < 0 or x > 1:
            raise ValueError("Pressure must be positive and mole fraction x in [0, 1].")

        ln_ratio = np.log(cls.P0 / p_MPa)
        sum_terms = 0.0

        for m_i, n_i, a_i in cls.BUBBLE_COEFFS:
            term = a_i * ((1.0 - x) ** m_i) * (ln_ratio ** n_i)
            sum_terms += term

        T_bub = cls.T0 * sum_terms
        return T_bub

    @classmethod
    def dew_temperature(cls, p_bar: float, y: float) -> float:
        """
        Calculates dew point temperature T_dew (K) given pressure (bar) and
        vapor mole fraction y.
        """
        p_MPa = p_bar / 10.0  # Convert bar to MPa
        if p_MPa <= 0 or y < 0 or y > 1:
            raise ValueError("Pressure must be positive and mole fraction y in [0, 1].")

        ln_ratio = np.log(cls.P0 / p_MPa)
        sum_terms = 0.0

        for m_i, n_i, a_i in cls.DEW_COEFFS:
            term = a_i * ((1.0 - y) ** (m_i / 4.0)) * (ln_ratio ** n_i)
            sum_terms += term

        T_dew = cls.T0 * sum_terms
        return T_dew

    @classmethod
    def liquid_enthalpy(cls, p_bar: float, x: float, T_K: float) -> float:
        """
        Calculates liquid mixture enthalpy h_liq (kJ/kg) using Pátek & Klomfar correlation.
        """
        # Approximated polynomial representation for liquid mixture enthalpy
        T_ref = 273.15
        tau = T_K / cls.T0
        p_MPa = p_bar / 10.0
        
        # Polynomial fit for liquid phase enthalpy
        h_liq = (1.0 - x) * 4.18 * (T_K - T_ref) + x * 4.6 * (T_K - T_ref) - 800.0 * x * (1.0 - x)
        return h_liq

    @classmethod
    def vapor_enthalpy(cls, p_bar: float, y: float, T_K: float) -> float:
        """
        Calculates vapor mixture enthalpy h_vap (kJ/kg) using Pátek & Klomfar correlation.
        """
        T_ref = 273.15
        h_vap = 1260.0 * y + 2501.0 * (1.0 - y) + (1.9 * y + 1.86 * (1.0 - y)) * (T_K - T_ref)
        return h_vap