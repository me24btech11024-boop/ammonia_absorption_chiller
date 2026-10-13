"""
Core thermodynamic property evaluation modules.

Contains self-implemented thermodynamic property functions:
- moist_air: Ideal-gas moist air psychrometrics and wet-bulb solver.
- pure_ammonia: Real-gas EOS departure functions & Lee-Kesler/Pitzer vapor pressure.
- patek_klomfar: Pátek & Klomfar (1995) correlations for NH3-H2O mixture properties.
"""

from .moist_air import MoistAir
from .pure_ammonia import PureAmmoniaEOS
from .patek_klomfar import PatekKlomfar

__all__ = ["MoistAir", "PureAmmoniaEOS", "PatekKlomfar"]