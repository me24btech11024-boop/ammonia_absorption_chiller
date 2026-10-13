"""
Water Decanter (Condensate Collection System) Model.

Calculates condensate mass and volume rates and sizes a central
holding tank from the governing daily condensate-generation rate.
"""


class DecanterModel:
    """Condensate recovery and central holding-tank sizing."""

    WATER_DENSITY_KG_M3 = 1000.0

    @staticmethod
    def _validate_nonnegative_finite(value, name):
        """Reject negative, non-finite, or non-numeric input values."""
        import math

        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"{name} must be a finite, nonnegative number.")

        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f"{name} must be a finite, nonnegative number.")

    @classmethod
    def calculate_condensate_rates(cls, m_condensate_district_kg_s: float):
        """
        Convert district condensate mass flow to hourly and daily rates.

        Args:
            m_condensate_district_kg_s: District condensate mass flow (kg/s).

        Returns:
            Dictionary containing kg/h, m³/h, and m³/day.
        """
        cls._validate_nonnegative_finite(
            m_condensate_district_kg_s,
            "District condensate mass flow",
        )

        mass_rate_kg_h = m_condensate_district_kg_s * 3600.0
        vol_rate_m3_h = mass_rate_kg_h / cls.WATER_DENSITY_KG_M3
        daily_vol_m3_day = vol_rate_m3_h * 24.0

        return {
            "mass_rate_kg_h": mass_rate_kg_h,
            "vol_rate_m3_h": vol_rate_m3_h,
            "daily_vol_m3_day": daily_vol_m3_day,
        }

    @classmethod
    def size_decanter_tank(
        cls,
        max_daily_vol_m3_day: float,
        retention_hours: float = 12.0,
    ):
        """
        Size a holding tank using a constant design-rate assumption.

        Args:
            max_daily_vol_m3_day: Governing daily condensate volume (m³/day).
            retention_hours: Required storage retention period (hours).

        Returns:
            Dictionary containing the required tank volume and assumptions.
        """
        cls._validate_nonnegative_finite(
            max_daily_vol_m3_day,
            "Maximum daily condensate volume",
        )
        cls._validate_nonnegative_finite(
            retention_hours,
            "Retention period",
        )

        if retention_hours > 24.0:
            raise ValueError(
                "Retention period must not exceed 24 hours "
                "for this daily-rate sizing model."
            )

        hourly_vol_m3_h = max_daily_vol_m3_day / 24.0
        required_tank_volume_m3 = hourly_vol_m3_h * retention_hours

        return {
            "required_tank_volume_m3": required_tank_volume_m3,
            "retention_hours": retention_hours,
            "design_daily_volume_m3_day": max_daily_vol_m3_day,
            "design_hourly_volume_m3_h": hourly_vol_m3_h,
            "water_density_kg_m3": cls.WATER_DENSITY_KG_M3,
            "model_status": (
                "Idealized sizing: constant condensate rate, "
                "no freeboard or operating reserve included."
            ),
        }
