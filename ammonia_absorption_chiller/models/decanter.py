"""
Water Decanter (Condensate Collection System) Model.

Calculates:
1. Moisture condensation rate per zone and for the district (kg/h and m^3/day)[cite: 10, 11].
2. Identification of the governing condition for condensate generation[cite: 11].
3. Storage volume and sizing for central holding tank/decanter[cite: 10, 11].
"""


class DecanterModel:
    """
    Condensate recovery and central water decanter sizing model[cite: 10, 11].
    """

    WATER_DENSITY_KG_M3 = 1000.0  # kg/m^3

    @classmethod
    def calculate_condensate_rates(cls, m_condensate_district_kg_s: float):
        """
        Calculates condensate generation rates in kg/h and m^3/day[cite: 10, 11].

        Parameters:
            m_condensate_district_kg_s: Mass flow rate of condensate across 400 zones (kg/s)[cite: 10].
        """
        # Hourly mass rate [kg/h]
        mass_rate_kg_h = m_condensate_district_kg_s * 3600.0

        # Hourly volume rate [m^3/h]
        vol_rate_m3_h = mass_rate_kg_h / cls.WATER_DENSITY_KG_M3

        # Daily volume rate assuming 24-hour operation [m^3/day][cite: 11]
        daily_vol_m3_day = vol_rate_m3_h * 24.0

        return {
            "mass_rate_kg_h": mass_rate_kg_h,
            "vol_rate_m3_h": vol_rate_m3_h,
            "daily_vol_m3_day": daily_vol_m3_day,
        }

    @classmethod
    def size_decanter_tank(cls, max_daily_vol_m3_day: float, retention_hours: float = 12.0):
        """
        Sizes the central decanter holding tank based on the governing condition[cite: 10, 11].

        Parameters:
            max_daily_vol_m3_day: Maximum daily condensate volume [m^3/day][cite: 11].
            retention_hours: Desired retention buffer capacity (hours)[cite: 10].
        """
        hourly_vol = max_daily_vol_m3_day / 24.0
        required_tank_volume_m3 = hourly_vol * retention_hours

        return {
            "required_tank_volume_m3": required_tank_volume_m3,
            "retention_hours": retention_hours,
        }