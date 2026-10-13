"""
Zone Psychrometrics and ADP/Bypass-Factor Cooling Coil Model.
"""

from core.moist_air import MoistAir


class ZonePsychrometrics:
    """
    Solves psychrometric states, ADP/BF coil performance, and district loads.
    """

    N_ZONES = 400
    AIR_DENSITY = 1.204  # kg/m^3
    CP_AIR = 1.006       # kJ/kg*K

    @classmethod
    def solve_condition(cls, T_db_out: float, rh_out: float, solar_on: float, fresh_air_frac: float = 0.30):
        # Room Design Parameters
        T_room = 22.0
        rh_room_target = 75.0

        # Mixing Box
        T_mix = (1.0 - fresh_air_frac) * T_room + fresh_air_frac * T_db_out
        w_room = MoistAir.humidity_ratio(T_room, rh_room_target)
        w_out = MoistAir.humidity_ratio(T_db_out, rh_out)
        w_mix = (1.0 - fresh_air_frac) * w_room + fresh_air_frac * w_out

        # Zone Sensible & Latent Heat Gains (per zone)
        Q_sens_zone = 12.0 if solar_on else 5.0  # kW
        Q_lat_zone = 3.0  # kW

        # Supply Air State Determination
        m_air_zone = 1.5  # kg/s per zone
        dT_supply = Q_sens_zone / (m_air_zone * cls.CP_AIR)
        T_supply = T_room - dT_supply

        w_supply = w_room - (Q_lat_zone / (m_air_zone * 2501.0))

        # ADP and Bypass Factor Model
        ADP = T_supply - 2.0
        BF = 0.15

        # Chiller Evaporator Load Calculation
        h_mix = MoistAir.enthalpy(T_mix, w_mix)
        h_supply = MoistAir.enthalpy(T_supply, w_supply)

        # Single Zone Load
        m_air_total_zone = m_air_zone / (1.0 - BF)
        Q_evap_zone_kW = m_air_total_zone * (h_mix - h_supply)

        # District Load (400 Zones)
        Q_evap_chiller_kW = Q_evap_zone_kW * cls.N_ZONES

        # Clamp cooling load to 0 if building is in heating/free-cooling mode
        if Q_evap_chiller_kW < 0:
            Q_evap_chiller_kW = 0.0

        # Condensate Mass Flow Rates
        w_ADP = MoistAir.humidity_ratio(ADP, 100.0)
        m_cond_zone = m_air_total_zone * max(w_mix - w_ADP, 0.0)
        m_cond_district = m_cond_zone * cls.N_ZONES

        # Room RH Check
        rh_room_actual = MoistAir.relative_humidity(T_room, w_room)
        rh_violated = rh_room_actual > 75.0

        return {
            "T_mix": T_mix,
            "w_mix": w_mix,
            "T_supply": T_supply,
            "w_supply": w_supply,
            "ADP": ADP,
            "rh_room": rh_room_actual,
            "rh_violated": rh_violated,
            "Q_evap_chiller_kW": Q_evap_chiller_kW,
            "m_condensate_zone_kg_s": m_cond_zone,
            "m_condensate_district_kg_s": m_cond_district,
        }