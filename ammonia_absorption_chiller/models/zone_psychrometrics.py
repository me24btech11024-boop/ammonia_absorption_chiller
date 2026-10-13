"""
Zone psychrometrics and cooling-coil model for Project 2.

Basis:
- 400 identical zones
- 0.40 kg/s dry air per zone
- 70% recirculated room air and 30% outdoor air by dry-air mass
- Room dry-bulb setpoint: 22 C
- Coil bypass factor: 0.15
- Atmospheric pressure: 101.325 kPa

This module uses the ideal-gas moist-air relations in core.moist_air.
It iterates room humidity using the zone moisture balance and calculates
coil leaving humidity from the ADP/bypass-factor relation.
"""

from __future__ import annotations

from typing import Dict, Any

from core.moist_air import MoistAir


class ZonePsychrometrics:
    """Solve a zone/AHU condition and scale the load to the district."""

    N_ZONES = 400
    T_ROOM_C = 22.0
    RH_LIMIT_PERCENT = 75.0
    M_DA_ZONE_KG_S = 0.40
    BF = 0.15
    FRESH_AIR_DEFAULT = 0.30
    PRESSURE_KPA = 101.325

    # Assignment data
    UA_ENVELOPE_W_K = 80.0
    SOLAR_GAIN_W = 1200.0
    OCCUPANTS = 4
    OCCUPANT_SENSIBLE_W = 75.0
    OCCUPANT_LATENT_W = 55.0
    LIGHTING_EQUIPMENT_W_M2 = 20.0
    FLOOR_AREA_M2 = 40.0

    @classmethod
    def _temperature_from_enthalpy(cls, h_kj_kg_da: float, w: float) -> float:
        """Invert h = 1.006*T + w*(2501 + 1.86*T), with T in C."""
        return (h_kj_kg_da - 2501.0 * w) / (1.006 + 1.86 * w)

    @classmethod
    def _sensible_zone_load_kw(cls, T_out_C: float, solar_on: bool) -> float:
        """Net zone sensible load in kW; negative means net envelope cooling."""
        envelope_W = cls.UA_ENVELOPE_W_K * (T_out_C - cls.T_ROOM_C)
        solar_W = cls.SOLAR_GAIN_W if solar_on else 0.0
        occupant_W = cls.OCCUPANTS * cls.OCCUPANT_SENSIBLE_W
        lighting_W = cls.LIGHTING_EQUIPMENT_W_M2 * cls.FLOOR_AREA_M2
        return (envelope_W + solar_W + occupant_W + lighting_W) / 1000.0

    @classmethod
    def solve_condition(
        cls,
        T_db_out: float,
        rh_out: float,
        solar_on: bool,
        fresh_air_frac: float = 0.30,
    ) -> Dict[str, Any]:
        """
        Solve one ambient condition.

        Args:
            T_db_out: Outdoor dry-bulb temperature, deg C.
            rh_out: Outdoor relative humidity, percent (e.g. 50, not 0.50).
            solar_on: True during the specified day conditions.
            fresh_air_frac: Outdoor-air fraction by dry-air mass.

        Returns:
            A dictionary of zone/district states, loads, condensate, and flags.
            Cooling duty is zero-clamped for reporting when net coil duty is
            negative; the signed raw duty is also returned.
        """
        if not 0.0 <= fresh_air_frac <= 1.0:
            raise ValueError("fresh_air_frac must be between 0 and 1.")
        if not 0.0 <= rh_out <= 100.0:
            raise ValueError("rh_out must be in percent, between 0 and 100.")

        P = cls.PRESSURE_KPA
        T_room = cls.T_ROOM_C
        m_da = cls.M_DA_ZONE_KG_S
        bf = cls.BF

        # Outdoor state
        w_out = MoistAir.humidity_ratio(T_db_out, rh_out, P)
        h_out = MoistAir.enthalpy(T_db_out, w_out)

        # Fixed zone loads from the assignment
        Q_sensible_zone_kW = cls._sensible_zone_load_kw(T_db_out, solar_on)
        Q_latent_zone_kW = (
            cls.OCCUPANTS * cls.OCCUPANT_LATENT_W / 1000.0
        )

        # Moisture generation from occupants, kg water/s per zone.
        # Latent heat of moist air at room temperature is used consistently
        # with h_v(T) = 2501 + 1.86*T, kJ/kg.
        h_fg_room_kJ_kg = 2501.0 + 1.86 * T_room
        m_water_generation = Q_latent_zone_kW / h_fg_room_kJ_kg

        # Supply temperature from the zone sensible balance.
        # This assumes the AHU supply stream handles the net sensible load.
        T_supply = T_room - Q_sensible_zone_kW / (
            m_da * (1.006 + 1.86 * 0.010)
        )

        # The bypass-factor relation gives:
        # T_supply = ADP + BF*(T_mix - ADP)
        # Thus ADP = (T_supply - BF*T_mix)/(1-BF).
        # T_mix depends on room humidity through enthalpy mixing, so solve
        # the room moisture balance iteratively.
        w_room = MoistAir.humidity_ratio(T_room, 50.0, P)
        converged = False
        max_iterations = 500
        tolerance = 1e-10

        for iteration in range(1, max_iterations + 1):
            h_room = MoistAir.enthalpy(T_room, w_room)

            # Adiabatic mixing on a dry-air mass basis.
            w_mix = (
                (1.0 - fresh_air_frac) * w_room
                + fresh_air_frac * w_out
            )
            h_mix = (
                (1.0 - fresh_air_frac) * h_room
                + fresh_air_frac * h_out
            )
            T_mix = cls._temperature_from_enthalpy(h_mix, w_mix)

            ADP = (T_supply - bf * T_mix) / (1.0 - bf)

            # Cooling-coil model requires a physically meaningful ADP below
            # the entering mixed-air temperature.
            if ADP >= T_mix:
                return {
                    "feasible": False,
                    "status": "No cooling-coil solution: ADP is not below mixed-air temperature.",
                    "T_mix": T_mix,
                    "w_mix": w_mix,
                    "h_mix": h_mix,
                    "T_supply": T_supply,
                    "ADP": ADP,
                    "Q_sensible_zone_kW": Q_sensible_zone_kW,
                    "Q_latent_zone_kW": Q_latent_zone_kW,
                    "Q_evap_chiller_kW": 0.0,
                    "rh_room": MoistAir.relative_humidity(T_room, w_room, P),
                }

            w_adp = MoistAir.humidity_ratio(ADP, 100.0, P)
            w_supply = w_adp + bf * (w_mix - w_adp)

            # Zone moisture balance:
            # m_da*w_room = m_da*w_supply + water generated in zone.
            w_room_new = w_supply + m_water_generation / m_da

            if w_room_new < 0.0:
                raise ValueError("Calculated room humidity ratio is negative.")

            if abs(w_room_new - w_room) < tolerance:
                w_room = w_room_new
                converged = True
                break

            # Damped iteration improves stability for high outdoor humidity.
            w_room = 0.5 * w_room + 0.5 * w_room_new

        # Recompute all states with the final room humidity ratio.
        h_room = MoistAir.enthalpy(T_room, w_room)
        w_mix = (
            (1.0 - fresh_air_frac) * w_room
            + fresh_air_frac * w_out
        )
        h_mix = (
            (1.0 - fresh_air_frac) * h_room
            + fresh_air_frac * h_out
        )
        T_mix = cls._temperature_from_enthalpy(h_mix, w_mix)
        ADP = (T_supply - bf * T_mix) / (1.0 - bf)

        if ADP >= T_mix:
            return {
                "feasible": False,
                "status": "No cooling-coil solution: ADP is not below mixed-air temperature.",
                "T_mix": T_mix,
                "w_mix": w_mix,
                "h_mix": h_mix,
                "T_supply": T_supply,
                "ADP": ADP,
                "Q_sensible_zone_kW": Q_sensible_zone_kW,
                "Q_latent_zone_kW": Q_latent_zone_kW,
                "Q_evap_chiller_kW": 0.0,
                "rh_room": MoistAir.relative_humidity(T_room, w_room, P),
            }

        w_adp = MoistAir.humidity_ratio(ADP, 100.0, P)
        w_supply = w_adp + bf * (w_mix - w_adp)
        h_supply = MoistAir.enthalpy(T_supply, w_supply)

        # Coil duty and condensate are based on the actual entering/leaving
        # air states. The enthalpy difference includes the condensate effect
        # through the moist-air enthalpy formulation used by the assignment.
        m_cond_zone_kg_s = max(m_da * (w_mix - w_supply), 0.0)
        h_condensate_kJ_kg = 4.186 * ADP  # liquid-water enthalpy, kJ/kg, 0 C reference
        Q_coil_zone_kW_signed = (
            m_da * (h_mix - h_supply)
            - m_cond_zone_kg_s * h_condensate_kJ_kg
        )
        Q_coil_zone_kW = max(Q_coil_zone_kW_signed, 0.0)

        rh_room_actual = MoistAir.relative_humidity(T_room, w_room, P)
        rh_violated = rh_room_actual > cls.RH_LIMIT_PERCENT + 1e-6

        # District heat gains: assignment specifies 5% distribution gain and
        # 2% chilled-water pump heat, applied to delivered cooling.
        Q_coil_district_kW = Q_coil_zone_kW * cls.N_ZONES
        Q_evap_chiller_kW = Q_coil_district_kW * 1.05 * 1.02

        # Convert condensate rates to convenient units.
        m_cond_district_kg_s = m_cond_zone_kg_s * cls.N_ZONES

        return {
            "feasible": converged,
            "status": "converged" if converged else "Room humidity iteration did not converge",
            "iterations": iteration,
            "T_mix": T_mix,
            "w_out": w_out,
            "h_out": h_out,
            "w_room": w_room,
            "h_room": h_room,
            "w_mix": w_mix,
            "h_mix": h_mix,
            "T_supply": T_supply,
            "w_supply": w_supply,
            "h_supply": h_supply,
            "ADP": ADP,
            "w_ADP": w_adp,
            "rh_room": rh_room_actual,
            "rh_violated": rh_violated,
            "Q_sensible_zone_kW": Q_sensible_zone_kW,
            "Q_latent_zone_kW": Q_latent_zone_kW,
            "Q_coil_zone_kW_signed": Q_coil_zone_kW_signed,
            "Q_coil_zone_kW": Q_coil_zone_kW,
            "Q_coil_district_kW": Q_coil_district_kW,
            "Q_evap_chiller_kW": Q_evap_chiller_kW,
            "m_condensate_zone_kg_s": m_cond_zone_kg_s,
            "m_condensate_district_kg_s": m_cond_district_kg_s,
            "m_condensate_zone_kg_h": m_cond_zone_kg_s * 3600.0,
            "m_condensate_district_kg_h": m_cond_district_kg_s * 3600.0,
            "m_condensate_district_m3_day": m_cond_district_kg_s * 86400.0 / 1000.0,
            "fresh_air_frac": fresh_air_frac,
            "rh_violated": rh_violated,
        }
