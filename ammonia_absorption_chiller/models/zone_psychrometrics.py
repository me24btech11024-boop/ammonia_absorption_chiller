"""
Zone psychrometrics and cooling-coil model for Project 2.

Basis:
- 400 identical zones
- 0.40 kg/s dry air per zone
- Outdoor-air fraction specified by fresh_air_frac (baseline 0.30)
- Room dry-bulb setpoint: 22 C
- Coil bypass factor: 0.15
- Atmospheric pressure: 101.325 kPa

The model calculates mixed-air states, supply temperature, coil ADP,
room humidity, cooling duty, and condensate generation.

Important:
- A returned feasible=False means the calculated operating point must not
  be treated as a valid design result.
- This is a simplified steady-state zone and coil model.
"""

from __future__ import annotations

import math
from typing import Any, Dict

from core.moist_air import MoistAir


class ZonePsychrometrics:
    """Solve one zone/AHU condition and scale results to the district."""

    N_ZONES = 400
    T_ROOM_C = 22.0
    RH_LIMIT_PERCENT = 75.0
    M_DA_ZONE_KG_S = 0.40
    BF = 0.15
    FRESH_AIR_DEFAULT = 0.30
    PRESSURE_KPA = 101.325

    # Assignment zone-load data
    UA_ENVELOPE_W_K = 80.0
    SOLAR_GAIN_W = 1200.0
    OCCUPANTS = 4
    OCCUPANT_SENSIBLE_W = 75.0
    OCCUPANT_LATENT_W = 55.0
    LIGHTING_EQUIPMENT_W_M2 = 20.0
    FLOOR_AREA_M2 = 40.0

    MAX_ITERATIONS = 500
    HUMIDITY_TOLERANCE = 1e-10

    @classmethod
    def _temperature_from_enthalpy(cls, h_kJ_kg_da: float, w: float) -> float:
        """Invert the moist-air enthalpy relation; temperature in °C."""
        return (h_kJ_kg_da - 2501.0 * w) / (1.006 + 1.86 * w)

    @classmethod
    def _sensible_zone_load_kw(cls, T_out_C: float, solar_on: bool) -> float:
        """Calculate the net zone sensible load in kW."""
        envelope_W = cls.UA_ENVELOPE_W_K * (T_out_C - cls.T_ROOM_C)
        solar_W = cls.SOLAR_GAIN_W if solar_on else 0.0
        occupants_W = cls.OCCUPANTS * cls.OCCUPANT_SENSIBLE_W
        lighting_W = cls.LIGHTING_EQUIPMENT_W_M2 * cls.FLOOR_AREA_M2

        return (
            envelope_W + solar_W + occupants_W + lighting_W
        ) / 1000.0

    @classmethod
    def _base_result(
        cls,
        T_db_out: float,
        rh_out: float,
        solar_on: bool,
        fresh_air_frac: float,
    ) -> Dict[str, Any]:
        """Return common output fields so all results have consistent keys."""
        return {
            "feasible": False,
            "status": "not solved",
            "iterations": 0,
            "T_db_out": T_db_out,
            "rh_out": rh_out,
            "solar_on": solar_on,
            "fresh_air_frac": fresh_air_frac,
            "T_mix": float("nan"),
            "w_out": float("nan"),
            "h_out": float("nan"),
            "w_room": float("nan"),
            "h_room": float("nan"),
            "w_mix": float("nan"),
            "h_mix": float("nan"),
            "T_supply": float("nan"),
            "w_supply": float("nan"),
            "h_supply": float("nan"),
            "ADP": float("nan"),
            "w_ADP": float("nan"),
            "rh_room": float("nan"),
            "room_rh_percent": float("nan"),
            "rh_violated": False,
            "Q_sensible_zone_kW": float("nan"),
            "Q_latent_zone_kW": float("nan"),
            "Q_coil_zone_kW_signed": float("nan"),
            "Q_coil_zone_kW": 0.0,
            "Q_coil_district_kW": 0.0,
            "Q_evap_chiller_kW": 0.0,
            "m_condensate_zone_kg_s": 0.0,
            "m_condensate_district_kg_s": 0.0,
            "m_condensate_zone_kg_h": 0.0,
            "m_condensate_district_kg_h": 0.0,
            "m_condensate_district_m3_day": 0.0,
        }

    @classmethod
    def solve_condition(
        cls,
        T_db_out: float,
        rh_out: float,
        solar_on: bool,
        fresh_air_frac: float = 0.30,
    ) -> Dict[str, Any]:
        """Solve the zone humidity balance and cooling-coil condition."""

        numeric_inputs = {
            "T_db_out": T_db_out,
            "rh_out": rh_out,
            "fresh_air_frac": fresh_air_frac,
        }
        for name, value in numeric_inputs.items():
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
            ):
                raise ValueError(f"{name} must be a finite number.")

        if not isinstance(solar_on, bool):
            raise ValueError("solar_on must be True or False.")
        if not 0.0 <= fresh_air_frac <= 1.0:
            raise ValueError("fresh_air_frac must be between 0 and 1.")
        if not 0.0 <= rh_out <= 100.0:
            raise ValueError("rh_out must be between 0 and 100 percent.")

        result = cls._base_result(
            T_db_out, rh_out, solar_on, fresh_air_frac
        )

        P = cls.PRESSURE_KPA
        T_room = cls.T_ROOM_C
        m_da = cls.M_DA_ZONE_KG_S
        bf = cls.BF

        # Outdoor-air state
        w_out = MoistAir.humidity_ratio(T_db_out, rh_out, P)
        h_out = MoistAir.enthalpy(T_db_out, w_out)

        # Zone sensible and latent loads
        Q_sensible_zone_kW = cls._sensible_zone_load_kw(
            T_db_out, solar_on
        )
        Q_latent_zone_kW = (
            cls.OCCUPANTS * cls.OCCUPANT_LATENT_W / 1000.0
        )

        # Occupant moisture generation based on the specified latent load
        h_fg_room_kJ_kg = 2501.0 + 1.86 * T_room
        m_water_generation_kg_s = (
            Q_latent_zone_kW / h_fg_room_kJ_kg
        )

        # Supply temperature from the zone sensible heat balance.
        # The humidity-ratio approximation is retained from the original model.
        w_for_sensible_balance = 0.010
        cp_moist_air = 1.006 + 1.86 * w_for_sensible_balance
        T_supply = T_room - Q_sensible_zone_kW / (
            m_da * cp_moist_air
        )

        # Initial guess for room humidity
        w_room = MoistAir.humidity_ratio(T_room, 50.0, P)
        converged = False

        for iteration in range(1, cls.MAX_ITERATIONS + 1):
            h_room = MoistAir.enthalpy(T_room, w_room)

            # Adiabatic mixing on a dry-air mass basis
            w_mix = (
                (1.0 - fresh_air_frac) * w_room
                + fresh_air_frac * w_out
            )
            h_mix = (
                (1.0 - fresh_air_frac) * h_room
                + fresh_air_frac * h_out
            )
            T_mix = cls._temperature_from_enthalpy(h_mix, w_mix)

            # BF relation: T_supply = ADP + BF*(T_mix - ADP)
            ADP = (T_supply - bf * T_mix) / (1.0 - bf)

            if not math.isfinite(ADP) or ADP >= T_mix:
                result.update({
                    "status": (
                        "Infeasible coil state: ADP must be below "
                        "the mixed-air temperature."
                    ),
                    "iterations": iteration,
                    "T_mix": T_mix,
                    "w_out": w_out,
                    "h_out": h_out,
                    "w_room": w_room,
                    "h_room": h_room,
                    "w_mix": w_mix,
                    "h_mix": h_mix,
                    "T_supply": T_supply,
                    "ADP": ADP,
                    "Q_sensible_zone_kW": Q_sensible_zone_kW,
                    "Q_latent_zone_kW": Q_latent_zone_kW,
                })
                return result

            w_adp = MoistAir.humidity_ratio(ADP, 100.0, P)
            w_supply = w_adp + bf * (w_mix - w_adp)

            # Steady zone moisture balance
            w_room_new = (
                w_supply + m_water_generation_kg_s / m_da
            )

            if not math.isfinite(w_room_new) or w_room_new < 0.0:
                result.update({
                    "status": "Invalid room humidity ratio during iteration.",
                    "iterations": iteration,
                    "T_mix": T_mix,
                    "w_mix": w_mix,
                    "h_mix": h_mix,
                    "T_supply": T_supply,
                    "ADP": ADP,
                    "Q_sensible_zone_kW": Q_sensible_zone_kW,
                    "Q_latent_zone_kW": Q_latent_zone_kW,
                })
                return result

            if abs(w_room_new - w_room) < cls.HUMIDITY_TOLERANCE:
                w_room = w_room_new
                converged = True
                break

            # Damped update for numerical stability
            w_room = 0.5 * w_room + 0.5 * w_room_new

        if not converged:
            result.update({
                "status": "Room humidity iteration did not converge.",
                "iterations": cls.MAX_ITERATIONS,
                "Q_sensible_zone_kW": Q_sensible_zone_kW,
                "Q_latent_zone_kW": Q_latent_zone_kW,
            })
            return result

        # Recompute final states with the converged humidity ratio
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

        if not math.isfinite(ADP) or ADP >= T_mix:
            result.update({
                "status": "Infeasible final coil state.",
                "iterations": iteration,
                "T_mix": T_mix,
                "w_out": w_out,
                "h_out": h_out,
                "w_room": w_room,
                "h_room": h_room,
                "w_mix": w_mix,
                "h_mix": h_mix,
                "T_supply": T_supply,
                "ADP": ADP,
                "Q_sensible_zone_kW": Q_sensible_zone_kW,
                "Q_latent_zone_kW": Q_latent_zone_kW,
            })
            return result

        w_adp = MoistAir.humidity_ratio(ADP, 100.0, P)
        w_supply = w_adp + bf * (w_mix - w_adp)
        h_supply = MoistAir.enthalpy(T_supply, w_supply)

        # Condensate rate from the air-side moisture balance
        m_cond_zone_kg_s = max(
            m_da * (w_mix - w_supply), 0.0
        )

        # Liquid-water enthalpy relative to 0 °C
        h_condensate_kJ_kg = 4.186 * ADP

        # Signed coil duty including condensate enthalpy
        Q_coil_zone_kW_signed = (
            m_da * (h_mix - h_supply)
            - m_cond_zone_kg_s * h_condensate_kJ_kg
        )
        Q_coil_zone_kW = max(Q_coil_zone_kW_signed, 0.0)

        rh_room_actual = MoistAir.relative_humidity(
            T_room, w_room, P
        )
        rh_violated = (
            rh_room_actual > cls.RH_LIMIT_PERCENT + 1e-6
        )

        # District loads: 5% distribution gain and 2% pump heat
        Q_coil_district_kW = Q_coil_zone_kW * cls.N_ZONES
        Q_evap_chiller_kW = Q_coil_district_kW * 1.05 * 1.02

        # District condensate rate
        m_cond_district_kg_s = (
            m_cond_zone_kg_s * cls.N_ZONES
        )

        result.update({
            "feasible": True,
            "status": (
                "converged; room RH exceeds limit"
                if rh_violated
                else "converged"
            ),
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
            "room_rh_percent": rh_room_actual,
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
            "m_condensate_district_kg_h": (
                m_cond_district_kg_s * 3600.0
            ),
            "m_condensate_district_m3_day": (
                m_cond_district_kg_s * 86400.0
                / cls.WATER_DENSITY_KG_M3
                if hasattr(cls, "WATER_DENSITY_KG_M3")
                else m_cond_district_kg_s * 86400.0 / 1000.0
            ),
        })

        return result
