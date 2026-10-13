"""
Automated validation tests for the ammonia-water absorption chiller.

Checks:
1. Psychrometric mixing
2. Pure-ammonia property estimates
3. Absorption-cycle sanity and mass-balance diagnostics
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import pytest

from core.moist_air import MoistAir
from core.pure_ammonia import PureAmmoniaEOS
from models.absorption_cycle import AbsorptionCycle


def test_check_1_psychrometric_mixing():
    """Check 1: Mixing humidity ratio, enthalpy, and temperature."""

    w_room = MoistAir.humidity_ratio(22.0, 75.0)
    h_room = MoistAir.enthalpy(22.0, w_room)

    w_out = MoistAir.humidity_ratio(40.0, 50.0)
    h_out = MoistAir.enthalpy(40.0, w_out)

    w_mix_calc = 0.7 * w_room + 0.3 * w_out
    h_mix_calc = 0.7 * h_room + 0.3 * h_out

    # Derive mixed temperature from the moist-air enthalpy relation.
    T_mix_calc = (
        h_mix_calc - 2501.0 * w_mix_calc
    ) / (
        1.006 + 1.86 * w_mix_calc
    )

    w_mix_hand = 0.01575
    h_mix_hand = 67.84
    T_mix_hand = 27.48

    diff_w_pct = abs(w_mix_calc - w_mix_hand) / w_mix_hand * 100.0
    diff_h_pct = abs(h_mix_calc - h_mix_hand) / h_mix_hand * 100.0
    diff_T_K = abs(T_mix_calc - T_mix_hand)

    assert diff_w_pct <= 0.5, (
        f"Humidity-ratio deviation {diff_w_pct:.3f}% exceeds 0.5%"
    )
    assert diff_h_pct <= 0.5, (
        f"Enthalpy deviation {diff_h_pct:.3f}% exceeds 0.5%"
    )
    assert diff_T_K <= 0.1, (
        f"Mixed-temperature deviation {diff_T_K:.3f} K exceeds 0.1 K"
    )


def test_check_2_ammonia_real_gas_eos():
    """Check 2: Ammonia saturation pressure and latent heat."""

    P_sat_0_calc = PureAmmoniaEOS.p_sat_lee_kesler(273.15) * 100.0
    h_fg_0_calc = PureAmmoniaEOS.h_fg(273.15)

    P_sat_40_calc = PureAmmoniaEOS.p_sat_lee_kesler(313.15) * 100.0
    h_fg_40_calc = PureAmmoniaEOS.h_fg(313.15)

    P_sat_0_hand = 429.6
    h_fg_0_hand = 1261.8
    P_sat_40_hand = 1554.9
    h_fg_40_hand = 1098.8

    diff_P0_pct = abs(P_sat_0_calc - P_sat_0_hand) / P_sat_0_hand * 100.0
    diff_hfg0_pct = abs(h_fg_0_calc - h_fg_0_hand) / h_fg_0_hand * 100.0
    diff_P40_pct = abs(P_sat_40_calc - P_sat_40_hand) / P_sat_40_hand * 100.0
    diff_hfg40_pct = abs(h_fg_40_calc - h_fg_40_hand) / h_fg_40_hand * 100.0

    assert diff_P0_pct <= 5.0, (
        f"Ammonia saturation-pressure deviation at 0 C: {diff_P0_pct:.3f}%"
    )
    assert diff_hfg0_pct <= 8.0, (
        f"Ammonia latent-heat deviation at 0 C: {diff_hfg0_pct:.3f}%"
    )
    assert diff_P40_pct <= 5.0, (
        f"Ammonia saturation-pressure deviation at 40 C: {diff_P40_pct:.3f}%"
    )
    assert diff_hfg40_pct <= 8.0, (
        f"Ammonia latent-heat deviation at 40 C: {diff_hfg40_pct:.3f}%"
    )


def test_check_3_absorption_cycle_balances():
    """
    Check basic absorption-cycle output and reported mass-balance residuals.

    These assertions do not independently validate the thermodynamic model
    or prove that its assumed solution concentrations are physically valid.
    """

    T_evap_K = 275.15
    T_cond_K = 318.15
    T_gen_K = 403.15
    T_abs_K = T_cond_K

    COP_rev = (
        T_evap_K / (T_cond_K - T_evap_K)
    ) * (
        (T_gen_K - T_abs_K) / T_gen_K
    )

    cycle = AbsorptionCycle.solve_cycle(
        Q_evap_kW=5000.0,
        T_cond_C=45.0,
        T_evap_C=2.0,
    )

    assert cycle["COP"] > 0, "COP must be positive"
    assert cycle["COP"] < COP_rev, (
        f"COP {cycle['COP']:.3f} exceeds the stated reversible bound "
        f"{COP_rev:.3f}"
    )
    assert 0.45 <= cycle["COP"] <= 0.65, (
        f"COP {cycle['COP']:.3f} is outside the expected range [0.45, 0.65]"
    )

    assert cycle["check_p_high_pass"] is True
    assert cycle["check_flue_pass"] is True

    assert cycle["Q_gen_kW"] > 0
    assert cycle["W_pump_kW"] >= 0
    assert cycle["Q_evap_kW"] == pytest.approx(5000.0)

    assert abs(cycle["mass_balance_solution_residual_kg_s"]) < 1e-8
    assert abs(cycle["mass_balance_ammonia_residual_kg_s"]) < 1e-8
