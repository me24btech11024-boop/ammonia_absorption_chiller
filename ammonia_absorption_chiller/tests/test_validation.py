"""
Automated Validation Test Suite asserting exact numerical tolerances (Checks 1, 2, and 3).
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
    """Check 1: Adiabatic mixing tolerances (w, h <= 0.5%, T <= 0.1 K)."""
    w_room = MoistAir.humidity_ratio(22.0, 75.0)
    h_room = MoistAir.enthalpy(22.0, w_room)

    w_out = MoistAir.humidity_ratio(40.0, 50.0)
    h_out = MoistAir.enthalpy(40.0, w_out)

    w_mix_calc = 0.7 * w_room + 0.3 * w_out
    h_mix_calc = 0.7 * h_room + 0.3 * h_out
    
    w_mix_calc = 0.7 * w_room + 0.3 * w_out
    h_mix_calc = 0.7 * h_room + 0.3 * h_out

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

    assert diff_w_pct <= 0.5, f"w mismatch: {diff_w_pct:.3f}%"
    assert diff_h_pct <= 0.5, f"h mismatch: {diff_h_pct:.3f}%"
    assert diff_T_K <= 0.1, f"T mismatch: {diff_T_K:.3f} K"


def test_check_2_ammonia_real_gas_eos():
    """Check 2: Saturation properties vs Table B.2.1 within cubic EOS bounds."""
    P_sat_0_calc = PureAmmoniaEOS.p_sat_lee_kesler(273.15) * 100.0
    h_fg_0_calc = PureAmmoniaEOS.h_fg(273.15)

    P_sat_40_calc = PureAmmoniaEOS.p_sat_lee_kesler(313.15) * 100.0
    h_fg_40_calc = PureAmmoniaEOS.h_fg(313.15)

    assert abs(P_sat_0_calc - 429.6) / 429.6 * 100.0 <= 5.0
    assert abs(h_fg_0_calc - 1261.8) / 1261.8 * 100.0 <= 8.0
    assert abs(P_sat_40_calc - 1554.9) / 1554.9 * 100.0 <= 5.0
    assert abs(h_fg_40_calc - 1098.8) / 1098.8 * 100.0 <= 8.0


def test_check_3_absorption_cycle_balances():
    """Check 3: First law energy conservation and Carnot COP limit."""
    T_evap_K = 2.0 + 273.15
    T_cond_K = 45.0 + 273.15
    T_gen_K = 130.0 + 273.15
    T_abs_K = T_cond_K

    # Correct single-effect Carnot COP upper bound
    COP_rev = (T_evap_K / (T_cond_K - T_evap_K)) * ((T_gen_K - T_abs_K) / T_gen_K)

    cycle = AbsorptionCycle.solve_cycle(Q_evap_kW=5000.0, T_cond_C=45.0, T_evap_C=2.0)

    assert cycle["COP"] < COP_rev, f"COP ({cycle['COP']:.3f}) exceeds Carnot bound ({COP_rev:.3f})"
    assert 0.45 <= cycle["COP"] <= 0.65, f"COP ({cycle['COP']:.3f}) out of expected range [0.45, 0.65]"
    assert bool(cycle["check_p_high_pass"]) is True
        # Total solution mass balance
    assert abs(
        cycle["mass_balance_solution_residual_kg_s"]
    ) < 1e-8

    # Ammonia component mass balance
    assert abs(
        cycle["mass_balance_ammonia_residual_kg_s"]
    ) < 1e-8

    # Generator duty and COP must be positive
    assert cycle["Q_gen_kW"] > 0
    assert cycle["COP"] > 0

```python
    # Basic physical and numerical sanity checks
    assert cycle["Q_evap_kW"] == pytest.approx(5000.0)
    assert cycle["Q_gen_kW"] > 0
    assert cycle["Q_cond_kW"] > 0
    assert cycle["Q_abs_kW"] > 0
    assert cycle["COP"] > 0
    assert cycle["W_pump_kW"] >= 0

    # The reported heat rejection must satisfy the overall first-law balance.
    expected_rejection = (
        cycle["Q_gen_kW"]
        + cycle["Q_evap_kW"]
        + cycle["W_pump_kW"]
    )
    assert cycle["Q_reject_required_kW"] == pytest.approx(
        expected_rejection, rel=1e-8
    )
```
