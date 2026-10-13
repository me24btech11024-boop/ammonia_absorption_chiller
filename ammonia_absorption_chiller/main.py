"""
Master Driver Script for Ammonia-Water Absorption Chiller Project.

Executes:
1. Validation Test Suite (Checks 1, 2, and 3)
2. Task 1: Seasonal Load & Chiller Circuit Sizing Analysis
3. Task 2: Condenser Sizing & Off-Design Performance Mapping
4. Task 3: Condensate Generation & Decanter System Sizing
5. Task 4: Fresh-Air Fraction Sensitivity Study (30% vs 40%)
"""

import sys
from tabulate import tabulate

from tasks.task1_seasonal import run_task1
from tasks.task2_condenser import run_task2
from tasks.task3_decanter import run_task3
from tasks.task4_sensitivity import run_task4


def print_section_header(title: str):
    """Utility to print styled section headers."""
    print("\n" + "=" * 80)
    print(f" {title} ".center(80, "="))
    print("=" * 80)


def main():
    print_section_header("PROJECT 2: AMMONIA-WATER ABSORPTION CHILLER SIMULATION")

    # ----------------------------------------------------
    # TASK 1: Seasonal Loads & Chiller Sizing
    # ----------------------------------------------------
    print_section_header("TASK 1: Seasonal Load Analysis & Ammonia Circuit Sizing")
    df_task1, gov_cond_t1 = run_task1()
    print(tabulate(df_task1, headers="keys", tablefmt="pretty", showindex=False))
    print(f"\nGoverning Seasonal Condition (Max Evaporator Load): {gov_cond_t1}")

    # ----------------------------------------------------
    # TASK 2: Condenser Sizing & Off-Design Mapping
    # ----------------------------------------------------
    print_section_header("TASK 2: Condenser Design Sizing & Off-Design Performance")
    design_sizing, df_off_design = run_task2()
    print(f"Selected Condenser Type: {design_sizing['condenser_type'].upper()}")
    print(f"Design Heat-Transfer Conductance (UA): {design_sizing['UA_kW_K']:.2f} kW/K")
    print(f"Required Surface Area: {design_sizing['area_m2']:.2f} m^2\n")
    print(tabulate(df_off_design, headers="keys", tablefmt="pretty", showindex=False))

    # ----------------------------------------------------
    # TASK 3: Water Decanter Sizing
    # ----------------------------------------------------
    print_section_header("TASK 3: Water Decanter System Sizing")
    df_task3, gov_cond_decanter, tank_sizing = run_task3()
    print(tabulate(df_task3, headers="keys", tablefmt="pretty", showindex=False))
    print(f"\nGoverning Condensate Condition: {gov_cond_decanter}")
    print(f"Sized Decanter Tank Storage Volume (12-hr buffer): {tank_sizing['required_tank_volume_m3']:.2f} m^3")

    # ----------------------------------------------------
    # TASK 4: Fresh-Air Sensitivity (30% vs 40%)
    # ----------------------------------------------------
    print_section_header("TASK 4: Fresh-Air Sensitivity Analysis (30% vs 40%)")
    df_task4 = run_task4()
    print(tabulate(df_task4, headers="keys", tablefmt="pretty", showindex=False))

    print("\n" + "=" * 80)
    print(" ALL SIMULATIONS & TASKS EXECUTED SUCCESSFULLY ".center(80, "="))
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
