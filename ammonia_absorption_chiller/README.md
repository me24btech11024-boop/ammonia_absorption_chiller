# Waste-Heat-Driven Ammonia-Water Absorption Chiller for District Cooling

A physics-based thermodynamic simulation, design, and validation framework for a waste-heat-driven single-effect $\text{NH}_3\text{-H}_2\text{O}$ absorption refrigeration unit paired with an outdoor Air Handling Unit (AHU) psychrometric model and a 400-zone district cooling network.

---

## Technical Overview

This project models the complete energy chain:
1. **Flue Gas Waste Heat Source**: Process flue gas ($300^\circ\text{C}$, up to $80\text{ kg/s}$) driving the generator ($\epsilon_{gen}=0.75$, stack floor $\ge 140^\circ\text{C}$).
2. **Ammonia-Water Absorption Chiller**: 
   - Real-gas Peng-Robinson / SRK EOS with Lee-Kesler/Pitzer vapor pressure for pure ammonia.
   - Pátek & Klomfar (1995) thermodynamic correlations for $\text{NH}_3\text{-H}_2\text{O}$ mixtures ($T_{bub}, T_{dew}, h_{liq}, h_{vap}$).
   - Components: Generator with rectifier ($99.5\%\text{ NH}_3$), Condenser, Evaporator ($T_{evap}=2^\circ\text{C}$), Absorber, Solution Pump ($\eta=0.70$), Solution Heat Exchanger ($\epsilon_{SHX}=0.70$), and Weak-Solution Throttle.
3. **District Cooling Network & Psychrometrics**:
   - 400 identical zones ($40\text{ m}^2 \times 3\text{ m}$) served by $7^\circ\text{C}$ supply / $12^\circ\text{C}$ return chilled water.
   - AHU adiabatic mixing (70% return / 30% fresh air standard; 60/40 sensitivity analysis).
   - Iterative room humidity ratio convergence to enforce $\le 75\%\text{ RH}$ ceiling.
4. **Condenser Selection & Off-Design Mapping**:
   - Air-Cooled ($U = 30\text{ W/m}^2\text{K}$) vs. Evaporative ($U = 600\text{ W/m}^2\text{K}$).
   - 6-season performance simulation & feasibility verification ($P_{high} \le 20\text{ bar}$, $T_{gen}$ limit, $T_{cond} \ge 25^\circ\text{C}$).

---

## Repository Structure

```text
ammonia_absorption_chiller/
├── requirements.txt           # Environment dependencies
├── README.md                  # Project documentation
├── main.py                    # Driver script executing Tasks 1-4
├── core/                      # Thermodynamic property solvers
│   ├── moist_air.py           # Ideal-gas moist air & wet-bulb solver
│   ├── pure_ammonia.py        # Real-gas EOS & departure functions
│   └── patek_klomfar.py       # Pátek & Klomfar (1995) correlations
├── models/                    # Subsystem physics models
│   ├── zone_psychrometrics.py # Zone psychrometrics & humidity iteration
│   ├── absorption_cycle.py    # Cycle mass & energy balance solver
│   ├── condenser.py           # Off-design condenser solver
│   └── decanter.py            # District condensate collection model
├── tasks/                     # Task execution handlers
│   ├── task1_seasonal.py      # Seasonal load & chiller sizing
│   ├── task2_condenser.py     # Condenser sizing & off-design mapping
│   ├── task3_decanter.py      # Water decanter sizing
│   └── task4_sensitivity.py   # 30% vs 40% fresh air sensitivity study
├── tests/                     # Validation suite
│   ├── test_validation.py     # Checks 1, 2, and 3 automated testing
│   ├── test_patek_klomfar.py  # Standalone Pátek-Klomfar unit tests
│   └── test_eos.py            # EOS departure function unit tests
└── outputs/                   # Exported tables and figures