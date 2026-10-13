"""
Task Execution Handlers package.

Modules:
- task1_seasonal: Task 1 seasonal load analysis, chiller sizing, turndown/flue-gas checks.
- task2_condenser: Task 2 condenser type selection, design sizing, off-design performance mapping.
- task3_decanter: Task 3 seasonal condensate tracking and decanter tank sizing.
- task4_sensitivity: Task 4 fresh-air fraction sensitivity analysis (30% vs 40%).
"""

from .task1_seasonal import run_task1
from .task2_condenser import run_task2
from .task3_decanter import run_task3
from .task4_sensitivity import run_task4

__all__ = [
    "run_task1",
    "run_task2",
    "run_task3",
    "run_task4",
]