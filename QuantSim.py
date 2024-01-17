import sympy as sp
import warnings

from Scenarios.FiniteWell import finite_well
from Scenarios.InfiniteWell import infinite_well
from Scenarios.QuantumOscillator import quantum_oscillator

# Suppress RuntimeWarnings for prod
warnings.filterwarnings('ignore', category=RuntimeWarning)


def main():
    choice = int(input("Please choose a scenario number:\nFinite Well (E < U_o) [1]\nInfinite Well [2]\nHarmonic "
                       "Oscillator [3] Quantum Oscillator\n"))

    match choice:
        case 1:
            finite_well()
        case 2:
            print("Infinite Well Scenario")
            infinite_well()
        case 3:
            print("Quantum Oscillator Scenario")
            quantum_oscillator()
        case _:
            print("Please enter a valid scenario input (1, 2, or 3).")
            main()


main()
