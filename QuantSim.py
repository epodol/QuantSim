import sympy as sp
import warnings

from Scenarios.FiniteWell import finite_well
from Scenarios.QuantumOscillator import quantum_oscillator

# Suppress RuntimeWarnings for prod
warnings.filterwarnings('ignore', category=RuntimeWarning)

x = sp.symbols('x')
choice = int(input("Please choose a scenario number:\nFinite Well (E < U_o) [1]\nInfinite Well [2]\nHarmonic "
                   "Oscillator [3] Quantum Oscillator\n"))
print("Next, you will enter parameters. A couple of tips while choosing these parameters:\n"
      "1. Potential energy should be in the eV range (1 eV = 1.602*10^(-19) Joules).\n"
      "2. Mass should be on a quantum scale (similar to the mass of an electron).")

match choice:
    case 1:
        finite_well()
    case 2:
        print("Infinite Well Scenario")
        # TODO: Implement Infinite Well
    case 3:
        print("Quantum Oscillator Scenario")
        quantum_oscillator()
    case _:
        print("Please enter a valid scenario input (1, 2, or 3).")
        quit()
