"""QuantSim provides three scenarios that the user
can choose from to simulate quantum mechanical behavior."""
import warnings

from Scenarios.FiniteWell import finite_well
from Scenarios.InfiniteWell import infinite_well
from Scenarios.QuantumOscillator import quantum_oscillator

# Suppress RuntimeWarnings for prod
warnings.filterwarnings('ignore', category=RuntimeWarning)


def main():
    """Main function that prompts the user to choose a scenario."""
    choice = int(input("Please choose a scenario number.\n[1] Finite Well (E < U_o)\n"
                       "[2] Infinite Well\n[3] Harmonic/Quantum Oscillator\nScenario #: "))

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
