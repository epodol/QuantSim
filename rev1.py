import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import sympy as sp

x = sp.symbols('x')
choice = int(input("Please choose a scenario number:\nFinite Well (E < U_o) [1]\nInfinite Well [2]\nHarmonic Oscillator [3]\n"))
print("Next, you will enter parameters. A couple of tips while choosing these parameters:\n"
      "1. Potential energy should be in the eV range (1 eV = 1.602*10^(-19) Joules).\n"
      "2. Mass should be on a quantum scale (similar to the mass of an electron).")

if choice not in [1, 2, 3]:
    print("Program terminated. Please enter a valid scenario input (1, 2, or 3).")
    quit()

mass = sp.symbols('mass')
dimensions = sp.symbols('dimensions')
domainL = sp.symbols('domainL')
domainU = sp.symbols('domainU')
potential_energy = sp.symbols('potential_energy')

mass_val = float(sp.sympify(input("Enter mass: ")))
dimensions_val = float(sp.sympify(input("Enter dimension: ")))
domainL_val = float(sp.sympify(input("Enter lower domain: ")))
domainU_val = float(sp.sympify(input("Enter upper domain: ")))
potential_energy_val = float(sp.sympify(input("Enter potential energy: ")))

# Define the equations
equation1 = sp.sqrt((potential_energy - x) / x) - sp.tan(
    sp.sqrt((2 * mass * (dimensions ** 2)) / (4 * ((6.941*(10**(-50)))))))


equation2 = -sp.sqrt((potential_energy - x) / x) - sp.cot(
    sp.sqrt((2 * mass * (dimensions ** 2)) / (4 * ((6.941* (10 ** (-50)))))))

# Substitute numerical values into the equations
equation1 = equation1.subs({mass: mass_val, dimensions: dimensions_val, potential_energy: potential_energy_val})
equation2 = equation2.subs({mass: mass_val, dimensions: dimensions_val, potential_energy: potential_energy_val})

# Convert symbolic equations to NumPy functions
equation1_np = sp.lambdify(x, equation1, 'numpy')
equation2_np = sp.lambdify(x, equation2, 'numpy')

# Use fsolve with NumPy functions
initial_guess = [float(domainL_val + (domainU_val - domainL_val) * i / 50) for i in range(50)]
solution1 = fsolve(equation1_np, initial_guess)
solution2 = fsolve(equation2_np, initial_guess)

# Filter solutions within the specified domain
valid_solutions_1 = [point for point in solution1 if domainL_val <= point <= domainU_val]
valid_solutions_2 = [point for point in solution2 if domainL_val <= point <= domainU_val]

print("Solutions 1:", set([round(point, 3) for point in valid_solutions_1]))
print("Solutions 2:", set([round(point, 3) for point in valid_solutions_2]))
