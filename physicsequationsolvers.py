import numpy as np
from scipy.optimize import brentq
from scipy.constants import physical_constants
from sympy.parsing.sympy_parser import parse_expr
import matplotlib.pyplot as plt
import warnings

# ex: length: b_r*2, eV: 40

# Constants
m_e = physical_constants['electron mass'][0]
hbar_eVs = physical_constants['reduced Planck constant in eV s'][0]
ch = physical_constants['atomic unit of charge'][0]
b_r = physical_constants['Bohr radius'][0]


# Suppress RuntimeWarnings for prod
# warnings.filterwarnings('ignore', category=RuntimeWarning)

# Transcendental equations for even and odd states
def trans_eq_even(x):
    return np.sqrt(x / (V_0 - x)) - 1 / (np.tan(np.sqrt(m_e * x * L ** 2 / (2 * hbar_eVs ** 2 * ch))))


def trans_eq_odd(x):
    return np.sqrt(x / (V_0 - x)) + (np.tan(np.sqrt(m_e * x * L ** 2 / (2 * hbar_eVs ** 2 * ch))))


# Finding roots of the transcendental equations, using Brent's method
def find_roots(func, min_val, max_val):
    if func(min_val) == 0:
        min_val = min_val + 1e-10
    # define array of floats empty, explicitly a float
    roots = []

    for val in np.linspace(min_val, max_val, 10000):
        if func(val) * func(val + (max_val - min_val) / 10000) <= 0:
            try:
                roots.append(brentq(func, val, val + (max_val - min_val) / 10000))
            except Exception as e:
                print(f"Error finding root: {e}")
    # Filter out near-duplicate roots
    return np.unique(np.around(roots, 6)).tolist()


def filter_energy_levels(energies, parity, length, depth):
    filtered_energies = []
    tolerance = 5  # Adjust the tolerance to a reasonable value for your calculations
    for energy in energies:
        ko = np.sqrt(2 * m_e * (depth - energy)) / np.sqrt(hbar_eVs ** 2 * ch)
        ki = np.sqrt(2 * m_e * energy) / np.sqrt(hbar_eVs ** 2 * ch)
        D = np.sqrt(ko / ((0.5 * length * ko) + 1))
        # Check the continuity condition at the boundary for ψ and dψ/dx
        if parity == 'even':
            B = np.sqrt(ko / ((0.5 * length * ko) + 1)) * (np.cos((0.5 * length * ki))) * np.exp(0.5 * length * ko)
            # For even parity, match the derivative of the wavefunction
            right_derivative = -ko * B * np.exp(-ko * (length / 2))
            left_derivative = -ki * D * np.sin(ki * (length / 2))
        elif parity == 'odd':
            B = -1 * np.sqrt(ko / ((0.5 * length * ko) + 1)) * (np.sin((0.5 * length * ki))) * np.exp(0.5 * length * ko)
            # For odd parity, match the derivative of the wavefunction
            right_derivative = -ko * B * np.exp(-ko * (length / 2))
            left_derivative = ki * D * np.cos(ki * (length / 2))
        else:
            print("Invalid parity. Please try again.")
            quit()
        # Check if the derivatives and wave function values match within the specified tolerance

        if np.isclose(float(abs(right_derivative) / abs(left_derivative)), 1, tolerance):
            filtered_energies.append(energy)
    return filtered_energies


# Define the potential well parameters

# Length of the well
L = float(parse_expr(str(input("Enter the length of the well: ")), local_dict={'b_r': b_r}))
# add checks for L input to be a float
if L <= 0:
    print("Length must be greater than 0. Please try again.")
    quit()

# Depth of the well
V_0 = float(input("Enter the well height (in eV): "))
# add checks for V_0 input to be a float
if V_0 <= 0:
    print("Depth must be greater than 0. Please try again.")
    quit()

# Find roots (energy levels) for even and odd states
even_energies = list(find_roots(trans_eq_even, 0, V_0))
odd_energies = list(find_roots(trans_eq_odd, 0, V_0))

# Filter out the energy levels
even_energies_filtered = filter_energy_levels(even_energies, 'even', L, V_0)
odd_energies_filtered = filter_energy_levels(odd_energies, 'odd', L, V_0)

# Print energy levels

if len(even_energies_filtered) == 0 and len(odd_energies_filtered) == 0:
    print("No energy levels found. Please adjust input parameters and try again.")
    quit()

if len(even_energies_filtered) == 0:
    print("No even energy levels found.")
else:
    print("Even energy levels (eV):", even_energies_filtered)

if len(odd_energies_filtered) == 0:
    print("No odd energy levels found.")
else:
    print("Odd energy levels (eV):", odd_energies_filtered)


# Calculate the wave functions
def calculate_wave_function(x, E, L, parity):
    ko = np.sqrt(2 * m_e * (V_0 - E)) / np.sqrt(hbar_eVs ** 2 * ch)
    ki = np.sqrt(2 * m_e * E) / np.sqrt(hbar_eVs ** 2 * ch)
    B = (-1 if parity == 'odd' else 1) * np.sqrt(ko / ((0.5 * L * ko) + 1)) * (
        np.sin((0.5 * L * ki)) if parity == 'odd' else np.cos((0.5 * L * ki))) * np.exp(0.5 * L * ko)
    D = np.sqrt(ko / ((0.5 * L * ko) + 1))
    result = np.zeros_like(x)
    if parity == 'even':
        result[x <= (-L / 2)] = B * np.exp(ko * x[x <= (-L / 2)])
        result[np.logical_and((-L / 2) < x, x < (L / 2))] = D * np.cos(
            ki * x[np.logical_and((-L / 2) < x, x < (L / 2))])
        result[x >= (L / 2)] = B * np.exp(-ko * x[x >= (L / 2)])
    elif parity == 'odd':
        result[x <= (-L / 2)] = B * np.exp(ko * x[x <= (-L / 2)])
        result[np.logical_and((-L / 2) < x, x < (L / 2))] = D * np.sin(
            ki * x[np.logical_and((-L / 2) < x, x < (L / 2))])
        result[x >= (L / 2)] = -B * np.exp(-ko * x[x >= (L / 2)])
    return result


# Plotting the wavefunctions
x_values = np.linspace(-3 * L, 3 * L, 10000)
plt.figure(figsize=(12, 6))


# Function to plot a limited number of wavefunctions
def plot_wavefunctions(energies, parity, max_plots=5):
    for i, E in enumerate(energies[:max_plots]):
        psi = calculate_wave_function(x_values, E, L, parity)
        label = f'{parity.capitalize()}, E = {E:.3f} eV'
        plt.plot(x_values, psi, label=label, linestyle='-' if parity == 'even' else '--')


plot_wavefunctions(even_energies_filtered, 'even')
plot_wavefunctions(odd_energies_filtered, 'odd')
plt.title('Wavefunctions by Energy Level')
plt.xlabel('Position (x)')
plt.yticks(color='w')
plt.legend()
plt.grid(True)
plt.show()
