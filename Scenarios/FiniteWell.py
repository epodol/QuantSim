# Test Values:
# length: b_r * 4
# depth: 40
# mass: m_e

import numpy as np
from scipy.optimize import brentq
from scipy.constants import physical_constants
from sympy.parsing.sympy_parser import parse_expr
import matplotlib.pyplot as plt
from scipy.integrate import quad

# Constants
electron_mass = physical_constants['electron mass'][0]
proton_mass = physical_constants['proton mass'][0]
neutron_mass = physical_constants['neutron mass'][0]
alpha_particle_mass = physical_constants['alpha particle mass'][0]
hbar_eVs = physical_constants['reduced Planck constant in eV s'][0]
atomic_charge = physical_constants['atomic unit of charge'][0]
bohr_radius = physical_constants['Bohr radius'][0]
a0 = bohr_radius  # Bohr radius in meters
hbar_si = physical_constants['Planck constant over 2 pi'][0]  # hbar in SI units (Joule seconds)


# Transcendental equations for even and odd states
def trans_eq_even(x, depth, mass, length):
    return np.sqrt(x / (depth - x)) - 1 / (
        np.tan(np.sqrt(mass * x * length ** 2 / (2 * hbar_eVs ** 2 * atomic_charge))))


def trans_eq_odd(x, depth, mass, length):
    return np.sqrt(x / (depth - x)) + (
        np.tan(np.sqrt(mass * x * length ** 2 / (2 * hbar_eVs ** 2 * atomic_charge))))


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


def filter_energy_levels(energies, parity, length, depth, mass):
    filtered_energies = []
    tolerance = 5  # Adjust the tolerance to a reasonable value for your calculations
    for energy in energies:
        ko = np.sqrt(2 * mass * (depth - energy)) / np.sqrt(hbar_eVs ** 2 * atomic_charge)
        ki = np.sqrt(2 * mass * energy) / np.sqrt(hbar_eVs ** 2 * atomic_charge)
        D = np.sqrt(ko / ((0.5 * length * ko) + 1))
        # Check the continuity condition at the boundary for ψ and dψ/dx
        if parity == 'even':
            B = np.sqrt(ko / ((0.5 * length * ko) + 1)) * (np.cos((0.5 * length * ki))) * np.exp(0.5 * length * ko)
            # For even parity, match the derivative of the wave function
            right_derivative = -ko * B * np.exp(-ko * (length / 2))
            left_derivative = -ki * D * np.sin(ki * (length / 2))
        elif parity == 'odd':
            B = -1 * np.sqrt(ko / ((0.5 * length * ko) + 1)) * (np.sin((0.5 * length * ki))) * np.exp(0.5 * length * ko)
            # For odd parity, match the derivative of the wave function
            right_derivative = -ko * B * np.exp(-ko * (length / 2))
            left_derivative = ki * D * np.cos(ki * (length / 2))
        else:
            print("Invalid parity. Please try again.")
            quit()
        # Check if the derivatives and wave function values match within the specified tolerance

        if np.isclose(float(abs(right_derivative) / abs(left_derivative)), 1, tolerance):
            filtered_energies.append(energy)
    return filtered_energies


# New Function: Normalize Wave Function
def normalize_wave_function(psi, x_values):
    norm, _ = quad(lambda x: abs(psi(x)) ** 2, min(x_values), max(x_values))
    return lambda x: psi(x) / np.sqrt(norm)


# New Function: Expectation Value of Position
def expectation_value_position(psi, x_values):
    integrand = lambda x: np.conj(psi(x)) * x * psi(x)
    return quad(integrand, min(x_values), max(x_values))[0]


# New Function: Expectation Value of Momentum

def expectation_value_momentum(psi, x_values, hbar):
    dx = x_values[1] - x_values[0]  # Assuming uniform spacing

    def dpsi_dx(x):
        # Find the closest index to x in x_values
        index = np.searchsorted(x_values, x) - 1
        # Use central difference for interior points
        if 0 < index < len(x_values) - 1:
            return (psi(x_values[index + 1]) - psi(x_values[index - 1])) / (2 * dx)
        # Use forward/backward difference for endpoints
        elif index == 0:
            return (psi(x_values[index + 1]) - psi(x_values[index])) / dx
        else:
            return (psi(x_values[index]) - psi(x_values[index - 1])) / dx

    integrand = lambda x: np.conj(psi(x)) * (-1j * hbar * dpsi_dx(x))
    return quad(integrand, min(x_values), max(x_values))[0]


# New Function: RMS Momentum
def rms_momentum(psi, x_values, hbar):
    dx = x_values[1] - x_values[0]  # Assuming uniform spacing

    def d2psi_dx2(x):
        # Find the closest index to x in x_values
        index = np.searchsorted(x_values, x) - 1
        # Use central difference for interior points
        if 1 < index < len(x_values) - 2:
            return (psi(x_values[index + 1]) - 2 * psi(x_values[index]) + psi(x_values[index - 1])) / (dx ** 2)
        # Use forward/backward difference for endpoints and near-endpoints
        elif index == 0 or index == 1:
            return (psi(x_values[index + 2]) - 2 * psi(x_values[index + 1]) + psi(x_values[index])) / (dx ** 2)
        else:
            return (psi(x_values[index]) - 2 * psi(x_values[index - 1]) + psi(x_values[index - 2])) / (dx ** 2)

    return np.sqrt(quad(lambda x: np.conj(psi(x)) * (-hbar ** 2 * d2psi_dx2(x)), min(x_values), max(x_values))[0])


def calculate_wave_function(energy, length, parity, mass, depth):
    ko = np.sqrt(2 * mass * (depth - energy)) / np.sqrt(hbar_eVs ** 2 * atomic_charge)
    ki = np.sqrt(2 * mass * energy) / np.sqrt(hbar_eVs ** 2 * atomic_charge)
    B = (-1 if parity == 'odd' else 1) * np.sqrt(ko / ((0.5 * length * ko) + 1)) * (
        np.sin((0.5 * length * ki)) if parity == 'odd' else np.cos((0.5 * length * ki))) * np.exp(0.5 * length * ko)
    D = np.sqrt(ko / ((0.5 * length * ko) + 1))

    def psi(x):
        if np.isscalar(x):
            x_vals = np.array([x])
        else:
            x_vals = x
        result = np.zeros_like(x_vals)
        if parity == 'even':
            result[x_vals <= (-length / 2)] = B * np.exp(ko * x_vals[x_vals <= (-length / 2)])
            result[np.logical_and((-length / 2) < x_vals, x_vals < (length / 2))] = D * np.cos(
                ki * x_vals[np.logical_and((-length / 2) < x_vals, x_vals < (length / 2))])
            result[x_vals >= (length / 2)] = B * np.exp(-ko * x_vals[x_vals >= (length / 2)])
        elif parity == 'odd':
            result[x_vals <= (-length / 2)] = B * np.exp(ko * x_vals[x_vals <= (-length / 2)])
            result[np.logical_and((-length / 2) < x_vals, x_vals < (length / 2))] = D * np.sin(
                ki * x_vals[np.logical_and((-length / 2) < x_vals, x_vals < (length / 2))])
            result[x_vals >= (length / 2)] = -B * np.exp(-ko * x_vals[x_vals >= (length / 2)])
        return result[0] if np.isscalar(x) else result

    return psi


def probability_density_at_x(x, energy, length, mass, depth, parity, x_values):
    psi = calculate_wave_function(energy, parity, length, mass, depth)
    normalized_psi = normalize_wave_function(psi, x_values)
    probability_density = abs(normalized_psi(x)) ** 2
    return probability_density


# Function to plot a limited number of wave functions
def plot_wave_functions(energies, parity, x_values, length, mass, depth, max_plots=5):
    for i, energy in enumerate(energies[:max_plots]):
        psi_function = calculate_wave_function(energy, length, parity, mass, depth)
        psi_values = psi_function(x_values)
        label = f'{parity.capitalize()}, E = {energy:.3f} eV'
        plt.plot(x_values, psi_values, label=label, linestyle='-' if parity == 'even' else '--')


def finite_well():
    print("Finite Well Scenario")

    # Define the potential well parameters

    # Length of the well
    user_length = float(parse_expr(str(input("Enter the length of the well: ")), {'b_r': bohr_radius}))
    # add checks for L input to be a float
    if user_length <= 0:
        print("Length must be greater than 0. Please try again.")
        quit()

    # Depth of the well
    user_depth = float(input("Enter the well height (in eV): "))
    # add checks for V_0 input to be a float
    if user_depth <= 0:
        print("Depth must be greater than 0. Please try again.")
        quit()

    user_mass = float(parse_expr(str(input("Enter the mass of the trapped particle: ")),
                                 {'m_e': electron_mass, 'm_p': proton_mass, 'm_n': neutron_mass,
                                  'm_a': alpha_particle_mass}))
    # add checks for mass input to be greater than 0
    if user_mass <= 0:
        print("Mass must be greater than 0. Please try again.")
        quit()

    # Find roots (energy levels) for even and odd states
    even_energies = list(find_roots(lambda x: trans_eq_even(x, user_depth, user_mass, user_length), 0, user_depth))
    odd_energies = list(find_roots(lambda x: trans_eq_odd(x, user_depth, user_mass, user_length), 0, user_depth))

    # Filter out the energy levels
    even_energies_filtered = filter_energy_levels(even_energies, 'even', user_length, user_depth, user_mass)
    odd_energies_filtered = filter_energy_levels(odd_energies, 'odd', user_length, user_depth, user_mass)

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

    # Plotting the wave functions
    x_values = np.linspace(-3 * user_length, 3 * user_length, 10000)
    plt.figure(1, (12, 6))

    # Calculate and Print Expectation Values
    print("\nExpectation Values for position, momentum, and rms momentum at each energy level:")
    for energy in even_energies_filtered + odd_energies_filtered:
        parity = 'even' if energy in even_energies_filtered else 'odd'
        psi = calculate_wave_function(energy, user_length, parity, user_mass, user_depth)
        normalized_psi = normalize_wave_function(psi, x_values)
        x_expectation = expectation_value_position(normalized_psi,
                                                   x_values) / user_length  # Convert to multiples of well length
        p_expectation = expectation_value_momentum(normalized_psi, x_values, hbar_eVs) * (
                hbar_si / a0)  # Convert to kg*m/s
        rms_p = rms_momentum(normalized_psi, x_values, hbar_eVs) * (hbar_si / a0)  # Convert to kg*m/s
        print(f"\nEnergy: {energy:.3f} eV, Parity: {parity}")
        print(f"  <x>: {x_expectation:.2f}L, <p>: {p_expectation:.2f} kg*m/s, RMS p: {rms_p:.4e} kg*m/s")

    plot_wave_functions(even_energies_filtered, 'even', x_values, user_length, user_mass, user_depth)
    plot_wave_functions(odd_energies_filtered, 'odd', x_values, user_length, user_mass, user_depth)
    plt.title('Wave Functions by Energy Level')
    plt.xlabel('Position (x)')
    plt.yticks(color='w')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Ask the user for a specific point x
    x_point = float(input("\nEnter a value for x (in multiples of well length L) to find the probability density: "))
    x_point_actual = x_point * user_length  # Convert to actual length
    # Display the probability density for each energy level

    print("\nProbability Density at x = {:.2f} L:".format(x_point))

    # For counter with i=0 to i=number of energy levels
    for i in range(len(even_energies_filtered) + len(odd_energies_filtered)):
        if i < len(even_energies_filtered):
            energy = even_energies_filtered[i]
            parity = 'even'
        else:
            energy = odd_energies_filtered[i - len(even_energies_filtered)]
            parity = 'odd'
        prob_density = probability_density_at_x(x_point_actual, energy, user_length, user_mass, user_depth, parity,
                                                x_values)
        print(
            f"At E_{i}: {energy:.3f} eV, the probability of finding particle at {x_point:.2f} L is {prob_density:.3e}")
