import sympy as sp
import numpy as np  
import matplotlib.pyplot as plt 

def infinite_well():
    width = float(sp.sympify((input("Enter the width of the well: "))))  # Width of the well
    energy_levels = int(input("Enter the energy level of the well: "))  # Number of energy levels to plot
    x_values = np.linspace(0, width, 10000)  # Values of x within the well

    # Create subplots
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(8, 6))

    # Plot wave functions
    for n in range(1, energy_levels + 1):
        wave_function = np.sqrt(2 / width) * np.sin(n * np.pi * x_values / width)
        axs[0].plot(x_values, wave_function, label=f'n={n}')

    axs[0].set_ylabel('Wave Function')
    axs[0].set_title('Wave Functions for an Infinite Square Well')
    axs[0].legend()

    # Plot probability density functions
    for n in range(1, energy_levels + 1):
        wave_function = np.sqrt(2 / width) * np.sin(n * np.pi * x_values / width)
        probability_density = np.abs(wave_function) ** 2
        axs[1].plot(x_values, probability_density, label=f'P(x) for n={n}', linestyle='dashed')

    axs[1].set_xlabel('Position (x)')
    axs[1].set_ylabel('Probability Density')
    axs[1].set_title('Probability Density for an Infinite Square Well')
    axs[1].legend()

    plt.show()

infinite_well()
