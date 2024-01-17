import sympy as sp
import numpy as np  
import matplotlib.pyplot as plt 


def infinite_well():
    width = float(sp.sympify((input("Enter the width of the well: "))))  # Width of the well
    energy_levels = int(input("Enter the energy level of the well: "))  # Number of energy levels to plot
    x_values = np.linspace(0, width, 10000)  # Values of x within the well

    # Define the wave function for an infinite well
    for n in range(1, energy_levels + 1):
        plt.plot(x_values, np.sqrt(2 / width) * np.sin(n * np.pi * x_values / width), label=f'n={n}')

    plt.xlabel('Position (x)')
    plt.ylabel('Wave Function')
    plt.title('Wave Functions for an Infinite Square Well')

    plt.legend()
    plt.show()
