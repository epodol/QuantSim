import numpy as np
import matplotlib.pyplot as plt

# Define constants
L = float(input("Enter the width of the well: "))  # Width of the well
N = float(input("Enter the energy level of the well: "))  # Number of energy levels to plot
x_values = np.linspace(0, L, 1000)  # Values of x within the well

# Define the wave function for an infinite well
def psi_n(n, x):
    return np.sqrt(2 / L) * np.sin(n * np.pi * x / L)

# Plot the wave functions for the first N energy levels
for n in range(1, N + 1):
    plt.plot(x_values, psi_n(n, x_values), label=f'n={n}')

# Set plot labels and title
plt.xlabel('Position (x)')
plt.ylabel('Wave Function')
plt.title('Wave Functions for an Infinite Well')

# Show the plot with a legend
plt.legend()
plt.show()
