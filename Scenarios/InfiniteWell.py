import sympy as sp
import numpy as np  
import matplotlib.pyplot as plt 
# Define constants
L = float(sp.sympify((input("Enter the width of the well: "))))  # Width of the well
N = int(input("Enter the energy level of the well: "))  # Number of energy levels to plot
x_values = np.linspace(0, L, 10000)  # Values of x within the well

# Define the wave function for an infinite well
def psi_n(n, x):
    return np.sqrt(2 / L) * np.sin(n * np.pi * x / L)

for n in range(1, N + 1):
    plt.plot(x_values, psi_n(n, x_values), label=f'n={n}')

plt.xlabel('Position (x)')
plt.ylabel('Wave Function')
plt.title('Wave Functions for an Infinite Square Well')

plt.legend()
plt.show()
