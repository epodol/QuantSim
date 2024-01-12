import sympy as sp
import numpy as np
from sympy import diff
import math
import matplotlib.pyplot as plt
from scipy.integrate import quad


def quantum_oscillator():
    print(
        "Welcome. You will be entering values for w (measured in radians per second and n (an integer)). This will "
        "return the energy state for the quantum harmonic oscillator. \nInput your value for w in the 10^7 and 10^8 "
        "range.")
    h_bar = 1.054571817 * 10 ** (-34)
    x = sp.symbols('x')
    y = sp.symbols('y')
    n = int(input("Enter the value for n (integer): "))
    w = float(sp.sympify(input("Enter the angular frequency: ")))
    m = float(sp.sympify(input("Enter the mass: ")))
    x1 = float(sp.sympify(input("Enter the lower bound: ")))
    x2 = float(sp.sympify(input("Enter the upper bound: ")))
    energy_state = (1 / 2 + n) * h_bar * w
    print("Energy state is ", energy_state, " joules")
    a = m * w / h_bar
    hermite = ((-1) ** n) * (np.e ** y ** 2) * (diff(np.e ** (-y ** 2), y, n))
    print(hermite)
    hermite = hermite.subs(y, sp.sqrt(a) * x)
    print(hermite)
    wave_fncn = (((a / np.pi) ** 0.25) * (1 / np.sqrt((2 ** n) * math.factorial(n))) * hermite * (
                np.e ** ((-(x ** 2)) / 2)))
    print(wave_fncn)
    fncn = sp.lambdify(x, wave_fncn, 'numpy')
    x_vals = np.linspace(-10000, 10000, 10000)

    def integrand(x_values):
        return (wave_fncn.subs(x, x_values)) ** 2

    probability = quad(integrand, x1, x2)
    pdf_fncn = sp.lambdify(x, integrand(x), 'numpy')
    print("Probability is ", probability)
    plt.plot(x_vals, fncn(x_vals), label='g(x)')
    plt.plot(x_vals, pdf_fncn(x_vals), label='PDF')
    plt.legend()
    plt.show()
