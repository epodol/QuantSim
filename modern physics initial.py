import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import sympy as sp

equation_str = input("Welcome. Please enter the equation in terms of x: ")
x = sp.symbols('x')
equation1 = sp.sympify(equation_str.split('=')[0])
equation2 = sp.sympify(equation_str.split('=')[1])
f = sp.lambdify(x, equation1, 'numpy')
g = sp.lambdify(x, equation2, 'numpy')

lower = float(input("Enter the lower limit of the domain: "))
upper = float(input("Enter the upper limit of the domain: "))

try:
    initial_guesses = [lower + (upper - lower) * i / 50 for i in range(50)]
    intersection_points = [fsolve(lambda x: f(x) - g(x), guess)[0] for guess in initial_guesses]

    valid_intersection_points = [point for point in intersection_points if lower <= point <= upper]

    print("Intersection Points:", set([str(round(point, 3)) for point in valid_intersection_points]))



    x_values = np.linspace(lower, upper)

    plt.plot(x_values, f(x_values), label='Equation 1')
    plt.plot(x_values, g(x_values), label='Equation 2')

    plt.scatter(valid_intersection_points, [f(point) for point in valid_intersection_points], color='red', zorder=5,
                label='Intersection Points')
    plt.legend()
    plt.show()

except:
    print("No intersection points found in the given domain. Make sure constants (pi, e, etc.) are not used in the "
          "equation, rather explicitly write out their values.")
