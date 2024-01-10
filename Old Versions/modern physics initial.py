import numpy as np                  # python number library
import matplotlib.pyplot as plt     # graphing 
from scipy.optimize import fsolve   # solving the actual equation 
import sympy as sp                  # converting user input to python expressions 
equation_str = input("Welcome. Please enter the transcendental equation in terms of x. Make sure to write constants (e, pi) as their values (2.71, 3.14): ") # user input
x = sp.symbols('x') # recognizes x as a variable/symbol
equation1 = sp.sympify(equation_str.split('=')[0]) # splits equation into two expressions
equation2 = sp.sympify(equation_str.split('=')[1]) 
f = sp.lambdify(x, equation1, 'numpy') # translates sympy expressions into python functions, using numpy to convert to the specified number library
g = sp.lambdify(x, equation2, 'numpy')

lower = float(input("Enter the lower limit of the domain: ")) # domain input limits
upper = float(input("Enter the upper limit of the domain: "))


try: # execute the generation of code in a try condition to catch potential errors
    initial_guesses = [lower + (upper - lower) * i / 50 for i in range(50)] # initial approximations for the solution; increase if domain is high
    intersection_points = [fsolve(lambda x: f(x) - g(x), guess)[0] for guess in initial_guesses] # looks at where the sum changes signs, and recognizes that as an intersection point. the number of guesses may affect this
    valid_intersection_points = [point for point in intersection_points if lower <= point <= upper] # validates this based on domain limits
    print("Intersection Points:", set([str(round(point, 3)) for point in valid_intersection_points])) # prints the points, rounding to 3 decimal places (refer to graph to approxmimate intersection; exact values listed in array)
    x_values = np.linspace(lower, upper, 1000) # linear space, used for the graph
    plt.plot(x_values, f(x_values), label='Equation 1') # graphs equation 1
    plt.plot(x_values, g(x_values), label='Equation 2') # graphs equation 2 
    plt.scatter(valid_intersection_points, [f(point) for point in valid_intersection_points], color='red', zorder=5, label='Intersection Points') # highlights intersection points with red dots 
    plt.legend() # graph legend
    plt.show() # shows graph in a seperate window
except: # exceptions; if there's a problem with the user input (entering eulers number as is), returns this statement
    print("There was a problem finding intersection points found in the given domain, or graphing those intersection points. Make sure constants (pi, e, etc.) are not used in the "
          "equation, rather explicitly write out their values.")
