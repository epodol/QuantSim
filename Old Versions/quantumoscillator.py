import sympy as sp
import numpy as np  
from sympy import diff 
import math
print("Welcome. You will be entering values for w (measured in radians per second and n (an integer)). This will return the energy state for the quantum harmonic oscillator. \nInput your value for w in the 10^7 and 10^8 range.")
h_bar = 1.054571817 * 10**(-34) 
x = sp.symbols('x')
n = int(input("Enter the value for n (integer): "))
w = float(sp.sympify(input("Enter the angular frequency: ")))
m = float(sp.sympify(input("Enter the mass: ")))
energy_state = (1/2 + n) * h_bar * w
print("Energy state is ", energy_state, " joules")
a = m*w/h_bar
y = np.sqrt(a) * x
hermite = (-1)**n*(np.e**x**2)*(diff(np.e**(-x**2), x, n))
wave_fncn = ((a/np.pi)**0.25)*(1/np.sqrt((2**n)*math.factorial(n))) * hermite * (np.e*(((-y)**2)/2))
print(wave_fncn)
