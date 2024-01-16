import sympy as sp
import numpy as np  
from sympy import diff 
import math
import matplotlib.pyplot as plt 
from scipy.integrate import quad 
print("Welcome. You will be entering values for w (measured in radians per second and n (an integer)). This will return the energy state, the hermite polynomial, the spacial function, graphs of the probability density function and g(x), and the probability for the quantum harmonic oscillator. \nInput your value for w in the 10^7 and 10^8 range.")
h_bar = 1.054571817 * 10**(-34) 
x = sp.symbols('x')
y = sp.symbols('y')
n = int(input("Enter the value for n (integer): "))
w = float(sp.sympify(input("Enter the angular frequency: ")))
m = float(sp.sympify(input("Enter the mass: ")))
x1 = float(sp.sympify(input("Enter the lower bound: ")))
x2 = float(sp.sympify(input("Enter the upper bound: ")))
energy_state = (1/2 + n) * h_bar * w
print("Energy state is: ", energy_state, "Joules")
a = m*w/h_bar
hermite = ((-1)**n)*(np.e**y**2)*(diff(np.e**(-y**2), y, n))
print("Hermite before substitution: " , hermite)  
hermite = hermite.subs(y, (sp.sqrt(a) * x))
print("Hermite after substitution: ", hermite)
wave_fncn = (((a/np.pi)**0.25)*(1/(np.sqrt((2**n)*math.factorial(n)))) * hermite * (np.e**((-(y**2))/2)))
wave_fncn = wave_fncn.subs(y, sp.sqrt(a) * x)
print("Spacial function g(x): ", wave_fncn)
fncn = sp.lambdify(x, wave_fncn, 'scipy')
x_vals = np.linspace(x1, x2, 1000000) 
def integrand(x_values) :
        return ((wave_fncn.subs(x, x_values))**2)
print("Probability Density Function (PDF): ", integrand(x))
probability, err= quad(integrand, x1, x2)
pdf_fncn = sp.lambdify(x, integrand(x), 'scipy') 
print("Probability: ", probability)
print("Error: ", err)
def expected_position(x_values) : 
        return ((integrand(x) * x).subs(x, x_values))
expected_position_integrated, error = quad(expected_position, x1, x2)
print("Expected position <x>: ", expected_position_integrated)
print("Error: ", error)
plt.plot(x_vals, pdf_fncn(x_vals), label = 'PDF')
plt.legend()
plt.show()
plt.plot(x_vals, fncn(x_vals), label = 'g(x)')
plt.legend() 
plt.show()
