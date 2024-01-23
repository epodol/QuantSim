import sympy as sp
import numpy as np  
from sympy import diff 
import cmath
import matplotlib.pyplot as plt 
from scipy.integrate import quad

def infinite_well():
h_bar =  1.054571817*10**(-34)
m = float(sp.sympify(input("Enter mass in AMU: ")))
m = m * 1.66054*10**(-27)
w = float(sp.sympify(input("Enter width (w) of well: ")))
print("For your reference later on, -w/2 is ", (-w/2), "and w/2 is", (w/2))
n = int(input("Enter energy level: "))
x = sp.symbols('x')
k = n * np.pi * 2 / w
wave_fncn = sp.sqrt(2/w) * sp.sin((n*np.pi/w)*x) 
pdf = (2/w)*((sp.sin((n*np.pi/w)*x))**2)
print("Spacial function g(x): ", wave_fncn)
print("Probability Density Function (PDF): ", pdf)
wave_fncn_g = sp.lambdify(x, wave_fncn, 'scipy')
pdf_g = sp.lambdify(x, pdf, 'scipy')
# find probability using a rectangle, user inputs x (between -L/2 and L/2), a small dx
x_val = (float(sp.sympify(input("Now, you will select a value x between (-w/2 and w/2) to find the probability of a particle being at that point: " ))))
pdf = lambda x: (2/w)*((sp.sin((n*np.pi/w)*x))**2)
probability, err = quad(pdf, x_val, (x_val + (w*0.001)))
print("Probability", probability)
print("Error: ", err)
expected_position_integrand = lambda x: (2/w)*((sp.sin((n*np.pi/w)*x))**2) * x
expected_position, err_x = quad(expected_position_integrand, 0, w)
print("Expeected position <x>: ", expected_position)
print("Error: ", err_x)
expected_momentum_integrand = lambda x: ((sp.sqrt(2/w) * sp.sin((n*np.pi/w)*x) ) * ((n * np.pi * x / w) * sp.sqrt(2/w) * sp.cos(n * np.pi * x / w)))
expected_momentum, err_p = quad(expected_momentum_integrand, x_val, (x_val + (w*0.001))) # bounds?
expected_momentum = (expected_momentum * h_bar * cmath.sqrt(-1)).imag
print("Expected momentum <p>: ", expected_momentum)
print("Error: ", err_p)
x_vals = np.linspace(-w/2, w/2, 1000000) 
plt.plot(x_vals, pdf_g(x_vals), label = 'PDF')
plt.legend()
plt.show()
plt.plot(x_vals, wave_fncn_g(x_vals), label = 'Wave Function')
plt.legend()
plt.show()
infinite_well()
