import sympy as sp
import numpy as np
from sympy import diff
import math
import matplotlib.pyplot as plt
from scipy.integrate import quad
import warnings
from scipy.integrate import IntegrationWarning

warnings.filterwarnings('ignore', category=IntegrationWarning)

def quantum_oscillator():
    # complex number i is "j" in py, use this for momentum ?

    print(
        "Welcome. You will be entering the sigma and epsilon values, mass (in AMU). "
        "This will return the energy state, the hermite polynomial, the spacial function, graphs of the probability density function and g(x), the probability for the quantum harmonic oscillator, the expected position, the expected momentum, and momentum RMS."
    )
    h_bar = 1.054571817 * 10**(-34)  # in joule-seconds
    x = sp.symbols('x')
    y = sp.symbols('y')
    sigma_value = float(sp.sympify(input("Enter the sigma value: ")))
    epsilon_value = float(sp.sympify(input("Enter the epsilon value: ")))
    m = float(sp.sympify(input("Enter the mass in AMU: ")))
    x1 = float(sp.sympify(input("Enter the lower bound: ")))
    x2 = float(sp.sympify(input("Enter the upper bound: ")))

    m = m * 1.66054*10**(-27)
    mu = m/2
    k = 2*epsilon_value*((6/sigma_value)**2)  # n = 6
    w = np.sqrt(k/mu)
    n = (epsilon_value/(h_bar*w))-0.5
    n = round(n)

    print("Energy state integer: ", n)
    energy_state = (1/2 + n) * h_bar * w
    print("Energy state is: ", energy_state, "Joules")
    a = mu*w/h_bar

    hermite = ((-1)**n)*(np.e**y**2)*(diff(np.e**(-y**2), y, n))
    print("Hermite Polynomial for n =", n ,": ", hermite)
    hermite = hermite.subs(y, (sp.sqrt(a) * x))
    print("Hermite after substitution: ", hermite)

    wave_fncn = (((a/np.pi)**0.25)*(1/(np.sqrt((2**n)*math.factorial(n)))) * hermite * (np.e**((-(y**2))/2)))
    wave_fncn = wave_fncn.subs(y, sp.sqrt(a) * x)
    print("Spacial Function g(x): ", wave_fncn)
    fncn = sp.lambdify(x, wave_fncn, 'scipy')
    x_vals = np.linspace(x1, x2, 1000000)

    def integrand(x_values):
        return (wave_fncn.subs(x, x_values)) ** 2

    print("Probability Density Function (PDF): ", integrand(x))
    probability, err = quad(integrand, x1, x2)
    pdf_fncn = sp.lambdify(x, integrand(x), 'scipy')
    print("Probability: ", probability)
    print("Error: ", err)

    def expected_position(x_values):
        return (integrand(x) * x).subs(x, x_values)

    expected_position_integrated, error = quad(expected_position, x1, x2)
    print("Expected position <x>: ", expected_position_integrated)
    print("Error: ", error)
     def expected_momentum(x_values) : 
        wave_fncn_conjugate = (((a/np.pi)**0.25)*(1/(np.sqrt((2**n)*math.factorial(n)))) * hermite * (np.e**(((-y**2))/2)))
        wave_fncn_conjugate = wave_fncn_conjugate.subs(y, sp.sqrt(a)*x)
        momentum_integrand = wave_fncn_conjugate*wave_fncn
        integrand_derivated = diff(momentum_integrand, x)
        return (integrand_derivated.subs(x, x_values))
expected_momentum_integrated, error_m = quad(expected_momentum, x1, x2) 
expected_momentum_integrated = (expected_momentum_integrated * h_bar*cmath.sqrt(-1)*(-1)).imag

print("Expected momentum (in terms of i) <p>: ", expected_momentum_integrated)
print("Error: ", error_m)

def momentum_rms(x_values) : 
        wave_fncn_conjugated = (((a/np.pi)**0.25)*(1/(np.sqrt((2**n)*math.factorial(n)))) * hermite * (np.e**(((-y**2))/2))).subs(y, sp.sqrt(a)*x)
        integrand_rms = diff(wave_fncn, x, 2)
        momentum_rms_integrand = (h_bar**2) * wave_fncn_conjugated * integrand_rms
        return ((momentum_rms_integrand.subs(x, x_values)))
momentum_rms_integrated, err_rms = quad(momentum_rms, x1, x2)
momentum_rooted = sp.sqrt(abs(momentum_rms_integrated))
print("Momentum RMS <p_rms>: ", momentum_rooted)
print("Error: ", err_rms)

    plt.plot(x_vals, pdf_fncn(x_vals), label='PDF')
    plt.legend()
    plt.show()
    plt.plot(x_vals, fncn(x_vals), label='g(x)')
    plt.legend()
    plt.show()

# Call the quantum_oscillator function
quantum_oscillator()
