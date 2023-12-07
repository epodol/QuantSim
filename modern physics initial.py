import numpy as np 
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

print("Enter your equation") 
test = input()
print("")
def fcn(x): 
    return (np.cos(x)-x)

print(fsolve(fcn, 0))
xpoints = np.arange(1,2, 8)
def RHS(y):
    return np.log(4*y)
# np.array([3, 10])

# plt.xlim() and plt.ylim() can be used for domain intervals? just a suggestion
# idx = np.argwhere(np.diff(np.sign(f - g))).flatten()  // returns points where two graphs (f and g) cross; print(idx) returns all the x-value points of intersection


plt.plot(np.linspace(0, 4*np.pi, 100),np.sin(np.linspace(0, 4*np.pi, 100)))
plt.show()

--------------------------------------------------------------------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import sympy
x = np.arange(0, 1000)
f = np.arange(0, 1000)

g = np.cos(np.arange(0, 10, 0.01) * 2) * 1000 #why does it return value error? 
lower = int(input("Enter a domain lower: "))
upper = int(input("Enter a domain upper: "))

plt.plot(x, f, '-')
plt.plot(x, g, '-')
plt.xlim(lower, upper)
idx = np.argwhere(np.diff(np.sign(f - g))).flatten()
plt.plot(x[idx], f[idx], 'ro')
plt.show()
idx = [i for i in idx if lower <= i <= upper]

print(idx)
# updates: added domain input [working] 
# goal: remove values higher than upper domain limit, return (print) updated array; recognize user input
# updates: removed values higher than upper domain limit (or lower than lower domain limit) 
------------------------------------------------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import sympy as sp

equation1_str = input("Welcome. To solve the transcendental equation, move the negative (if possible) to the other side (i.e x - sin(x)---> x = sin(x). Equation one will be the left hand side, equation two wil be the right hand side. Please enter the first equation in terms of x: ")
equation2_str = input("Enter the second equation in terms of x: ")

x = sp.symbols('x')
equation1 = sp.sympify(equation1_str)
equation2 = sp.sympify(equation2_str)

f = sp.lambdify(x, equation1, 'numpy')
g = sp.lambdify(x, equation2, 'numpy')

lower = float(input("Enter the lower limit of the domain: "))
upper = float(input("Enter the upper limit of the domain: "))

initial_guesses = [lower + (upper - lower) * i / 1000 for i in range(1000)]
intersection_points = [fsolve(lambda x: f(x) - g(x), guess)[0] for guess in initial_guesses]

valid_intersection_points = [point for point in intersection_points if lower <= point <= upper]
print("Intersection Points:", valid_intersection_points)

x_values = np.linspace(lower, upper, 1000)

plt.plot(x_values, f(x_values), label='Equation 1')
plt.plot(x_values, g(x_values), label='Equation 2')

plt.scatter(valid_intersection_points, [f(point) for point in valid_intersection_points], color='red', zorder=5, label='Intersection Points')
plt.legend()
plt.show()
#updates: recognizes user input, works for trig graphs 
# goals: TBD; array of intersection points returns BS values 
