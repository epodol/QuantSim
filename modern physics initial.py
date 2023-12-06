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
for nIdx in idx :
    if (upper < idx[nIdx]) :
     idx.remove(nIdx) 
print(idx)
# updates: added domain input [working] 
# goal: remove values higher than upper domain limit, return (print) updated array; recognize user input
