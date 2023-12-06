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
