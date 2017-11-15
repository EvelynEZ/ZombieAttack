import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

y0 = [0.5, 0.3] # the first entry is h and second entry is z
tau = np.linspace(0, 70, num=1000)

# the original params
a = 0.5 # human birth rate
b = 0.1 # human natural death rate 
delta = 0.001 # human death rate caused by zombie killing
alpha = 0.01 # human death rate caused by zombie infection
k = 100.0 # human carrying capacity
beta = 0.005 # zombie death rate caused by human killing
c = 0.5 # zombie natural death rate

# dimensionless variables
gamma = b / a # 0.2
mu = alpha * k / a # 2
nu = delta * k / a # 0.2
epsilon = beta * k / a # 1
omega = c / a # 1

params = [gamma, mu, nu, epsilon, omega]

def f(y, t0, args):
    h = y[0]
    z = y[1]
    gamma = args[0]
    mu = args[1]
    nu = args[2]
    epsilon = args[3]
    omega = args[4]
    return [h * (1 - h) - gamma * h - mu * h * z - nu * h * z, mu * h * z - epsilon * h * z - omega * z]


y = odeint(f, y0, tau, args=(params,))
line = plt.plot(tau, y)
plt.xlabel("tau")
plt.ylabel("[h, z]")
plt.show()
