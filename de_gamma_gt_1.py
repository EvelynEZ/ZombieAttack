import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

y0 = [0.5, 0.3] # the first entry is h and second entry is z
tau = np.linspace(0, 20, num=1000)

# the original params
a = 0.29 # human birth rate
b = 0.3 # human natural death rate
delta = 0.002 # human death rate caused by zombie killing
alpha = 0.002 # human death rate caused by zombie infection
k = 900.0 # human carrying capacity
beta = 0.01 # zombie death rate caused by human killing
c = 0.01 # zombie natural death rate

# dimensionless variables
gamma = b / a
mu = alpha * k / a
nu = delta * k / a
epsilon = beta * k / a
omega = c / a

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


def f2(y, t0, args):
    H = y[0]
    Z = y[1]
    return [a * H * (1 - H/ k) - b * H - alpha * H * Z - delta * H * Z, alpha * H * Z - beta * H * Z - c * Z]

y = odeint(f, y0, tau, args=(params,))
line = plt.plot(tau, y)
plt.xlabel("tau")
plt.ylabel("[h, z]")
plt.show()
