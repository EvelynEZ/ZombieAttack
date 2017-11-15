import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import itertools

ep = 8e-4
y00 = np.linspace(0, 1, num=150)
y01 = np.linspace(0, 1, num=150)
allys = [] # record all ys such that stalemate happens

tau = np.linspace(0, 200, num=1000)

# the original params
a = 0.29 # human birth rate
b = 0.3 # human natural death rate 
delta = 0.01 # human death rate caused by zombie killing
alpha = 0.01 # human death rate caused by zombie infection
k = 900.0 # human carrying capacity
beta = 0.015 # zombie death rate caused by human killing
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

def allZero(y):
    return cmpz(y[0], 0, ep) and cmpz(y[1], 0, ep)


def cmpz(num1, num2, ep):
    return abs(num1 - num2) < ep

for y00, y01 in itertools.product(y00, y01):
    y0 = [y00, y01] # the first entry is h and second entry is z


    y = odeint(f, y0, tau, args=(params,))
    for i in range(len(y)):
        ysol = y[i]
        if (i != 0) and cmpz(ysol[0], ysol[1], ep) and cmpz(ysol[1], 0, ep):
            allys.append([y00, y01])
            print "stalemate is h={}, z={}".format(y00, y01)
            break
        elif cmpz(ysol[0], 0, ep) or cmpz(ysol[1], 0, ep):
            break

y = odeint(f, [allys[len(allys)/2][0], allys[len(allys)/2][1]], tau, args=(params,)) # choose one stalemate to plot
# find first point where h and z are both zero
stalemateIdx = [x for x in y if allZero(x)][0]

# find first tau when h and z are both zero
for i in range(len(y)):
    if y[i][0] == stalemateIdx[0] and y[i][1] == stalemateIdx[1]:
        print "stalemate tau is {}".format(tau[i])
plt.plot(tau, y)
plt.xlabel("tau")
plt.ylabel("[h, z]")
plt.show()
