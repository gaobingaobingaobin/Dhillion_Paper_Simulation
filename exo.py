__author__ = 'qsong'


import matplotlib.pyplot as plt
import math
from scipy import integrate


line1, = plt.plot([1,2,3], label="Line 1", linestyle='--')
line2, = plt.plot([3,2,1], label="Line 2", linewidth=4)

# Create a legend for the first line.
first_legend = plt.legend(handles=[line1], loc=1)

# Add the legend manually to the current Axes.
ax = plt.gca().add_artist(first_legend)

# Create another legend for the second line.
plt.legend(handles=[line2], loc=4)


def f(x, mu, delta):
    return 1.0/(1+x)*1.0/(math.sqrt(2*math.pi*delta**2))*math.exp(-(x-mu)**2/(2*delta**2))



I = integrate.quad(f, -0.25, 0.25, args=(0, 1.6))

print I


def g(x, a, b):
    return a*x+b


I = integrate.quad(g, 0, 1, args=(2,1))

print I

def h(x, mu, delta):
    return 1.0/(math.sqrt(2*math.pi*delta**2))*math.exp(-(x-mu)**2/(2*delta**2))*(1.0/(1+x))

I = integrate.quad(h, -0.75, 0.75, args=(0, 0.5))
print I



