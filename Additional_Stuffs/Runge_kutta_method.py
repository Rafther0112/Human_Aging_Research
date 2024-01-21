#%%
import numpy as np
import matplotlib.pyplot as plt
import math 

def f(x,t):
    return -x**3 + math.sin(t)

a = 0.0
b = 10.0
N = [10,20,50,100]

x_points = []
temporal_points = []

for n in N:
    x = 0.0
    h = (b-a)/n
    temporal_individuals_points = np.arange(a,b,h)
    x_individual_points = []
    for t in temporal_individuals_points:
        x_individual_points.append(x)
        k1 = h*f(x,t)
        k2 = h*f(x+0.5*k1, t +0.5*h)
        x += k2
    x_points.append(x_individual_points)
    temporal_points.append(temporal_individuals_points)
#%%
plt.figure()
for i in range(len(temporal_points)):
    plt.plot(temporal_points[i], x_points[i], label = f"N: {N[i]}")
plt.xlabel("t")
plt.ylabel("X(t)")
plt.legend()
plt.show()
