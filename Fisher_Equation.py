import numpy as np
import scipy as sp
from scipy import sparse
from scipy.sparse import linalg
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import animation

x_0 = 0
x_N = 20
steps = 40
h = (x_N - x_0)/steps

identity = np.identity(steps)
e = np.ones(steps)
diagonals = [e, -2*e, e]
offset = [-1, 0, 1]
L_m =  sp.sparse.spdiags(diagonals, offset, steps, steps)
k = 0.2*h**2
A_m = identity + (k/h**2)*(L_m)
time_steps = int(42/k)

time_sol = np.zeros((time_steps,steps))

initial_value = np.zeros(steps)
initial_value[0:5] = 1

for i in range(time_steps):
    if i == 0:
        B = initial_value
        time_sol[i] =  A_m.dot(initial_value) +  k*(initial_value*(1-initial_value))
    else:
        B = time_sol[0]
        time_sol[i] =  A_m.dot(time_sol[i-1]) + k*(time_sol[i-1]*(1- time_sol[i]))
        time_sol[i,0] = 1
        time_sol[i,-1] = 0

x_ax = np.linspace(0,20,steps)
plt.plot(x_ax,time_sol[200])
plt.show()

fig = plt.figure()
ims = []
for i in range(time_steps):
    im = plt.plot(x_ax, time_sol[i], animated = True, color = 'red')
    ims.append(im)
ani = animation.ArtistAnimation(fig,ims,interval = (1000*k), blit = True)
plt.show()