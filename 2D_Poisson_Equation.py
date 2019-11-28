import numpy as np
import scipy as sp
from scipy import sparse
from scipy.sparse import linalg
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d

x_0 = 0
x_N = 1
steps = 10
h = (x_N - x_0)/10

e = np.ones(steps)
diagonals = [e, -2*e, e]
offset = [-1, 0, 1]
L_m =  sp.sparse.spdiags(diagonals, offset, steps, steps)

identity = np.identity(steps)
X = sparse.kron(identity, L_m)
Y = sparse.kron(L_m,identity)
L2D = X + Y


A = L2D * (1/h**2)
B = -1* np.ones(steps**2)
sol = sp.sparse.linalg.spsolve(A,B)

sol_matrix = np.reshape(sol, (steps,steps))
x = np.linspace(-1,1, steps)
y = np.linspace(-1,1, steps)
x, y = np.meshgrid(x, y)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(x, y, sol_matrix)
plt.show()
