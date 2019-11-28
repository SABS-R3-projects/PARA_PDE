import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt
import time
from matplotlib.animation import FuncAnimation

def construct_laplace_matrix_1d(N: int, h: float):
    '''Method to construct a sparse laplace matrix using scipy

    param N: dimensions of the laplace matirx
    param h: step size 
    '''
    e = np.ones(N)
    diagonals = [e, -2*e, e]
    offsets = [-1, 0, 1]
    L = scipy.sparse.spdiags(diagonals, offsets, N, N) / h**2
    return L

def solve_fitzhuge_nagumo():
    '''Iterative method of solving the Fitzhuge-Nagumo system when episolon is very small as in most
    neuroscience applications know as the Nagumo equation:

                dv/dt=d^2V/dx^2 + v(1-v)(v-a), where t > 0 and X exists in the reals
    '''
    eps = 1.0
    h = 0.05
    k = 0.2 * h**2 / eps
    a= 0.13
    Tf = 42.0
    x = np.arange(0+h, 20-h, h)

    N = len(x)
    L = construct_laplace_matrix_1d(N, h)
    bc = np.concatenate(([0], np.zeros(N-2), [1]))/h**2 # boundary conditions

    def initial_conditions(x: list, alpha: float = 0.2):
        '''Method to set the initial conditions of the Nagumo equation for iterative solving 
        using the analytical solution

        param x: list of position coordinate
        param alpha: alpha is a constant of the equation which should obey 0 < alpha < 0.5 
                     the default value 0.2 is used if not specified
        return: initial condition with given parameters
        '''

        size = len(x) 
        u = np.zeros(size)
        for i in range(size):
            u[i] = 0.5*(1 + alpha) + 0.5*(1- alpha)*(np.tanh((np.sqrt(2)*(1-alpha)*x[i])/4))

        return u

    u = initial_conditions(x,0.2)
    initial = u
    plt.plot(x,u)
    plt.show()

    fig, ax = plt.subplots()
    ln, = plt.plot(x, u)
    ax.set_ylim(0, 1.1)

    numsteps = int(np.ceil(Tf/k)/10) # Based on the final timestep and the step size K, it works out how many frames we have
    Tf = numsteps*k
    #solv = np.zeros()

    def update(frame):
        for i in range(10):
            u_new = u + k*( eps*(L@u + bc) + (u**2 - u**3 - a*u + a*u**2) )
            u[:] = u_new

        ln.set_data(x, u)
        ax.set_title('t = {}'.format(10*frame*k))
        #print(frame)
        return ln,

    ani = FuncAnimation(fig, update, frames=numsteps, interval=30, blit=False, repeat=False)
    plt.show()

    plt.plot(x,u)
    plt.plot(x,initial)
    plt.show()

    return u

solution = solve_fitzhuge_nagumo()
print(solution.shape)