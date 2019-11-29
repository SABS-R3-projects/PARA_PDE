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

def solve_fitzhuge_nagumo(alpha: float=0.2, beta:float = 1.0, gamma:float = 1.0, eps: float = 1.0, h: float = 0.05, a: float = 0.13, Tf: float = 42.0):
    '''Iterative method of solving the Fitzhuge-Nagumo system when episolon is very small as in most
    neuroscience applications know as the Nagumo equation:
                dv/dt=d^2V/dx^2 + v(1-v)(v-a), where t > 0 and X exists in the reals
    '''
    
    k = 0.2 * h**2 / eps
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
    u_output = u
    print("u:",type(u))
    print("u_output:",type(u_output))
    initial = u
    plt.plot(x,u)
    plt.show()

    fig, ax = plt.subplots()
    ln, = plt.plot(x, u)
    ax.set_ylim(0, 1.1)

    numsteps = int(np.ceil(Tf/k)/10) # Based on the final timestep and the step size K, it works out how many frames we have
    Tf = numsteps*k
    
    def update(frame):
        u_new = u[:]
        for i in range(10):
            u_new = u_new + k*( eps*(L@u_new + bc) + (u_new**2 - u_new**3 - alpha*u_new + alpha*u_new**2) )
            u_new[:] = u_new
        u[:] = u_new

        ln.set_data(x, u)
        ax.set_title('t = {}'.format(10*frame*k))
        return ln,

    t = np.linspace(0,1,numsteps)

    def solver(u_out, t):
        for t_iter,i in enumerate(t):
            if i==0:
                pass
            else:
                u_new = u[:]
                for i in range(10):
                    u_new = u_new + k*( eps*(L@u_new + bc) + (u_new**2 - u_new**3 - a*u_new + a*u_new**2) )
                    u_new[:] = u_new
                u[:] = u_new
                u_out = np.append(u_out, u_new, axis=0)
        return u_out

    u_output = solver(u_output, t)
    u_output = np.reshape(u_output,np.meshgrid(t,u)[0].shape)

    ani = FuncAnimation(fig, update, frames=numsteps, interval=30, blit=False, repeat=False)
    plt.show()

    plt.plot(x,u)
    plt.plot(x,initial)
    plt.show()

    return u_output

solve_fitzhuge_nagumo()