import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt
import time
from matplotlib.animation import FuncAnimation

class FitzHugh_Nagumo_solver(object):
    '''
    Class to create and solve Fitzhugh-Nagumo equation

    '''

    def __init__ (self):
        pass
    
    def _fitzhugh_nagumo_initial_conditions(self):
        '''
        Method to set the initial conditions of the Nagumo equation for iterative solving 
        using the analytical solution

        return: initial condition with given parameters
        '''
        u = np.zeros(self.N_dim)

        for i in range(self.N_dim):
            u[i] = 0.5*(1 + self.alpha) + 0.5*(1- self.alpha)*(np.tanh((np.sqrt(2)*(1-self.alpha)*self.x_range[i])/4))
        
        return u


    def __laplace_matrix(self):
        '''
        Defines Laplace Matrix with dimensions of X
        Returns: The laplace matrix for a second order differential 
        '''
        e = np.ones(self.N_dim)
        diagonals = [e, -2*e, e]
        offsets = [-1, 0, 1]
        self.L = scipy.sparse.spdiags(diagonals, offsets, self.N_dim, self.N_dim) /  self.h **2

        return self.L

    def FN_solver(self, x_0 : int, x_n: int, boundary_conditions: tuple = [0,1], step_size: float = 0.05,
                 time_steps: int = 8000, alpha: float = 0.13):
        '''Iterative method of solving the Fitzhuge-Nagumo system when episolon is very small as in most
            neuroscience applications know as the Nagumo equation:

                dv/dt=d^2V/dx^2 + v(1-v)(v-a), where t > 0 and X exists in the reals
            param x_0: lower bound of positional dimensions x
            param x_N: upper bound of positional dimensions x
            boundary_conditions
            param N: dimensions of the laplace matirx
            param h: step size 
            param alpha: alpha is a constant of the equation which should obey 0 < alpha < 0.5 
            the default value 0.2 is used if not specified
        '''

        self.h = step_size
        self.x_range = np.arange(x_0+self.h, x_n-self.h, self.h)
        self.alpha = alpha
        self.N_dim = len(self.x_range)
        self.k_N = time_steps
        lower, upper = boundary_conditions[:]
        self.bc = bc = np.concatenate(([lower], np.zeros(self.N_dim-2), [upper]))/self.h**2

        u = np.empty((self.N_dim, self.k_N))
        u[:,0] = self._fitzhugh_nagumo_initial_conditions()
        k = 0.2 * self.h**2
        L = self.__laplace_matrix()
        #int(np.ceil(42/k)/10) DO WE NEED THIS?!

        for i in range(1, self.k_N):
            u[:,i] = u[:,i-1] + k*( (L@u[:,i-1] + self.bc) + (u[:,i-1]**2 - u[:,i-1]**3 - self.alpha*u[:,i-1] + self.alpha*u[:,i-1]**2) )
        
        return u

# Tf: float = 42.0):
    def plot(self, parameter_list):
        pass

    def __update(self, frame):
        u =_fitzhugh_nagumo_initial_conditions()
        for i in range(10):
            u_new = u + k*( eps*(L@u + bc) + (u**2 - u**3 - a*u + a*u**2) )
            u[:] = u_new

        ln.set_data(x, u)
        ax.set_title('t = {}'.format(10*frame*k))
    
    def animation(self):
    pass
