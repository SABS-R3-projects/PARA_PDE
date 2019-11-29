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
        Defines Laplace Matrix with dimensions of X as a sparse matrix

        Returns: The laplace matrix for a second order differential 
        '''

        e = np.ones(self.N_dim)
        diagonals = [e, -2*e, e]
        offsets = [-1, 0, 1]
        self.L = self.beta*(scipy.sparse.spdiags(diagonals, offsets, self.N_dim, self.N_dim)) /  self.h **2

        return self.L


    def FN_solver(self, x_0 : int, x_n: int, boundary_conditions: tuple = [0,1], step_size: float = 0.05,
                 time_steps: int = 8000, alpha: float = 0.2, beta: float = 1.0, gamma: float = 1.0):
        '''
            Iterative method of solving the Fitzhuge-Nagumo system when episolon is very small as in most
            neuroscience applications know as the Nagumo equation:

                dv/dt= beta * d^2V/dx^2 + gamma * v(1-v)(v-alpha), where t > 0 and X exists in the reals

            param x_0: lower bound of positional dimensions x, this must be specified.
            param x_N: upper bound of positional dimensions x, this must be specified.
            param boundary_conditions: The condition of the system at the lower and upper bound respectively, the default value is [0,1].
            param step_size: Step Size of the positional dimension, the default value is 1.
            param time_steps: Number of steps in time to be used, the default value is 8000.
            param alpha: A constant of the equation which should obey 0 < alpha < 0.5, the default value is 0.2.
            param Beta: A constant of the equation, the default value is 1.
            param gamma: A constant of the equation, the default value is 1.
            return: A matrix containing columns of V values for each position x at a given t (each column is a given t and each row is a given x)
        '''

        #defining class variables
        self.h = step_size
        self.x_range = np.arange(x_0+self.h, x_n-self.h, self.h)
        self.alpha = alpha
        self.beta= beta
        self.gamma = gamma
        self.N_dim = len(self.x_range)
        self.k_N = time_steps
        lower, upper = boundary_conditions[:]
        self.bc = bc = np.concatenate(([lower], np.zeros(self.N_dim-2), [upper]))/self.h**2

        #initialising an empty matrix to contain the calculated solutions 
        u = np.empty((self.N_dim, self.k_N))
        #placing initial conditions into the first column
        u[:,0] = self._fitzhugh_nagumo_initial_conditions()

        #defining th time step size k
        k = 0.2 * self.h**2
        #calculating the time the simulation runs for 
        self.end_time = k*self.k_N 

        #calculating the laplace matrix
        L = self.__laplace_matrix()

        #iterative finite difference method
        for i in range(1, self.k_N):
            u[:,i] = u[:,i-1] +  k*( (L@u[:,i-1] + self.bc) 
            + self.gamma*(u[:,i-1]**2 - u[:,i-1]**3 - self.alpha*u[:,i-1] + self.alpha*u[:,i-1]**2) )
        
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
