import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt
import time
import matplotlib.animation as animation
from animation_and_polting import animate


class FitzHugh_Nagumo_solver(object):
    '''
    Class to create and solve Fitzhugh-Nagumo equation

    '''

    def __init__ (self):
        pass
    
    def _fitzhugh_nagumo(self, x, time):

        '''
        Method to set the initial conditions of the Nagumo equation for iterative solving 
        using the analytical solution

        return: initial condition with given parameters
        '''

        u = np.zeros(self.N_dim)
        if type(x) == int:
            u = 1.0 / 2.0 * (1 + self.alpha) + 1.0/2.0 * (1 - self.alpha) * np.tanh( np.sqrt(2.0)/4.0 * (1-self.alpha)*x + (1-np.power(self.alpha,2)) * time)
        else:
            size = len(x)
            for i in range(size):
                u[i] = 1.0 / 2.0 * (1 + self.alpha) + 1.0/2.0 * (1 - self.alpha) * np.tanh( np.sqrt(2.0)/4.0 * (1-self.alpha)*x[i] + (1-np.power(self.alpha,2)) * time)
                
        return u


    def __laplace_matrix(self):
        '''
        Defines Laplace Matrix with dimensions of X
        Returns: The laplace matrix for a second order differential

        '''

        e = np.ones(self.N_dim)
        diagonals = [e, -2*e, e]
        offsets = [-1, 0, 1]
        self.L = self.beta*(scipy.sparse.spdiags(diagonals, offsets, self.N_dim, self.N_dim)) /  self.h **2

        return self.L

    def FN_solver(self, x_0 : int, x_n: int, step_size: float = 0.05,
                 time_steps: int = 80, alpha: float = 0.13):
        '''Iterative method of solving the Fitzhuge-Nagumo system when episolon is very small as in most

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

        self.h = x_n/time_steps # positions step size
        self.x_range = np.arange(x_0, x_n+self.h, self.h) # numpy array of x poistions used 
        self.alpha = alpha
        self.beta= 1.0 # beta vairable may add later
        self.gamma = 1.0 # gamma vairable may add later
        self.N_dim = len(self.x_range) # dimensions of x poistions varaible (number of X poistions)
        self.k_N = time_steps # number of time steps

        #initialising an empty matrix to contain the calculated solutions 
        u = np.empty((self.N_dim, self.k_N))
        u[:,0] = self._fitzhugh_nagumo(x = self.x_range,time=0)
        self.k = 0.2 * self.h**2
        #print(self.k)

        L = self.__laplace_matrix()

        #iterative finite difference method
        for i in range(1, self.k_N):
            lower = self._fitzhugh_nagumo(x_0,time=i*self.k)
            #print(lower)
            #print(i*k)
            upper = self._fitzhugh_nagumo(x_n,time=i*self.k)
            bc = np.concatenate(([lower], np.zeros(self.N_dim-2), [upper]))/self.h**2
            u[:,i] = u[:,i-1] + self.k*( (L@u[:,i-1] + bc) + (u[:,i-1]**2 - u[:,i-1]**3 - self.alpha*u[:,i-1] + self.alpha*u[:,i-1]**2) )
            #print(u.shape)
        return u


if __name__ == '__main__':
    trial = FitzHugh_Nagumo_solver()
    u = trial.FN_solver(0,20)
    plt.plot(u[:,-1])
    plt.show()


    animate(u, trial.k_N, trial.x_range, trial.k)

