import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt
import time
import matplotlib.animation as animation

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
        self.L = scipy.sparse.spdiags(diagonals, offsets, self.N_dim, self.N_dim) /  self.h **2

        return self.L

    def FN_solver(self, x_0 : int, x_n: int, step_size: float = 0.05,
                 time_steps: int = 80, alpha: float = 0.13):
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

        self.h = x_n/time_steps
        self.x_range = np.arange(x_0+self.h, x_n-self.h, self.h)
        self.alpha = alpha
        self.N_dim = len(self.x_range)
        self.k_N = time_steps

        u = np.empty((self.N_dim, self.k_N))
        u[:,0] = self._fitzhugh_nagumo(x = self.x_range,time=0)
        k = 0.2 * self.h**2
        print(k)
        L = self.__laplace_matrix()
        #int(np.ceil(42/k)/10) DO WE NEED THIS?!

        for i in range(1, self.k_N):
            lower = self._fitzhugh_nagumo(x = 0,time=i*k)
            #print(lower)
            #print(i*k)
            upper = self._fitzhugh_nagumo(x = 20,time=i*k)
            bc = np.concatenate(([lower], np.zeros(self.N_dim-2), [upper]))/self.h**2
            u[:,i] = u[:,i-1] + k*( (L@u[:,i-1] + bc) + (u[:,i-1]**2 - u[:,i-1]**3 - self.alpha*u[:,i-1] + self.alpha*u[:,i-1]**2) )
            print(u.shape)
        return u





trial = FitzHugh_Nagumo_solver()
u = trial.FN_solver(0,20)
print(u.shape)
plt.plot(u[:,-1])
plt.show()


timeSteps = trial.k_N
postionSteps = 398
fig= plt.figure()
ims = []
for i in range(timeSteps):
    im = plt.plot(u[:,i] , animated = True, color = 'red')
    ims.append(im)

ani = animation.ArtistAnimation(fig, ims, interval = (10), blit = True)
plt.show()

