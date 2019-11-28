'''
Class to create and solve Fitzhugh-Nagumo equation

'''

import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt
import time
from matplotlib.animation import FuncAnimation

class FitzHugh_Nagumo_solver(object):
    def __init__ (self):
        pass 
    
    def fitzhugh_nagumo(self, x_0, a_0: float):

        self.x_0 = x_0
        self.a_0 = a_0
        alpha = self.a_0
        size = len(self.x_0)
        u = np.zeros(size)

        for i in range(size):
            u[i] = 0.5*(1 + alpha) + 0.5*(1- alpha)*(np.tanh((np.sqrt(2)*(1-alpha)*x[i])/4))
            return u


    def __laplace_matrix(self):

        N = len(self.x_0)
        e = np.ones(N)
        diagonals = [e, -2*e, e]
        offsets = [-1, 0, 1]
        self.L = scipy.sparse.spdiags(diagonals, offsets, N, N) / h**2

        return self.L

    def FN_solver(self, boundary_conditions, step_size):
        self.h = step_size
        self.bc = boundary_conditions

        

    def plot(self, parameter_list):
        pass

