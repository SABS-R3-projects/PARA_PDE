import numpy as np
from scipy.optimize import minimize, shgo
# from PDE_solver import FitzHugh_Nagumo_solver
from AnalyticalSolution import AnalyticalSolution
from autograd import grad
from solver_using_fourier_transforms import fitzhugh_nagumo_solver

import matplotlib.pyplot as plt

def ModelFit():

    '''Function which performs PDE Model Fitting using the Solver the Analytical Solutions of the FitzHugh-Nagumo Operation.

    The function loops through several optimizers, which depending on hardware might make the runtime long. 

    Two solvers are provided for this function. One is the PDE_solver, which represents a simple approach, while the other is the fitzhugh_nagumo_solver, which employs a accurate space derivative method (ASD)

    As this code is still in development, the PDE_solver does not yet work, so the user is encouraged to select the fitzhugh_nagumo_solver. Instructions are in code. 
    
    '''
    var = 0.0 # scale of random normal noise
    np.random.seed(100) #setting the seed to ensure cosnsitency in random number generation

    ##### To select the PDE_solver, un-cooment this code #####

    # solver = FitzHugh_Nagumo_solver() # Loading the solver
    # u = solver.FN_solver(x_0 = 0, x_n = 20)
    # t = np.linspace(0,solver.end_time,solver.k_N)
    # x = solver.x_range

    ##### To select the fitzhugh_nagumo_solverl, un-comment this code #####

    u, max_t, max_x = fitzhugh_nagumo_solver()
    num_timesteps = np.size(u[:,0])
    num_x_steps = np.size(u[0,:])
    t = np.linspace(0, max_t, num_timesteps)
    x = np.linspace(0, max_x, num_x_steps)

    ##### END #####

    T,X = np.meshgrid(t, x)


    def cost_function(x):

        '''Function which returns the L2 squared error between solver data and noisy data obtained using the analytical solution.

        Args:
            x (float): variable(s) to be optimized 

        Returns:
            (float): L2 squared error
        '''

        ##### To select the PDE_solver, un-cooment this code #####

        ### Due to the PDE_solver currently working correctly, the y_data is generated using itself plus added noise, to test the correctness of the optimizer methods

        # y_data = AnalyticalSolution(input_data[0], input_data[1]) + var * np.random.randn(len(solver.x_range), len(t))

        # u = solver.FN_solver(x_0 = 0, x_n = 20, alpha=0.2)
        # y_solver = solver.FN_solver(x_0 = 0, x_n = 20, alpha=x)
        # y_data = u + var * np.random.randn(len(solver.x_range), len(t))


        ##### To select the fitzhugh_nagumo_solverl, un-comment this code #####

        y_data = AnalyticalSolution(X, T, alpha=0.2) + var * np.random.randn(len(x), len(t))

        y_solver, _, _ = fitzhugh_nagumo_solver(alpha=x, beta=1, gamma=1, max_x=1, max_t=10)

        y_solver = y_solver.T

        ##### END #####

        return np.sum(np.power(y_solver - y_data, 2))

    # Randomly initializing the model parameters.
    # As stated in the original paper, alpha needs to be within (0, 0.5)

    alpha0 = np.random.normal(1) * var
    while alpha0 < 0 or alpha0 > 0.5 or alpha0 == 0.2:
        alpha0 = np.random.normal(1) * var*10

    # beta0 = np.random.normal(1) * var 
    # gama0 = np.random.normal(1) * var

    # original_parameters = [alpha0, beta0, gama0]
    # bounds = [(0.0, 0.5), (-10.0, 10.0), (-10.0, 10.0)]

    original_parameters=np.array([alpha0])
    bounds = [(0.0, 0.5)]

    desired_parameters= np.array( [0.2, 1.0, 1.0] )

    def optimize(cost_function, method, autodiff, original_parameters, bounds):

        '''Function which runs the different optimizers for to determine parameter fitting. 

        Args:
            cost_function (function): variable(s) to be optimized 
            method (string): string indicating the optimisation method to be used
            autodiff (boolean): boolean indicating if auto differentiation is to be employed with the considered method
            original_parameters (floats): vector of floats containing the parameters wich initialize the optimizers
            bounds (floats): vector containing the bounds for the shgo optimizer

        Returns:
            res.x (float): a vector contains the predicted float values obtained after optimisation
        '''

        if autodiff:
            jac = grad(cost_function)
        else:
            jac = None
        
        if method == 'shgo':
            res = shgo(cost_function, bounds, options={'disp': True})
        else:
            res = minimize(cost_function, original_parameters, method=method, jac = jac, options={'disp': True})

        return res.x

    for method in ['shgo','nelder-mead','cg','bfgs']:
        for autodiff in [False, True]:
            if method == 'newton-cg' and autodiff == False:
                continue
            if method == 'shgo' and autodiff == True:
                continue
            if method == 'cg' or method == 'bfgs' and autodiff == True:
                continue
            
            if autodiff==True:
                txt = " with auto differentiation "
            else:
                txt = " without auto differentiation"

            print("Running method ", method, txt)
            estimated_params = optimize(cost_function, method, autodiff, original_parameters, bounds)
            error = original_parameters - estimated_params
            print("The estimated [alpha, beta, gamma] parameters are: ", estimated_params, "and the original parameters are:", desired_parameters)
    

if __name__ == '__main__':
    ModelFit()