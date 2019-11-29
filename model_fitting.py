import numpy as np
from scipy.optimize import minimize, shgo
from PDE_solver import FitzHugh_Nagumo_solver
from AnalyticalSolution import AnalyticalSolution
from autograd import grad

import matplotlib.pyplot as plt

def ModelFit():

    '''Function which performs PDE Model Fitting using the Solver the Analytical Solutions of the FitzHugh-Nagumo Operation.

    The function loops through several optimizers, which depending on hardware might make the runtime long. 
    
    '''

    solver = FitzHugh_Nagumo_solver() # Loading the solver

    u = solver.FN_solver(x_0 = 0, x_n = 20)
    

    var = 0.0 # scale of random normal noise

    np.random.seed(100) #setting the seed to ensure cosnsitency in random number generation


    t = np.linspace(0,solver.end_time,solver.k_N)

    T,X = np.meshgrid(t, solver.x_range)

    # Generating noisy data using the Analytical Solution and randomly normal distributed noise

    y_data = AnalyticalSolution(X, T, alpha=solver.alpha) + var * np.random.randn(len(solver.x_range), len(t))
    y_solver = solver.FN_solver(x_0 = 0, x_n = 20)



    plt.figure()
    plt.plot(solver.x_range, y_data[:,-1], '*r',label='analytical')
    plt.plot(solver.x_range, y_solver[:,-1],label='solver')
    plt.legend()
    plt.show()


    def cost_function(x):

        '''Function which returns the L2 squared error between solver data and noisy data obtained using the analytical solution.

        Args:
            x (float): variable(s) to be optimized 

        Returns:
            (float): L2 squared error
        '''

        # y_data = AnalyticalSolution(input_data[0], input_data[1]) + var * np.random.randn(len(solver.x_range), len(t))
        u = solver.FN_solver(x_0 = 0, x_n = 20, alpha=0.3)
        y_solver = solver.FN_solver(x_0 = 0, x_n = 20, alpha=x)
        y_data = u + var * np.random.randn(len(solver.x_range), len(t))

        # print(np.sum(np.power(y_solver - y_data, 2)))

        return np.sum(np.power(y_solver - y_data, 2))

    # Randomly initializing the model parameters.
    # As stated in the original paper, alpha needs to be within (0, 0.5)

    alpha0 = np.random.normal(1) * var
    while alpha0 < 0 or alpha0 > 0.5 or alpha0 == 0.2:
        alpha0 = np.random.normal(1) * var*10

    # beta0 = np.random.normal(1) * var 
    # gama0 = np.random.normal(1) * var

    # original_parameters = [alpha0, beta0, gama0]
    original_parameters=np.array([alpha0])
    bounds=[(0.0,0.5)]

    desired_parameters= np.array( [0.3, 1.0, 1.0] )


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