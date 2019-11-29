import numpy as np
from scipy.optimize import minimize, shgo
from PDE_solver import FitzHugh_Nagumo_solver
from AnalyticalSolution import AnalyticalSolution

import matplotlib.pyplot as plt

def ModelFit():

    '''
    ADD DESCRIPTION!
    '''

    solver = FitzHugh_Nagumo_solver()

    u = solver.FN_solver(x_0 = 0, x_n = 20)
    
    # eps = 1.0
    # h = 0.05
    # k = 0.2 * h**2 / eps
    # Tf = 42.0
    var = 0.0 # scale of random normal noise

    np.random.seed(100) #setting the seed to ensure cosnsitency in random number generation

    # x = np.arange(0+h, 20-h, h)

    # numsteps = int(np.ceil(Tf/k)/10)

    t = np.linspace(0,1,solver.k_N)
    

    input_data = np.meshgrid(t, solver.x_range)

    # Generating noisy data using the Analytical Solution and randomly normal distributed noise

    # y_data = AnalyticalSolution(input_data[0], input_data[1]) + var * np.random.randn(len(solver.x_range), len(t),)

    # # Define the model parameters to be fitted

    # alpha = 0.2
    # beta = 1.0 
    # gama = 1.0

    y_data = AnalyticalSolution(input_data[0], input_data[1], alpha=solver.alpha) + var * np.random.randn(len(solver.x_range), len(t))
    y_solver = solver.FN_solver(x_0 = 0, x_n = 20)

    plt.figure()
    plt.plot(solver.x_range, y_data[:,0],label='analytical')
    plt.plot(solver.x_range, y_solver[:,-1],label='solver')
    plt.legend()
    plt.show()


    # def cost_function(parameters):

    #     '''
    #     ADD DESCRIPTION!
    #     '''

    #     # y_solver = solve_fitzhuge_nagumo(alpha=parameters[0],beta=parameters[1],gamma=parameters[2])
    #     y_data = AnalyticalSolution(input_data[0], input_data[1]) + var * np.random.randn(len(solver.x_range), len(t))
    #     y_solver = solver.FN_solver(x_0 = 0, x_n = 20,alpha=parameters[0])

    #     return np.sum(np.power(y_solver - y_data, 2))

    # # Randomly initializing the model parameters.
    # # As stated in the original paper, alpha needs to be within (0, 0.5)

    # alpha0 = np.random.normal(1) * var
    # while alpha0 < 0 or alpha0 > 0.5:
    #     alpha0 = np.random.normal(1) * var
    # # beta0 = np.random.normal(1) * var 
    # # gama0 = np.random.normal(1) * var

    # # original_parameters = [alpha0, beta0, gama0]
    # original_parameters=np.array([alpha0])

    # for method in ['nelder-mead', 'bfgs']:
    #     res = minimize(cost_function, original_parameters, method=method, options={'disp': True})
    #     print(res.x)
        
    # # res2 = shgo(cost_function, bounds=[(0.0,0.5),(-10,10),(-10,10)], options={'disp': True})
    # res2 = shgo(cost_function, bounds=[(0.0,0.5)], options={'disp': True})
    # print(res2.x)
    

if __name__ == '__main__':
    ModelFit()