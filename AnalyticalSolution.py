import numpy as np

def AnalyticalSolution(x, t, alpha=0.2):
    '''
    Function which returns the Analytical Solution of the FitzHugh-Nagumo Operation.
    Args:
        x (float): space variable
        t (float): time variable
        alpha (float): model constant, defined by default to be 0.2 based on the original paper. 

    Returns:
        u (float): potential of membrane

    Source of mathematical model is below:
    
    @article{parand2017numerical,
    title={A numerical method to solve the 1D and the 2D reaction diffusion equation based on Bessel functions and Jacobian free Newton-Krylov subspace methods},
    author={Parand, K and Nikarya, M},
    journal={The European Physical Journal Plus},
    volume={132},
    number={11},
    pages={496},
    year={2017},
    publisher={Springer}
    }
    '''

    u = 1.0 / 2.0 * (1 + alpha) + 1.0/2.0 * (1 - alpha) * np.tanh( np.sqrt(2.0)/4.0 * (1-alpha)*x + (1-np.power(alpha,2)) * t ) 

    return u