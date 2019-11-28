import numpy as np

def AnalyticalSolution(t, t0, kappa, tau, delta):
    '''
    Function which returns the Analytical Solution of the FitzHugh-Nagumo Operation.

    Args:
        t (float): model constant
        t0 (float): model constant
        kappa (float): model constant
        tau (float): model constant
        delta (float): model constant

    Returns:
        V (float): potential of membrane
    
    Source of mathematical model is below, with equation being (4.20):

    @article{kudryashov2018asymptotic,
    title={Asymptotic and exact solutions of the Fitzhugh--Nagumo model},
    author={Kudryashov, Nikolay A},
    journal={Regular and Chaotic dynamics},
    volume={23},
    number={2},
    pages={152--160},
    year={2018},
    publisher={Springer}
    }

    '''

    V = np.sqrt( 1.0 / 2.0 + kappa / (2 * tau) + np.sqrt( delta ) / 2 * np.tanh ( np.sqrt( delta ) * ( t - t0 ) ) )

    return V