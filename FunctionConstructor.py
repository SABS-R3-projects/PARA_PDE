import numpy as np

def Constructor(v,a):
    '''
    Function which constructs the righthand side of the FitzHugh-Nagumo Operation.

    Args:
        v (float): transmembrane potential
        a (float): positive parameter e [0, 0.5]

    Returns:
        (float): difference between dvdt and d^2v/dx^2 
    
    Source of mathematical model is below:

    @article{feng2015finite,
    title={A finite difference method for the FitzHugh-Nagumo equations},
    author={Feng, Hongsong and Lin, Runchang},
    journal={Dynamics of Continuous, Discrete and Impulsive Systems. Series B: Applications \& Algorithms},
    volume={22},
    pages={401--412},
    year={2015}
    }

    '''

    return  v * (1 - v) * (v - a)