import sys
import math
import numpy as np
from Iterative_PDE_solver import FitzHugh_Nagumo_solver

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

import unittest
from typing import List

class test_Iterative_PDE_solver(unittest.TestCase):
    """This test suite is about testing the iterative  FitzHugh_Nagumo_solver """


    def test_solver_output_exists(self):
        """Tests the PDE solver outputs a argument that is not of null value
        """
        trial = FitzHugh_Nagumo_solver()
        u = trial.FN_solver(x_0=0, x_n=20, step_size  = 0.05,
                 time_steps = 80, alpha = 0.13)
        self.assertIsNotNone(u)
    
    def test_solver_output_dimensions(self):
        """Tests the PDE solver outputs a argument of the correct dimensions
        """
        trial = FitzHugh_Nagumo_solver()
        u = trial.FN_solver(x_0=0, x_n=20, step_size  = 0.05,
                 time_steps = 80, alpha = 0.13)
        Dimensions = u.shape
        self.assertEqual((trial.N_dim, trial.k_N), Dimensions)
if __name__ == '__main__':
    unittest.main()