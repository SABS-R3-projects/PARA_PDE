from setuptools import setup, find_packages

setup(
    name='FitzHugh_Nagumo_PDE_Solver',
    version='0.0.1',
    description='Solver for the reaction-diffusion FitzHugh Nagumo PDE',
    license='BSD 3-clause license',
    maintainer='AndreiRoibu, PavanChaggar, Rebecca-Rumney, alisterde',
    maintainer_email='andrei-claudiu.roibu@dtc.ox.ac.uk, rebecca.rumney@wolfson.ox.ac.uk, alister.dale-evans@linacre.ox.ac.uk, pavanjit.chaggar@exeter.ox.ac.uk',
    install_requires=[
        'numpy',
        'matplotlib',
	    'scipy',
	    'ipython',
        'pandas',
        'sympy'
    ],
)
