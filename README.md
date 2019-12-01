# FitzHugh-Nagumo Reaction-Diffusion Equation Solver

## Motivation
This projects aims at providing a solver for the FitzHugh-Nagumo Reaction-Diffusion nonlinear PDE, as well as code which can extract the PDE's parameters.

This code was created as an assignment for the _Modelling & Scientific Computing_ module, at the _EPSRC CDT in Sustainable Approached to Biomedical Sciences: Responsible and Reproducible Research - SABS R3_

## Installation
In order to install the packaged, the user will require the presence of Python3 and the [pip3](https://pip.pypa.io/en/stable/) installer. 

For installation on Linux or OSX, use the following commands. This will create an environment and automatically install all the requirements.

```bash
python3 -m venv env
source env/bin/activate
pip install -e .
```

## Usage
In order to run the solver, type the following commands int the activated python environment. For a **simple solver** type in the first command, for a **fourier transform based solver** ued the second command, while for the **moldel fitting function** type in the second command. 

```python
python PDE_solver.py
```

```python
python solver_using_fourier_transforms.py
```

```python
python model_fitting.py
```

The user should know that, in it's current format, the codes are not completely operational. The _PDE-solver_ is currently fails to capture the analytical solution. The _solver-using-fourier-transforms_ currently is operational, however it still has certain issues. The _model-fitting_ script is operational, however it still suffers from dimensionality issues due to incompatibilities with the previous two solvers. In previous tests, however, it proved to be working.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[BSD 3-Clause License](https://opensource.org/licenses/BSD-3-Clause) © [Alister Dale-Evans](https://github.com/alisterde), [Andrei Roibu](https://github.com/AndreiRoibu), [Pavan Chaggar](https://github.com/PavanChaggar), [Rebecca Rumney](https://github.com/Rebecca-Rumney)

## References
In the creation of this code, material was used from the following papers:

```
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

@article{10.2307/3212689,
    ISSN = {00219002},
    URL = {http://www.jstor.org/stable/3212689},
    author = {Jenö Gazdag and José Canosa},
    journal = {Journal of Applied Probability},
    number = {3},
    pages = {445--457},
    publisher = {Applied Probability Trust},
    title = {Numerical Solution of Fisher's Equation},
    volume = {11},
    year = {1974}
    }
```

## Credits
The team would like to express their gratitude to [Martin Robinson](https://github.com/martinjrobins) and [Aleksandra Ardaseva](https://www.maths.ox.ac.uk/people/aleksandra.ardaseva) for their help and support during the development of this code.

## Build status
Due to problems and issues emerging during the development of the code, no test files were yet created. These will be added later. 




