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
In order to run the solver, type the following commands int the activated python environment. For the **solver** type in the first command, while for the **moldel fitting function** type in the second command. 

```python
for solver
```

```python
for model_fitting
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[BSD 3-Clause License](https://opensource.org/licenses/BSD-3-Clause) © MIT © [Alister Dale-Evans](https://github.com/alisterde), [Andrei Roibu](https://github.com/AndreiRoibu), [Pavan Chaggar](https://github.com/PavanChaggar), [Rebecca Rumney](https://github.com/Rebecca-Rumney)

## References
In the creation of this code, material was used from the following paper:

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
```

## Credits
The team would like to express their gratitude to [Martin Robinson](https://github.com/martinjrobins) and [Aleksandra Ardaseva](https://www.maths.ox.ac.uk/people/aleksandra.ardaseva) for their help and support during the development of this code.

## Build status
Describe and show how to run the tests with code examples.

Build status of continus integration i.e. travis, appveyor etc. Ex. - 

[![Build Status](https://travis-ci.org/akashnimare/foco.svg?branch=master)](https://travis-ci.org/akashnimare/foco)
[![Windows Build Status](https://ci.appveyor.com/api/projects/status/github/akashnimare/foco?branch=master&svg=true)](https://ci.appveyor.com/project/akashnimare/foco/branch/master)




