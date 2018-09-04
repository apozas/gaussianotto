[![DOI](https://zenodo.org/badge/100956997.svg)](https://zenodo.org/badge/latestdoi/100956997)

## Code to accompany *[A quantum Otto engine with finite heat baths: energy, correlations and degradation](https://dx.doi.org/10.1088/1367-2630/aaba02)*
#### Alejandro Pozas-Kerstjens, Karen V. Hovhannisyan, and Eric G. Brown

This is a repository for all code which was written for the article "*A quantum Otto engine with finite heat baths: energy, correlations and degradation*. Alejandro Pozas-Kerstjens, Karen V. Hovhannisyan, and Eric G. Brown. [New Journal of Physics 20, 043034 (2018)](https://dx.doi.org/10.1088/1367-2630/aaba02)."

All code is written both in MATLAB and Python. The MATLAB version requires no additional toolboxes, and the Python version requires NumPy, SciPy and Matplotlib.

The code is separated in two parts: the functions required to compute the simulations and the simulations themselves.

- Functions required: 
The functions found in MATLAB's [functions](https://github.com/apozas/gaussianotto/tree/master/MATLAB/functions) folder (or Python's [functions.py](https://github.com/apozas/gaussianotto/tree/master/Python/functions.py) file) are those needed for creating the Hamiltonian matrices, the symplectic time evolution operators, and the initial states of the baths, among others.

  - [Energy](https://github.com/apozas/gaussianotto/blob/master/MATLAB/functions/Energy.m): calculates the energy of a machine-bath(s) system.
  
  - [Entropy](https://github.com/apozas/gaussianotto/blob/master/MATLAB/functions/Entropy.m): calculates the entropy of a state characterized by a covariance matrix.
  
  - [FreeRing](https://github.com/apozas/gaussianotto/blob/master/MATLAB/functions/FreeRing.m): creates the Hamiltonian matrix of a translation-invariant chain of harmonic oscillators with nearest-neighbour interactions.
  
  - [Initialize](https://github.com/apozas/gaussianotto/blob/master/MATLAB/functions/Initialize.m): creates the covariance matrix of the initial (thermal) state of a bath.
  
  - [MakeInt](https://github.com/apozas/gaussianotto/blob/master/MATLAB/functions/MakeInt.m): constructs the machine-bath(s) interaction Hamiltonian matrix.
  
  - [MakeS](https://github.com/apozas/gaussianotto/blob/master/MATLAB/functions/MakeS.m): calculates an array of symplectic time evolution operators of the system for a series of time steps.
  
  - [MakeStimeIndep](https://github.com/apozas/gaussianotto/blob/master/MATLAB/functions/MakeStimeIndep.m): calculates the symplectic time evolution operator of the machine-bath(s) system for one full interaction between the machine and one bath.
  
  - [Ordering](https://github.com/apozas/gaussianotto/blob/master/MATLAB/functions/Ordering.m): creates the matrix that transforms between the mode-mode ordering and the $q$-$p$ ordering of the covariance matrices.
  
  - [Switching](https://github.com/apozas/gaussianotto/blob/master/MATLAB/functions/Switching.m): creates the switching function $\lambda(t)$.
  
- Simulations:
The simulations we run are those that appear in the figures and discussions of the manuscript. The code made available here is that used for creating the figures, and can be used for further considerations.
  
  - correlations: models an Otto cycle between a hot and a cold bath, and computes the machine-baths and bath-bath correlations during the joint evolution over a number of interaction cycles. This code is used for creating Figs. 5 and 6 in the manuscript, as well as the animation of the bath-bath correlations.
  
  - ottocycle: models an Otto cycle between a hot and a cold bath, computing the energetics of the cycle (work output, heat loss and efficiency). This code is used for creating Fig. 4 in the manuscript.
  
  - singlebath: models the interaction of the machine with a single heat bath, and computes the temperature of the machine, its athermality, and different mutual informations during the joint evolution. This code is used for creating Figs. 1 and 7 in the manuscript.

  - temperature: models the interaction of the machine with a single heat bath, and computes the temperature of the machine. This code is used for creating Fig. 2 in the manuscript.