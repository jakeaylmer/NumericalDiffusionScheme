# NumericalDiffusionScheme

## Overview
This repository contains some Python (2.7) code which may be used to solve the one-dimensional diffusion equation with variable diffusivity *k*(*x*) and non-zero source term *S*(*x*, *t*), subject to Neumann boundary conditions:

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;q}{\partial&space;t}-\frac{\partial}{\partial&space;x}\left[k(x)\frac{\partial&space;q}{\partial&space;x}&space;\right&space;]=S(x,t)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;q}{\partial&space;t}-\frac{\partial}{\partial&space;x}\left[k(x)\frac{\partial&space;q}{\partial&space;x}&space;\right&space;]=S(x,t)." title="\frac{\partial q}{\partial t}-\frac{\partial}{\partial x}\left[k(x)\frac{\partial q}{\partial x} \right ]=S(x,t)." /></a>

with

<a href="https://www.codecogs.com/eqnedit.php?latex=\left.k(x)\frac{\partial&space;q}{\partial&space;x}\right|_{x=0}=\left.k(x)\frac{\partial&space;q}{\partial&space;x}\right|_{x=L}=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left.k(x)\frac{\partial&space;q}{\partial&space;x}\right|_{x=0}=\left.k(x)\frac{\partial&space;q}{\partial&space;x}\right|_{x=L}=0" title="\left.k(x)\frac{\partial q}{\partial x}\right|_{x=0}=\left.k(x)\frac{\partial q}{\partial x}\right|_{x=L}=0" /></a>

where 0 < _x_ < _L_.

The Documentation sub-directory contains information on the numerical scheme and how it is derived and set up. The code has been written to match the symbols in the document "Implementation_Notes.pdf" as closely as possible.

## Usage
In the code sub-directory are three python files:
- `diffusion_scheme.py`
- `plotting.py`
- `test_me.py`

The actual numerical scheme is contained within two functions in `diffusion_scheme.py`. `test_me.py` can be run directly from the command line set in the top (repository) directory using

`python -i Code\test_me.py`

which will generate a figure (using the sub-routines in `plotting.py`) showing the results of the numerical integration on a specific test case (constant *k*, zero *S* - the reason for this is so that an analytic solution can be used for comparison - see the documentation for more details). It can be run again with different parameters by then using the command `main(L=4.0, nt=15)`, for example. Type `help(main)` to see which other parameters may be changed.

Also shown is a plot which demonstrates that the numerical scheme is conserving *q* - since Neumann boundary conditions are used and the scheme is a finite volume method, the total integral quantity of *q*(*x*) over the domain should remain constant as it is integrated forwards in time. The scheme has been verified to conserve *q* for (i) constant diffusivity *k* and (ii) position-dependent diffusivity *k*(*x*) (18/04/2018). Unfortunately this numerical scheme does not appear to conserve *q* when the source *S* is time dependent (although the solutions retain stability (and probably accuracy depending on the nature of the time dependency)); it is possible this may be due to the back-end matrix inversion methods (i.e. potentially solvable) but this is not being investigated further at this time.

## Dependencies
  * Python 2.7.14
  * MatPlotLib 2.2.2
  * NumPy 1.14.3
