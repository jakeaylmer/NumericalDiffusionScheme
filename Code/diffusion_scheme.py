### NumericalDiffusionScheme
### Jake Aylmer
###
### This code contains two sub-routines needed to solve the diffusion equation
### with variable diffusivity and non-zero source term. It does not include
### advection terms (dq/dx). SolveDiffusionEquation() integrates forward by one
### time-step using SchemeMatrix() to calculate the diffusion operator, A.
### 
### See the repository documentation for further details.
### ---------------------------------------------------------------------------

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def SchemeMatrix(N, k, L=1.0):
    """Calculate the matrix A which expresses the diffusion operator in the
    numerical scheme for solving the diffusion equation with spatially-variable
    diffusivity k(x) and source term S:
    
        [d/dt - A]q(x,t) - S(x,t) = 0
    
    where q is the discretised prognostic variable, and Neumann boundary
    conditions (dq/dx = 0 at the boundaries 0, L) are used.
    
    See http://www.csc.kth.se/utbildning/kth/kurser/DN2255/ndiff13/Lecture3.pdf
    for details of how the matrix in this scheme is derived, or the repository
    documentation for a brief summary.
    
    --Args--
    N   : integer; number of grid cells.
    k   : function; of x, which should return the diffusivity at x.
    (L) : float, upper limit of spatial domain (i.e. 0 < x < L), default L=1.0.
    """
    h = L / N
    A = np.zeros( (N, N) )
    
    # Inner elements (away from boundaries) - all rows except first and last:
    for i in xrange(1, N-1):
        for j in xrange(0, N):
            A[i][j] = (i==(j+1)) * k((j+1)*h) + \
                      (i==(j-1)) * k(j*h) + \
                      ( i == j ) * -(k((j+1)*h) + k(j*h))
    
    # Neumann boundaries - first and last rows:
    A[0][0] = -k(h)
    A[0][1] = k(h)
    A[N-1][N-2] = k((N-1)*h)
    A[N-1][N-1] = -k((N-1)*h)
    
    A /= h**2 
    
    return A


def SolveDiffusionEquation(q_old, S_old, S_new, k, dt, L=1.0, theta=1.0):
    """Solves the diffusion equation with spatially-variable diffusivity and
    variable source term:
    
       dq/dt - d/dx[k(x)dq/dx] = S(x,t)
    
    using finite volume discretisation (method of lines). The x domain is
    defined on 0 < x < L, and is discretised into N grid cells of with h = L/N
    and q is given at q_j=q(x_j) where x_j = h/2 + j*h. Neumann boundary
    conditions are assumed. Returns the profile of q at the next time level.
    
    See http://www.csc.kth.se/utbildning/kth/kurser/DN2255/ndiff13/Lecture3.pdf
    for details of how this scheme is derived, or the repository documentation
    for a brief summary.
    
    --Args--
    q_old   : NumPy array of length N, q at the current time level.
    S_old   : NumPy array of length N, S(x) at the current time step.
    S_new   : NumPy array of length N, S(x) at the next time step.
    k       : function of x, which should return the diffusivity at x.
    dt      : float, time step.
    (L)     : float, upper limit of spatial domain (i.e. 0 < x < L), default
              L=1.0.
    (theta) : float, between 0 and 1, specifies which scheme is used (0 is
              forward-Euler, 0.5 is Crank-Nicholson, 1 is backward-Euler).
              Default theta=1.
    """
    
    N = len(q_old)
    
    A = SchemeMatrix(N, k, L)
    
    M1 = np.linalg.inv( np.eye(N) - theta*dt*A )
    M2 = np.dot(np.eye(N) + (1-theta)*dt*A, q_old) + \
        dt*(theta*S_new + (1-theta)*S_old)
    
    q_new = np.dot(M1, M2) 
    
    return q_new
