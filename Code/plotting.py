### NumericalDiffusionScheme
### Jake Aylmer
###
### Some plotting functions to display results of test runs on the diffusion
### scheme.
### ---------------------------------------------------------------------------

from __future__ import division
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


def MakePlots(x, initial, analytic, CrankNicolson, BackwardEuler, nt, dt, k0,
    xlim=[0,1], ylim=[0,4]):
    """Plot the sample run shown in the documentation.
    
    --Args--
    x             : NumPy array, coordinates of grid-cell centres.
    initial       : NumPy array, initial conditions
    analytic      : NumPy array, analytic solution at time nt*dt.
    CrankNicolson : NumPy array, final q using Crank Nicolson scheme.
    BackwardEuler : NumPy array, final q using Backward Euler scheme.
    nt            : integer, number of time steps used.
    dt            : float, time step used (s).
    k0            : float, constant diffusivity used (m^2 s^-1).
    """
    fig, ax = plt.subplots()
    
    ax.plot(x, initial, color='k', linestyle='--',label='Initial')
    ax.plot(x, analytic, color='k', label='Analytic')
    ax.plot(x, CrankNicolson, color='g', label='Crank-Nicolson')
    ax.plot(x, BackwardEuler, color='b', label='Backward Euler')
    
    title = r'$N=%i$, $\Delta t=%.1f$ s, $n_t=%i$' % (len(x), dt, nt)
    title += r', $k_0 = %.1f\times 10^{-3}$ m$^2$s$^{-1}$' % (1000*k0)
    ax.set_title(title)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel(r'$x$ (m)')
    ax.set_ylabel(r'$q(x, t=%.1f$ s$)$' % (nt*dt) )
    ax.legend()
    fig.tight_layout()
    
    return fig, ax


def EnergyPlot(t, E_CN, E_BE, xlim=[0,1], ylim=[5,6]):
    """"""
    
    fig, ax = plt.subplots()
    
    ax.plot(t, E_CN, color='g', label='Crank-Nicolson')
    ax.plot(t, E_BE, color='b', label='Backward Euler')
    
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel(r'Time, $t$ (s)')
    ax.set_ylabel(r'$\int_{0}^{{\ }L} q(x,t)\mathrm{d}x$ ([q] m)')
    ax.legend()
    fig.tight_layout()
    
    return fig, ax


def SetRCParams():
    """Set default MatPlotLib formatting styles (rcParams) which will be set
    automatically for any plotting method.
    """
    # FONTS (NOTE: SOME OF THESE ARE SET-ORDER DEPENDENT):
    mpl.rcParams['font.sans-serif'] = 'Calibri' #Set font for sans-serif style
    mpl.rcParams['font.family'] = 'sans-serif' #Choose sans-serif font style
    mpl.rcParams['mathtext.fontset'] = 'custom' #Allow customising maths fonts
    mpl.rcParams['mathtext.rm'] = 'sans' #Maths roman font in sans-serif format
    mpl.rcParams['mathtext.it'] = 'sans:italic' #Maths italic font
    mpl.rcParams['mathtext.default'] = 'it' #Maths in italic by default
    
    # PLOT ELEMENT PROPERTIES:
    mpl.rcParams['lines.linewidth'] = 1.5 #Default plot linewidth (thickness)
    mpl.rcParams['lines.markersize'] = 4 #Default marker size (pts)
    mpl.rcParams['lines.markeredgewidth'] = 0 #Default marker edge width (pts)
    
    # LABEL PROPERTIES:
    mpl.rcParams['axes.titlesize'] = 20 #Title font size (pts)
    mpl.rcParams['axes.labelsize'] = 19 #Axis label font sizes (pts)
    mpl.rcParams['xtick.labelsize'] = 18 #X-tick label font size (pts)
    mpl.rcParams['ytick.labelsize'] = 18 #Y-tick label font size (pts)
    
    # GRID PROPERTIES:
    mpl.rcParams['axes.grid'] = True #Major grid on by default
    mpl.rcParams['grid.color'] = 'bfbfbf' #Grid line color
    mpl.rcParams['xtick.minor.visible'] = True #X-minor ticks on by default
    mpl.rcParams['ytick.minor.visible'] = True #Y-minor ticks on by default
    mpl.rcParams['xtick.major.pad'] = 8 #X-major tick padding
    mpl.rcParams['ytick.major.pad'] = 8 #Y-major tick padding
    mpl.rcParams['axes.axisbelow'] = True
    
    # LEGEND PROPERTIES:
    mpl.rcParams['legend.fancybox'] = False #Whether to use a rounded box
    mpl.rcParams['legend.fontsize'] = 16 #Legend label font size (pts)
    mpl.rcParams['legend.framealpha'] = 1 #Legend alpha (transparency)
    mpl.rcParams['legend.edgecolor'] = '#000000' #
    
    # GENERAL FIGURE PROPERTIES
    mpl.rcParams['figure.figsize'] = 8, 6 #Figure window size (inches)
    mpl.rcParams['savefig.format'] = 'pdf' #Default format to save to
    
    pass
