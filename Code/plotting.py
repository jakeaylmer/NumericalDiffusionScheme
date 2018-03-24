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


def MakePlots(x, initial, analytic, CrankNicolson, BackwardEuler, nt, dt, k0):
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
    ax.set_title(title, y=1.02)
    ax.set_xlabel(r'$x$ (m)')
    ax.set_ylabel(r'$q(x, t=%.1f$ s$)$' % (nt*dt) )
    ax.legend(loc=0, fontsize=17)
    
    return FormatAxis(fig, ax)


def PlotDefaults():
    """Set the global default formatting styles for some of the plot elements.
    This should be called before any other plotting functions.
    """
    mpl.rcParams['font.sans-serif'] = 'Calibri' #font for sans-serif style
    mpl.rcParams['font.family'] = 'sans-serif' #sans-serif font style
    mpl.rcParams['mathtext.fontset'] = 'custom' #allow customising maths fonts
    mpl.rcParams['mathtext.rm'] = 'sans' #maths roman font in sans-serif
    mpl.rcParams['mathtext.it'] = 'sans:italic' #maths italic font
    mpl.rcParams['mathtext.default'] = 'it' #maths in italic by default
    mpl.rcParams['axes.titlesize'] = 20 #default plot title font size
    mpl.rcParams['axes.labelsize'] = 18 #default axis label font size
    mpl.rcParams['lines.linewidth'] = 1.5 #default plot line width
    pass


def FormatAxis(fig, ax, gridon=True):
    """Set the layout and formatting of the plot on axis ax belonging to figure
    object fig.
    
    --Args--
    fig, ax    : MatPlotLib figure and axis objects respectively.
    (gridon)   : boolean; whether to set the grid on or not.
    """
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', direction='out') #outward ticks
    ax.tick_params(axis='both', which='major', labelsize=18, pad=8)
    if gridon:
        ax.grid(which='major', linestyle='-', color=[.75, .75, .75])
    ax.set_axisbelow(True)
    fig.tight_layout()
    return fig, ax
