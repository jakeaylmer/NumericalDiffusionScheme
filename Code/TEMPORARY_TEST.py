### TEMPORARY_TEST
###
### Testing for the case with:
###
### dq/dt - d/dx[c1*(1-x^2)d/dx]q = S(x,t)
###
### ---------------------------------------------------------------------------

from __future__ import division
import matplotlib as mpl, matplotlib.pyplot as plt
import numpy as np
import diffusion_scheme as dif, plotting as pl


def S(x, t, S0=341.75):
    #return S0*(1.229 - 0.542*np.cos(2*np.pi*t) - 3*0.229*x**2)
    return S0*(1.229 - 0.0*np.cos(2*np.pi*t) - 3*0.229*x**2)


def Diffusivity(x, c1=0.0684):
    return c1*(1-x**2)


def main(L=1.0, N=30, t_tot=5, nt=1000):
    
    h = L/N # width of grid cells in m
    dt = t_tot/nt # time step in yr

    x = np.linspace(h/2, L-h/2, N) # grid-cell centres (m)
    q_init = 30*x*(L - x) # initial conditions
    q_old_CN = q_init.copy(); q_old_BE = q_init.copy()
    
    diag_q_CN = np.zeros(nt); diag_q_BE = np.zeros(nt)
    
    for i in xrange(nt):
        
        print "Calculating %i of %i...\r" % (i, nt),
        t = i*dt # current time (at "old")
        
        q_new_CN = dif.SolveDiffusionEquation(q_old_CN, S(x, t), S(x, t+dt),
            Diffusivity, dt, L, 0.5) # Crank-Nicolson (CN)
        q_new_BE = dif.SolveDiffusionEquation(q_old_BE, S(x, t), S(x, t+dt),
            Diffusivity, dt, L, 1.0) # Backward-Euler (BE)
        
        ddt_int_q_CN = ( h*np.sum(q_new_CN) - h*np.sum(q_old_CN) ) / dt
        ddt_int_q_BE = ( h*np.sum(q_new_BE) - h*np.sum(q_old_BE) ) / dt
        diag_q_CN[i] = ddt_int_q_CN - h*np.sum(S(x,t))
        diag_q_BE[i] = ddt_int_q_BE - h*np.sum(S(x,t))
        
        q_old_CN = q_new_CN.copy()
        q_old_BE = q_new_BE.copy()
    
    fig1, ax1 = plt.subplots()
    ax1.plot(x, q_new_BE, color='b', label=r'Backward-Euler')
    ax1.plot(x, q_new_CN, color='g', label=r'Crank-Nicolson')
    ax1.legend(loc='best')
    fig1, ax1= pl.FormatAxis(fig1, ax1)
    
    fig2, ax2 = plt.subplots()
    ax2.plot(np.linspace(0, t_tot, nt), diag_q_BE, color='b', label=r'Backward-Euler')
    ax2.plot(np.linspace(0, t_tot, nt), diag_q_CN, color='g', label=r'Crank-Nicolson')
    ax2.legend(loc='best')
    fig2, ax2 = pl.FormatAxis(fig2, ax2)
    
    fig1.show()
    fig2.show()
    pass


if __name__ == '__main__':
    pl.PlotDefaults()
    main()
