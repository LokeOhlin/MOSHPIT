import CGS as cgs
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import argparse


muHe = 4.002602
muC  = 12.011
muO  = 15.9994
muSi = 28.0855

abundO = 3e-4
abundC = 1e-4
abundSi = 1.5e-5
abundHe = 0.1


mf_scale = 1.0 + abundC*muC
abar     = 1.0 + abundHe*muHe + abundC*muC + abundO*muO + abundSi*muSi

AvConversion = 6.289e-22

def plotFile(filename, axes):
    time = 0
    with open(filename, 'r') as fi:
        header = True
        while(header):
            line = fi.readline()
            if line[0] == '#':
                line = line.split(' = ')
                if(line[0] == '# time'):
                    time = float(line[1])
            else:
                header = False

    data = np.loadtxt(filename)

    # get xdistance and calulate dx
    xs = data[:,0]
    dx = np.zeros(xs.shape)
    xR = np.zeros(xs.shape)
    xL = np.zeros(xs.shape)
    xR[:-1] = (xs[:-1] + xs[1:]) * 0.5
    xL[1:] = xR[:-1]
    xR[-1] = xs[-1] + (xs[-1] - xL[-1])
    xL[0]  = 0

    dx = xR - xL
    
    density = data[:,1]
    pres    = data[:,3]
    
    xH0 = data[:,4]*mf_scale/density
    xH2 = data[:,5]*mf_scale/density/2.0
    xHp = data[:,6]*mf_scale/density
    xCO = data[:,7]*mf_scale/muC/density
    xCp = data[:,8]*mf_scale/muC/density
    ctot = xCO + xCp
    numH = density/(cgs.mH*abar)

    numH0 = xH0 * numH
    numH2 = xH2 * numH
    numHp = xHp * numH
    numCO = xCO * numH
    numCp = xCp * numH

    numTot = numH*(1 - xH2 + xHp + abundHe) 
    temp   = pres/(numTot*cgs.kb)
    print(np.min(ctot), np.max(ctot), np.mean(ctot))

    colmH0 = np.cumsum(numH0*dx)
    colmH2 = np.cumsum(numH2*dx)
    colmHp = np.cumsum(numHp*dx)
    colmCO = np.cumsum(numCO*dx)
    colmCp = np.cumsum(numCp*dx)
    colmHtot = colmH0+colmH2*2 + colmHp
    
    Av = colmHtot * AvConversion
    
    axes[0].plot(Av, temp)
    axes[1].plot(Av, colmH0/xR)
    axes[1].plot(Av, colmH2/xR)
    axes[2].plot(Av, colmCO/xR)
    axes[3].plot(Av, colmCp/xR)
    
adi = 5./3

idxs = [5, 10, 20]
dirs = ['V1', 'V2', 'V3', 'V4']
xlims = [[1e-4, 2e1], [1e-3,2e1], [1e-7, 1e-2], [1e-3, 2e1]]
lscale = [[False, True, True, True],
          [False, True, True, True],
          [False, True, True, True],
          [False, True, True, True]]
ylims =[[[0, 220], [1e-2, 1e4], [1e-9, 2], [1e-3, 2]],
        [[0, 3e3], [1e-7, 1e5], [1e-21,1e1], [1e-2, 2]],
        [[0, 9e2], [1e3, 1e6], [1e-2, 1e3], [1e-1, 1e3]],
        [[0, 1.25e4], [1e-1, 1e6], [1e-10, 1e3], [1e0, 1e3]]] 

for idx in idxs:
    print(idx)
    fig, axes = plt.subplots(4,4,sharex = 'col')
    for i, di in enumerate(dirs):
        filename = di+'/output_{:04d}'.format(idx)
        plotFile(filename, axes[:,i])
        axes[0,i].set_xlim(xlims[i])
        axes[3,i].set_xlabel(r'$A_v$')
        for j, ax in enumerate(axes[:,i]):
            ax.set_xscale('log')
            ax.set_ylim(ylims[i][j])
            if(lscale[i][j]):
                ax.set_yscale('log')
    axes[0,0].set_ylabel(r'$T$ [K]')
    axes[1,0].set_ylabel(r'$\Sigma_{H,H_2}/z$ [cm$^{-3}$]')
    axes[2,0].set_ylabel(r'$\Sigma_{CO}/z$ [cm$^{-3}$]')
    axes[3,0].set_ylabel(r'$\Sigma_{C_p}/z$[cm$^{-3}$]')
    plt.show()
