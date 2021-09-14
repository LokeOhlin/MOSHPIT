import numpy
import numpy as np
import matplotlib as mpl
import matplotlib.colors as mcol
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse
import sys
import CGS as cgs
import h5py as hp
ap=argparse.ArgumentParser()

#---------------outputs-----------------------------
ap.add_argument('f', nargs='+')
ap.add_argument('--log', action = 'store_true')
args=ap.parse_args()

for f in args.f:
    print(f)
    data = hp.File(f,'r')
    fig = plt.figure()
    ax1 = plt.subplot(231)
    coords = data['coordinates'][:]
    dens   = data['density'][:]
    temp   = data['temperature'][:]
    pres   = data['pressure'][:]
    velx   = data['velocity'][:]
    
    if('chemicalAbundances' in data.keys()):
        xH0 = data['chemicalAbundances'][:,0]
        xH2 = data['chemicalAbundances'][:,1]
        xHp = data['chemicalAbundances'][:,2]
        xCO = data['chemicalAbundances'][:,3]
        xCp = data['chemicalAbundances'][:,4]
    
    
    ax1.plot(coords, dens)
    if(args.log):
        ax1.set_yscale('log')
        ax1.set_xscale('log')
#   ax1.set_ylim(0,1.1e5)
    
    ax1 = plt.subplot(232)
    ax1.plot(coords, velx/1e5)
    if(args.log):
        ax1.set_xscale('log')
    #if(args.log):
    #    ax1.set_yscale('log')
#    ax1.set_ylim(0,None)
    
    ax1 = plt.subplot(233)
    #ax1.plot(data[:,0], data[:,3])
    #print(np.max(data[:,3]), np.min(data[:,3]))
    ax1.plot(coords, temp)
    if(args.log):
        ax1.set_yscale('log')
        ax1.set_xscale('log')
    
    ax1 = plt.subplot(234)
    ax1.plot(coords, pres)
    if(args.log):
        ax1.set_yscale('log')
        ax1.set_xscale('log')
    

    if('chemicalAbundances' in data.keys()):
        ax1 = plt.subplot(235)
        ax1.plot(coords, xH0, ls = '-', c= 'k' , label = r"$x_{H}$")
        ax1.plot(coords, xH2, ls = '--', c= 'k', label = r"$x_{H_2}$")
        ax1.plot(coords, xHp, ls = ':', c= 'k' , label = r"$x_{H^+}$")
        ax1.legend()
        if(args.log):
            ax1.set_yscale('log')
            ax1.set_ylim(1e-5, 1.5)
            ax1.set_xscale('log')
        
        ax1 = plt.subplot(236)
        ax1.plot(coords, xCO, ls = '-', c= 'k' , label = r"$x_{CO}$")
        ax1.plot(coords, xCp, ls = '--', c= 'k', label = r"$x_{C^+}$")
        ax1.legend() 
        if(args.log):
            ax1.set_yscale('log')
            ax1.set_ylim(1e-6, 5e-4)
            ax1.set_xscale('log')
    
    plt.savefig("hydro"+f[-4:]+'.png')
    plt.close(fig)
