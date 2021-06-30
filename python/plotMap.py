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
    data = hp.File(f,'r')
    fig = plt.figure()
    ax1 = plt.subplot(231)
    coords = data['coordinates'][:]
    dens   = data['density'][:]
    temp   = data['temperature'][:]
    pres   = data['pressure'][:]
    velx   = data['velocity'][:]
    ax1.plot(coords, dens)
    if(args.log):
        ax1.set_yscale('log')
        ax1.set_xscale('log')
#   ax1.set_ylim(0,1.1e5)
    
    ax1 = plt.subplot(232)
    ax1.plot(coords, velx)
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
    plt.savefig("hydro"+f[-4:]+'.png')
    plt.close(fig)

