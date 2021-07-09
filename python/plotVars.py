import numpy as np
import h5py as hp
import matplotlib.pyplot as plt
import argparse
import cgs

ap=argparse.ArgumentParser()
#---------------outputs-----------------------------
ap.add_argument('f', nargs='+', 
        help = 'HDF5 output to plot')
#---------------Variables---------------------------
ap.add_argument('--vars', nargs = '+', 
        help = 'Variables to plot')
ap.add_argument('--ncols', type = int, default = 3, 
        help = 'Number of columns in plot grid (ignored not enough vaiable)')
args = ap.parse_args()

def plotVar(coords, var, ax, xlim = None, ylim = None, label = None, xlabel = None, ylabel = None, color = 'navy', ls = '-'):
    ax.plot(coords, var, c = color, ls = ls, label = label)
    if not xlim is None:
        ax.set_xlim(xlim)
    if not ylim is None:
        ax.set_ylim(ylim)
    if not xlabel is None:
        ax.set_xlabel(xlabel)
    if not ylabel is None:
        ax.set_ylabel(ylabel)

numVars = len(args.vars)
ncols = min(numVars, args.ncols)
nrows = int(np.ceil(numVars/ncols))
fscale = 3
rfig = 1.0
rw   = 0.1
rh   = 0.1
fsize = (ncols * fscale * (rfig + rw), nrows * fscale * (rfig +rh)) 
for f in args.f:
    fig, axes = plt.subplots(nrows = nrows, ncols = ncols, figsize = fsize, sharex = True)
    axes = axes.flatten()
    for i, varname in enumerate(args.vars):
        h5file = hp.File(f, 'r')
        print([varname]) 
        coords = h5file['coordinates'][:]/cgs.pc
        # determine if special variable or synonym name, otherwise just grab straight from file
        if varname == 'numd':
            abar = h5file['/Headers'].attrs.get('abar')*cgs.mH
            var  = h5file['density'][:]/abar                 
        elif varname == 'dens':
            var  = h5file['density'][:]
        elif varname == 'temp':
            var  = h5file['temperature'][:]
        elif varname == 'velx':
            var  = h5file['velocity'][:]
        elif varname == 'xHI':
            var  = h5file['chemicalAbundances'][:,0]
        elif varname == 'xH2':
            var  = h5file['chemicalAbundances'][:,1]
        elif varname == 'xHp':
            var  = h5file['chemicalAbundances'][:,2]
        elif varname == 'xCO':
            var  = h5file['chemicalAbundances'][:,3]
        elif varname == 'xCp':
            var  = h5file['chemicalAbundances'][:,4]
        else:
            var  = h5file[varname]

        if( i >= numVars - ncols ):
            plotVar(coords, var, axes[i], xlabel = r'$x$ [pc]', ylabel = varname)
        else:
            plotVar(coords, var, axes[i], ylabel = varname)
    
    plt.subplots_adjust(bottom=0.15,top=0.95,left=0.05,right=0.95,hspace=0.1,wspace=0.3)
    plt.show()
