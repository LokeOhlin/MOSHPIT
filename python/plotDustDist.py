import CGS as cgs
import numpy as np
import h5py as hp
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm
import argparse
ap=argparse.ArgumentParser()
ap.add_argument('f', nargs='+')
ap.add_argument('--linRadbins', action = 'store_true')
args=ap.parse_args()


cmap = matplotlib.cm.get_cmap('plasma')
for idx, f in enumerate(args.f):
    
    hfile = hp.File(f, 'r')
    print(f, hfile['Headers'].attrs.get('time')/cgs.yr) 
    rad = hfile['coordinates'][:]
    abin_c = hfile['agrain'][:]
    isilicone = hfile['Headers'].attrs.get('isilicone')
    num_c = hfile['number'][:,:isilicone]
    num_s = hfile['number'][:,isilicone:]

    idxs = [0, int(len(rad)/4), int(len(rad)/2), int(3*len(rad)/4), -1]
    if args.linRadbins:
        norm = Normalize(vmin = np.min(rad), vmax = np.max(rad))
    else:
        norm = LogNorm(vmin = np.min(rad), vmax = np.max(rad))
    fig, ax = plt.subplots(1,1)
    for i in range(len(idxs)):
        color = cmap(norm(rad[idxs[i]]))
        ax.plot(abin_c, num_c[idxs[i],:], c = color, label = r'$R$ = {:.2e} pc'.format(rad[idxs[i]]/cgs.pc))
        ax.plot(abin_c, num_s[idxs[i],:], c = color, ls = '--')
        
    
    ax.set_ylim(1e-15, 1e-5)
    ax.plot([],[], c='k', label = 'graphite')
    ax.plot([],[], c='k', ls = '--', label = 'silicate')
    ax.legend()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$a$ [cm]')
    ax.set_ylabel(r'$f(a)$')
    
    plt.show()
    
