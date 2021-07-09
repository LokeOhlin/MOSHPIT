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
    
    rad = hfile['coordinates'][:]
    ephots = hfile['photonBins'][:][:]
    isilicone = hfile['Headers'].attrs.get('isilicone')
    opacity = np.sum(hfile['opticalDepth'][:,:,:], axis = 1)

    idxs = [0, int(len(rad)/4), int(len(rad)/2), int(3*len(rad)/4), -1]
    if args.linRadbins:
        norm = Normalize(vmin = np.min(rad), vmax = np.max(rad))
    else:
        norm = LogNorm(vmin = np.min(rad), vmax = np.max(rad))
    fig, ax = plt.subplots(1,1)
    for i in range(len(idxs)):
        color = cmap(norm(rad[idxs[i]]))
        ax.plot(ephots/cgs.eV, opacity[idxs[i],:], c = color, label = r'$R$ = {:.2e} pc'.format(rad[idxs[i]]/cgs.pc))
        
    
    ax.legend()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$E_\gamma$ [eV]')
    ax.set_ylabel(r'$\tau_\mathrm{dust}/\mathrm{d}r$ [cm$^{-1}$]')
    
    plt.show()
    
