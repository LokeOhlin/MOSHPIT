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


def numdist(a, Ni, Si, da, ac):
    return Ni/da + Si * (a-ac)
def mass_to_number(Mi, Si, ac, ae, aep):
    NfactP = (aep**4)/(4*(aep-ae))
    Nfact  = (ae**4 )/(4*(aep-ae))
    SfactP = (aep**5)/5. -ac*(aep**4)/4.
    Sfact  = (ae**5 )/5. -ac*(ae**4 )/4.

    number = (Mi - Si*(SfactP - Sfact))/(NfactP - Nfact) 
    number[Mi <= 0] = 0
    return number

mnorm_c = 4*np.pi/3 * 3 * 2.2
mnorm_s = 4*np.pi/3 * 3 * 3.5 

cmin = 1.3335214321633213e-06
cmax = 7.498942093324532e-06
dadt = 1.5e-6

cmap = matplotlib.cm.get_cmap('plasma')
for idx, f in enumerate(args.f):
    
    hfile = hp.File(f, 'r')
    
    rad = hfile['radius'][:]
    abin_c = hfile['agrain'][:]
    isilicone = hfile['Header'].attrs.get('isilicone')
    num_c = hfile['number'][:,:isilicone]
    num_s = hfile['number'][:,isilicone:]
    dadt_c = hfile['dadt'][:, :isilicone]
    dadt_s = hfile['dadt'][:, isilicone:]
    idxs = [0, int(len(rad)/4), int(len(rad)/2), int(3*len(rad)/4), -1]
    if args.linRadbins:
        norm = Normalize(vmin = np.min(rad), vmax = np.max(rad))
    else:
        norm = LogNorm(vmin = np.min(rad), vmax = np.max(rad))
    fig, ax = plt.subplots(1,1)
    mint = 1e13
    for i in range(len(idxs)):
        color = cmap(norm(rad[idxs[i]]))

        ij = np.where(np.abs(dadt_c[idxs[i],:])>0)[0]
        timescale = -abin_c[ij]/dadt_c[idxs[i],ij]/cgs.yr
        mint = min(mint, np.min(timescale[timescale > 0])) 
        ax.plot(abin_c[ij], timescale, c = color, label = r'$R$ = {:.2e} pc'.format(rad[idxs[i]]/cgs.pc))
        
        
        ij = np.where(np.abs(dadt_s[idxs[i],:])>0)[0]
        timescale = -abin_c[ij]/dadt_s[idxs[i],ij]/cgs.yr
        mint = min(mint, np.min(timescale[timescale > 0])) 
        ax.plot(abin_c[ij], timescale, c = color, ls = '--')

        
    
    ax.plot([],[], c='k', label = 'graphite')
    ax.plot([],[], c='k', ls = '--', label = 'silicate')
    ax.legend()
    ax.set_ylim(mint, 1e11)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$a$ [cm]')
    ax.set_ylabel(r'$\tau_\mathrm{evap} [yr]$')
    
    plt.show()
    
