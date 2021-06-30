import numpy
import numpy as np
import matplotlib as mpl
import matplotlib.colors as mcol
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse
import sys
import CGS as cgs
ap=argparse.ArgumentParser()

#---------------outputs-----------------------------
ap.add_argument('f', nargs='+')
args=ap.parse_args()

scale_l = 1#cgs.pc
scale_t = 1#cgs.yr*1e6
scale_m = 1#cgs.Msun

scale_v = scale_l/scale_t
scale_e = scale_m*scale_v**2
scale_d = scale_m/scale_l**3


n0   = 100
Etot = 1e51/scale_e
xi0  = 1.1515

dens0 = n0*1.2*cgs.mH/scale_d


#  Load semi analytic solution
ST_semi = np.loadtxt('ST_semiAnalytic.dat')

for f in args.f:
    time = 0
    with open(f, 'r') as fi:
        header = True
        while(header):
            line = fi.readline()
            if line[0] == '#':
                line = line.split(' = ')
                if(line[0] == '# time'):
                    time = float(line[1])
            else:
                header = False
    
    # Determine scaling factors
    U_sh = 2./5. * xi0 *(Etot/(dens0*time**3))**(1./5.)
    r_sh = xi0   * (Etot*time**2/dens0)**(1./5.)

    Velnorm = 2*U_sh/(5./3.+1)
    Prenorm = 2*dens0*U_sh**2/(5./3.+1)

    data = np.loadtxt(f)
    fig = plt.figure(figsize = (15,6))
    ax1 = plt.subplot(131)
    ax1.plot(data[:,0], data[:,1], c = 'b')
    ax1.plot(ST_semi[:,0]*r_sh, ST_semi[:,1]*dens0*4, c = 'r')
    ax1.set_ylabel(r"$\rho$")

    ax2 = plt.subplot(132)
    ax2.plot(data[:,0], data[:,2])
    ax2.plot(ST_semi[:,0]*r_sh, ST_semi[:,2]*Velnorm, c = 'r')
    ax2.set_ylabel(r"$v_r$")
    ax2.set_xlabel(r"$r$ [pc]")
    ax3 = plt.subplot(133)
    ax3.plot(data[:,0], data[:,3])
    ax3.plot(ST_semi[:,0]*r_sh, ST_semi[:,3]*Prenorm, c = 'r')
    ax3.set_ylabel(r"$P$")

    plt.savefig("stCompare_"+f[-4:]+'.png')
    plt.close(fig)

