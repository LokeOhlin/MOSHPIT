import CGS as cgs
import numpy as np
import matplotlib.pyplot as plt
import argparse

adi = 5./3

ap=argparse.ArgumentParser()
ap.add_argument('f', nargs='+')
args=ap.parse_args()


def Nphots3(T, R, emin, emax):
    
    ll = emin*cgs.eV/(T*cgs.kb)
    ul = emax*cgs.eV/(T*cgs.kb)

    prefac2 = T**3*2*cgs.kb**3/(cgs.c**2*cgs.h**3)*np.pi
    its =1000
    tmp = 0
    
    ratell = 0
    rateul = 0

    for i in range(its):
        ii = i+1
        ratell += np.exp(-ii*ll)*(2/ii**3+ 2*ll/ii**2+ll**2/ii)
    if ul > 0 :
        for i in range(its):
            ii = i+1
            ratelu += np.exp(-ii*ul)*(2/ii**3+ 2*ul/ii**2+ul**2/ii)
    return prefac2*(ratell-rateul)*4*np.pi*R**2

def Stromgren(Nphot,alpha,numdens):
    return (3*Nphot/(4*np.pi*alpha))**(1/3)*numdens**(-2/3)

Lstar = 352316*cgs.Lsun
Teff  = 45258.5
Rstar = np.sqrt(Lstar/(4*np.pi*Teff**4*cgs.sigm_sb)) 

nH = 100

Nion = Nphots3(Teff, Rstar, 13.6, -1)
RST = Stromgren(Nion, 2.5e-13, nH)

print(RST/cgs.pc)


for f in args.f:
    print(f)
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
    data = np.loadtxt(f)

    rad = data[:,0]/cgs.pc
    dens = data[:,1]
    vel  = data[:,2]
    Pres = data[:,3]
    xH0  = data[:,4]/dens
    xH2  = data[:,5]/dens
    xHp  = data[:,6]/dens
    abar = 1 # We assume no other species here
    numd = dens/(abar*cgs.mH)*(1+xHp-xH2)
    T = Pres / (cgs.kb*numd)
    ionReg = xHp > 0.1
    if(sum(ionReg) == 0):
        continue
    rHII = np.max(rad[ionReg])
 
    fig = plt.figure(figsize = (15,9))
    ax1 = plt.subplot(151)
    ax1.plot(rad, dens, c = 'b')
    ax1.axvline(rHII, c = 'r', ls = ':')
    ax1.set_yscale('log')

    ax2 = plt.subplot(152)
    ax2.plot(rad, T)
    ax2.axvline(rHII, c = 'r', ls = ':')
    ax2.set_yscale('log')
    
    ax3 = plt.subplot(153)
    ax3.plot(rad, xH0)
    ax3.axvline(rHII, c = 'r', ls = ':')
    ax3.set_yscale('log')
    
    ax4 = plt.subplot(154)
    ax4.plot(rad, xH2)
    ax4.axvline(rHII, c = 'r', ls = ':')
    ax4.set_yscale('log')
    
    ax5 = plt.subplot(155)
    ax5.plot(rad, Pres)
    ax5.axvline(rHII, c = 'r', ls = ':')
    ax5.set_yscale('log')
    plt.savefig("plot_"+f[-4:]+'.png')
    plt.close()

