import CGS as cgs
import numpy as np
import matplotlib.pyplot as plt
import argparse

def numdist(a, Ni, Si, da, ac):
    return Ni/da + Si * (a-ac)
def mass_to_number(Mi, Si, ac, ae, aep):
    NfactP = (aep**4)/(4*(aep-ae))
    Nfact  = (ae**4 )/(4*(aep-ae))
    SfactP = (aep**5)/5. -ac*(aep**4)/4.
    Sfact  = (ae**5 )/5. -ac*(ae**4 )/4.

    return (Mi - Si*(SfactP - Sfact))/(NfactP - Nfact) 


def distribution(a, time, amin, amax, plaw, dadt, rho_grain, mtot):
    norm = (plaw + 1)/(amax**(plaw+1) - amin**(plaw+1))
    mfirst = 4*np.pi *  rho_grain/3 * norm * (amax**(plaw+4)-amin**(plaw+4))/(plaw+4)
    print(norm, mfirst)
    amax = amax - dadt * (np.cos(2*np.pi*time)-1)/(2*np.pi)
    amin = amin - dadt * (np.cos(2*np.pi*time)-1)/(2*np.pi)
    norm = (plaw + 1)/(amax**(plaw+1) - amin**(plaw+1))
    norm = norm * mtot / mfirst

    dist = np.zeros(a.shape)
    dist[ (a<amin) ] = 0
    dist[ (a>amax) ] = 0
    mask = (a>amin)*(a<amax) 
    dist[mask] = norm*a[mask]**plaw
    return dist
    

mnorm = 4*np.pi/3 * 3
 
cmin = 1.3335214321633213e-06
cmax = 7.498942093324532e-06
dadt = 1.5e-6

fig, axes = plt.subplots(figsize = (11,1.4*5), nrows = 5)

files = ['output_0000', 'output_0001','output_0002', 'output_0004', 'output_0016']

for idx, f in enumerate(files):
    print(f)
    time = 0
    
    amin = 0
    amax = 0
    nabins = 0
    isilicone = 0
    idust_start = 0
    nvar = 0
    fSi  = 0
    with open(f, 'r') as fi:
        header = True
        while(header):
            line = fi.readline()
            if line[0] == '#':
                line = line.split(' = ')
                if(line[0] == '# time'):
                    time = float(line[1])
                if(line[0] == '# amin'):
                    amin = float(line[1])
                if(line[0] == '# amax'):
                    amax = float(line[1])
                if(line[0] == '# fSi'):
                    fSi = float(line[1])
                if(line[0] == '# Nabins'):
                    nabins = int(line[1])
                if(line[0] == '# nvar'):
                    nvar = int(line[1])
                if(line[0] == '# isilicone'):
                    isilicone = int(line[1])
                if(line[0] == '# IDUST_START'):
                    idust_start = int(line[1])
                   
            else:
                header = False

    abin_e = np.zeros(nabins+1)
    abin_c = np.zeros(nabins)
    da  = (np.log(amax) - np.log(amin))/nabins
    abin_e[0] = amin
    for i in range(nabins):
        abin_e[i+1] = np.exp(np.log(abin_e[i])+da)
    for i in range(nabins):
        abin_c[i] = (abin_e[i+1] + abin_e[i])/2.
    data = np.loadtxt(f)
   
    rho = data[1]
    print("rho", rho) 
    idx_mass_c = np.arange(0,2*nabins,2) + idust_start  + 1
    idx_slope_c = np.arange(0,2*nabins,2) + idust_start + 2
    mass_c  = data[idx_mass_c]
    slope_c = data[idx_slope_c]
    num_c = mass_to_number(mass_c/mnorm, slope_c, abin_c, abin_e[:-1], abin_e[1:])
    num_c[num_c<1e-20] = 0 
    idx_mass_s  = np.arange(0,2*nabins,2) + idust_start + 2*isilicone + 1
    idx_slope_s = np.arange(0,2*nabins,2) + idust_start + 2*isilicone + 2
    mass_s  = data[idx_mass_s]
    slope_s = data[idx_slope_s]
    num_s = mass_to_number(mass_s/mnorm, slope_s, abin_c, abin_e[:-1], abin_e[1:])
   
    
    aplot = np.zeros(3)
    cplot = np.zeros(3)
    splot = np.zeros(3)
    pmax = 0
    pmin = 1e99
    for i in range(nabins):
        aplot[0] = abin_e[i]
        aplot[1] = abin_c[i]
        aplot[2] = abin_e[i+1]
        if num_c[i] <= 0:
            continue
        
        cplot = numdist(aplot, num_c[i], slope_c[i], abin_e[i+1]-abin_e[i], abin_c[i]) 
            
        splot = numdist(aplot, num_s[i], slope_s[i], abin_e[i+1]-abin_e[i], abin_c[i]) 
        pmax = max(pmax, np.max(cplot))
        pmin = min(pmin, np.min(cplot))
        axes[idx].plot(aplot/1e-4, cplot , c='red')
        axes[idx].plot(aplot[1]/1e-4, cplot[1], ls = '', marker = 'x', c='red')
        
        axes[idx].plot(aplot/1e-4, splot , c='blue')
        axes[idx].plot(aplot[1]/1e-4, splot[1], ls = '', marker = 'x', c='blue')
    
    aplot = np.logspace(np.log10(amin),np.log10(amax),1000)
    axes[idx].plot(aplot/1e-4, distribution(aplot, time, cmin, cmax, 0, dadt, 3.0, rho*(1-fSi)), c = 'red', alpha = 0.4, ls = '--')
    axes[idx].plot(aplot/1e-4, distribution(aplot, time, cmin, cmax, 0, dadt, 3.0, rho*fSi), c = 'blue', alpha = 0.4, ls = '--')
    axes[idx].set_xscale("log")
    axes[idx].set_xlim(5e-3,0.2)
    #plt.yscale("log")
    axes[idx].set_ylim(0, 0.3)
    axes[idx].text(0.1,0.2,'$t = {0:.2f}T$'.format(time), size = 12)

axes[-1].set_xlabel(r'$a$ [$\mu$m]')
axes[2].set_ylabel(r'$\mathrm{d}n/\mathrm{d}a$ [cm$^{-3}\,$cm$^{-1}$]')
axes[0].plot([],[],c='red',label = 'Species 1')
axes[0].plot([],[],c='blue',label = 'Species 2')
axes[0].plot([],[],c='k',marker = 'x',label = 'Numerical')
axes[0].plot([],[],c='k', alpha = 0.4, label = 'Analytic')
axes[0].legend(loc='upper left')
plt.savefig('dust_singrowth.png')
plt.close(fig)
