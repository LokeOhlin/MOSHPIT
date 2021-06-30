import CGS as cgs
import numpy as np
import matplotlib.pyplot as plt
import argparse
ap=argparse.ArgumentParser()
ap.add_argument('f', nargs='+')
args=ap.parse_args()


def numdist(a, Ni, Si, da, ac):
    return Ni/da + Si * (a-ac)
def mass_to_number(Mi, Si, ac, ae, aep):
    NfactP = (aep**4)/(4*(aep-ae))
    Nfact  = (ae**4 )/(4*(aep-ae))
    SfactP = (aep**5)/5. -ac*(aep**4)/4.
    Sfact  = (ae**5 )/5. -ac*(ae**4 )/4.

    return (Mi - Si*(SfactP - Sfact))/(NfactP - Nfact) 



mnorm_c = 4*np.pi/3 * 3 * 2.2
mnorm_s = 4*np.pi/3 * 3 * 3.5 

cmin = 1.3335214321633213e-06
cmax = 7.498942093324532e-06
dadt = 1.5e-6


for idx, f in enumerate(args.f):
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
   
    rho = data[:,1]
    idx_mass_c = np.arange(0,2*nabins,2) + idust_start + 1
    idx_slope_c = np.arange(0,2*nabins,2) + idust_start + 2
    mass_c  = data[:,idx_mass_c]
    slope_c = data[:,idx_slope_c]
    num_c = mass_to_number(mass_c/mnorm_c, slope_c, abin_c, abin_e[:-1], abin_e[1:])
    aveSize = np.sum(num_c*abin_c, axis = 1)/np.sum(num_c, axis = 1)
    plt.plot(data[:,0], aveSize)

    idx_mass_s = np.arange(0,2*nabins,2) + idust_start + 2*isilicone + 1
    idx_slope_s = np.arange(0,2*nabins,2) + idust_start + 2*isilicone + 2
    mass_s  = data[:,idx_mass_s]
    slope_s = data[:,idx_slope_s]
    num_s = mass_to_number(mass_s/mnorm_s, slope_s, abin_c, abin_e[:-1], abin_e[1:])
  
    aveSize = np.sum(num_s*abin_c, axis = 1)/np.sum(num_s, axis = 1)
    plt.plot(data[:,0], aveSize)
    #plt.xscale('log')
    plt.yscale('log')
    plt.show()
    
    
