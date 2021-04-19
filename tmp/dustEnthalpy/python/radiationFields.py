import CGS as cgs
import numpy as np
import matplotlib.pyplot as plt
def bbody_wl(wavelength, Teff):
    exp = cgs.h*cgs.c/(cgs.kb*Teff*wavelength)
    if(exp > 500):
        return 0
    return 2*cgs.h*cgs.c**2/(wavelength**5)/(np.exp(exp)-1)

def bbody_fr(frequency, Teff):
    exp = cgs.h*frequency/(cgs.kb*Teff)
    if(exp > 500):
        return 0
    return 2*cgs.h*frequency**3/(cgs.c**2)/(np.exp(exp)-1)
# 8 - 1000 micrometers from Mathis et al 1983. No formula found, so do linear interp
Math_L = np.array([8, 10, 12 , 20, 25, 30, 40, 50, 60 ,70 ,80, 90, 100, 150, 200, 300, 400, 600,800,1000])*1e-4
Math_E = np.array([4.33e-5, 3.26e-5, 2.21e-5, 7.12e-6, 1.14e-5, 1.7e-5, 3.06e-5, 4.08e-5, 4.76e-5, 4.51e-5, 4.2e-5, 3.62e-5, 3.16e-5, 1.36e-5, 5.61e-6, 1.39e-6, 5.61e-7, 1.32e-7, 4.2e-8, 2.14e-8])/(1e-4)

def interpUlam(wavelength):
    idx = np.abs(wavelength-Math_L).argmin()
    if(idx == len(Math_L)-1):
        return Math_E[-1]
    if (Math_L[idx]>wavelength):
        idx = idx - 1
    E0 = Math_E[idx]
    E1 = Math_E[idx+1]
    L0 = Math_L[idx]
    L1 = Math_L[idx+1]

    return E0 + (wavelength -L0)*(E1-E0)/(L1-L0)


def galacticRF(wavelength):
    Teffs = [2.9, 7500, 4000, 3000]
    Weights = [1, 1e-14, 1e-13, 4e-13]
    
    # Mathis 1983
    u_lam = 0
    if wavelength < 0.0912e-4:
        return 0
    elif wavelength < 0.11e-4:
        u_lam += 38.57 * (wavelength/1e-4)**3.4172 /(1e-4)
    elif wavelength < 0.134e-4:
        u_lam += 2.045e-2/1e-4
    elif wavelength < 0.246e-4:
        u_lam += 7.115e-4 * (wavelength/1e-4)**(-1.6678)  / (1e-4)
    elif(wavelength < 8e-4):
        for i in range(1,4):
            u_lam = u_lam + Weights[i] * 4 * np.pi * bbody_wl(wavelength,Teffs[i])
    elif(wavelength < 0.1):
        u_lam += interpUlam(wavelength)
        
    
    u_lam = u_lam + 4* np.pi*bbody_wl(wavelength,Teffs[0])
    u_lam = u_lam/cgs.c
    return u_lam





def testGalaxy():
    lams = np.logspace(np.log10(700*cgs.A), -1,1000)
    edens = np.zeros(lams.shape)
    for i in range(len(lams)):
        edens[i] = galacticRF(lams[i])*lams[i]
    plt.plot(lams, edens)

    for i in range(len(lams)):
        edens[i] = lams[i]*4*np.pi*bbody_wl(lams[i], 2.9)/cgs.c
    plt.plot(lams, edens)
    
    for i in range(len(lams)):
        edens[i] = lams[i]*(1e-14)*4*np.pi*bbody_wl(lams[i], 7500)/cgs.c
    plt.plot(lams, edens)
    
    for i in range(len(lams)):
        edens[i] = lams[i]*(1e-13)*4*np.pi*bbody_wl(lams[i], 4000)/cgs.c
    plt.plot(lams, edens)
    
    for i in range(len(lams)):
        edens[i] = lams[i]*(4e-13)*4*np.pi*bbody_wl(lams[i], 3000)/cgs.c
    plt.plot(lams, edens)
   
    edens = edens/(1e-4*cgs.c)
    plt.plot(lams, edens*lams)

    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(500*cgs.A, 1e-3)
    plt.ylim(1e-15, 3e-10)
    plt.show()

if __name__ == '__main__':
    testGalaxy()
