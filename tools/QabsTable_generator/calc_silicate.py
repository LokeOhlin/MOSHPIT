import numpy as np
import CGS as cgs
import matplotlib.pyplot as plt
from bhmie import *

def linInterp(xin, xtab, ytab):
    if xin > np.max(xtab):
        return 0
    if xin < np.min(xtab):
        return 0

    idm = np.abs(xin-xtab).argmin()
    if(xin < xtab[idm]):
        idm = idm -1
    idp = idm+1

    xm = xtab[idm]
    xp = xtab[idp]
    ym = ytab[idm]
    yp = ytab[idp]

    return ym + (xin - xm) * (yp - ym)/(xp - xm)



def delef(omega, a, T, Pars):
    omegaP  = Pars[0]
    tauBulk = Pars[1]
    Tf = Pars[2]
    vf = Pars[3]

    veff = vf * np.sqrt(1+T/Tf)

    tau = 1/(tauBulk**-1 + veff/a)


    return -(omegaP*tau)**2 /((omega*tau)**2 + 1j * omega*tau)
def getQabs(omegas, agrain, eps):
    vol = 4*np.pi*agrain**3 / 3
    
    Ceabs = omegas*vol/(cgs.c*(1/3)**2) * eps.imag/((eps.real+3-1)**2+eps.imag**2)
    
    ypar = np.sqrt(eps*(omegas*agrain/cgs.c)**2)
    am = - 0.5*agrain**3*(1+(3*np.tan(ypar)**-1)/ypar - 3/ypar**2)
    Cmabs = 4*np.pi*omegas/cgs.c * am.imag
    
    Qabs = (Ceabs + Cmabs)/(np.pi*agrain**2)
    return Qabs

SiDat = np.loadtxt('eps_Sil', skiprows=6 )
wavs = SiDat[:,0]*1e-4
omegas = 2*np.pi*cgs.c/wavs

eSi = (SiDat[:,1]+1) + 1j*SiDat[:,2]
nSi = np.sqrt((np.abs(eSi)+eSi.real)/2)
kSi = np.sqrt((np.abs(eSi)-eSi.real)/2)
refract = nSi + 1j*kSi

plt.plot(cgs.hbar*omegas/cgs.eV, eSi.imag)
plt.plot(cgs.hbar*omegas/cgs.eV, eSi.real)
plt.xscale('log')
plt.yscale('log')
plt.ylim(0.05,12)
plt.xlim(1e-3, 30)
plt.show()

#agrains = np.array([0.003, 0.01, 0.03, 0.1, 0.3 , 1])*1e-4
agrains = np.logspace(np.log10(3*cgs.A), np.log10(1e-3), 30)

freqs = np.logspace(np.log10(0.002), np.log10(1000),380)*cgs.eV/cgs.h
wavelengths = cgs.c/freqs
Qabs = np.zeros(wavelengths.shape)
angles = np.linspace(0,180,180)
with open('dustAborption_Silicone.txt', 'w') as f:
    for agrain in agrains:
        print(agrain)
        for i in range(len (wavelengths)):
            x = 2*np.pi*agrain/wavelengths[i]
            refPer = linInterp(wavelengths[i], wavs, refract)
            S1, S2, Qext, Qabs[i], Qsca, Qback, gsca = bhmie(x, refPer, angles)
            f.write('{} {} {} {} {}\n'.format(agrain, freqs[i], Qabs[i], Qsca, gsca))
        plt.plot(cgs.h*freqs/cgs.eV, Qabs)
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.01,1e4)
plt.ylim(1e-3, 10)
plt.show()




