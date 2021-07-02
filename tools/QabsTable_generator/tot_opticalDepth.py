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

PerDat = np.loadtxt('callindex.out_CpeD03_0.01', skiprows = 5)
ParDat = np.loadtxt('callindex.out_CpaD03_0.01', skiprows = 5)
wavsPer = PerDat[:,0]*1e-4
omegasPer = 2*np.pi*cgs.c/wavsPer

e1Per    = PerDat[:,1]+1
e2Per    = PerDat[:,2]
ePer001  = e1Per+1j*e2Per 

Temp = 20
omegaP  = 4.33*np.sqrt(1-6.24e-3*Temp + 3.66e-5*Temp*Temp)
tauBulk = 4.2e-11/(1+0.322*Temp + 0.00130*Temp*Temp)
Teff = 255
vf = 4.5e7
ParsPer = [omegaP, tauBulk, Teff, vf]

delefPer = delef(omegasPer, 0.01*1e-4, Temp, ParsPer)
ePer = ePer001-delefPer

wavsPar = ParDat[:,0]*1e-4
omegasPar = 2*np.pi*cgs.c/wavsPar

e1Par    = ParDat[:,1]+1
e2Par    = ParDat[:,2]
ePar001  = e1Par+1j*e2Par 

Temp = 20
omegaP  = 1.53
tauBulk = 1.4e-14
Teff = 255
vf = 3.7e6
ParsPar = [omegaP, tauBulk, Teff, vf]

delefPar = delef(omegasPar, 0.01*1e-4, Temp, ParsPar)
ePar = ePar001-delefPar

SiDat = np.loadtxt('eps_Sil', skiprows=6 )
wavsSi = SiDat[:,0]*1e-4
omegaSi = 2*np.pi*cgs.c/wavsSi

eSi = (SiDat[:,1]+1) + 1j*SiDat[:,2]
nSi = np.sqrt((np.abs(eSi)+eSi.real)/2)
kSi = np.sqrt((np.abs(eSi)-eSi.real)/2)
refractSi = nSi + 1j*kSi

agrains = np.logspace(np.log10(5e-7), np.log10(2.5e-5), 200)
dagrains = np.zeros(agrains.shape)
for ia in range(len(agrains)-1):
    amid = (agrains[ia]+agrains[ia+1])/2
    dagrains[ia]   += amid
    dagrains[ia+1] -= amid

dagrains[0] -= agrains[0]
print(dagrains[-1], agrains[-1])
dagrains[-1] += agrains[-1]
Agraph = 10**(-25.16)
num_graph = Agraph * agrains**-3.5 * dagrains
Asili = 10**(-25.11)
num_sili  = Asili  * agrains**-3.5 * dagrains


freqs = np.logspace(-2, 3, 500)*cgs.eV/cgs.h
wavelengths = cgs.c/freqs
tauExt_tot = np.zeros(freqs.shape)
tauExt_gra = np.zeros(freqs.shape)
tauExt_sil = np.zeros(freqs.shape)
angles = np.linspace(0,180,180)
for i, wavs in enumerate(wavelengths):
    for ia, agrain in enumerate(agrains):
        x = 2*np.pi*agrain/wavs
        
        #graphite
        delefPer = delef(omegasPer, agrain, Temp, ParsPer)
        ePer_grain = ePer + delefPer
        nPer = np.sqrt((np.abs(ePer_grain)+ePer_grain.real)/2)
        kPer = np.sqrt((np.abs(ePer_grain)-ePer_grain.real)/2)
        refractPer = nPer + 1j*kPer
        
        delefPar = delef(omegasPar, agrain, Temp, ParsPar)
        ePar_grain = ePar + delefPar
        nPar = np.sqrt((np.abs(ePar_grain)+ePar_grain.real)/2)
        kPar = np.sqrt((np.abs(ePar_grain)-ePar_grain.real)/2)
        refractPar = nPar + 1j*kPar
        
        
        refPer = linInterp(wavs, wavsPer, refractPer)
        S1, S2, QextPer, QabsPer, QscatPer, Qback, gscaPer = bhmie(x, refPer, angles)
        
        refPar = linInterp(wavs, wavsPar, refractPar)
        S1, S2, QextPar, QabsPar, QscatPar, Qback, gscaPar = bhmie(x, refPar, angles)
        Qabs_graph = (2*QabsPer + QabsPar)/3
        Qsca_graph = (2*QabsPer + QabsPar)/3
        gsca_graph = (2*gscaPer + gscaPar)/3
        
        Qext_graph = (2*QextPer + QextPar)/3
        
        #silicaone
        refSi = linInterp(wavs, wavsSi, refractSi)
        S1, S2, Qext_sili, Qabs_sili, Qsca_sili, Qback_sili, gsca = bhmie(x, refSi, angles)
        
        tauExt_sil[i] += np.pi*agrain**2 * Qext_sili*num_sili[ia]
        tauExt_gra[i] += np.pi*agrain**2 * Qext_graph*num_graph[ia] 
        tauExt_tot[i] += np.pi*agrain**2 * (Qext_graph*num_graph[ia] +  Qext_sili*num_sili[ia])


np.savetxt("averageExtinction.txt", [freqs, tauExt_tot], delimiter = ', ')
plt.plot(freqs/cgs.c * 1e-4, tauExt_tot/1e-21)
plt.plot(freqs/cgs.c * 1e-4, tauExt_gra/1e-21)
plt.plot(freqs/cgs.c * 1e-4, tauExt_sil/1e-21)
plt.xscale('log')
plt.yscale('log')
#plt.xlim(0, 10)
#plt.ylim(0,2.5)
plt.show()
