import numpy as np
import matplotlib.pyplot as plt
import CGS as cgs




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

Mstar=np.array([10, 20, 50, 100])
Lstar = np.array([5572, 40728.5, 352316,1.24e6])*cgs.Lsun
Teff  = np.array([25320.5, 34850, 45258.5,50893])

print(Mstar)
print(Lstar/cgs.Lsun)
print(Teff)

nH = 100
Rstar = np.sqrt(Lstar/(4*np.pi*Teff**4*cgs.sigm_sb)) 
nphot3 = np.array([Nphots3(Teff[i],Rstar[i],13.6,-1) for i in range(len(Rstar))])
Rs3 = Stromgren(nphot3,2.59e-13,nH)

print('numd = ', nH)
print('trec = {:.4e}'.format( 1/(2.59e-13*nH)))
print(Rs3)
print(10*Rs3)

nH = 1000
Rstar = np.sqrt(Lstar/(4*np.pi*Teff**4*cgs.sigm_sb)) 
nphot3 = np.array([Nphots3(Teff[i],Rstar[i],13.6,-1) for i in range(len(Rstar))])
Rs3 = Stromgren(nphot3,2.59e-13,nH)

print('numd = ', nH)
print('trec = {:.4e}'.format( 1/(2.59e-13*nH)))
print(Rs3)
print(10*Rs3)

nH = 10000
Rstar = np.sqrt(Lstar/(4*np.pi*Teff**4*cgs.sigm_sb)) 
nphot3 = np.array([Nphots3(Teff[i],Rstar[i],13.6,-1) for i in range(len(Rstar))])
Rs3 = Stromgren(nphot3,2.59e-13,nH)

print('numd = ', nH)
print('trec = {:.4e}'.format( 1/(2.59e-13*nH)))
print(Rs3)
print(10*Rs3)
