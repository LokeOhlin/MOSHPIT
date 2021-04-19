import CGS as cgs
import numpy as np
import matplotlib.pyplot as plt
from Tdist import getDist, get_UfromT
from scipy.special import gamma as gammaFunc
from scipy.special import loggamma as loggammaFunc
def U_graphite(T, N):
    Uatom = 4.15*1e-22*T**(3.3) / (1+6.51*1e-3*T+1.5*1e-6*T*T+8.3*1e-7*T**(2.3))
    return (1-2/N)*N*Uatom

def U_silicate(T, N, a):
    if (T < 50):
        U = 1.4*1e3/3 * T**3
    else:
        U = 1.4*1e3/3* (50)**3
        if (T < 150): 
            U += 2.2*1e4 / (2.3) * ( T**(2.3) - 50**(2.3))
        else: 
            U += 2.2*1e4 / (2.3) * ( 150**(2.3) - 50**(2.3))
            if (T < 500):
                U += 4.8*1e5/(1.68) * (T**(1.68) - 150**(1.68))
            else:
                U += 4.8*1e5/(1.68) * (500**(1.68) - 150**(1.68))
                U += 3.41*1e7*(T - 500)

    return U * (1-2/N) * 4/3*np.pi*a**3 


def U_silicate_arr(T, N):
    pass

rho_c   = 2.2
mmol_c  = cgs.mH*12.011


a = 3.7* cgs.A
N_c = 4*np.pi*rho_c*a**3 /(3*mmol_c)

U_c = get_UfromT(a,2800)
epm_c = 0.75*cgs.kb*420
m = U_c/(epm_c)
f = 3*N_c - 6
g = m/f
b = (81200 - 20000*N_c**(-1./3.))*cgs.kb/(epm_c)
SN = np.exp((np.log(1+g)-np.log(g))*b + loggammaFunc(g*f+1) + loggammaFunc(g*f+f-b) - loggammaFunc(g*f+1-b) - loggammaFunc(g*f+f))
print(a, N_c, U_c, f, b*epm_c/cgs.kb, b, SN)


Ts, delT, Tdist = getDist(a, 500, get_UfromT(a, 7.5e3))
Tdist = Tdist/np.log((Ts+delT/2)/(Ts-delT/2))

U_c = U_graphite(Ts,N_c) 
plt.plot(Ts,Tdist)
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-12,10)
plt.xlim(5,5e3)
plt.show()
plt.plot(U_c/cgs.eV,Tdist)
plt.axvline(13.6)
plt.xscale('log')
plt.yscale('log')
plt.show()
epm_c = 0.75*cgs.kb*420
m = U_c/(epm_c)
f = 3*N_c - 6
g = m/f
b = (81200 - 20000*N_c**(-1./3.))*cgs.kb/(epm_c)
SN = np.zeros(U_c.shape)
SN = np.exp((np.log(1+g)-np.log(g))*b + loggammaFunc(g*f+1) + loggammaFunc(g*f+f-b) - loggammaFunc(g*f+1-b) - loggammaFunc(g*f+f))

SN[g*f+1-b <= 0]=0


An = 4*np.pi*a*a
Rn = 4.6e29*np.exp(-(81200-20000*N_c**(-1/3))/Ts)
#print(np.max(U_c/cgs.eV), np.min(U_c/cgs.eV))
dR = -Rn*An*SN
#print(An, SN, Rn, -dR)
da = Rn*An*SN/(4*np.pi*rho_c/mmol_c *a*a)


plt.plot(U_c/cgs.eV,-dR*Tdist)
plt.xlim(8,27)
plt.yscale('log')
plt.ylim(7e-20, 1e-10)
plt.show()
