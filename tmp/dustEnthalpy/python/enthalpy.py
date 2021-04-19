import CGS as cgs
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma as gammaFunc
from scipy.special import loggamma as loggammaFunc
def U_graphite(T, N):
    Uatom = 4.15e-22*T**3.3 / (1+6.51e-3*T+1.5*1e-6*T*T+8.3*1e-7*T**(2.3))
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



numa = 6
a = np.logspace(np.log10(3.75*cgs.A), np.log10(1e-4), num = numa)
labs = ['{:.2e} $\mu$m'.format(ai/1e-4) for ai in a]
a = a.reshape((1,numa))


rho_s   = 3.5
mmol_s  = cgs.mH*(28.086+2*24.305+4*15.999)
N_s = 4*np.pi*a**3/3 * rho_s /mmol_s 


rho_c   = 2.2
mmol_c  = cgs.mH*12.011
N_c = 4*np.pi*a**3/3 * rho_c /mmol_c 

nTs = 100
T = np.logspace(0, np.log10(4000), num = nTs)
T = T.reshape((nTs, 1))
U_c = U_graphite(T, N_c)
U_s = np.zeros(U_c.shape)

for i in range(len(T)):
    U_s[i,:] = U_silicate(T[i],N_s, a)

for ia in range(numa):
    P = plt.plot(T, U_c[:,ia])
    plt.plot(T, U_s[:,ia], ls = '--', c=P[0].get_color())
plt.xscale('log')
plt.yscale('log')
plt.show()

epm_c = 0.75*cgs.kb*420
m = U_c/(epm_c)
f = 3*N_c-6
g = m/(3*N_c - 6)
b = (81200 - 20000*N_c**(-1./3.))*cgs.kb/(epm_c)
SN = np.zeros(U_c.shape)
SN = np.exp(np.log((1+g)/g)*b + loggammaFunc(g*f+1) + loggammaFunc(g*f+f-b) - loggammaFunc(g*f+1-b) - loggammaFunc(g*f+f))

SN[g*f+1-b <= 0]=0
SN[g*f+f-b <= 0]=0

fig = plt.figure()
for ia in range(numa):
    print(a[0,ia],SN[:,ia])
    P = plt.plot(T, SN[:,ia], label=labs[ia])


plt.legend()
plt.xlabel(r'$T$ [k]')
plt.ylabel(r'Suppresion factor')
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-8,1)
plt.savefig('suppresion_graphite.png')
plt.close(fig)


fig = plt.figure()
for ia in range(numa):
    print(a[0,ia],SN[:,ia])
    P = plt.plot(T*a[0,ia]*1e4, SN[:,ia], label=labs[ia])
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-1,1e1)
plt.ylim(1e-8,1)
plt.show()


An = 4*np.pi*a*a
Rn = 4.6e29*np.exp(-(81200-20000*N_c**(-1/3))/T)

da = -Rn*An*SN/(4*np.pi*rho_c/mmol_c *a*a)*cgs.yr/1e-4
print(np.min(da))
fig = plt.figure()
for ia in range(numa):
    P = plt.plot(T, -da[:,ia], label=labs[ia])
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$T$ [K]')
plt.ylabel(r'$\frac{\mathrm{d}a}{\mathrm{d}t}$ [$\mu$m yr$^{-1}$]')
plt.ylim(1e-8,1e13)
plt.xlim(1e3, 4e3)
plt.savefig('dadt_graphite.png')
plt.close(fig)
