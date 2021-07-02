import CGS as cgs
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode as ode



def finite_grain_correction(E, agrain, Zp, graphite):
    if graphite:
        p1 = 4.9
        p2 = 0.55
        p3 = 0.77
        p4 = 4.7
        p5 = 3.0
        p6 = 1.2
        ap = -4.37
    else:
        p1 = 1.0
        p2 = 0.5
        p3 = 1.0
        p4 = 1.8
        p5 = 2.1
        p6 = 0.76
        ap = -3.34


    rp = 10**(2.8*Zp**-0.21 + ap) * (E/cgs.eV) * 1e-7

    xs = agrain/(0.7*rp)

    fx = 1 + p1*np.exp(-np.log(xs/p2)**2/(2*p3**2)) - p4*np.exp(-(p5*xs - p6)**2)
    return max(fx,0)

def stopping_sigma(E, Mp, Zp, Md, Zd):
    a0 = 0.529*cgs.A
    asc = 0.885*a0*(Zp**(2/3) + Zd**(2/3))**(-1/2)

    ered = Md/(Mp + Md) * asc * E/(Zp*Zd*cgs.qe_ESU**2)

    si = 3.441*np.sqrt(ered) * np.log(ered + 2.718) / (1+6.35*np.sqrt(ered) + ered*(-1.708 + 6.882*np.sqrt(ered)))


    return 4*np.pi*asc*Zp*Zd*cgs.qe_ESU**2 * Mp/(Mp+Md) * si



def sputteringYield(E, Mp, Zp, Md, Zd, U0, K):
    mui = Md/Mp
    gi  = 4*Mp*Md/(Mp + Md)**2
    
    if 1/mui <= 0.3:
        Ethresh = U0/(gi*(1-gi))*cgs.eV
    else:
        Ethresh = 8*U0*mui**(-1/3)*cgs.eV
    if E < Ethresh:
        return 0 

    if mui <= 0.5:
        alphai =  0.2
    elif mui < 1:
        alphai =  0.1/mui + 0.25*(mui - 0.5)**2
    else:
        alphai =  0.3*(mui - 0.6)**(2/3)


    Si = stopping_sigma(E, Mp, Zp, Md, Zd)

    yld = 4.2e14 * (Si/cgs.eV)/U0 * alphai/(K*mui +1) * (1 - (Ethresh/E)**(2/3))*(1 - Ethresh/E)**2
    return yld

def Maxwellian(Ener, kbTgas, Mpart):
    return (Mpart/(2*np.pi*kbTgas))**(1/2)*4*Ener/(kbTgas)*np.exp(-Ener/(kbTgas))

def maxwellSputteringYield(vel, integral, Mp, Zp, Md, Zd, U0, K, kbTgas, agrain, graphite):
    Ener = 0.5*Mp*cgs.mH*vel**2

    yld = sputteringYield(Ener, Mp, Zp, Md, Zd, U0, K)
    Maxwell = Maxwellian(Ener, kbTgas, Mp*cgs.mH)
    fx =  finite_grain_correction(Ener, agrain, Zp, graphite)
    if(yld*Maxwell*vel*fx < 0):
        print(vel, yld, Maxwell, fx)
    return yld*Maxwell*vel*fx

def getAverageSputtering(Mp, Zp, Md, Zd, U0, K, kbTgas, agrain, graphite):
    mui = Md/Mp
    gi  = 4*Mp*Md/(Mp + Md)**2
    
    if 1/mui <= 0.3:
        Ethresh = U0/(gi*(1-gi))*cgs.eV
    else:
        Ethresh = 8*U0*mui**(-1/3)*cgs.eV
    
    integrator = ode(maxwellSputteringYield).set_integrator('dopri5')
    integrator.set_initial_value(0, np.sqrt(2*Ethresh/(Mp*cgs.mH))).set_f_params(Mp, Zp, Md, Zd, U0, K, kbTgas, agrain, graphite)
    y = integrator.integrate(np.sqrt(100*kbTgas/(Mp*cgs.mH)))
    return max(y[0],0)


def getTotalLoss(Md, Zd, U0, K, kbTgas, proj_mass, proj_atom, proj_abund, agrain, graphite):

    totSputter = 0
    for pi in range(len(proj_mass)):
        totSputter += getAverageSputtering(proj_mass[pi], proj_atom[pi], Md, Zd, U0, K, kbTgas, agrain, graphite)*proj_abund[pi]

    return totSputter

# projectile properties
proj_mass  = np.array([1, 4  , 12  , 14  , 16  ])
proj_atom  = np.array([1, 2  , 6   , 7   , 8   ])
proj_abund = np.array([1, 0.1, 1e-4, 1e-5, 3e-4]) 

# carbon
Md_c = 12
Zd_c = 6
K_c  = 0.61
U0_c = 4.0

#silicate
Md_s = 20
Zd_s = 10
K_s  = 0.1
U0_s = 5.8
numTemps   = 50
numGrains  = 30
Tgas    = np.logspace(3,10, numTemps)
totLoss_graphite = np.zeros(Tgas.shape)
totLoss_silicone = np.zeros(Tgas.shape)
agrains = np.logspace(np.log10(3*cgs.A), np.log10(1e-3), numGrains)
with open('sputtering_yield.dat', 'w') as f:
    f.write('# ngrain = {}\n'.format(numGrains))
    f.write('# ntemps = {}\n'.format(numTemps))
    for agrain in agrains:
        print(agrain)
        for i in range(len(Tgas)):
            totLoss_graphite[i] = getTotalLoss(Md_c, Zd_c, U0_c, K_c, cgs.kb*Tgas[i], proj_mass, proj_atom, proj_abund, agrain, 1)
            totLoss_silicone[i] = getTotalLoss(Md_s, Zd_s, U0_s, K_s, cgs.kb*Tgas[i], proj_mass, proj_atom, proj_abund, agrain, 0)
            f.write('{} {} {} {}\n'.format(agrain, Tgas[i], totLoss_graphite[i], totLoss_silicone[i]))
        totLoss_graphite = totLoss_graphite * (Md_c*cgs.mH/(2*2.25))
        totLoss_silicone = totLoss_silicone * (Md_s*cgs.mH/(2*3.2))
        P = plt.plot(Tgas, totLoss_graphite/(1e-4/cgs.yr))
        plt.plot(Tgas, totLoss_silicone/(1e-4/cgs.yr), c = P[0].get_color(), ls = '--')
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e4, 1e10)
plt.ylim(1e-9, 1e-4)
plt.show()
'''
proj_mass  = np.array([1, 2, 4, 20, 40])
proj_atom  = np.array([1, 1, 2, 10, 18])
proj_cols = ['k', 'r', 'g', 'orange', 'b']

Eners = np.logspace(1, 6) * cgs.eV
yelds = np.zeros(Eners.shape)

fig, ax = plt.subplots(1,1)
for pi in range(len(proj_mass)):
    for i in range(len(Eners)):
        yelds[i] = sputteringYield(Eners[i], proj_mass[pi], proj_atom[pi], Md, Zd, U0, K)

    ax.plot(Eners/cgs.eV, yelds, c = proj_cols[pi])
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(10, 1e6)
ax.set_ylim(1e-4, 10)
plt.show()
'''





