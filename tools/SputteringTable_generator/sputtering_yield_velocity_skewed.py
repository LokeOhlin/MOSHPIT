from moshpit_utils.units import  cgs
from integrator import dopri5_integrator
from sputtering_utils import *

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode as ode
from multiprocessing import Pool


def getAverageSputtering(Mp, Zp, zp, Md, Zd, U0, K, kbTgas, vrel, agrain, graphite):
    maxnorm    = np.sqrt(Mp*cgs.mH/(2*np.pi*kbTgas) )
    maxexpnorm = Mp*cgs.mH/(2*kbTgas)
    Enorm      = 0.5*Mp*cgs.mH
    vexpmin = np.sqrt(7e2/maxexpnorm)
    mui = Md/Mp
    gi  = 4*Mp*Md/(Mp + Md)**2
    
    if 1/mui <= 0.3:
        Ethresh = U0/(gi*(1-gi))*cgs.eV
    else:
        Ethresh = 8*U0*mui**(-1/3)*cgs.eV
    Emax_finiteGrain = getEofRp(10*agrain, Zp, graphite)
    vmin = np.sqrt(Ethresh/Enorm)
    vmax = np.sqrt(Emax_finiteGrain/Enorm)
    if vrel > 0:
        #if kbTgas / cgs.kb < 1e5:
        #    Ener = Enorm*vrel*vrel
        #    yld = sputteringYield(Ener, Mp, Zp, Md, Zd, U0, K, Ethresh)
        #    fx =  finite_grain_correction(Ener, agrain, Zp, graphite)
        #    return yld*vrel*fx
        #odep = ode(maxwellUpperYield).set_integrator('dopri5', nsteps = 1000000000)#, atol = 1e-3, rtol = 1e-2)
        #odep.set_f_params(Mp, Zp, Md, Zd, U0, K, kbTgas, vrel, agrain, graphite, Ethresh, Enorm, maxnorm, maxexpnorm)
        #odep.set_initial_value(0, vmin)
        odep = dopri5_integrator(maxwellUpperYield, Mp, Zp, zp, Md, Zd, U0, K, kbTgas, vrel, agrain, graphite, Ethresh, Enorm, maxnorm, maxexpnorm, min_stepsize = 0, y_independent = True, rtol = 1e-8, atol = 1e-8)
        vmin = min(max(-vexpmin + vrel,0), vmax)
        odep.set_ic(vmin, 0)
        
        #odem = ode(maxwellLowerYield).set_integrator('dopri5', nsteps = 1000000000)#, atol = 1e-3, rtol = 1e-2)
        #odem.set_f_params(Mp, Zp, Md, Zd, U0, K, kbTgas, vrel, agrain, graphite, Ethresh, Enorm, maxnorm, maxexpnorm)
        #odem.set_initial_value(0, vmin)
        odem = dopri5_integrator(maxwellLowerYield, Mp, Zp, zp, Md, Zd, U0, K, kbTgas, vrel, agrain, graphite, Ethresh, Enorm, maxnorm, maxexpnorm, min_stepsize = 0, y_independent = True, rtol = 1e-8, atol = 1e-8)
        vmin = min(max(-vexpmin + vrel,0), vmax)
        vmin = min(max(-vexpmin - vrel,0), vmax)
        odem.set_ic(vmin,0)
   
        y  = odep.integrate(vmax, first_step = 1e-7)
        y -= odem.integrate(vmax, first_step = 1e-7)
        y =  y 
        #if not odep.successful() or not odem.successful():
    else:
        odep = ode(maxwellYield).set_integrator('dopri5', nsteps = 1000000000)
        odep.set_f_params(Mp, Zp, zp, Md, Zd, U0, K, kbTgas, agrain, graphite, Ethresh, Enorm, maxnorm, maxexpnorm)
        odep.set_initial_value(0, 0)

        y = odep.integrate(vmax*2)
    return max(y[0],0)


def getTotalLoss(Md, Zd, U0, K, kbTgas, vrel, proj_mass, proj_atom, proj_abund, proj_charge, agrain, graphite):

    totSputter = 0
    for pi in range(len(proj_mass)):
        totSputter += getAverageSputtering(proj_mass[pi], proj_atom[pi], proj_charge[pi], Md, Zd, U0, K, kbTgas, vrel, agrain, graphite)*proj_abund[pi]
    return totSputter


def pool_over_T(Md, Zd, U0, K, Tgas, vrel, proj_mass, proj_atom, proj_abund, proj_charge, agrain, graphite, nproc):
    args = [(Md, Zd, U0, K,  cgs.kb* T, vrel, proj_mass, proj_atom, proj_abund, proj_charge, agrain, graphite) for T in Tgas]

    pool = Pool(processes = nproc)
    return pool.starmap(getTotalLoss, args)


use_multiprocessing = True
nproc = 16

# projectile properties
proj_mass   = np.array([1, 4  , 12  , 14  , 16  ])
proj_atom   = np.array([1, 2  , 6   , 7   , 8   ])
proj_abund  = np.array([1, 0.1, 1e-4, 1e-5, 3e-4]) 
proj_charge = np.array([1, 1, 1, 1, 1]) 

# carbon
Md_c = 12
Zd_c = 6
K_c  = 0.61
U0_c = 4.0

#silicate
Md_s = 20
Zd_s = 10
K_s  = 0.1
U0_s = 6.0

# Table parameters
# Temperature
numTemps   = 30
Tmin = 1e3
Tmax = 1e10

# Grain size
numGrains  = 30
amin = 3*cgs.A
amax = 1e-3

# velocities 
numVels    = 30
minv_fact  = 1
rp_nagrain = 15
# velocites are generated from v(Ethresh) + minv to v(Ethresh) + maxv 
# these are in units of vscale = v(RP = rp_nagrain) - v(Ethresh)
minv = 1e-2
maxv = 1e0
vtilde = np.logspace(np.log10(minv), np.log10(maxv), numVels - 1)


#Space out gas temperatures between Tmax and Tmin
Tgas    = np.logspace(np.log10(Tmin), np.log10(Tmax), numTemps)

#Preallocate arrays
totLoss_graphite = np.zeros(Tgas.shape)
totLoss_silicone = np.zeros(Tgas.shape)

#agrains between amax and amin
agrains = np.logspace(np.log10(amin), np.log10(amax), numGrains)
#agrains = np.array([1e-4])#, 1e-3])

# space out velocities between vmin and vmax in units of maximum velocity 
# (defined as the velocity at which the penetration depth is = 10)
# and add the case for zero relative velocity
#vrels = np.array([5e-1])
# For each grain, both for carbonates and silicates we define the maximum and minimum relative velocity considered
minvels_c = np.zeros(agrains.shape)
maxvels_c = np.zeros(agrains.shape)

minvels_s = np.zeros(agrains.shape)
maxvels_s = np.zeros(agrains.shape)

for ia, agrain in enumerate(agrains):
    maxRP = rp_nagrain*agrain 
    minvels_c[ia], maxvels_c[ia] = getVlims(maxRP, agrain, proj_mass, proj_atom, Md_c, U0_c, 1, minv_fact=minv_fact)
    minvels_s[ia], maxvels_s[ia] = getVlims(maxRP, agrain, proj_mass, proj_atom, Md_s, U0_s, 0, minv_fact=minv_fact)
    
#acols = ['b', 'r', 'g', 'orange']


with open('sputtering_yield.dat', 'w') as f:
    f.write('# ngrain = {}\n'.format(numGrains))
    f.write('# nvels  = {}\n'.format(numTemps))
    f.write('# ntemps = {}\n'.format(numTemps))
    for ia, agrain in enumerate(agrains):
        #if ia != 0:
        #    continue
        # scale vtilde and add to vmin
        vrels_c = np.append(np.array([0]), minvels_c[ia]+vtilde*(maxvels_c[ia] - minvels_c[ia]))
        vrels_s = np.append(np.array([0]), minvels_s[ia]+vtilde*(maxvels_s[ia] - minvels_s[ia]))
        #vrels_s = np.append(np.array([0]), np.logspace(np.log10(minvels_s[ia]), np.log10(maxvels_s[ia]), numVels - 1))
        for iv in range(len(vrels_c)) :
            #if iv != 2:
            #    continue
            graphDone = False
            silicDone = False
            if(iv == 0):
                vtilde_i = 0
            else:
                vtilde_i = vtilde[iv-1]
            if use_multiprocessing:
                print(agrain, vrels_c[iv], vrels_s[iv])
                totLoss_graphite[:] = pool_over_T(Md_c, Zd_c, U0_c, K_c, Tgas, vrels_c[iv], proj_mass, proj_atom, proj_abund, proj_charge, agrain, 1, nproc)
                totLoss_silicone[:] = pool_over_T(Md_s, Zd_s, U0_s, K_s, Tgas, vrels_s[iv], proj_mass, proj_atom, proj_abund, proj_charge, agrain, 0, nproc)
                for i in range(len(Tgas)):
                    f.write('{} {} {} {} {} {} \n'.format(agrain, vrels_c[iv], vrels_s[iv], Tgas[i], totLoss_graphite[i], totLoss_silicone[i]))
            else:
                for i in range(len(Tgas)):
                    #print("\n ----------------- \n")
                    if not graphDone :
                        totLoss_graphite[i] = getTotalLoss(Md_c, Zd_c, U0_c, K_c, cgs.kb*Tgas[i], vrels_c[iv], proj_mass, proj_atom, proj_abund, agrain, 1)
                    if not silicDone :
                        totLoss_silicone[i] = getTotalLoss(Md_s, Zd_s, U0_s, K_s, cgs.kb*Tgas[i], vrels_s[iv], proj_mass, proj_atom, proj_abund, agrain, 0)

                    print(agrain, vrels_c[iv], vrels_s[iv], Tgas[i], totLoss_graphite[i], totLoss_silicone[i])
                    # save in format of agrain, vrel, Tgas and rate for both graphite and silicates
                    f.write('{} {} {} {} {} {} \n'.format(agrain, vrels_c[iv], vrels_s[iv], Tgas[i], totLoss_graphite[i], totLoss_silicone[i]))
                    # Past a givne temperature the sputtering rate only goes down as more of the particles will
                    # hit the dust grains with penetration depths much larger than the grain size
                    # if we hit zero we can stop
                    #if totLoss_graphite[i] <=0 and Tgas[i] > 1e6:
                    #    graphDone = True
                    #if totLoss_silicone[i] <=0 and Tgas[i] > 1e6:
                    #    silicDone = True
                    #if i > 13:
                    #    break
            totLoss_graphite = totLoss_graphite * (Md_c*cgs.mH/(2*2.25))
            totLoss_silicone = totLoss_silicone * (Md_s*cgs.mH/(2*3.2))
            P = plt.plot(Tgas, totLoss_graphite/(1e-4/cgs.yr))#, c = acols[iv])
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





