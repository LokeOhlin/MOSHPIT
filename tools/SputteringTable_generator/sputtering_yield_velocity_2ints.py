import CGS as cgs
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode as ode
import sys

def getEofRp(rp, Zp, graphite):
    if graphite:
        ap = -4.37
    else:
        ap = -3.34

    E = (rp * 1e7) * 10**(-2.8*Zp**-0.21 - ap) * cgs.eV
    return E

def getEthresh(Md, Mp, U0):
    mui = Md/Mp
    gi  = 4*Mp*Md/(Mp + Md)**2
    
    if 1/mui <= 0.3:
        Ethresh = U0/(gi*(1-gi))*cgs.eV
    else:
        Ethresh = 8*U0*mui**(-1/3)*cgs.eV
    return Ethresh
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
    a0 = 5.29e-9
    asc = 0.885*a0/np.sqrt(Zp**(2/3) + Zd**(2/3))

    ered = Md/(Mp + Md) * asc * E/(Zp*Zd*cgs.qe_ESU**2)
    sqrt_ered = np.sqrt(ered)
    si = 3.441*sqrt_ered * np.log(ered + 2.718) / (1+6.35*sqrt_ered + ered*(-1.708 + 6.882*sqrt_ered))


    return 4*np.pi*asc*Zp*Zd*cgs.qe_ESU**2 * Mp/(Mp+Md) * si



def sputteringYield(E, Mp, Zp, Md, Zd, U0, K, Ethresh):
    if E < Ethresh:
        return 0 
    
    mui = Md/Mp
    gi  = 4*Mp*Md/(Mp + Md)**2
    

    if mui <= 0.5:
        alphai =  0.2
    elif mui < 1:
        alphai =  0.1/mui + 0.25*(mui - 0.5)**2
    else:
        alphai =  0.3*(mui - 0.6)**(2/3)


    Si = stopping_sigma(E, Mp, Zp, Md, Zd)

    yld = 2.6217228464419475e+26 * Si/U0 * alphai/(K*mui +1) * (1 - (Ethresh/E)**(2/3))*(1 - Ethresh/E)**2
    return yld

def maxwell2DSputteringYield(vort, integral, vpar, usq, Mp, Zp, Md, Zd, U0, K, kbTgas, vrel, agrain, graphite, Ethresh, Enorm, maxnorm, maxexpnorm):
    vsqr = vort*vort + vpar*vpar
    Ener = Enorm*(vort*vort + usq)
    yld = sputteringYield(Ener, Mp, Zp, Md, Zd, U0, K, Ethresh)
    Maxwell = maxnorm * vort * np.exp(-maxexpnorm*vsqr)
    fx =  finite_grain_correction(Ener, agrain, Zp, graphite)
    vel = np.sqrt(vort*vort + usq) 
    
    return yld*Maxwell*vel*fx

def dvpar(vpar, integral, Mp, Zp, Md, Zd, U0, K, kbTgas, vrel, agrain, graphite, Ethresh, Emax_finiteGrain, Enorm, maxnorm, maxexpnorm, odeort):
    vsq=vpar*vpar
    dv = vpar-vrel
    usq=dv*dv
    # initial value for vort
    Emin = max(Ethresh - Enorm*usq,0)
    Emax = min(max(Emax_finiteGrain-Enorm*usq,0), np.sqrt(25/maxexpnorm))
    if Emax <= Emin:
        return 0

    vmin = np.sqrt(Emin/Enorm)
    vmax = np.sqrt(Emax/Enorm)
    odeort.set_f_params(vpar, usq, Mp, Zp, Md, Zd, U0, K, kbTgas, vrel, agrain, graphite, Ethresh, Enorm, maxnorm, maxexpnorm)
    odeort.set_initial_value(0,vmin)
    odeort.integrate(vmax)
    #if not odeort.successful():
    #    print(vpar, vrel, agrain, odeort.y[0])
    return odeort.y[0]

def getAverageSputtering(Mp, Zp, Md, Zd, U0, K, kbTgas, vrel, agrain, graphite):
    maxnorm    = np.sqrt((Mp*cgs.mH)**3/(kbTgas**3 * 2*np.pi)) 
    maxexpnorm = Mp*cgs.mH/(2*kbTgas)
    Enorm      = 0.5*Mp*cgs.mH

    mui = Md/Mp
    gi  = 4*Mp*Md/(Mp + Md)**2
    
    if 1/mui <= 0.3:
        Ethresh = U0/(gi*(1-gi))*cgs.eV
    else:
        Ethresh = 8*U0*mui**(-1/3)*cgs.eV
    Emax_finiteGrain = getEofRp(10*agrain, Zp, graphite)
    dvmax = np.sqrt(Emax_finiteGrain/Enorm)
    vmin = min(vrel - dvmax, -np.sqrt(25/maxexpnorm))
    vmax = max(vrel - dvmax,  np.sqrt(25/maxexpnorm))

    odepar = ode(dvpar).set_integrator('dopri5', nsteps = 1000000000 )
    odeort = ode(maxwell2DSputteringYield).set_integrator('dopri5', nsteps = 1000000000 )
    odepar.set_f_params(Mp, Zp, Md, Zd, U0, K, kbTgas, vrel, agrain, graphite, Ethresh, Emax_finiteGrain, Enorm, maxnorm, maxexpnorm, odeort)
    
    odepar.set_initial_value(0, vrel)
    y = odepar.integrate(vmax)

    odepar.set_initial_value(0, vrel)
    y -= odepar.integrate(vmin)
    
    return max(y[0],0)


def getTotalLoss(Md, Zd, U0, K, kbTgas, vrel, proj_mass, proj_atom, proj_abund, agrain, graphite):

    totSputter = 0
    for pi in range(len(proj_mass)):
        totSputter += getAverageSputtering(proj_mass[pi], proj_atom[pi], Md, Zd, U0, K, kbTgas, vrel, agrain, graphite)*proj_abund[pi]
    return totSputter


def getVlims(maxRP, agrain, proj_mass, proj_atom, Md, U0, graphite, minv_fact=0.1):
    print("a = {:.4e}".format(agrain))
    # maximum velocity from the highest penetration depth that we consider
    Emax    = np.array([getEofRp(maxRP, proj_atom[i], graphite) for i in range(len(proj_atom))])
    fg_vel  = np.sqrt(2*Emax/(proj_mass*cgs.mH))
    fg_corr = np.array([finite_grain_correction(Emax[i], agrain, proj_atom[i], graphite) for i in range(len(proj_atom))])
    maxvel = np.max(fg_vel)

    # minimum velocity as energy threshold times a factor
    Ethresh = np.array([getEthresh(Md, proj_mass[i], U0) for i in range(len(proj_mass))])
    minvel  = np.min(np.sqrt(2*Ethresh/(proj_mass*cgs.mH)))*minv_fact
    print("{} min = {:.4e}; max = {:.4e}".format(graphite, minvel, maxvel))
    #print("corr = "+ " ".join("{:.4e}".format(corr) for corr in fg_corr ))
    #print("velo = "+ " ".join("{:.4e}".format(velo) for velo in fg_vel ))

    return minvel, maxvel



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

# Table parameters
# Temperature
numTemps   = 1
Tmin = 1e3
Tmax = 1e6

# Grain size
numGrains  = 30
amin = 3*cgs.A
amax = 1e-3

# velocities 
numVels    = 4
minv_fact  = 1
rp_nagrain = 10
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
agrains = np.array([1e-7, 1e-6, 1e-4])#, 1e-3])

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
    
acols = ['b', 'r', 'g', 'orange']


with open('sputtering_yield.dat', 'w') as f:
    f.write('# ngrain = {}\n'.format(numGrains))
    f.write('# nvels  = {}\n'.format(numTemps))
    f.write('# ntemps = {}\n'.format(numTemps))
    for ia, agrain in enumerate(agrains):
        # scale vtilde and add to vmin
        vrels_c = np.append(np.array([0]), minvels_c[ia]+vtilde*(maxvels_c[ia] - minvels_c[ia]))
        vrels_s = np.append(np.array([0]), minvels_s[ia]+vtilde*(maxvels_s[ia] - minvels_s[ia]))
        #vrels_s = np.append(np.array([0]), np.logspace(np.log10(minvels_s[ia]), np.log10(maxvels_s[ia]), numVels - 1))
        for iv in range(len(vrels_c)) :
            graphDone = False
            silicDone = False
            if(iv == 0):
                vtilde_i = 0
            else:
                vtilde_i = vtilde[iv-1]

            for i in range(len(Tgas)):
                print(agrain, vrels_c[iv], vrels_s[iv], Tgas[i])
                if not graphDone :
                    totLoss_graphite[i] = getTotalLoss(Md_c, Zd_c, U0_c, K_c, cgs.kb*Tgas[i], vrels_c[iv], proj_mass, proj_atom, proj_abund, agrain, 1)
                if not silicDone :
                    totLoss_silicone[i] = getTotalLoss(Md_s, Zd_s, U0_s, K_s, cgs.kb*Tgas[i], vrels_s[iv], proj_mass, proj_atom, proj_abund, agrain, 0)

                # save in format of agrain, vrel, Tgas and rate for both graphite and silicates
                f.write('{} {} {} {} {} \n'.format(agrain, vtilde_i, Tgas[i], totLoss_graphite[i], totLoss_silicone[i]))
                # Past a givne temperature the sputtering rate only goes down as more of the particles will
                # hit the dust grains with penetration depths much larger than the grain size
                # if we hit zero we can stop
                if totLoss_graphite[i] <=0 and Tgas[i] > 1e6:
                    graphDone = True
                if totLoss_silicone[i] <=0 and Tgas[i] > 1e6:
                    silicDone = True

            totLoss_graphite = totLoss_graphite * (Md_c*cgs.mH/(2*2.25))
            totLoss_silicone = totLoss_silicone * (Md_s*cgs.mH/(2*3.2))
            P = plt.plot(Tgas, totLoss_graphite/(1e-4/cgs.yr), c = acols[iv])
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





