import numpy as np

from moshpit_utils.units import cgs



def getEofRp(rp, zp, graphite):
    if graphite:
        ap = -4.37
    else:
        ap = -3.34

    E = (rp * 1e7) * 10**(-2.8*zp**-0.21 - ap) * cgs.eV
    return E

def getEthresh(Md, Mp, U0):
    mui = Md/Mp
    gi  = 4*Mp*Md/(Mp + Md)**2
    
    if 1/mui <= 0.3:
        Ethresh = U0/(gi*(1-gi))*cgs.eV
    else:
        Ethresh = 8*U0*mui**(-1/3)*cgs.eV
    return Ethresh
def finite_grain_correction(E, agrain, zp, graphite):
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
    
    rp = 10**(2.8*zp**-0.21 + ap) * (E/cgs.eV) * 1e-7
    if(rp == 0):
        return 1
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


def getVlims(maxRP, agrain, proj_mass, proj_atom, Md, U0, graphite, minv_fact=0.1):
    # maximum velocity from the highest penetration depth that we consider
    Emax    = np.array([getEofRp(maxRP, proj_atom[i], graphite) for i in range(len(proj_atom))])
    fg_vel  = np.sqrt(2*Emax/(proj_mass*cgs.mH))
    fg_corr = np.array([finite_grain_correction(Emax[i], agrain, proj_atom[i], graphite) for i in range(len(proj_atom))])
    maxvel = np.max(fg_vel)

    # minimum velocity as energy threshold times a factor
    Ethresh = np.array([getEthresh(Md, proj_mass[i], U0) for i in range(len(proj_mass))])
    minvel  = np.min(np.sqrt(2*Ethresh/(proj_mass*cgs.mH)))*minv_fact

    return minvel, maxvel


def maxwellLowerYield(vel, integral, Mp, Zp, zp, Md, Zd, U0, K, kbTgas, vrel, agrain, graphite, Ethresh, Enorm, maxnorm, maxexpnorm):
    Ener = Enorm*vel*vel
    yld = sputteringYield(Ener, Mp, Zp, Md, Zd, U0, K, Ethresh)
    Maxwell = maxnorm * vel/vrel * np.exp(-maxexpnorm*(vel + vrel)**2)
    fx =  finite_grain_correction(Ener, agrain, zp, graphite)
    return yld*Maxwell*vel*fx

def maxwellUpperYield(vel, integral, Mp, Zp, zp, Md, Zd, U0, K, kbTgas, vrel, agrain, graphite, Ethresh, Enorm, maxnorm, maxexpnorm):
    Ener = Enorm*vel*vel
    yld = sputteringYield(Ener, Mp, Zp, Md, Zd, U0, K, Ethresh)
    Maxwell = maxnorm * vel/vrel * np.exp(-maxexpnorm*(vel - vrel)**2)
    fx =  finite_grain_correction(Ener, agrain, zp, graphite)
    return yld*Maxwell*vel*fx

def maxwellYield(vel, integral, Mp, Zp, zp, Md, Zd, U0, K, kbTgas, agrain, graphite, Ethresh, Enorm, maxnorm, maxexpnorm):
    Ener = Enorm*vel*vel
    yld = sputteringYield(Ener, Mp, Zp, Md, Zd, U0, K, Ethresh)
    Maxwell = maxnorm**3 * 4*np.pi* vel**2 * np.exp(-maxexpnorm*vel**2)
    fx =  finite_grain_correction(Ener, agrain, zp, graphite)
    return yld*Maxwell*vel*fx


