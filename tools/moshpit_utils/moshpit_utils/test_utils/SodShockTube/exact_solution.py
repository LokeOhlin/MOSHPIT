import numpy as np
#
# semi-analytic solution to the rieman problem based of code from bruce fryxell
#

def soundspeed_gas(Pstate, rpar):
    gamma = rpar[0]
    dens = Pstate[0]
    pres = Pstate[2]
    return np.sqrt(gamma * pres / dens)


def soundspeed_dust(Pstate, rpar):
    gamma = rpar[0]
    dust_to_gas = rpar[2]

    dens = Pstate[0]
    pres = Pstate[2]
    return np.sqrt(gamma * pres / dens / (1+dust_to_gas))


def sod_function(pres1, pres4, pres5, dens1, dens5, cs1, cs5, gamma):
    z = pres4/pres5 - 1
    gammam = gamma - 1
    gammap = gamma + 1
    fact = gammam / (2*gamma) * cs5/cs1 * z/np.sqrt(1 + gammap/(2*gamma)* z)
    fact = (1-fact)**(2*gamma/gammam)
    return pres1*fact - pres4

def solve_sodshock(PL, PR, rpar, time, xs, cs_method = soundspeed_gas, max_error = 1e-6, maxiter = 1000):
   
    gamma = rpar[0]
    xdisc = rpar[1]

    dens1 = PL[0]
    velo1 = PL[1]
    pres1 = PL[2]
    cs1   = cs_method(PL, rpar)

    dens5 = PR[0]
    velo5 = PR[1]
    pres5 = PR[2]
    cs5   = cs_method(PR, rpar)

    
    # initial guess
    pres4_0 = pres1
    pres4_1 = pres5
    pres4 = pres4_1
    f0 = sod_function(pres1, pres4_0, pres5, dens1, dens5, cs1, cs5,  gamma)
    i = 0 
    while(True):
        f1 = sod_function(pres1, pres4_1, pres5, dens1, dens5, cs1, cs5,  gamma)
        if f1 == f0:
            break
        
        pres4 = pres4_1 - (pres4_1 - pres4_0)*f1/(f1-f0)
        error = np.abs(pres4 - pres4_1)/pres4_1

        if error < max_error:
            break

        pres4_0 = pres4_1
        pres4_1 = pres4
        f0 = f1

        i = i+1
        if i > maxiter:
            sys.exit("SOD SHOCK TUBE ERROR: Solution to post-shock pressure failed to converge within maxiter = %d steps for max error = %.4e.\nReached error = %4e"%(maxiter, max_error, error))

    z = pres4/pres5 - 1

    gammam = gamma - 1
    gammap = gamma + 1
    gmfact1 = 0.5*gammam / gamma
    gmfact2 = 0.5*gammap / gamma
    
    tmp = np.sqrt(1 + gmfact2 * z)
    velo4 = cs5 * z  / (gamma * tmp)
    dens4 = dens5 * (1 + gmfact2 * z)/(1 + gmfact1 * z)
    #velo4 = velo5 + (pres4 - pres5)/np.sqrt(dens5/2 * (gammap*pres4 + gammam*pres5))
    #dens4 = dens5 * (pres4 + pres5*gammam/gammap) / (pres5 + pres4*gammam/gammap)

    
    veloS = cs5 * tmp
    
    pres3 = pres4
    velo3 = velo4
    dens3 = dens1 * (pres3/pres1)**(1/gamma)

    cs3 = cs_method([dens3, velo3, pres3], rpar)

    print((0.75-xdisc)/veloS)
    pos_shock = xdisc + veloS * time
    pos_cdisc = xdisc + velo3 * time
    pos_rfac1 = xdisc + (velo3 - cs3) * time
    pos_rfac0 = xdisc - cs1 * time

    dens = np.zeros(xs.shape)
    velo = np.zeros(xs.shape)
    pres = np.zeros(xs.shape)
    
    zoneI = xs < pos_rfac0
    if np.sum(zoneI) > 0:
        dens[zoneI] = dens1
        velo[zoneI] = velo1
        pres[zoneI] = pres1

    zoneII = (xs >= pos_rfac0) * (xs < pos_rfac1)
    if np.sum(zoneII) > 0:
        xII = xs[zoneII]
        velII = 2/gammap * (cs1 + (xII - xdisc)/time)
        tmp = 1 - 0.5*gammam * velII / cs1
        dens[zoneII] = dens1 * tmp**(2/gammam)
        velo[zoneII] = velII
        pres[zoneII] = pres1 * tmp**(2*gamma/gammam)


    zoneIII = (xs >= pos_rfac1) * (xs < pos_cdisc)
    if np.sum(zoneIII) > 0:
        dens[zoneIII] = dens3
        velo[zoneIII] = velo3
        pres[zoneIII] = pres3

    zoneIV = (xs >= pos_cdisc) * (xs < pos_shock)
    if np.sum(zoneIV) > 0:
        dens[zoneIV] = dens4
        velo[zoneIV] = velo4
        pres[zoneIV] = pres4

    zoneV = (xs >= pos_shock)
    if np.sum(zoneV) > 0:
        dens[zoneV] = dens5
        velo[zoneV] = velo5
        pres[zoneV] = pres5
    return dens, velo, pres
