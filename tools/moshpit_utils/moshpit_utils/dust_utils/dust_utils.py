import moshpit_utils.units.cgs as cgs
import numpy as np



def dustBinMass(number, slope, am, ac, ap):
    Nfact = (ap**4 - am**4)/(4*(ap-am)) 
    Sfact =   ap**4*(ap/5 - ac/4)
    Sfact -=  am**4*(am/5 - ac/4)

    return number*Nfact + slope*Sfact


def dustBinEdges(hfile):
    abin_c = hfile['agrain'][:]
    # get bin edges
    amin = hfile['Parameters'].attrs.get("dust_amin")
    amax = hfile['Parameters'].attrs.get("dust_amax")
    Nabin = len(abin_c)
    da = np.exp((np.log(amax)-np.log(amin))/Nabin)
    aes = np.zeros(Nabin + 1)
    for ibin in range(Nabin + 1):
        aes[ibin] = amin*da**(ibin)
    return aes

def dustDensity(hfile):
    isilicone = hfile['Headers'].attrs.get('isilicone')
    num_c = hfile['number'][:,:isilicone]
    num_s = hfile['number'][:,isilicone:]

    slope_c = hfile['slope'][:,:isilicone]
    slope_s = hfile['slope'][:,isilicone:]

    abin_c = hfile['agrain'][:]
    aes = hfile["agrain_binEdges"][:]

    if isilicone > 0:
        mass_c = dustBinMass(num_c, slope_c, aes[:-1][np.newaxis,:], abin_c[np.newaxis,:], aes[1:][np.newaxis,:])
    else:
        mass_c = np.array([[0]])
    if isilicone < len(abin_c):
        mass_s = dustBinMass(num_s, slope_s, aes[:-1][np.newaxis,:], abin_c[np.newaxis,:], aes[1:][np.newaxis,:])
    else:
        mass_s = np.array([[0]])

    print(np.sum(mass_c))
    mass_c = 4*np.pi/3 * 2.26 * mass_c
    mass_s = 4*np.pi/3 * 3.5  * mass_s


    totmas = np.sum(mass_c, axis = 1) + np.sum(mass_s, axis = 1)
    return totmas

def dustToGas(hfile):
    ddust = dustDensity(hfile)
    dgas = hfile["density"][:]
    return ddust/dgas
