import CGS as cgs
import numpy as np
import matplotlib.pyplot as plt
import argparse
ap=argparse.ArgumentParser()
ap.add_argument('--Emin', type = float, default = 5.6*cgs.eV)
ap.add_argument('--numBinsSubIon', type = int, default = 1)
ap.add_argument('--Emax', type = float, default = 100*cgs.eV)
ap.add_argument('--numBinsFullIon', type = int, default = 1)
ap.add_argument('--dxOrt', type = float, default = 2.44*1e12)
ap.add_argument('--strength', type = float, default = 1.0)
args=ap.parse_args()



def draineField(emin, emax):
    lammin = cgs.h*cgs.c/emax
    lammax = cgs.h*cgs.c/emin
    lams   = np.linspace(lammin, lammax)
    

    lams_3 = lams/(1000*cgs.A)
    lamUlam = 6.8e-14 * (31.016*lams_3**-3 - 49.913*lams_3**-4 + 19.897*lams_3**-5)
    dlams = lams[1:] - lams[:-1]

    Num  = np.sum(lamUlam[:-1]/cgs.h * dlams)
    Ener = np.sum(lamUlam[:-1]*cgs.c/lams[:-1] * dlams)/Num

    return Num/2, Ener




nbins = args.numBinsSubIon + 2 + args.numBinsFullIon

EbinEdges = np.zeros(nbins + 1)


#below 11.2
delE = np.exp((np.log(11.2*cgs.eV) - np.log(args.Emin))/args.numBinsSubIon)
for i in range(args.numBinsSubIon):
    EbinEdges[i] = args.Emin * delE**i

EbinEdges[args.numBinsSubIon]     = 11.2*cgs.eV
EbinEdges[args.numBinsSubIon +1 ] = 13.6*cgs.eV


delE = np.exp((np.log(args.Emax) - np.log(15.2*cgs.eV))/args.numBinsFullIon)

for i in range(args.numBinsFullIon):
    EbinEdges[args.numBinsSubIon + 2 + i] = 15.2*cgs.eV  *delE**i
EbinEdges[-1] = args.Emax

with open('sed.dat', 'w') as f:
    for i in range(nbins):  
        if(i <= args.numBinsSubIon):
            Nphot, aveE = draineField(EbinEdges[i], EbinEdges[i+1])
        else:
            Nphot = 0
            aveE  = (EbinEdges[i] + EbinEdges[i+1])/2
        Nphot = Nphot * args.dxOrt**2 * args.strength
        f.write('{:.4e}, '.format(Nphot) + '{:.4e}\n'.format(aveE))




    

