import CGS as cgs
import numpy as np
import matplotlib.pyplot as plt
import argparse

adi = 5./3

ap=argparse.ArgumentParser()
ap.add_argument('f', nargs='+')
args=ap.parse_args()



def Stromgren(Nphot,alpha,numdens):
    return (3*Nphot/(4*np.pi*alpha))**(1/3)*numdens**(-2/3)


nH = 1000
alpha = 2.6e-13

trecomb = 1/(nH*alpha)

Nion = 1e49
RST = Stromgren(Nion, alpha, nH)
RST = 4.49904612e+17 
print(Nion)
print(RST, RST/cgs.pc)

times = np.zeros(len(args.f))
radIF = np.zeros(len(args.f))
for i, f in enumerate(args.f):
    print(f)
    time = 0
    with open(f, 'r') as fi:
        header = True
        while(header):
            line = fi.readline()
            if line[0] == '#':
                line = line.split(' = ')
                if(line[0] == '# time'):
                    time = float(line[1])
            else:
                header = False
    data = np.loadtxt(f)

    rad = data[:,0]
    dens = data[:,1]
    xH0  = data[:,4]/dens
    xH2  = data[:,5]/dens
    xHp  = data[:,6]/dens 

    if(np.sum(xHp > 5e-1)==0):
        continue
    idx = np.argmin(np.abs(xHp-0.5))
    radIF[i] = rad[idx]
    times[i]  = time


print(radIF/(RST*(1-np.exp(-times/trecomb))**(1./3.)))

tnorm = np.linspace(0,2,1000)
plt.plot(times/trecomb,radIF/RST)
plt.plot(tnorm, (1-np.exp(-tnorm))**(1./3.), c='r')
plt.xlim(0,1.4)
plt.xlabel(r"$t/t_\mathrm{rec}$")
plt.ylabel(r"$R_\mathrm{ion}/R_\mathrm{St}$")
plt.ylim(0,1)
plt.show()


