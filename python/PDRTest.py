import CGS as cgs
import numpy as np
import matplotlib.pyplot as plt
import argparse
adi = 5./3

muHe = 4.002602
muC  = 12.011
muO  = 15.9994
muSi = 28.0855
abundc = 1e-4
abundo = 3e-4
abundhe = 1e-1

abar = 1 + abundhe*muHe + abundc*muC + abundo*muO
mf_scale = 1 + abundc*muC

ap=argparse.ArgumentParser()
ap.add_argument('f', nargs='+')
args=ap.parse_args()

nH = 1000
alpha = 2.59e-13

trecomb = 1/(nH*alpha)



radius = []
abH  = []
abH2 = []
abHp = []

Temp  = []
Dens  = []
Pres  = []

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
    print(time, time/trecomb)
    data = np.loadtxt(f)

    rad = data[:,0]
    dens = data[:,1]
    pres = data[:,3]
    xH0  = data[:,4]*mf_scale
    xH2  = data[:,5]*mf_scale
    xHp  = data[:,6]*mf_scale
    
  
    numd = dens/(abar*cgs.mH)*(1+xHp-xH2/2)
    temp = pres / (cgs.kb*numd)
    

    radius.append(rad)
    abH.append(xH0)
    abH2.append(xH2)
    abHp.append(xHp)
    Pres.append(pres)
    Dens.append(dens)
    Temp.append(temp)



fig = plt.figure()
pid = [131,132,133]
xmin = [0.2, 0.5, 0.9]
xmax = [1.2, 1.5, 1.9]
for i in range( len(args.f)):
    ax = plt.subplot(pid[i])
    ax.plot(radius[i]/cgs.pc, abH[i], c = 'k', ls = ':')
    ax.plot(radius[i]/cgs.pc, abHp[i], c = 'k', ls = '--')
    ax.plot(radius[i]/cgs.pc, abH2[i], c = 'k')
    ax.set_xlim(xmin[i],xmax[i])
    ax.set_ylim(0,1.05)
plt.show()

fig = plt.figure()
ymin = [0., 0., 0.]
ymax = [12000, 2.5e-20, 3.5e-9]
axT = plt.subplot(311)
axD = plt.subplot(312)
axP = plt.subplot(313)
for i in range( len(args.f)):
    axT.plot(radius[i]/cgs.pc, Temp[i], c = 'k')
    axD.plot(radius[i]/cgs.pc, Dens[i], c = 'k')
    axP.plot(radius[i]/cgs.pc, Pres[i], c = 'k')
axT.set_ylim(ymin[0],ymax[0])
axD.set_ylim(ymin[1],ymax[1])
axP.set_ylim(ymin[2],ymax[2])

axT.set_xlim(0,2)
axD.set_xlim(0,2)
axP.set_xlim(0,2)
plt.show()
