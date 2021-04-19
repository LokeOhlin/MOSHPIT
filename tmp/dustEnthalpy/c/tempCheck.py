import numpy as np
import matplotlib.pyplot as plt
import argparse
ap=argparse.ArgumentParser()

#---------------outputs-----------------------------
ap.add_argument('f', nargs='+')

args=ap.parse_args()

for fil in args.f:
    TSS = 50
    with open(fil, 'r') as f:
        while True:
            line = f.readline()
            if(line[0] != '#'):
                break
            else:
                line = line.split('=')
                if line[0] == '#TSS':
                    TSS = float(line[1])
    data = np.loadtxt(fil)
    
    Ts = data[:,0]
    Tl = (Ts[1:] + Ts[:-1])/2.0
    dT = np.zeros(Ts.shape)
    dT[1:-1] = Tl[1:]-Tl[:-1]
    dT[-1] = 2*(Ts[-1] - Tl[-1])
    dT[0] = 2*(Tl[0] - Ts[0])
    print(fil, np.min(Ts), dT[0], np.max(Ts), np.sum(Ts*data[:,1]), TSS)
    P = plt.plot(Ts, data[:,1]/(np.log((Ts+dT*0.5)/(Ts-dT*0.5))))
    plt.axvline(TSS, ls = '--', c = P[0].get_color())
plt.xscale('log')
plt.yscale('log')
plt.xlim(5,5e3)
plt.ylim(1e-12,10)
plt.show()
