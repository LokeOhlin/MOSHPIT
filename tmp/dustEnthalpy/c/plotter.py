import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt('testDist_3.5e-8.dat')
Ts = data[:,0]
Tl = (Ts[1:] + Ts[:-1])/2.0
dT = np.zeros(Ts.shape)
dT[1:-1] = Tl[1:]-Tl[:-1]
dT[-1] = 2*(Ts[-1] - Tl[-1])
dT[0] = 2*(Tl[0] - Ts[0])


plt.plot(Ts, data[:,1]/(np.log((Ts+dT*0.5)/(Ts-dT*0.5))))

data = np.loadtxt('testDist_1.75e-6.dat')
Ts = data[:,0]
Tl = (Ts[1:] + Ts[:-1])/2.0
dT = np.zeros(Ts.shape)
dT[1:-1] = Tl[1:]-Tl[:-1]
dT[-1] = 2*(Ts[-1] - Tl[-1])
dT[0] = 2*(Tl[0] - Ts[0])

plt.plot(Ts, data[:,1]/(np.log((Ts+dT*0.5)/(Ts-dT*0.5))))
plt.xscale('log')
plt.yscale('log')
plt.xlim(5,1e4)
plt.ylim(1e-12,10)
plt.show()
