import dustywaves
import numpy as np
import matplotlib.pyplot as plt
time = 0.5
ampl = 1 
cs = 1
Kdragin = 1
lam = 1
x0 = 0
rhogeq = 1
rhodeq = 1
xplot = np.linspace(0,1)
vgaso = np.zeros(xplot.shape)
vdusto = np.zeros(xplot.shape)
rhogaso = np.zeros(xplot.shape)
rhodusto = np.zeros(xplot.shape)
ierr = 0
vgaso,vdusto,rhogaso,rhodusto,ierr = dustywaves.dustywaves.exact_dustywave(time,ampl,cs,Kdragin,lam,x0,rhogeq,rhodeq, xplot)
plt.plot(xplot, vgaso)
plt.plot(xplot, vdusto)
plt.ylim(-1,1)
plt.xlim(0,1)
plt.show()
plt.plot(xplot, rhogaso)
plt.plot(xplot, rhodusto)
plt.ylim(0,2)
plt.xlim(0,1)
plt.show()


