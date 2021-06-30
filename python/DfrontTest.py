import CGS as cgs
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import argparse

adi = 5./3

ap=argparse.ArgumentParser()
ap.add_argument('f', nargs='+')
args=ap.parse_args()

RST = 0
def spitzer(time, Temp):
    cs = np.sqrt(5/3 * cgs.kb * Temp / cgs.mH)

    return RST * (1 + 7/4 * cs * time / RST)**(4/7)
def rfront(tnorm):
    return (1 - np.exp(-tnorm))**(1/3)

def Stromgren(Nphot,alpha,numdens):
    return (3*Nphot/(4*np.pi*alpha))**(1/3)*numdens**(-2/3)


nH = 1e-3
alpha = 2.6e-13

trec = 1/(nH*alpha)

Nion = 1e49
RST = Stromgren(Nion, alpha, nH)

print(Nion)
print(RST/cgs.pc,RST, "{:.4e}".format(trec))

times = np.zeros(len(args.f))
dts   = np.zeros(len(args.f))
radIF = np.zeros(len(args.f))
Tmin  = np.zeros(len(args.f))
Tavg  = np.zeros(len(args.f))
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
    pres = data[:,3]
    xH0  = data[:,4]/dens
    xH2  = data[:,5]/dens
    xHp  = data[:,6]/dens 
    
    ionReg = xHp > 5e-1
    if(np.sum(ionReg)==0):
        continue
  
    abar = 1 # We assume no other species here
    numd = dens/(abar*cgs.mH)
    temp = pres / (cgs.kb*numd*(1-xH2+xHp))
    Tmin   = np.min (temp[ionReg])
    Tmax   = np.max (temp[ionReg])
    Tmean  = np.mean(temp[ionReg])
    print(Tmin, Tmean, Tmax)
    dr = rad[0] # constant dr, r[0] = dr
    rm = rad - dr
    rp = rad + dr

    vol = 4*np.pi*(rp**3-rm**3)/3

    Tavg[i] = np.average(temp[ionReg], weights = vol[ionReg]*dens[ionReg])

    idx = np.abs(xHp - 0.5).argmin()
    radIF[i] = rad[idx]
    #plt.plot(rad, xHp)
    #plt.plot(rad, temp/np.max(temp))
    #plt.axvline(radIF[i])
    #plt.axvline(RST)
    #plt.show()
    times[i]  = time
    if(i==0):
        dts[i] = time
    else :
        dts[i] = time-dts[i-1]


idx = np.abs(radIF - RST).argmin()
if(radIF[idx] > RST):
    idx = idx - 1
tm = times[idx]
tp = times[idx + 1]
rm = radIF[idx]
rp = radIF[idx + 1]

tSP = tm + (tp - tm) * (RST - rm)/(rp - rm)
print(idx, tm, tp, rm, rp, tSP)


tfit = times[times > tSP] - tSP
Rfit = radIF[times > tSP]

pars = curve_fit(spitzer, tfit, Rfit, 5e3)
Tfit = pars[0][0]

plt.plot(times/trec,radIF/RST)
tplot = np.linspace(tSP, np.max(times))
plt.plot(tplot/trec, spitzer(tplot-tSP, Tfit)/RST , ls = '--')
plt.plot(tplot/trec, spitzer(tplot-tSP, 10000)/RST, ls = '--')
plt.plot(tplot/trec, spitzer(tplot-tSP, 20000)/RST, ls = '--')
tplot = np.linspace(0, np.max(times)/trec, 1000)
plt.plot(tplot, rfront(tplot), ls = ':')
plt.xlim(0, None)
plt.ylim(0, None)
plt.show()


