import numpy as np
import CGS as cgs
import matplotlib.pyplot as plt
import argparse 
ap=argparse.ArgumentParser()

#---------------outputs-----------------------------
ap.add_argument('f', nargs='+')
args=ap.parse_args()


MdotWind = 1e21
M6 = MdotWind/(1e-6*cgs.Msun/cgs.yr)

vWind    = 1e8
v8       = vWind /1e8

Lwind = 0.5* MdotWind * vWind**2

rho      = 2.35418624e-22

RFE = np.sqrt(Lwind/(2*np.pi*rho*vWind**3))
tFE = RFE /  vWind

def rads_Adi(times):
    return 0.88*(0.5*MdotWind*vWind**2/rho)**(1/5) * times**(3/5)

def radFunc(Rs):
    return 0.91*Rs/(1+0.3*Rs/RFE)**(1/3) * (1+0.49*(Rs/RFE))
def getRadsAdi(time):
    if(time == 0):
        return 0, 0
    Rm = rads_Adi(time)
    R0 = Rm/10
    R1 = Rm*10

    while(True):
        Rm = (R0 + R1)/2
        F = radFunc(Rm)
        if(np.abs(F-vWind*time)/(vWind*time) < 1e-4):
            return 0.91*Rm/(1+0.3*Rm/RFE)**(1/3), Rm
        if(F > vWind*time):
            R1 = Rm
        else:
            R0 = Rm
    



times = np.zeros(len(args.f))
radS  = np.zeros(len(args.f))
radSW = np.zeros(len(args.f))
idxmax = len(args.f)
for i, f in enumerate(args.f):
    time = 0
    with open(f, 'r') as fil :
        for lines in fil:
            if lines[0] == "#":
                line = lines.split('=')
                if line[0][:6] == "# time":
                    time = float(line[1])
            else:
                break

    data = np.loadtxt(f)
    rads = data[:,0]
    velx = data[:,2]
    Temp = data[:,3]/(cgs.kb*data[:,1]/cgs.mH)

    times[i] = time
    if np.sum(velx> 1.0) > 0 : 
        radS[i] = np.max(rads[velx > 1.0])
        if(np.max(rads[velx>1.0]) == rads[-1]):
            idxmax = i
            break
        radSW[i] = np.min(rads[velx < vWind*0.9])
    print(i,f,times[i], radSW[i], radS[i])

tana = np.linspace(0, np.max(times), 1000)
RSa = np.zeros(tana.shape)
RSWa = np.zeros(tana.shape)
for i, t in enumerate(tana):
    RSWa[i], RSa[i] = getRadsAdi(t)
plt.plot(times[:idxmax], radS[:idxmax], c = 'navy')
plt.plot(tana, RSa, ls = '--', c = 'navy')


plt.plot(times[:idxmax], radSW[:idxmax], c = 'orange')
plt.plot(tana, RSWa, ls = '--', c = 'orange')
plt.yscale('log')
plt.show()



