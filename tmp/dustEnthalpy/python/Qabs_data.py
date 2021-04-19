import sys
import numpy as np
import CGS as cgs
import matplotlib.pyplot as plt
from scipy.integrate import ode
from radiationFields import bbody_fr
import os
from mpmath import polylog

GraphData = np.loadtxt('dustAborption_graphite.txt')
GraphData =  GraphData.reshape(30, 380, 5)
GraphFreq =  GraphData[0,:,1]
GraphSize =  GraphData[:,0,0]
GraphQabs =  GraphData[:,:,2]

SilicData =  np.loadtxt('dustAborption_silicone.txt')
SilicData =  SilicData.reshape(30, 380, 5)
SilicFreq =  SilicData[0,:,1]
SilicSize =  SilicData[:,0,0]
SilicQabs =  SilicData[:,:,2]



def getQabs_nu(ida, nu, graphite):
    if(graphite):
        freq = GraphFreq
        qabs = GraphQabs
    else:
        freq = SilicFreq
        qabs = SilicQabs

    if(nu < freq[0]):
        if graphite:
            return 0.01*(cgs.c/nu/1e-2)**(-2)*GraphSize[ida]*1e4
        else:
            return 0.014*(cgs.c/nu/1e-2)**(-2)*SilicSize[ida]*1e4

    elif(nu >= freq[-1]):
        idn = -5
        idnp = -1
        
        f0 = freq[idn]
        f1 = freq[idnp]
        Q0 = qabs[ida, idn]
        Q1 = qabs[ida, idnp]
        power = (np.log(Q1)-np.log(Q0))/(np.log(f1) - np.log(f0))
        return Q1 * (nu/f1)**power
    else:   
        idn = np.abs(freq-nu).argmin()
        if(freq[idn] > nu):
            idn = idn-1
        idnp = idn + 1
    
    dQ  = qabs[ida,idnp]-qabs[ida,idn]
    dnu =  freq[idnp] - freq[idn]
    return qabs[ida, idn] + (nu - freq[idn])*dQ/dnu 

def getQabs(agrain, nu, graphite):
    if(graphite) :
        size = GraphSize
    else:
        size = SilicSize
    if(agrain < size[0]) :
        Qp = getQabs_nu(1, nu, graphite)
        Qm = getQabs_nu(0, nu, graphite)
        da = size[1] - size[0]
        a0 = size[0]
    
    elif(agrain >= size[-1]):
        Qp = getQabs_nu(-1, nu, graphite)
        Qm = getQabs_nu(-2, nu, graphite)
        da = size[-1] - size[-2]
        a0 = size[-2]
    
    else:
        ida = np.abs(size-agrain).argmin()
        if(size[ida] > agrain):
            ida = ida - 1
        Qp = getQabs_nu(ida+1, nu, graphite)
        Qm = getQabs_nu(ida, nu  , graphite)
        da = size[ida+1] - size[ida]
        a0 = size[ida]

    return Qm + (agrain-a0)*(Qp-Qm)/da

def qemintegral(freq,Qem, T, agrain, graphite, norm):
    return getQabs(agrain, freq*1e15, graphite)*bbody_fr(freq*1e15, T)*norm

def solPowerlaw(freq0, freq1, T, agrain, graphite, norm):
    expfact = cgs.h/cgs.kb/T
    exp0 = expfact*freq0
    exp1 = expfact*freq1

    expm0 = np.exp(-exp0)
    expm1 = np.exp(-exp1)
    
    
    int0 = (exp0**5 * np.log(1-expm0) - 5*exp0**4*polylog(2, expm0) - 20*exp0**3 * polylog(3, expm0) - 60*exp0**2*polylog(4,expm0) - 120 * exp0 * polylog(5, expm0) - 120 * polylog(6, expm0))/expfact**6
   
    int1 = (exp1**5 * np.log(1-expm1) - 5*exp1**4*polylog(2, expm1) - 20*exp1**3 * polylog(3, expm1) - 60*exp1**2*polylog(4,expm1) - 120 * exp1 * polylog(5, expm1) - 120 * polylog(6, expm1))/expfact**6
    if graphite:
        inte =  0.01*1e4*(cgs.c/0.01)**-2 * (int1 - int0)*2*cgs.h/cgs.c**2 * norm * agrain
        return inte
    else:
        return 0.014*1e4*(cgs.c/0.01)**-2 * (int1 - int0)*2*cgs.h/cgs.c**2 * norm * agrain

def calcQemAve(T, agrain, fmin, fmax, fpeak, graphite):
    norm =1/agrain/(cgs.sigm_sb*T**4/np.pi)
    start = solPowerlaw(1, fmin, T, agrain, graphite, norm)
    if(start < 0):
        start = 0
    #integrator = ode(qemintegral).set_integrator('dopri5', max_step = 1000000000000000)
    integrator = ode(qemintegral).set_integrator('dop853', nsteps = 10000000000000, atol = 1e-10, rtol = 1e-10)
    #integrator = ode(qemintegral).set_integrator('lsoda', max_step = 1000000000000000)
    integrator.set_initial_value(start/1e15, fmin/1e15).set_f_params(T,agrain, graphite, norm)
    if graphite:
        freqs = GraphFreq
    else:
        freqs = SilicFreq
    for i in range(1, len(freqs)):
        Qem_ave = integrator.integrate(freqs[i]/1e15)*agrain*1e15
        if(not integrator.successful()):
            print(T, fmin, integrator.t*1e15, Qem_ave, integrator.successful())
            sys.exit(0)
    return Qem_ave[0]

def tabulateQemAve():
    Temps = np.logspace(0,5,1000)
    with open('QemAve.dat', 'w') as f:
        for agrain in GraphSize:
            print(agrain)
            for Teff in Temps:
                fpeak = Teff * 5.879e10
                Qem_ave_graph = calcQemAve(Teff, agrain, GraphFreq[0], max(fpeak, GraphFreq[0])*100, fpeak, 1)
                Qem_ave_silic = calcQemAve(Teff, agrain, SilicFreq[0], max(fpeak, GraphFreq[0])*100, fpeak, 0)
                f.write('{} {} {} {}\n'.format(agrain, Teff, Qem_ave_graph, Qem_ave_silic))
if ( not os.path.isfile('QemAve.dat')):
    print("Qalculating <Qem> for all dust grains. This will take time")
    tabulateQemAve()

QemAveData = np.loadtxt('QemAve.dat')
QemAveData = QemAveData.reshape(30,1000,4)
QemAveSize = QemAveData[:,0,0]
QemAveTemp = QemAveData[0,:,1]
QemAveQemGraph  = QemAveData[:,:,2]
QemAveQemSilic  = QemAveData[:,:,3]

def get_QemAve_T(T, ida, graphite):
    if graphite:
        AveQem = QemAveQemGraph  
    else:
        AveQem = QemAveQemSilic
    if T > QemAveTemp[-1]:
        return AveQem[ida, -1]
    idT = np.abs(QemAveTemp-T).argmin()
    Tm = QemAveTemp[idT]
    if(Tm >= T):
        idT = idT-1
        Tm = QemAveTemp[idT]
    Tp = QemAveTemp[idT+1]
    dT =  Tp-Tm
    Qm = AveQem[ida, idT]
    Qp = AveQem[ida, idT+1]
    dQ  = Qp-Qm
    return Qm + (T-Tm)*dQ/dT 

    

def get_QemAve(T, agrain, graphite):
    if(agrain < QemAveSize[0]) :
        Qp = get_QemAve_T(T, 1, graphite)
        Qm = get_QemAve_T(T, 0, graphite)
        dQ = Qp - Qm
        a0 = QemAveSize[1]
        a1 = QemAveSize[0]
        da  = a1 - a0  
    
    elif(agrain >= QemAveSize[-1]):
        Qp = get_QemAve_T(T, -1, graphite)
        Qm = get_QemAve_T(T, -2, graphite)
        dQ = Qp - Qm
        a0 = QemAveSize[-2]
        a1 = QemAveSize[-1]
        da  = a1 - a0  
    else:
        ida = np.abs(QemAveSize-agrain).argmin()
        if(QemAveSize[ida] > agrain):
            ida = ida - 1
        Qp = get_QemAve_T(T, ida+1, graphite)
        Qm = get_QemAve_T(T, ida,   graphite)
        dQ = Qp - Qm
        a0 = QemAveSize[ida]
        a1 = QemAveSize[ida+1]
        da  = a1 - a0  

    return Qm + (agrain-a0) * dQ/da
if __name__ == '__main__':
    lams = np.logspace(np.log10(0.001*cgs.A), 1, 10000)
    Qabsa = np.zeros(lams.shape)
    agrains = [0.001*1e-4, 0.01e-4, 0.1e-4, 1e-4]
    for grain in agrains:
        for i in range(len(lams)):
            Qabsa[i] = getQabs(grain,cgs.c/lams[i], 1)/1e-1
        plt.plot(cgs.h*cgs.c/lams/cgs.eV, Qabsa)
    plt.axvline(GraphFreq[0]*cgs.h/cgs.eV)
    plt.axvline(GraphFreq[-1]*cgs.h/cgs.eV)
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    
    for grain in agrains:
        for i in range(len(lams)):
            Qabsa[i] = getQabs(grain,cgs.c/lams[i], 0)/1e-1
        plt.plot(cgs.h*cgs.c/lams/cgs.eV, Qabsa)
    plt.axvline(GraphFreq[0]*cgs.h/cgs.eV)
    plt.axvline(GraphFreq[-1]*cgs.h/cgs.eV)
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    

    agrains = [0.001*1e-4, 0.01e-4, 0.1e-4, 1e-4]
    temps = np.logspace(1, 5, 10000)
    Qave = np.zeros(temps.shape)
    for ia in range(len(agrains)):
        for it in range(len(temps)):
            Qave[it] = get_QemAve(temps[it], agrains[ia], 1)
        
        plt.plot(temps, Qave/(agrains[ia]/1e-4))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(10,1e3)
    plt.ylim(1e-4,1)
    plt.show()

    for ia in range(len(agrains)):
        for it in range(len(temps)):
            Qave[it] = get_QemAve(temps[it], agrains[ia], 0)
        
        plt.plot(temps, Qave/(agrains[ia]/1e-4))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(10,1e3)
    plt.ylim(1e-4,1)
    plt.show()



