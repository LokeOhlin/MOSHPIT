import numpy as np
import matplotlib.pyplot as plt
import CGS as cgs
import Qabs_data
import radiationFields

def print2DMatr(Matr):
    shape = Matr.shape
    Ni = shape[0]
    Nj = shape[1]
    for i in range(Ni):
        st = ''
        for j in range(Nj):
            st = st + '{:.2e}, '.format(Matr[Ni-1-i, Nj-1-j])
        print(st)


def get_dUdt(agrain, T_bin, dEdt_abs, graphite):
    #dEdt_em  =  8e-7*(agrain/(1e-5))*4*np.pi*agrain**2*cgs.sigm_sb*T_bin**6
    dEdt_em = 4*np.pi*agrain**2 * Qabs_data.get_QemAve(T_bin, agrain, graphite) * cgs.sigm_sb * T_bin**4
    return dEdt_abs - dEdt_em

def get_TSS(agrain, dEdt_abs, graphite):
    if graphite:
        dEdt_em  = 8e-7*(agrain/(1e-5))*4*np.pi*agrain**2*cgs.sigm_sb
    else: 
        dEdt_em  = 1.3e-6*(agrain/(1e-5))*4*np.pi*agrain**2*cgs.sigm_sb
    return (dEdt_abs/dEdt_em)**(1/6)

rho_c   = [3.5, 2.2]
mmol_c  = [cgs.mH*(28.086+2*24.305+4*15.999)/7,cgs.mH*12.011]
def get_UfromT(agrain, T, graphite):
    N =   1/3*4*np.pi*agrain**3 * rho_c[graphite]/mmol_c[graphite]
    if graphite:
        Uatom = 4.15e-22*T**3.3 / (1+6.51e-3*T+1.5*1e-6*T*T+8.3*1e-7*T**(2.3))
        return (1-2/N)*N*Uatom
    else:
        if (T < 50):
            U = 1.4*1e3/3 * T**3
        else:
            U = 1.4*1e3/3* (50)**3
            if (T < 150): 
                U += 2.2*1e4 / (2.3) * ( T**(2.3) - 50**(2.3))
            else: 
                U += 2.2*1e4 / (2.3) * ( 150**(2.3) - 50**(2.3))
                if (T < 500):
                    U += 4.8*1e5/(1.68) * (T**(1.68) - 150**(1.68))
                else:
                    U += 4.8*1e5/(1.68) * (500**(1.68) - 150**(1.68))
                    U += 3.41*1e7*(T - 500)
        return U * (1-2/N) * 4/3*np.pi*agrain**3 
            

def get_TfromU(agrain, U, graphite, T0 = 0, T1 = 1e4, tolT = 1e-4, tolU = 1e-8):
    # Midpoint search
    Tm = (T0+T1)/2
    dT = T1 - T0
    # accuracy 
    while True:
        Um = get_UfromT(agrain,Tm, graphite)
        if(Um == U):
            return Tm
        if(Um > U):
            T1 = Tm
        else:
            T0 = Tm
    
        Tm = (T0+T1)/2
        dT = T1 - T0
        if(abs(dT/Tm) < tolT or abs(Um-U)/U < tolU):
            return Tm

def get_dUdt_abs(agrain, lam_min, lam_max, graphite):
    wavelength = np.logspace(np.log10(lam_min), np.log10(lam_max), 100)
    dwlen = wavelength[1:] - wavelength[:-1]
    dEdt_abs = 0
    for i in range(1,len(wavelength)):
        freq = cgs.c/wavelength[i]
        dEdt_abs += Qabs_data.getQabs(agrain, freq, graphite)*radiationFields.galacticRF(wavelength[i])*dwlen[i-1]
    dEdt_abs = dEdt_abs *np.pi*agrain**2 * cgs.c
    return dEdt_abs 
def get_dE_sporadic(agrain, lam_min, lam_max, graphite):
    wavelength = np.logspace(np.log10(lam_min), np.log10(lam_max), 100)
    dwlen = wavelength[1:] - wavelength[:-1]
    dEdt_abs = 0
    for i in range(1,len(wavelength)):
        freq = cgs.c/wavelength[i]
        dEdt_abs += wavelength[i]*Qabs_data.getQabs(agrain, freq, graphite)*radiationFields.galacticRF(wavelength[i])*dwlen[i-1]
    dEdt_abs = dEdt_abs *np.pi*agrain**2 /cgs.h
    return dEdt_abs 

def getDist(agrain, N, graphite, Emax_init):
    lam_trunc = 0.1 # above this wavelenght we assume that absorption is continous
    lam_min = 800*cgs.A # max photon energy (for integration purposes)
    dUdt_abs = get_dUdt_abs(agrain, lam_trunc , 1000, graphite) # continous absorption
    A_matr = np.zeros((N,N))
    B_matr = np.zeros((N,N))
    Pi = np.zeros(N)
    Xi = np.zeros(N) 
    Emax = Emax_init
    #Emin = get_UfromT(agrain,2, graphite)
    Emin = get_UfromT(agrain,get_TSS(agrain, dUdt_abs, graphite), graphite)
    while True:
        Ubins = np.logspace(np.log10(Emin), np.log10(Emax), N)
        delU  = np.zeros(Ubins.shape)
        # define midpoint
        Ul = (Ubins[1:]+Ubins[:-1])/2
        # bin widths
        delU[1:-1] = Ul[1:] - Ul[:-1]
        delU[0]  = 2*(Ul[0]-Ubins[0])
        delU[-1] = 2*(Ubins[-1]-Ul[-1])
        Tbins = np.zeros(Ubins.shape)
        delT  = np.zeros(Ubins.shape)
        for i in range(N):
            Tbins[i] = get_TfromU(agrain, Ubins[i], graphite)
            if( (i != 0) and (i != N-1)):
                delT[i]  = get_TfromU(agrain, Ul[i], graphite)
                delT[i] -= get_TfromU(agrain, Ul[i-1], graphite)
        delT[0]   = get_TfromU(agrain, Ubins[0]+delU[0]/2, graphite)
        delT[0]  -= get_TfromU(agrain, Ubins[0]-delU[0]/2, graphite)
        delT[-1]  = get_TfromU(agrain, Ubins[-1]+delU[-1]/2, graphite)
        delT[-1] -= get_TfromU(agrain, Ubins[-1]-delU[-1]/2, graphite)
        dUdt = np.zeros(Tbins.shape)
        for i in range(N):
            dUdt[i] = get_dUdt(agrain, Tbins[i], dUdt_abs, graphite)
        dUdt[0] = 0
        # populate of diagonal A matrix
        for i in range(N-1):
            A_matr[i,i+1] = -dUdt[i+1]/delU[i+1]
            #if(dUdt[i] > 0):
            #    A_matr[i+1,i] = dUdt[i+1]/delU[i+1]
            #else:
            #    A_matr[i-1,i] = -dUdt[i+1]/delU[i+1]

            for j in range(i):
                lam = cgs.h*cgs.c/(Ubins[i] - Ubins[j])
                if lam > lam_trunc:
                    continue
                freq = cgs.c/lam
                Qabs = Qabs_data.getQabs(agrain, freq, graphite)
                ul = radiationFields.galacticRF(lam)
                A_matr[i,j] =  Qabs * np.pi*agrain**2*ul*delU[i]*lam**3/(cgs.h**2*cgs.c)

            #add last contribution to last bin

            lam = cgs.h*cgs.c/(Ubins[-1] - Ubins[i])

            if lam < lam_trunc:
                freq = cgs.c/lam
                Qabs = Qabs_data.getQabs(agrain, freq, graphite)
                ul = radiationFields.galacticRF(lam)
                A_matr[-1,i] = Qabs*np.pi*agrain**2*ul*delU[-1]*lam**3/(cgs.h**2*cgs.c)
            
            lam_i = cgs.h*cgs.c/(Ubins[-1] + delU[-1]/2 - Ubins[i])
            A_matr[-1,i] += get_dE_sporadic(agrain, lam_min, lam_i, graphite)

        # now do diagonal elements
        for i in range(N):
            A_matr[i,i] = 0
            A_matr[i,i] = - np.sum(A_matr[:,i])
        
        # Calculate B matrix
        B_matr[-1,:] = A_matr[-1,:]
        for i in range(N-1):
            B_matr[N-2-i,:] = B_matr[N-1-i,:] + A_matr[N-2-i,:]
        
        # Debug, print A and B
        #print('A')
        #print2DMatr(A_matr)
        #print('B')
        #print2DMatr(B_matr)


        # Calculate Pi/P0   
        xi = np.zeros(Ubins.shape)
        xi[0] = 1
        for i in range(1,N):
            xi[i] = (1/A_matr[i-1,i])*np.sum(B_matr[i,:i]*xi[:i])
        norm = np.sum(xi)
        print('norm = ', norm)
        xi = xi/norm
        if(xi[-1] < 1e-13):
            return Tbins, delT, xi
            break
        
        Emax = Emax * 1.5
        print("redo", Emax)



if __name__ == '__main__':
    abins = np.logspace(np.log10(3.5e-8), np.log10(1e-4), 3)
    abins[0] = 3.5e-8
    abins[1] = 3.5e-8
    abins[2] = 1.75e-6
    #abins[1] = 15*cgs.A
    #abins[2] = 50*cgs.A
    for i in range(len(abins)):
        print(abins[i], get_TfromU(abins[i], 13.6*cgs.eV,0))
        #Tbins, delT, dist = getDist(abins[i],100, 0,get_UfromT(abins[i], 1e4,0))
        #plt.plot(Tbins, dist/np.log((Tbins+delT/2)/(Tbins-delT/2)))
        Tbins, delT, dist = getDist(abins[i],100, 1,get_UfromT(abins[i], 5000, 1))
        plt.plot(Tbins, dist/np.log((Tbins+delT/2)/(Tbins-delT/2)),ls = '--')
    for i in range(len(Tbins)):
        print(Tbins[i], dist[i]) 
    plt.xscale('log')
    plt.xlim(5,1e4)
    plt.ylim(1e-12,10)
    plt.yscale('log')
    plt.show()
