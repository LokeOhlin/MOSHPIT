#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "dustRadiation.h"
#include "cgeneral.h"
#include "dust.h"

//
int NUbins;

// Energy bins and their widths
double *Ubins  = NULL;
double *dUbins = NULL;
double *Ul     = NULL;
// Transisiton matrices A and B
double *Amatr = NULL;
double *Bmatr = NULL;

double sigma_sb = 5.6704e-5;

double get_Em(double agrain, double piaa, int ida, double temp, int graphite){ 
    if(getQemAve(agrain, temp, graphite, ida, -1) < 0){
        printf("!!!!!! %.4e \n", getQemAve(agrain, temp, graphite, ida, -1));
    }
    return 4 * piaa * getQemAve(agrain, temp, graphite, ida, -1) * sigma_sb * pow(temp, 4);
}

// Silicate U(T) is based of a segmented function of the heat capacity, integrated to T
// Store the values of the integral up to segmented values 
double Ulim1 = 5.833e7; // U(50 K)
double Ulim2 = 9.486e8; // U(150 K)
double Ulim3 = 9.432e9; // U(500 K)


double getUfromT(double Natom, double volgrain, double temp, int graphite){
    double Uatom;
    // Graphite grains : functional expression
    if(graphite){
        Uatom = 4.15e-22 * pow(temp, 3.3) / (1 + 6.51e-3*temp + 1.5e-6*temp*temp+8.3e-7*pow(temp,2.3));
        return (Natom-2)*Uatom;
        //return (1-2/Natom)*Natom*Uatom;
    } 
    // Silicate grains : Piecewise function here of heat capacity integrated

    if(temp < 50){
        Uatom = 4.6667e2 * temp*temp*temp;
    } else if(temp < 150){
        Uatom = Ulim1 + 9.5652e3 * (pow(temp, 2.3) - 8.0841e3);
    } else if(temp < 500){
        Uatom = Ulim2 + 2.8571e5 * (pow(temp, 1.68) - 4.5272e3);
    } else {
        Uatom = Ulim3 + 3.41e7 * (temp - 500);
    } 
    return Uatom * (1-2/Natom) * volgrain;
}

double getTfromU(double Natom, double volgrain, double Ugrain, int graphite){
    if(graphite){
        // The graphite function appears non-invertable. use midpoint rootfinder
        double tolT = 1e-4;
        double tolU = 1e-8;
        double Umid;
    double Tmax = 1e5;
        double Tmin = 0.;
        double Tmid = (Tmax + Tmin)*0.5;
        double dT = Tmax - Tmin;
        
        //remove all energy independent factors
        double Ured = Ugrain/((Natom-2)*4.15e-22);
        while(1){
            Umid = pow(Tmid, 3.3) / (1 + 6.51e-3*Tmid + 1.5e-6*Tmid*Tmid + 8.3e-7*pow(Tmid, 2.3));
            if(Umid == Ured){
                return Tmid;
            }    
            if(Umid > Ured){
                Tmax = Tmid;
            } else {
                Tmin = Tmid;
            }
            Tmid = (Tmax + Tmin)*0.5;
            dT = Tmax - Tmin;
            if((fabs(dT/Tmid) < tolT) || (fabs(Umid-Ured)/Ured < tolU)){
                return Tmid;
            }
        }
    }   
    //Silicates we can do exactly
    //Take out grain dependent part
    double Ured = Ugrain/((1-2/Natom)*volgrain);
    
    
    if(Ured < Ulim1){
        return pow(Ured*2.1429e-3, 1./3.);
    } else if(Ured < Ulim2){
        return pow((Ured-Ulim1)*1.045455e-4 +  8.0841e3, 1./2.3);
    } else if(Ured < Ulim3){
        return pow((Ured-Ulim2)*3.5e-6 +  4.5272e3, 1./1.68);
    } else {
        return (Ured - Ulim3)*2.9325e-8 + 500;
    }
}




// get stable state temperature based on continous absorbption as emission of photons
// TODO: include gas temperature and heating by collisions
double StableState_temperature(double agrain, double piaa, int ida, int graphite, double continous_absorption, double Tmin, double Tmax){
    double continous_emission;
    double dUmid;
    //Mid point search TODO: change to better method ( newton raphson ? )
    double Tmid = (Tmax + Tmin)/2.;
    double dT = Tmax - Tmin;

    //Tolerance TODO: make changeable by user
    double tolT_SS = 1e-3;
    
    // Loop until accuracy has been found, or TODO: max steps have been reached
    while(1){
        // Get emmission of grain for mid temperature
        continous_emission = get_Em(agrain, piaa, ida, Tmid, graphite);
        // dU for the mid temperature
        dUmid = continous_absorption - continous_emission;
        if(dUmid == 0){
            return Tmid;
        }
        // if less than zero, then cooling to strong. decrease temperature
        if(dUmid < 0){
            Tmax = Tmid;
        } else {
            Tmin = Tmid;
        }
        Tmid = (Tmax + Tmin)*0.5;
        dT   =  Tmax - Tmin;
        // 
        if(fabs(dT) < tolT_SS){
            return Tmid;
        }
    }
}


int initTemperatureDist(int Nbins){
    NUbins = Nbins;
    Ubins  = (double *) malloc(Nbins * sizeof(double));
    dUbins = (double *) malloc(Nbins * sizeof(double));
    Ul     = (double *) malloc(Nbins * sizeof(double));
    Amatr  = (double *) malloc(Nbins * Nbins * sizeof(double));
    Bmatr  = (double *) malloc(Nbins * Nbins * sizeof(double));
}



int getTemperatureDist(int ibin, int graphite, double Tmax, double *Ts, long double *Ps){
    //  Main call for calculating the temperature distribution of a given dust grain size
    //  In :
    //          double agrain  -  (spherical) size of the dustgrain in cm
    //          int graphite   -  flag for type of dust (1 = graphite grains, 0 = silicate grains)
    //          double Tmax    -  Maximum temperature to distribution. Must be set high enough  
    //                            to not truncate distribution
//  out: 
    //          double *Ts     -  Array of temperature values to be written to. 
    //          double *Ps     -  Array of fraction of dust grains in each bin to be written to
    //
    //  NOTE: both out arrays need to allocated beforehand
    
    // Start by setting energy bins. We need maximum and minimum energy (Emax and Emin), where
    // Emax is user set and should be large enough to not truncate distrubtion too excessively
    // and Emin is set such that the cooling and heating from frequent low energy photons,
    // eg. the processes we consider to be continous, is in equilibrium
    int iUbin, fUbin, fUbinp, jUbin, idx, idxp;
    double AiiSum, freq_if, freq_inf, Qabs, radEdens;
    double agrain = agrains[ibin];
    double volgrain = volgrains[ibin];
    double Natom = Natoms[ibin];
    double piaa = pi_asquare[ibin];
    int ida_Qabs = ida_tabQabs[ibin];
    int ida_Qem = ida_tabQem[ibin];
    
    
    
    double Emax = getUfromT(Natom, volgrain, Tmax, graphite);
    
    double nuCont = 2997924580.000;
    //double nuCont = 24180876579417.275;
    double ECont = nuCont*planck;
    // Absorbed energy rate from frequent photons (eg low energy photons with nu < nuCont)
    double continous_absorption = getAbsorption_range(ibin, graphite, 0, nuCont);     
    // Get "stable state" temperature for continous processes. initial range between 0.01 and 100... TODO: CHECKIF ALWAYS VALID!
    double Tmin = StableState_temperature(agrain, piaa, ida_Qem,  graphite, continous_absorption, 0.01, 100);
    double Emin = getUfromT(Natom, volgrain, Tmin, graphite);
    // populate Energy bins in log scale
    double delE = exp((log(Emax) - log(Emin))/NUbins);
    //double delE = (Emax - Emin)/NUbins;
    //if(delE < ECont){
    //    printf("Too small bins\n");
    //    return -1;
    //} 
    printf("Tmin = %.8e %.4e %.4e %.4e\n", Tmin, continous_absorption, delE/electronVolt, getRadEdens(delE/planck));
    
    
    Ubins[0] = Emin;
    Ts[0] = Tmin;
    
    for(iUbin = 1; iUbin < NUbins; iUbin++){
        Ubins[iUbin] = Emin*pow(delE,iUbin);
        //Ubins[iUbin] = Emin + delE*iUbin;
        Ts[iUbin] = getTfromU(Natom, volgrain, Ubins[iUbin], graphite);
        // left edge of bin
        Ul[iUbin] = (Ubins[iUbin] + Ubins[iUbin-1])*0.5;
        dUbins[iUbin]    = -Ul[iUbin];
        dUbins[iUbin-1] +=  Ul[iUbin];
        // reset bin
        Ps[iUbin] = 0;
    }
    // fix dUbins at upper and lower limit (remember dUbins[-1] = - (Ubins[-1] + Ubins[-2])/2 still)
    dUbins[0] = 2*(dUbins[0] - Ubins[0]);       
    dUbins[NUbins - 1] = 2*(Ubins[NUbins - 1] + dUbins[NUbins -1]);
    // set matrices to zero
    for(iUbin = 0; iUbin < NUbins*NUbins; iUbin++){
        Amatr[iUbin] = 0;
        Bmatr[iUbin] = 0;
    }
    
    // start populating transition A matrix
    // transition is going from ibin (initial) to fbin (final)
    for(jUbin = 0; jUbin < NUbins-1; jUbin++){
        // continous cooling (since dUdt for Ubins[0] is defined as 0)
        fUbin = jUbin;
        iUbin = jUbin+1; 
        idx = fUbin*NUbins + iUbin;
        Amatr[idx] = (get_Em(agrain, piaa, ida_Qem, Ts[iUbin], graphite) - continous_absorption)/(Ubins[iUbin] - Ubins[fUbin]); //dUbins[iUbin];
        
        printf("%d %.4e %.4e %.4e \n", iUbin, Ts[iUbin], Ubins[iUbin]*planckInv, Ts[iUbin] * 5.879e10);
        // Discreete heating processes from lower energy bins
        for(iUbin = 0; iUbin < fUbin; iUbin++){
            idx = fUbin*NUbins + iUbin;
            // frequency of photons needed to get from ibin to fbin
            freq_if = (Ubins[fUbin] - Ubins[iUbin])*planckInv;         
            //if(freq_if < nuCont){
            //    continue;
            //}
            Qabs = getQabs(agrain, freq_if, graphite, ida_Qabs, -1);
            radEdens = getRadEdens(freq_if);
            double old = Qabs * piaa * radEdens * dUbins[fUbin]*clght*planckInv*planckInv/(freq_if);
            Amatr[idx] =  getUpwardTransition(ibin, graphite, Ubins[iUbin], Ul[iUbin], Ul[iUbin+1], dUbins[iUbin],
                                                              Ubins[fUbin], Ul[fUbin], Ul[fUbin+1], dUbins[fUbin], ECont);
            //if(Amatr[idx] != Amatr[idx] || Amatr[idx] < 0){
            //    printf("%d %d | %.4e %.4e %.4e", fUbin, iUbin, freq_if, Qabs, radEdens);
            //}
        }
        // Add discreete heating from jbin to final bin.
        iUbin = jUbin;
        fUbin = NUbins - 1;
        idx = fUbin*NUbins + iUbin;
        
        // frequency of photons needed to get from iUbin to fUbin
        freq_if = (Ubins[fUbin] - Ubins[iUbin])*planckInv;
        //if(freq_if > nuCont){
        Qabs = getQabs(agrain, freq_if, graphite, ida_Qabs, -1);
        radEdens = getRadEdens(freq_if);
        Amatr[idx] = Qabs * piaa * radEdens * dUbins[fUbin]*clght*planckInv*planckInv/(freq_if);
        Amatr[idx] =  getUpwardTransition(ibin, graphite, Ubins[iUbin], Ul[iUbin], Ul[iUbin+1], dUbins[iUbin],
                                                          Ubins[fUbin], Ul[fUbin], Ul[fUbin+1], dUbins[fUbin], ECont);
        //} else {
        //    Amatr[idx] = 0;
        //}
        
        // Last bin technically goes to infintity. Take this into account
        freq_inf = (Ubins[fUbin] + dUbins[fUbin]*0.5 - Ubins[iUbin])*planckInv;
        Amatr[idx] += getAbsorptionNum_range(ibin, graphite, freq_inf, 1e99); 
        if(Amatr[idx] != Amatr[idx] || Amatr[idx] < 0){
            printf("%d %d | %.4e %.4e | %.4e %.4e %.4e", fUbin, iUbin, freq_inf, getAbsorption_range(ibin, graphite, freq_inf, 1e99), freq_if, Qabs, radEdens);
        }
    }
    // now do the diagonal elements ( eg. everything going out from bins)
    double check;
    for(iUbin = 0; iUbin < NUbins; iUbin++){
        AiiSum  = 0;
        for(fUbin = 0; fUbin < NUbins; fUbin++){
            idx = fUbin*NUbins + iUbin;
            AiiSum -= Amatr[idx];
            if(isfinite(AiiSum)==0){
                printf("%d %d | %.4e", fUbin, iUbin, AiiSum);
            }
        }
        idx = (iUbin-1)*NUbins + iUbin;
        check= fabs(AiiSum + Amatr[idx]);
        printf("%d %.4e %.4e %.4e %.4e\n", iUbin, Ts[iUbin], check, Amatr[idx], check - Amatr[idx]);   
        idx = iUbin*NUbins + iUbin;
        Amatr[idx] = AiiSum;
    }

    // Populate B matrix
    // Start with top row
    fUbin = NUbins - 1;
    for(iUbin = 0; iUbin < fUbin; iUbin++){
        idx = fUbin*NUbins + iUbin;
        Bmatr[idx] = Amatr[idx];
        //if(Bmatr[idx] < 0 || Bmatr[idx] != Bmatr[idx]){
        //    printf("%d %d, %.4e \n", fUbin, iUbin, Amatr[idx]);
        //    return -1;
        //}
    }
    // Sum up the rest
    for(jUbin = 0; jUbin < NUbins-1; jUbin++){
        fUbin  = NUbins - 2 - jUbin;
        fUbinp = fUbin +1;
        for(iUbin = 0; iUbin < fUbin; iUbin++){
            idx  = fUbin  * NUbins + iUbin;
            idxp = fUbinp * NUbins + iUbin;
            Bmatr[idx] = Bmatr[idxp] + Amatr[idx];
            if(Bmatr[idx] < 0 || isfinite(Bmatr[idx])==0){
                printf("%d %d||| %.4e %.4e \n", fUbin, iUbin, Bmatr[idxp], Amatr[idx]);
                return -1;
            }
        }

    }
    //Calculate probablilites
    Ps[0] = 1e-200;
    long double norm = Ps[0];
    for(fUbin = 1; fUbin < NUbins; fUbin++){
        Ps[fUbin] = 0;
        for(iUbin = 0; iUbin < fUbin; iUbin++){
            idx = fUbin*NUbins + iUbin;
            Ps[fUbin] += Bmatr[idx]*Ps[iUbin];
        }
        idx = (fUbin-1)*NUbins + fUbin;
        
        long double Psold = Ps[fUbin];
        Ps[fUbin] = Ps[fUbin]/Amatr[idx];
        if(isfinite(Ps[fUbin])==0 || Ps[fUbin] < 0){
            printf("%.4Le %.4e %.4Le %.4Le\n", Ps[fUbin], Amatr[idx], Psold, Ps[fUbin-1]);
            Psold = 0;
            for(iUbin = 0; iUbin < fUbin; iUbin++){
                idx = fUbin*NUbins + iUbin;
                Psold += Bmatr[idx]*Ps[iUbin];
                printf("%.4Le %.4e %.4Le \n", Psold, Bmatr[idx], Ps[iUbin]);
            }
            return -1;
        }
        norm += Ps[fUbin];
    }
    norm = 1/norm;
    for(fUbin = 0; fUbin < NUbins; fUbin++){
        Ps[fUbin] = Ps[fUbin] * norm;
    }

    if(Ps[NUbins-1] > 1e-13){
        // above truncation treshold. throw warning but continue 
        return -2;
    }
    return 1;
}

