#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "dustRadiation.h"
#include "cgeneral.h"
#include <time.h>
#include "dust.h"

double parsec = 3.086e+18;
double clght  = 29979245800.0; // speed of light
double planck = 6.62507e-27;  // planck constant
double planckInv = 1.509418013696459e+26; // 1/planck
double electronVolt = 1.602e-12;
double boltzmann = 1.38065e-16;

double rho_s = 3.5;
double rho_c = 2.2;


double *agrains    = NULL;
double *volgrains  = NULL;
double *Natoms     = NULL;
double *pi_asquare = NULL;
int *ida_tabQabs   = NULL;
int *ida_tabQem    = NULL;

double bbody_wl(double wavelength, double Teff){
    double expo = planck*clght/(boltzmann*Teff*wavelength);
    if(expo > 500){
        return 0;
    }
    return 2*planck*clght*clght/pow(wavelength,5)/(exp(expo)-1);
}
double Math8_1000L[20] = {0.0008, 0.001 , 0.0012, 0.002, 0.0025, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.015 ,0.02 ,0.03  ,0.04  , 0.06  , 0.08  , 0.1};

double Math8_1000E[20] = {4.33e-01, 3.26e-01, 2.21e-01, 7.12e-02, 1.14e-01, 1.70e-01, 3.06e-01, 4.08e-01, 4.76e-01, 4.51e-01, 4.20e-01, 3.62e-01, 3.16e-01, 1.36e-01, 5.61e-02, 1.39e-02, 5.61e-03, 1.32e-03, 4.20e-04, 2.14e-04};

double interp8_1000(double wavelength){
    int idx = sequentialSearch(wavelength, Math8_1000L, 20);
    if(idx == 19){
        return Math8_1000E[idx];
    }
    double E0 = Math8_1000E[idx];
    double E1 = Math8_1000E[idx+1];

    double L0 = Math8_1000L[idx];
    double L1 = Math8_1000L[idx+1];

    return E0 + (wavelength - L0)*(E1-E0)/(L1-L0);
}


double galacticRF(double wavelength){
    double u_lam = 0;
    if(wavelength < 0.0912e-4){
        return 0;
    } else if(wavelength < 0.11e-4) {
        u_lam += 38.57e4 * pow(wavelength/1e-4 , 3.4172);
    } else if(wavelength < 0.134e-4) {
        u_lam += 2.045e2;
    } else if(wavelength < 0.246e-4) {
        u_lam += 7.115 * pow(wavelength/1e-4, -1.6678);
    } else if(wavelength < 8e-4){
        u_lam += 1e-14 * 4 * M_PI * bbody_wl(wavelength,7500);
        u_lam += 1e-13 * 4 * M_PI * bbody_wl(wavelength,4000);
        u_lam += 4e-13 * 4 * M_PI * bbody_wl(wavelength,3000);
    
    } else if(wavelength < 0.1){
        u_lam += interp8_1000(wavelength);
    }
    
    u_lam = u_lam + 4* M_PI * bbody_wl(wavelength, 2.9);
    u_lam = u_lam/clght;
    return u_lam;
    
}

double starRF(double wavelength, double Lstar, double Tstar, double dist){
    double sigma_sb = 5.6704e-5;
    double Rstar_SQ = Lstar/(4*M_PI*sigma_sb*Tstar*Tstar*Tstar*Tstar);
    double B_lam = bbody_wl(wavelength, Tstar);
    double flux_lam = B_lam * Rstar_SQ/(dist*dist);
    return flux_lam/clght;
}


int writeTempdist(int ibin, double *Ts, long double *Ps, int graphite, int NTbins){
    char filename[20];
    int iEbin;
    if(graphite){
        sprintf(filename, "dustTemp_g_%.2e", agrains[ibin]);
    } else {
        sprintf(filename, "dustTemp_s_%.2e", agrains[ibin]);
    }
    FILE *fptr = fopen(filename,"w"); 
    // get statble state temperature
    double absorbed = getAbsorption(ibin, graphite);
    double TSS  =  StableState_temperature(agrains[ibin], pi_asquare[ibin], ida_tabQem[ibin], graphite, absorbed, 0.1, 10000);
    //printf("TSS : %.4e %.4e \n", absorbed, get_Em(agrains[ibin], pi_asquare[ibin], ida_tabQem[ibin], TSS, graphite));
    fprintf(fptr, "#TSS=%.8e \n",TSS); 
    for(iEbin = 0; iEbin < NTbins; iEbin++){
        fprintf(fptr, "%.8e %.8Le \n", Ts[iEbin], Ps[iEbin]);
    }
    fclose(fptr);
    return 1;
}

int main(){
    int ierr, iEbin, iint;
    ierr = loadDustRadiationTables();
    // number of radiation bins
    int Nbins = 1000;
    // edges 
    int NbinsE = Nbins +1;
    // 
    int Nbins_lower = NbinsE/2; // integer division;
    int Nbins_upper = NbinsE - Nbins_lower;
    double *Ebins, *radData;
    
    Ebins = (double *) malloc(NbinsE * sizeof(double));
    radData = (double *) malloc(2 * Nbins * sizeof(double));

    double Emax = 13.6*electronVolt;
    double Emin = 0.00001*electronVolt;
    double Emid = 0.5*electronVolt;

    //Half of the energy bins are set between Emid and Emax
    //the other between Emin and Emid
    // Upper half;
    double delE = exp((log(Emax) - log(Emid))/Nbins_upper);
    for(iEbin = Nbins_lower; iEbin < NbinsE; iEbin ++){
        Ebins[iEbin] = Emid * pow(delE, iEbin - Nbins_lower);
    }
    // lower half;
    delE = exp((log(Emid) - log(Emin))/Nbins_lower);
    for(iEbin = 0; iEbin < Nbins_lower; iEbin ++){
        Ebins[iEbin] = Emin * pow(delE, iEbin);
    }
    
    // initialize internal bins;
    ierr = initDustRadiation(Nbins, Ebins);
    
    // number of steps going into 
    int Nint = 100;
    double Ener, freq, ulam, ufreq, delf;
    double wavelength, Etot, Ntot;
    double Tstar, Lstar, dist;
    int galactic = 0;
    if(galactic == 0){
        Tstar = 25320.5;
        Lstar = 2.138829126e+37;
        dist  = 0.03 * parsec;
    }
    
    
    
    for(iEbin = 0; iEbin < Nbins; iEbin++){
        delE = (Ebins[iEbin+1] - Ebins[iEbin])/Nint;
        delf = delE*planckInv;
        Etot = 0;
        Ntot = 0;
        for(iint = 0; iint < Nint; iint++){
            Ener = Ebins[iEbin] + delE * (iint + 0.5);
            freq = Ener*planckInv;
            wavelength = clght/freq;
            if(galactic){
                ulam = galacticRF(wavelength);
            } else {
                ulam = starRF(wavelength, Lstar, Tstar, dist);
            }
            ufreq = ulam * clght/freq/freq;
            Etot = Etot + ufreq*delf;
            Ntot = Ntot + ufreq*delf/Ener;
        }
        
        radData[iEbin] = Ntot*clght;
        if(Ntot > 0){
            radData[Nbins + iEbin] = Etot/Ntot;
        } else {
            radData[Nbins + iEbin] = 0;
        }
    }
    ierr = setRadiationBins(radData, 1.0, 1.0);
    int NTbins = 500;
    
    
    int ngrains = 100;
    int isilicone = 50;
    int ibin, iabin;
    int graphite;
    agrains = (double *) malloc(ngrains*sizeof(double));
    Natoms  = (double *) malloc(ngrains*sizeof(double));
    volgrains  = (double *) malloc(ngrains*sizeof(double));
    pi_asquare  = (double *) malloc(ngrains*sizeof(double));
    ida_tabQabs  = (int *) malloc(ngrains*sizeof(int));
    ida_tabQem  = (int *) malloc(ngrains*sizeof(int));
    
    double amin = 3.5e-8;
    double amax = 1.0e-4;
    double da = exp((log(amax) - log(amin))/(0.5*ngrains));
    double rho, matom;
    for(ibin = 0; ibin < ngrains; ibin++){
        if(ibin < isilicone){
            iabin = ibin;
            graphite = 1;
            rho = rho_c;
            matom = 1.675e-24 * 12.011;
        }else{
            iabin = ibin - isilicone;
            graphite = 0;
            rho = rho_s;
            matom = 1.675e-24 * (28.086+2*24.305+4*15.999)/7.;
        }
        agrains[ibin] = amin*pow(da, iabin);
        volgrains[ibin] = 4.*M_PI*pow(agrains[ibin],3.)/3.;
        Natoms[ibin] = volgrains[ibin]*rho/matom; 
        pi_asquare[ibin] = M_PI*agrains[ibin]*agrains[ibin];
        ida_tabQabs[ibin] = getQabs_ida(agrains[ibin], graphite);
        ida_tabQem[ibin]  = getQemAve_ida(agrains[ibin], graphite);
    }
    
    ierr = initTemperatureDist(NTbins);
    double *Ts;
    long double *Ps;
    Ts = (double *) malloc(NTbins * sizeof(double));
    Ps = (long double *) malloc(NTbins * sizeof(long double));
    
    clock_t start = clock();
    double TSS;
    double abs;
    double Tmax;
    //for(ibin = 0; ibin < ngrains; ibin++){
    //    if(ibin > isilicone){
    //        graphite = 0;
    //    } else {
    //        graphite = 1;
    //    }
    //    abs = getAbsorption(ibin, graphite);
    //    TSS  =  StableState_temperature(agrains[ibin], pi_asquare[ibin], ida_tabQem[ibin], graphite, abs, 0.1, 10000);  
    //    Tmax = fmin(TSS*pow(10,fmax(fmin(-3-log10(agrains[ibin]),3),1)), 4e4);
    //    printf("%.4e %.4e %.4e \n", agrains[ibin], TSS, Tmax);
    //    ierr = getTemperatureDist(ibin, graphite, Tmax, Ts, Ps);
    //    printf("last = %.4e %.4e \n", Ts[NTbins -1], Ps[NTbins-1]);
    //    ierr = writeTempdist(ibin, Ts, Ps, graphite, NTbins);
    //}
    clock_t end = clock();
    double cpu_time_used = ((double) (end - start))/CLOCKS_PER_SEC;
    printf("time per call: %f \n", cpu_time_used/ngrains);
    printf("time for 1000 cells with 40 bins: %f \n", cpu_time_used*1e3*40/100);
    ibin = 40;
    graphite = 1; 
    abs = getAbsorption(ibin, graphite);
    TSS  =  StableState_temperature(agrains[ibin], pi_asquare[ibin], ida_tabQem[ibin], graphite, abs, 0.1, 10000);  
    Tmax = fmin(TSS*pow(10,fmax(fmin(-3-log10(agrains[ibin]),3),1)), 4e4);
    //Tmax = 100;
    printf("%.4e %.4e %.4e \n", agrains[ibin], TSS, Tmax);
    ierr = getTemperatureDist(ibin, graphite, Tmax, Ts, Ps);
    ierr = writeTempdist(ibin, Ts, Ps, graphite, NTbins);
    //start = clock();
    //int imax = 1000;
    //for(ibin = 0; ibin < imax; ibin++){ 
    //    ierr = getTemperatureDist(50, 0, 5000, Ts, Ps);
    //    ierr = getTemperatureDist(0, 1, 10000, Ts, Ps);
    //}
    //end = clock();
    //cpu_time_used = ((double) (end - start))/CLOCKS_PER_SEC;
    //printf("time per call: %f \n", cpu_time_used/(2*imax));
    //printf("time for 1000 cells with 40 bins: %f \n", cpu_time_used*1e3*40/(2*imax));
    
    return 0;
}

