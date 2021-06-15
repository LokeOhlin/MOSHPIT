#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <tgmath.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dustRadiation.h>
#include <cgeneral.h>
#include <dust.h>
#include <constantsAndUnits.h>
#ifdef useDust
//
int NUbins;

// Energy bins and their widths
double *Ubins  = NULL;
double *dUbins = NULL;
double *Ul     = NULL;
// Transisiton matrix A 
double *Amatr = NULL;
double *Adivs = NULL;

// Array of temperature distribution
long double *Ps;
double *Ts;

int matr_idxMax;

double get_Em(double agrain, double piaa, int ida, double temp, int graphite){ 
    return 4 * piaa * getQemAve(agrain, temp, graphite, ida, -1) * sigma_stefanBoltz * pow(temp, 4);
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
    Ul     = (double *) malloc((Nbins + 1) * sizeof(double));
    Amatr  = (double *) malloc(Nbins * Nbins * sizeof(double));
    
    Ps = (long double *) malloc(Nbins *sizeof(long double));
    Ts = (double *) malloc(Nbins *sizeof(double));
    //Bmatr  = (double *) malloc(Nbins * Nbins * sizeof(double));
    Adivs  = (double *) malloc(Nbins * sizeof(double));
    matr_idxMax = Nbins*Nbins;
    return 1;
}



int getTemperatureDist(int ibin, int iabin, int graphite, double Tmax){
    //  Main call for calculating the temperature distribution of a given dust grain size
    //  In :
    //          int ibin       -  index of the dust grain in the dust size array
    //          int graphite   -  flag for type of dust (1 = graphite grains, 0 = silicate grains)
    //          double Tmax    -  Maximum temperature to distribution. Must be set high enough  
    //                            to not truncate distribution
//  out: 
    //          double *Ts      -  Array of temperature values to be written to. 
    //          long double *Ps -  Array of fraction of dust grains in each bin to be written to
    //
    //  NOTE: both out arrays need to allocated beforehand
    
    // Start by setting energy bins. We need maximum and minimum energy (Emax and Emin), where
    // Emax is user set and should be large enough to not truncate distrubtion too excessively
    // and Emin is set such that the cooling and heating from frequent low energy photons,
    // eg. the processes we consider to be continous, is in equilibrium
    int iUbin, fUbin, fUbinp, jUbin, idx, idxp, imatr, imatrp;
    double AiiSum, freq_if, E_inf, Qabs, radEdens;
    double Uf, dUf, UfR, UfL;
    double agrain = abin_c[iabin];
    double volgrain = volgrains[iabin];
    double Natom = Natoms[ibin];
    double piaa = pi_asquare[iabin];
    int ida_Qabs = ida_tabQabs[ibin];
    int ida_Qem = ida_tabQem[ibin];
    
    
    
    double Emax = getUfromT(Natom, volgrain, Tmax, graphite);
    
    //double nuCont = 2997924580.000;
    double nuCont = 2415068821914.3345;
    double ECont = nuCont*planck;
    ECont = 0;
    // Absorbed energy rate from frequent photons (eg low energy photons with nu < nuCont)
    double continous_absorption = getAbsorption_range(ibin, graphite, 0.0, nuCont);     
    // Get "stable state" temperature for continous processes. initial range between 0.01 and 100... TODO: CHECKIF ALWAYS VALID! Would need a lot of low energy photons (lambda>10 cm) for this to not be the case.
    double Tmin = StableState_temperature(agrain, piaa, ida_Qem,  graphite, continous_absorption, 0.01, 100);
    double Emin = getUfromT(Natom, volgrain, Tmin, graphite);
    // populate Energy bins in log scale 
    // TODO: we dont need high samle rates at low temperatures.
    // Make more eduacated distribution
    double delE = exp((log(Emax) - log(Emin))/NUbins);
    
    //printf("%d %d %.4e %.4e %.4e \n", ibin, iabin, continous_absorption, Tmin, Emin); 
    Ubins[0] = Emin;
    Ts[0] = Tmin;
    Ps[0] = 0.0; 
    for(iUbin = 1; iUbin < NUbins; iUbin++){
        Ubins[iUbin] = Emin*pow(delE,iUbin);
        //Ubins[iUbin] = Emin + delE*iUbin;
        Ts[iUbin] = getTfromU(Natom, volgrain, Ubins[iUbin], graphite);
        // left edge of bin
        Ul[iUbin] = (Ubins[iUbin] + Ubins[iUbin-1])*0.5;
        dUbins[iUbin]    = -Ul[iUbin];
        dUbins[iUbin-1] +=  Ul[iUbin];
        // reset bin
        Ps[iUbin] = 0.0;
    }
    // fix dUbins at upper and lower limit (remember dUbins[-1] = - (Ubins[-1] + Ubins[-2])/2 still)
    dUbins[0] = 2*(dUbins[0] - Ubins[0]);
    Ul[0] = Ubins[0] - 0.5 * dUbins[0];

    dUbins[NUbins - 1] = 2*(Ubins[NUbins - 1] + dUbins[NUbins -1]);
    Ul[NUbins] = Ubins[NUbins - 1] + 0.5*dUbins[NUbins - 1];
    // set transition matrix to zero since its reused for all grains
    for(iUbin = 0; iUbin < matr_idxMax; iUbin++){
        Amatr[iUbin] = 0;
    }
    
    // start populating transition matrix
    // transition is going from iUbin (initial) to fUbin (final)
    for(jUbin = 0; jUbin < NUbins-1; jUbin++){
        // continous cooling (since dUdt for Ubins[0] is defined as 0)
        fUbin = jUbin;
        iUbin = jUbin+1; 
        idx = fUbin*NUbins + iUbin;
        Amatr[idx] = (get_Em(agrain, piaa, ida_Qem, Ts[iUbin], graphite) - continous_absorption)/(Ubins[iUbin] - Ubins[fUbin]); 
        
        // Discreete heating processes from lower energy bins
        fUbinp = fUbin + 1;
        imatr = fUbin*NUbins;
        
        // parameters for the upper bin (constant in next loop)
        Uf  = Ubins[fUbin];
        UfL = Ul[fUbin];
        UfR = Ul[fUbinp];
        dUf = dUbins[fUbin];

        // Loop over all lower bins and determine transition rate uppwards
        for(iUbin = 0; iUbin < fUbin; iUbin++){
            idx = imatr + iUbin;
            // Upward transition from iUbin to fUbin 
            Amatr[idx] =  getUpwardTransition(ibin, graphite, Ubins[iUbin], Ul[iUbin], Ul[iUbin+1], dUbins[iUbin], Uf, UfL, UfR, dUf, ECont);
        }
        // Account for intrabin transitions
        if(fUbin > 0){
            iUbin = fUbin - 1;
            idx    = imatr + iUbin;   
            Amatr[idx] = Amatr[idx] + intrabinUpwardTransition(ibin, graphite, dUbins[iUbin], Ubins[iUbin], Uf);
        }

        // Add discreete heating from jbin to final bin.
        iUbin = jUbin;
        fUbin = NUbins - 1;
        idx = fUbin*NUbins + iUbin;
        
        // frequency of photons needed to get from iUbin to fUbin
        Amatr[idx] =  getUpwardTransition(ibin, graphite, Ubins[iUbin], Ul[iUbin], Ul[iUbin+1], dUbins[iUbin],
                                                          Ubins[fUbin], Ul[fUbin], Ul[fUbin+1], dUbins[fUbin], ECont);
        if(iUbin == fUbin-1){
            Amatr[idx] = Amatr[idx] + intrabinUpwardTransition(ibin, graphite, dUbins[iUbin], Ubins[iUbin], Ubins[fUbin]);
        }
        
        // Last bin technically goes to infintity. Take this into account
        E_inf = (Ubins[fUbin] + dUbins[fUbin]*0.5 - Ubins[iUbin]);
        Amatr[idx] += getAbsorptionNum_range(ibin, graphite, E_inf, 1e99); 
    }

    // now do the diagonal elements ( eg. everything going out from bins)
    for(fUbin = 0; fUbin < NUbins; fUbin++){
        AiiSum  = 0;
        imatr = fUbin * NUbins;
        // can only transition into fUbin from iUbin < fUbin
        // or iUbin = fUbin + 1
        for(iUbin = 0; iUbin < fUbin; iUbin++){
            idx = imatr + iUbin;
            AiiSum -= Amatr[idx];
        }
        // check that we are not in the uppermost bin 
        if(fUbin < NUbins - 1){
            idx = imatr + fUbin + 1;
            AiiSum -= Amatr[idx];
        }
        idx = fUbin*NUbins + fUbin;
        Amatr[idx] = AiiSum;
        // Save element used for normalisation later
        if(fUbin >= 1){
            idx = (fUbin-1)*NUbins + fUbin;
            Adivs[fUbin] = 1/Amatr[idx];
        }
    }

    
    // We have what we need from the A matrix,
    // we now store B matrix (guhathakurta & draine 1989) 
    // We reuse A matrix to save some time & memory
    
    for(jUbin = 0; jUbin < NUbins-1; jUbin++){
        // imatr selects starting point in 1D array. 
        // multiply with NUbins to save some time.
        fUbin  = (NUbins - 2 - jUbin);
        
        imatr  = fUbin * NUbins;
        imatrp = (fUbin + 1)*NUbins; 
        for(iUbin = 0; iUbin < fUbin; iUbin++){
            idx  = imatr  + iUbin;
            idxp = imatrp + iUbin;
            Amatr[idx] = Amatr[idxp] + Amatr[idx];
        }

    }

    // Calculate probablilites
    // If heated enough, then probability of finding a grain 
    // at the lowest temperature is vanishing. Use long doubles 
    // and start with initial guess at a low number to allow 
    // for several orders of magnitude growth
    Ps[0] = 1e-200;
    long double norm = Ps[0];
    long double Pbin;
    for(fUbin = 1; fUbin < NUbins; fUbin++){
        Pbin = 0;
        idx = fUbin * NUbins;
        for(iUbin = 0; iUbin < fUbin; iUbin++){
            Pbin += Amatr[idx + iUbin]*Ps[iUbin];
        }
        idx = (fUbin-1)*NUbins + fUbin;
        
        Ps[fUbin] = Pbin * Adivs[fUbin];
        norm += Ps[fUbin];
    }
    if(isfinite(norm) == 0){
        printf("ERROR: DUST TEMPERATURE DISTRIBUTION NOT FINITE \n");
        for(iUbin = 0; iUbin < NUbins; iUbin++){
            printf("%.4e, %.4e, %.4Le \n", Ts[iUbin], Ubins[iUbin], Ps[iUbin]);
        }
        exit(0);
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



// constants for sublimation rate
double Asub_c = 4.6e29;
double Asub_s = 4.9e30;
double debyeEnergy_c = 4.3490475e-14;
double debyeEnergy_s = 4.86679125e-14;


double Asub;
double Bsub;
double fsub;
double debyeEnergy;
double bsub;

double get_dadt_sublimation_iUbin(int iUbin){
    // Start by getting the suppresion factor for low energies
    double meanQuanta_pdof = Ubins[iUbin]/(fsub*debyeEnergy);
    double arg, Sn;

    if(meanQuanta_pdof*fsub+1-bsub > 0 && meanQuanta_pdof*fsub+fsub-bsub >0){
        arg = meanQuanta_pdof * fsub + 1;
        double logGam1 = lgamma(arg);
        
        arg = meanQuanta_pdof * fsub + fsub - bsub;     double logGam2 = lgamma(arg);
        
        arg = meanQuanta_pdof * fsub + 1 - bsub;
        double logGam3 = lgamma(arg);
        
        arg = meanQuanta_pdof * fsub + fsub;
        double logGam4 = lgamma(arg);

        double exparg = (log((1+meanQuanta_pdof)/meanQuanta_pdof)*bsub + logGam1 + logGam2) - (logGam3 + logGam4);

        Sn = exp(exparg);
    } else {
        Sn = 0;
    }
    // default rate
    double Rn = Asub*exp(-Bsub/Ts[iUbin]);
    return -Sn*Rn*(double)Ps[iUbin];    
}

double get_dadt_sublimation(int ibin, int iabin, int graphite){
    int iUbin, ierr;
    double abs = getAbsorption(ibin, graphite);
    // calculate stable state temperature (should never be higher than 1e4...)
    double TSS = StableState_temperature(abin_c[iabin], pi_asquare[iabin], ida_tabQem[ibin], graphite, abs, 0.1, 1e4);
    // Take the maximum temperature in the distribution to be the minimum of 4e4 K or
    // Tmax = Tss * max((10^-3 cm/max(a, 1e-6)),10)
    double Tmax = fmin(TSS*pow(10,fmax(fmin(-3-log10(abin_c[iabin]),3),1)), 4e4);
    // calculate temperature distribution
    ierr = getTemperatureDist(ibin, iabin, graphite, Tmax);
    if(ierr == -2){
        printf("WARNING: final temperature bin has more than 1e-13. Temperature dist may not wide enough\n");
    } 

    // Set parameters that are constant for given grain.
    if(graphite){
        Asub = Asub_c*aveMatom_c/rho_c;
        Bsub = 81200 - 20000*pow(Natoms[ibin]-1, -1./3.);
        debyeEnergy = debyeEnergy_c;
    } else {
        Asub = Asub_s*aveMatom_s/rho_s;
        Bsub = 68100 - 20000*pow(Natoms[ibin]-1, -1./3.);
        debyeEnergy = debyeEnergy_s;
    }
    bsub = boltzmann*Bsub/debyeEnergy;
    fsub = 3*Natoms[ibin] - 6;
    double dadt_sub = 0;
    for(iUbin = 0; iUbin < NUbins; iUbin++){
        dadt_sub += get_dadt_sublimation_iUbin(iUbin);
        //if(iabin < 4 && (iUbin == 175 || iUbin == 199 )){
        //    printf("%d %.4e| %.4e %.4e\n", ibin, abin_c[iabin], Ts[iUbin], get_dadt_sublimation_iUbin(iUbin)/(double)Ps[iUbin]);
        //}
    }
    //if(graphite == 0 && iabin > 3){
    //    exit(0);
    //}
    return dadt_sub;
}
#endif
