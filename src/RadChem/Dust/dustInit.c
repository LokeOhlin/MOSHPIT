#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "hydro.h"
#include "radchem.h"
#include "cgeneral.h"
#include "dust.h"
int nrealDustPars = 12;
real_list_t *dustDPars = NULL;
int nintDustPars = 2;
int_list_t *dustIPars = NULL;
//int nstrPars = 10;
//str_list_t chemSPars[];



double amin, amax;
double *abin_e = NULL;
double *abin_c = NULL;
double fSi;

int isilicone;
int Nabins;
// options for initial distribution
int initDist; 
// options for grain growth/evaporation rate
int dadt_mode; // 0 standard based on physical processes
               // 1 constant rate
               // 2 rate proportional to grain size 
               // 3 rate oscillating with time
double dadt_c;
double dadt_tscale;

// scratch arrays for calculations in each cell
double *dadt   = NULL;
double *number = NULL;
double *slope  = NULL;

double *Mnew = NULL;
double *Nnew = NULL;
double *Snew = NULL;

// SET TO PROPER VALUES
double rho_s = 3.0;
double rho_c = 3.0;

// timestep limiter
double dadt_lim;

// set default parameter lists 
int setDustPars(){
    // Allocate spaces
    dustDPars = (real_list_t *) malloc(nrealDustPars * sizeof(real_list_t));
    dustIPars = (int_list_t *) malloc(nintDustPars * sizeof(int_list_t));

    // Real/Double parameters
    strcpy(dustDPars[0].name, "dust_amin");      dustDPars[0].value = 5e-7;
    strcpy(dustDPars[1].name, "dust_amax");      dustDPars[1].value = 1e-3;
    strcpy(dustDPars[2].name, "dust_fsilicone"); dustDPars[2].value = 0.5;
    strcpy(dustDPars[3].name, "dust_Tsput");     dustDPars[3].value = 2e6;
    strcpy(dustDPars[4].name, "dust_dadt_lim");   dustDPars[4].value = 0.5;


    strcpy(dustDPars[5].name, "dust_cmin");  dustDPars[5].value = 1e-8;
    strcpy(dustDPars[6].name, "dust_cmax");  dustDPars[6].value = 1e-4;
    strcpy(dustDPars[7].name, "dust_smin");  dustDPars[7].value = 1e-8;
    strcpy(dustDPars[8].name, "dust_smax");  dustDPars[8].value = 1e-4;
    strcpy(dustDPars[9].name, "dust_plaw");  dustDPars[9].value = -3.5;

    strcpy(dustDPars[10].name, "dust_dadt_c");  dustDPars[10].value = 1e-13;
    strcpy(dustDPars[11].name, "dust_dadt_tscale");  dustDPars[11].value = 0.0;
    
    // Integer parameters
    strcpy(dustIPars[0].name, "dust_initDist"); dustIPars[0].value = 0;
    strcpy(dustIPars[1].name, "dust_dadt_mode"); dustIPars[1].value = 0;
    return 1;
}

int checkDustPars(char *name, char *value){
    // Method to see if parameter matches any of the defined chemistry parameters
    int ipar;
    // Check reals
    for(ipar = 0; ipar < nrealDustPars; ipar ++){
        if(compStr(name, dustDPars[ipar].name, 80) > 0){
            dustDPars[ipar].value = atof(value);
            return 1;
        }
    }

    for(ipar = 0; ipar < nintDustPars; ipar ++){
        if(compStr(name, dustIPars[ipar].name, 80) > 0){
            dustIPars[ipar].value = atoi(value);
            return 1;
        }
    }
    return -1;
}

// Method to get value of a double(real) parameter
void getrealdustpar(char *name, double *value){
    int ipar;
    for(ipar = 0; ipar < nrealDustPars; ipar ++){
        if(compStr(name, dustDPars[ipar].name, 80) > 0){
            *value = dustDPars[ipar].value;
            return;
        }
    }
    printf("%s NOT FOUND\n", name);

}

// Method to get value of an integer parameter
void getintegerdustpar(char *name, int *value){
    int ipar;
    for(ipar = 0; ipar < nintDustPars; ipar ++){
        if(compStr(name, dustIPars[ipar].name, 80) > 0){
            *value = dustIPars[ipar].value;
            return;
        }
    }
    printf("%s NOT FOUND\n", name);
}



int initDustPowerLaw(){
    int ierr, ibin, iabin;
    double cmax, cmin, smax, smin, plaw;
    double snorm, cnorm;
    double amin_bin, amax_bin;
    double slope_min, slope_max;
    double dadn_e, dadn_c, dadn_ep;
    double Nj, Sj, mass;
    getrealdustpar("dust_cmin",&cmin);    
    getrealdustpar("dust_cmax",&cmax);    
    getrealdustpar("dust_smin",&smin);    
    getrealdustpar("dust_smax",&smax);    
    getrealdustpar("dust_plaw",&plaw);

    cnorm = (plaw+1)/(pow(cmax,plaw+1) - pow(cmin,plaw+1));
    snorm = (plaw+1)/(pow(smax,plaw+1) - pow(smin,plaw+1));
    printf("cmax %.4e, cmin %.4e, smax %.4e, smin %.4e, plaw %.4e \n", cmax, cmin, smax, smin, plaw); 
    for(ibin = 0; ibin < isilicone; ibin++){
        if(abin_e[ibin] < cmin){
            amin_bin = cmin;        
        } else {
            amin_bin = abin_e[ibin];
        }

        if(abin_e[ibin+1] > cmax){
            amax_bin = cmax;        
        } else {
            amax_bin = abin_e[ibin+1];
        }
        
        if(amin_bin >= amax_bin){
            number[ibin] = 0;
            slope[ibin]  = 0;
            continue;
        }
        // total number in bin is  integral_bin dnda da
        number[ibin] = cnorm*(pow(amax_bin,plaw+1)-pow(amin_bin,plaw+1))/(plaw+1);
        // slope to match max bin edge  (or maximum in distribution)
        slope_max  = (cnorm*pow(amax_bin,plaw) - number[ibin]/(abin_e[ibin+1]-abin_e[ibin]));
        slope_max  = slope_max/(amax_bin-abin_c[ibin]);
        // min bin edge (or minimum in distribution)
        if(abin_e[ibin] < cmin){
            slope_min  = number[ibin]/(abin_e[ibin+1]-abin_e[ibin])/(abin_c[ibin]-abin_e[ibin]);
        }else{
            slope_min  = (cnorm*pow(amin_bin,plaw) - number[ibin]/(abin_e[ibin+1]-abin_e[ibin]));
            slope_min  = slope_min/(amin_bin-abin_c[ibin]);
        } 
        // take the average
        // printf("bin %d : slope_min %.4e slope_max %.4e %.4e %.4e %.4e %.4e %.4e\n", ibin, slope_min, slope_max, amin_bin, amax_bin, abin_c[ibin], cnorm*pow(amax_bin,plaw), number[ibin]/(abin_e[ibin+1]-abin_e[ibin]));
        slope[ibin] = (slope_min+slope_max)/2.;
        // get mass (needed for slope limiter)
        mass = getMass(number[ibin], slope[ibin], ibin);
        // make sure no negative dnda
        printf(" %d %.4e %.4e %.4e %.4e\n",ibin, number[ibin], slope[ibin],slope_min, slope_max);
        limitSlope(&Nj, &Sj, number[ibin], slope[ibin], mass, ibin);
        number[ibin] = Nj;
        slope[ibin]  = Sj;
        
        dadn_e  = number[ibin]/(abin_e[ibin+1]-abin_e[ibin]) + slope[ibin]*(abin_e[ibin]-abin_c[ibin]) ;
        dadn_c  = number[ibin]/(abin_e[ibin+1]-abin_e[ibin]) ;
        dadn_ep = number[ibin]/(abin_e[ibin+1]-abin_e[ibin]) + slope[ibin]*(abin_e[ibin+1]-abin_c[ibin]) ;
        printf(" %d %.4e %.4e %.4e %.4e %.4e %.4e %.4e\n",ibin, number[ibin], slope[ibin],dadn_e, dadn_c, dadn_ep, slope_min, slope_max);

        if(dadn_ep < 0){
            slope[ibin] = slope_max;
        } 
    }
    
    for(ibin = isilicone; ibin < NdustBins; ibin++){
        iabin = ibin - isilicone;
        if(abin_e[iabin] < smin){
            amin_bin = smin;        
        } else {
            amin_bin = abin_e[iabin];
        }

        if(abin_e[iabin+1] > smax){
            amax_bin = smax;        
        } else {
            amax_bin = abin_e[iabin+1];
        }
        
        if(amin_bin >= amax_bin){
            number[ibin] = 0;
            slope[ibin]  = 0;
            continue;
        }
        number[ibin] = snorm*(pow(amax_bin,plaw+1)-pow(amin_bin,plaw+1))/(plaw+1);
        // slope to match max bin edge  (or maximum in distribution)
        slope_max  = (snorm*pow(amax_bin,plaw) - number[ibin]/(abin_e[iabin+1]-abin_e[iabin]));
        slope_max  = slope_max/(amax_bin-abin_c[iabin]);
        // min bin edge (or minimum in distribution)
        slope_min  = (snorm*pow(amin_bin,plaw) - number[ibin]/(abin_e[iabin+1]-abin_e[iabin]));
        slope_min  = slope_min/(amin_bin-abin_c[iabin]);
        
        // take the average
        // printf("bin %d : slope_min %.4e slope_max %.4e\n", ibin, slope_min, slope_max);
        slope[ibin] = (slope_min+slope_max)/2.;
    }
    
    return 1;
}






int initDust(){
    int ierr, ibin;
    double amin, amax, da;

    // Number of dust size bins per species
    getrealdustpar("dust_amin", &amin);
    getrealdustpar("dust_amax", &amax);
    
    // Mass ratio between silicates and carbon species
    getrealdustpar("dust_fsilicone", &fSi);

    // if we have silicates and carbonates (eg two species) we need to reduce the number of bins to accomodate both
    if (fSi > 0 && fSi < 1) {
        if(NdustBins % 2 != 0){
            printf("Need even number of dust bins in order to have both silicates and graphite grains\n");
            return -1;
        }
        Nabins = NdustBins/2;
        isilicone = Nabins;
        printf("Ndustbins = %d, nsizebins = %d, isilicone = %d\n", NdustBins, Nabins, isilicone); 
    } else if(fSi == 1.0){ //only silicates
        Nabins = NdustBins;
        isilicone = 0;
    } else {
        Nabins = NdustBins;
        isilicone = NdustBins;
    }
    
    if (amin > amax ) {
        printf("ERROR:  amin > amax \n");
        return -1;
    }

    da = (log(amax) - log(amin))/ Nabins;
    
    // add ghost cells, one below, one above
    //    Nabins = Nabins + 2;

    // allocate arrays
    abin_c = (double *) malloc(Nabins*sizeof(double));
    abin_e = (double *) malloc((Nabins+1)*sizeof(double));
    // populate the edges of the cells
    abin_e[0] = amin;

    for(ibin = 1; ibin < Nabins+1; ibin++){
        abin_e[ibin] = exp(log(abin_e[ibin-1])+da);
    }
    // Calculate cell centers (linear)
    for(ibin = 0; ibin < Nabins; ibin ++){
        abin_c[ibin] = (abin_e[ibin]+abin_e[ibin+1])/2.;
    }
    



    // allocate all arrays
    dadt   = (double *) malloc(NdustBins*sizeof(double));  
    number = (double *) malloc(NdustBins*sizeof(double)); 
    slope  = (double *) malloc(NdustBins*sizeof(double)); 
    Mnew   = (double *) malloc(NdustBins*sizeof(double)); 
    Nnew   = (double *) malloc(NdustBins*sizeof(double)); 
    Snew   = (double *) malloc(NdustBins*sizeof(double)); 


    // Set initial distribution
    getintegerdustpar("dust_initDist", &initDist);
    printf("INITDIST %d", initDist);
    if(initDist == 0){
        ierr = initDustPowerLaw();
    }

    // get parameters for dust growth
    getintegerdustpar("dust_dadt_mode", &dadt_mode);
    getrealdustpar("dust_dadt_c", &dadt_c);
    getrealdustpar("dust_dadt_tscale", &dadt_tscale);

    // timestep limiters
    getrealdustpar("dust_dadt_lim", &dadt_lim);
    return 1;
}


// Set initial cell dust data, normalised such that everything in iabin >= 1 has a total mass matching the dustMass
int setCellInit(int icell, double dustMass){
    int ibin, iabin, graphite;
    double mass;
    double Mtot_s, Mtot_c;
    double norm_s, norm_c, norm, Sj;

    Mtot_c = 0;
    Mtot_s = 0;

    for(ibin = 0; ibin < NdustBins; ibin++){
        if(ibin < isilicone){
            graphite = 1;
        } else {
            graphite = 0;
        }

        // graphite and silicone have the same size bins. Use half array 
        iabin = ibin - isilicone*(1-graphite);    

        // Conversion to mass density
        norm = 4*M_PI/3.;
        if(graphite){
            norm = norm * rho_c;
        } else {
            norm = norm * rho_s;
        }
        // Mass (density) in bin
        mass = getMass(number[ibin], slope[ibin], iabin)*norm;
        if(graphite){
            Mtot_c += mass;
        } else {
            Mtot_s += mass;
        }
        ustate[icell*nvar + IDUST_START + ibin*NdustVar]   = mass;
    }
    
    norm_s = fSi*dustMass/Mtot_s;
    norm_c = (1-fSi)*dustMass/Mtot_c;

    for(ibin = 0; ibin < NdustBins; ibin++){
        if(ibin < isilicone){
            graphite = 1;
        } else {
            graphite = 0;
        }

        // Mass (density) in bin
        mass = ustate[icell*nvar + IDUST_START + ibin*NdustVar];
        Sj   = slope[ibin];
        // Normalize to match mass density
        if(graphite){
            mass = mass * norm_c;
            Sj   = Sj * norm_c; 
        } else {
            mass = mass * norm_s;
            Sj   = Sj * norm_s;
        }
        ustate[icell*nvar + IDUST_START + ibin*NdustVar] = mass;
        ustate[icell*nvar + IDUST_START + ibin*NdustVar+1] = Sj;
    }

    return 1;
}

