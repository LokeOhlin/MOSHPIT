#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "hydro.h"
#include "radchem.h"
#include "cgeneral.h"
int nrealPars = 24;
real_list_t *dustDPars = NULL;
int nintPars = 11;
int_list_t *dustIPars = NULL;
//int nstrPars = 10;
//str_list_t chemSPars[];

// set default parameter lists 


int amin, amax;
// scratch arrays for calculations in each cell
double *dadt   = NULL;
double *number = NULL;
double *slope  = NULL;
double *mass   = NULL;

double *Mnew = NULL;
double *Nnew = NULL;
double *Snew = NULL;



int setDustPars(){
    // Allocate spaces
    dustDPars = (real_list_t *) malloc(nrealPars * sizeof(real_list_t));
    dustIPars = (int_list_t *) malloc(nintPars * sizeof(int_list_t));

    // Real/Double parameters
    strcpy(chemDPars[0].name, "dust_amin");  chemDPars[0].value = 5e-7;
    strcpy(chemDPars[1].name, "dust_amax");  chemDPars[1].value = 1e-3;
    strcpy(chemDPars[2].name, "dust_fsilicone");  chemDPars[2].value = 0.5;
    strcpy(chemDPars[3].name, "dust_Tsput");  chemDPars[3].value = 2e6;
    strcpy(chemDPars[3].name, "dust_dadtlim");  chemDPars[3].value = 1.0;


    strcpy(chemDPars[3].name, "dust_cmin");  chemDPars[3].value = 0.0;
    strcpy(chemDPars[3].name, "dust_cmax");  chemDPars[3].value = 0.0;
    strcpy(chemDPars[3].name, "dust_smin");  chemDPars[3].value = 0.0;
    strcpy(chemDPars[3].name, "dust_smax");  chemDPars[3].value = 0.0;
    strcpy(chemDPars[3].name, "dust_plaw");  chemDPars[3].value = 0.0;
    

    // Integer parameters
    strcpy(chemIPars[0].name, "dust_nbins"); chemIPars[0].value =  30;
    strcpy(chemIPars[0].name, "dust_intitDist"); chemIPars[0].value = 0;
    return 1;
}

int checkDustPars(char *name, char *value){
    // Method to see if parameter matches any of the defined chemistry parameters
    int ipar;
    // Check reals
    for(ipar = 0; ipar < nrealPars; ipar ++){
        if(compStr(name, chemDPars[ipar].name, 80) > 0){
            chemDPars[ipar].value = atof(value);
            return 1;
        }
    }

    for(ipar = 0; ipar < nintPars; ipar ++){
        if(compStr(name, chemIPars[ipar].name, 80) > 0){
            chemIPars[ipar].value = atoi(value);
            return 1;
        }
    }
    return -1;
}

// Method to get value of a double(real) parameter
void getrealdustpar(char *name, double *value){
    int ipar;
    for(ipar = 0; ipar < nrealPars; ipar ++){
        if(compStr(name, chemDPars[ipar].name, 80) > 0){
            *value = chemDPars[ipar].value;
            return;
        }
    }
    printf("%s NOT FOUND\n", name);

}

// Method to get value of an integer parameter
void getintegerdustpar(char *name, int *value){
    int ipar;
    for(ipar = 0; ipar < nintPars; ipar ++){
        if(compStr(name, chemIPars[ipar].name, 80) > 0){
            *value = chemIPars[ipar].value;
            return;
        }
    }
    printf("%s NOT FOUND\n", name);
}



int initDustPowerLaw(){
    int ierr, ibin;
    double cmax, cmin, smax, smin, plaw;
    double snorm, cnorm;
    double amin, amax;
    double smin, smax;

    ierr = getrealdustpar("dust_cmin",cmin);    
    ierr = getrealdustpar("dust_cmax",cmax);    
    ierr = getrealdustpar("dust_smin",smin);    
    ierr = getrealdustpar("dust_smax",smax);    
    ierr = getrealdustpar("dust_plaw",plaw);

    cnorm = (plaw+1)/(cmax**(plaw+1) - cmin**(plaw+1));
    snorm = (plaw+1)/(smax**(plaw+1) - smin**(plaw+1));
    
    for(ibin = 0; ibin < isilicate; ibin++){
        if(abin_e[ibin] < cmin){
            amin = cmin;        
        } else {
            amin = abin_e[ibin];
        }

        if(abin_e[ibin+1] > cmax){
            amax = cmax;        
        } else {
            amax = abin_e[ibin+1];
        }
        
        if(amin > amax){
            number[ibin] = 0;
            slope[ibin]  = 0;
            continue;
        }
        number[ibin] = cnorm*(amax**(plaw+1)-amin**(plaw+1))/(plaw+1);
        // slope to match max bin edge  (or maximum in distribution)
        smax  = (cnorm*amax**plaw - number[ibin]/(abin_e[ibin+1]-abin_e[ibin]));
        smax  = smax/(amax-abin_c[ibin]);
        // min bin edge (or minimum in distribution)
        smin  = (cnorm*amin**plaw - number[ibin]/(abin_e[ibin+1]-abin_e[ibin]));
        smin  = smin/(amin-abin_c[ibin]);
        
        // take the average
        slope[ibin] = (smin+smax)/2.;
    }
    
    for(ibin = isilicate; ibin < NdustBins; ibin++){
        iabin = iabin - isilicate;
        if(abin_e[iabin] < smin){
            amin = smin;        
        } else {
            amin = abin_e[iabin];
        }

        if(abin_e[iabin+1] > smax){
            amax = smax;        
        } else {
            amax = abin_e[iabin+1];
        }
        
        if(amin > amax){
            number[ibin] = 0;
            slope[ibin]  = 0;
            continue;
        }
        number[ibin] = snorm*(amax**(plaw+1)-amin**(plaw+1))/(plaw+1);
        // slope to match max bin edge  (or maximum in distribution)
        smax  = (snorm*amax**plaw - number[ibin]/(abin_e[iabin+1]-abin_e[iabin]));
        smax  = smax/(amax-abin_c[iabin]);
        // min bin edge (or minimum in distribution)
        smin  = (snorm*amin**plaw - number[ibin]/(abin_e[iabin+1]-abin_e[iabin]));
        smin  = smin/(amin-abin_c[iabin]);
        
        // take the average
        slope[ibin] = (smin+smax)/2.;
    }

    return 1;
}






int initDust(){
    int ierr, ibin;
    double amin, amax, da;
    // Number of dust size bins per species
    ierr = getintegerdustpar("dust_nbins", &NaBins);
    ierr = getrealdustpar("dust_amin", &amin);
    ierr = getrealdustpar("dust_amax", &amax);
    
    if (amin > amax ) {
        printf("ERROR:  amin > amax \n");
        return -1;
    }

    da = (log(amax) - log(amin))/ Nabins;
    
    // add ghost cells, two below, one above
    Nabins = Nabins + 3;

    // allocate arrays
    abin_c = (double *) malloc(Nabins*sizeof(double));
    abin_e = (double *) malloc((Nabins+1)*sizeof(double));
    // populate the edges of the cells
    abin_e[0] = exp(log(amin)-2*da);
    abin_e[1] = exp(log(amin)-da);
    abin_e[2] = amin;

    for(ibin = 3; ibin < Nabins +1, ibin++){
        abin_e[ibin] = exp(log(abin_e[ibin-1])+da);
    }
    // Calculate cell centers (linear)
    for(ibin = 0; ibin < Nabins; ibin ++){
        abin_c[ibin] = (abin_e[ibin+1]+abin_e[ibin+2])/2.
    }
    

    // Mass ratio between silicates and carbon species
    ierr = getrealdustpar("dust_fsilicone", &fSi);

    // if we have silicates and carbonates (eg two species) we need to increase the bins
    if (fSi > 0) {
        NdustBins = 2*Nabins;
        isilicone = Nabins;
        printf("Ndustbins = %d, nsizebins = %d, isilicone = %d", NdustBins, Nabins, isilicone); 
    } else if(fSi == 1.0){ //only silicates
        NdustBins = Nabins;
        isilicone = 0;
    } else {
        NdustBins = Nabins;
        isilicone = Nabins;
    }


    // allocate all arrays
    dadt   = (double *) malloc(NdustBins*sizeof(double));  
    number = (double *) malloc(NdustBins*sizeof(double)); 
    slope  = (double *) malloc(NdustBins*sizeof(double)); 
    mass   = (double *) malloc(NdustBins*sizeof(double));    
    Mnew   = (double *) malloc(NdustBins*sizeof(double)); 
    Nnew   = (double *) malloc(NdustBins*sizeof(double)); 
    Snew   = (double *) malloc(NdustBins*sizeof(double)); 


    // Set initial distribution
    ierr = getintegerdustpar("dust_initDist", &initDist);
    if(initDist == 0){
        ierr = initDustPowerLaw();
    }
    return 1;
}


int setCellInit(int icell){
    int ibin;
    for(ibin = 0; ibin < NdustBins; ibin++){
        ustate[IDUST_BEGIN + ibin*ndustVar]   = number[ibin];
        ustate[IDUST_BEGIN + ibin*ndustVar+1] = slope[ibin];
    }
    return 1;
}

