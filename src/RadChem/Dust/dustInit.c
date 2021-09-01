#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <hydro.h>
#include <radchem.h>
#include <cgeneral.h>
#include <dust.h>
#include <dustRadiation.h>
#include <constantsAndUnits.h>
#ifdef useDust
int nrealDustPars = 13;
real_list_t *dustDPars = NULL;
int nintDustPars = 9;
int_list_t *dustIPars = NULL;
//int nstrPars = 10;
//str_list_t chemSPars[];

double amin, amax;
double *abin_e = NULL;
double *abin_c = NULL;

// Factors om front of N and S after integrating over bin to get average size
double *NfactA = NULL;
double *SfactA = NULL;

// Factors in front of N and S after integrating over bin to get mass
double *NfactM = NULL;
double *SfactM = NULL;

// Geometric quantities 
double *pi_asquare = NULL;
double *volgrains  = NULL;

double fSi;

int isilicone;
int Nabins;
int dust_nbins;
// options for initial distribution
int initDist; 
// options for grain growth/evaporation rate
int dadt_mode; // 0 standard based on physical processes
               // 1 constant rate
               // 2 rate proportional to grain size 
               // 3 rate oscillating with time
double dadt_c;
double dadt_tscale;
double dadt_min;

// whether the lower/upper bounds of the distribution are to be pile up (eg. gathered) or outflow when rebinned
int dust_lowerBound_pileUp;
int dust_upperBound_pileUp;

// Number of atoms in a grain 
double *Natoms = NULL;
// scratch arrays for calculations in each cell
double *dadt         = NULL;
double *dadt_fixed   = NULL;
double *number = NULL;
double *slope  = NULL;

double *dust_vrel = NULL;

double *Mnew = NULL;
double *Nnew = NULL;
double *Snew = NULL;

// density and average atom mass of silicates 
double rho_s = 3.5;
double aveMatom_s = 3.3665848928571423e-23; //(~20 mH)
// density and average atom mass of graphites
double rho_c = 2.26;
double aveMatom_c = 2.0118425e-23; //(~12 mH)

// timestep limiter
double dadt_lim;

// Maximum number of substeps taken in integration of dust sizes
int dust_maxSubSteps;

// saved indexes in tables
int *ida_tabQabs = NULL;
int *ida_tabQem  = NULL;
int *ida_tabSput = NULL;

// options for physical processes
int dust_useRadiation, dust_useSublimation, dust_useSputtering;

//flag for writing dust output
int outputDust;
int outputNum = 0;

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

    strcpy(dustDPars[10].name, "dust_dadt_c");  dustDPars[10].value = 0.0;
    strcpy(dustDPars[11].name, "dust_dadt_tscale");  dustDPars[11].value = 0.0;
    
    strcpy(dustDPars[12].name, "dust_dadt_min");  dustDPars[12].value = 3.17e-17;
    
    
    // Integer parameters
    strcpy(dustIPars[0].name, "dust_initDist"); dustIPars[0].value = 0;
    strcpy(dustIPars[1].name, "dust_dadt_mode"); dustIPars[1].value = 0;
    strcpy(dustIPars[2].name, "dust_useRadiation"); dustIPars[2].value = 1;
    strcpy(dustIPars[3].name, "dust_useSublimation"); dustIPars[3].value = 1;
    strcpy(dustIPars[4].name, "dust_useSputtering"); dustIPars[4].value = 1;
    strcpy(dustIPars[5].name, "dust_NtempBins"); dustIPars[5].value = 200;
    strcpy(dustIPars[6].name, "dust_maxSubSteps"); dustIPars[6].value = 1000;
    strcpy(dustIPars[7].name, "dust_lowerBound_pileUp"); dustIPars[7].value = 0;
    strcpy(dustIPars[8].name, "dust_upperBound_pileUp"); dustIPars[8].value = 1;
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
    int ibin, iabin;
    double cmax, cmin, smax, smin, plaw;
    double snorm, cnorm;
    double amin_bin, amax_bin;
    double slope_min, slope_max;
    double Nj, Sj, mass;
    getrealdustpar("dust_cmin",&cmin);    
    getrealdustpar("dust_cmax",&cmax);    
    getrealdustpar("dust_smin",&smin);    
    getrealdustpar("dust_smax",&smax);    
    getrealdustpar("dust_plaw",&plaw);

    cnorm = (plaw+1)/(pow(cmax,plaw+1) - pow(cmin,plaw+1));
    snorm = (plaw+1)/(pow(smax,plaw+1) - pow(smin,plaw+1));
    printf("cmax %.4e, cmin %.4e, smax %.4e, smin %.4e, plaw %.4e \n", cmax, cmin, smax, smin, plaw);
    if(fSi < 1.0){ 
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

            // slope to match max bin edge, if outside distribution the it is zero
            if(abin_e[ibin+1] > cmax){
                slope_max  = -number[ibin]/(abin_e[ibin+1]-abin_e[ibin])/(abin_e[ibin+1]-abin_c[ibin]);
            }else{
                slope_max  = (cnorm*pow(amax_bin,plaw) - number[ibin]/(abin_e[ibin+1]-abin_e[ibin]));
                slope_max  = slope_max/(amax_bin-abin_c[ibin]);
            }
            // min bin edge 
            if(abin_e[ibin] < cmin){
                slope_min  = number[ibin]/(abin_e[ibin+1]-abin_e[ibin])/(abin_c[ibin]-abin_e[ibin]);
            }else{
                slope_min  = (cnorm*pow(amin_bin,plaw) - number[ibin]/(abin_e[ibin+1]-abin_e[ibin]));
                slope_min  = slope_min/(amin_bin-abin_c[ibin]);
            } 
            // take the average
            slope[ibin] = (slope_min+slope_max)/2.;
            // get mass (needed for slope limiter)
            mass = getMass(number[ibin], slope[ibin], ibin);
            // make sure no negative dnda
            
            limitSlope(&Nj, &Sj, number[ibin], slope[ibin], mass, ibin);
            number[ibin] = Nj;
            slope[ibin]  = Sj;
        }
    }
    if(fSi > 0.0){ 
        for(ibin = isilicone; ibin < isilicone + Nabins; ibin++){
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
            
            // slope to match max bin edge, if outside distribution the it is zero
            if(abin_e[iabin+1] > smax){
                slope_max  = -number[ibin]/(abin_e[iabin+1]-abin_e[iabin])/(abin_e[iabin+1]-abin_c[iabin]);
            }else{
                slope_max  = (snorm*pow(amax_bin,plaw) - number[ibin]/(abin_e[iabin+1]-abin_e[iabin]));
                slope_max  = slope_max/(amax_bin-abin_c[iabin]);
            }
            // min bin edge 
            if(abin_e[iabin] < smin){
                slope_min  = number[ibin]/(abin_e[iabin+1]-abin_e[iabin])/(abin_c[iabin]-abin_e[iabin]);
            }else{
                slope_min  = (snorm*pow(amin_bin,plaw) - number[ibin]/(abin_e[iabin+1]-abin_e[iabin]));
                slope_min  = slope_min/(amin_bin-abin_c[iabin]);
            } 
            // take the average
            slope[ibin] = (slope_min+slope_max)/2.;
            // get mass (needed for slope limiter)
            mass = getMass(number[ibin], slope[ibin], iabin);
            // make sure no negative dnda
            
            limitSlope(&Nj, &Sj, number[ibin], slope[ibin], mass, iabin);
            number[ibin] = Nj;
            slope[ibin]  = Sj;
        }
    } 
    return 1;
}






int initDust(){
    int ierr, ibin, iabin, idx, graphite, NtempBins;
    double amin, amax, da, norm;
    double ae, aep, ac;
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
    } else {
        Nabins = NdustBins;
    }
    
    if (amin > amax ) {
        printf("ERROR:  amin > amax \n");
        return -1;
    }

    da = exp((log(amax) - log(amin))/ Nabins);
    
    // add ghost cells, one below, one above
    Nabins = Nabins + 2*dust_nghost;

    // allocate arrays
    abin_c     = (double *) malloc(Nabins*sizeof(double));
    pi_asquare = (double *) malloc(Nabins*sizeof(double));
    volgrains  = (double *) malloc(Nabins*sizeof(double));
    abin_e     = (double *) malloc((Nabins+1)*sizeof(double));

    NfactA  = (double *) malloc(Nabins*sizeof(double));
    NfactM  = (double *) malloc(Nabins*sizeof(double));
    SfactA  = (double *) malloc(Nabins*sizeof(double));
    SfactM  = (double *) malloc(Nabins*sizeof(double));
    // populate the edges of the cells
    abin_e[0] = 0;

    for(ibin = 1; ibin < Nabins+1; ibin++){
        abin_e[ibin] = amin*pow(da,ibin-1);
    }
    // Calculate cell centers (linear)
    for(ibin = 0; ibin < Nabins; ibin ++){
        abin_c[ibin] = (abin_e[ibin]+abin_e[ibin+1])/2.;
        pi_asquare[ibin] = M_PI*pow(abin_c[ibin], 2.0);
        volgrains [ibin] = 4*M_PI*pow(abin_c[ibin], 3.0)/3.0;
    }

    // pre calculate factors
    for(ibin = 0; ibin < Nabins; ibin ++){
        aep = abin_e[ibin+1];
        ae  = abin_e[ibin];
        ac  = abin_c[ibin];
        NfactA[ibin]  = 0.5*(aep*aep-ae*ae)/(aep-ae);
        SfactA[ibin]  = pow(aep, 2)*(aep/3. - ac/2.);
        SfactA[ibin] -= pow(ae , 2)*(ae/3.  - ac/2.);

        NfactM[ibin]  = (pow(aep, 4) - pow(ae,4))/(4*(aep-ae));
        SfactM[ibin]  = pow(aep,4)*(aep/5. - ac/4.);
        SfactM[ibin] -= pow(ae ,4)*(ae/5.  - ac/4.);
    }
    
    // allocate all arrays, these also have ghost cells
    // If we have both silicates and graphite, we need to add two sets of ghost cells
    if( fSi>0 && fSi<1.0){
        dust_nbins = NdustBins + 4*dust_nghost;
    } else {
        dust_nbins = NdustBins + 2*dust_nghost;
    }

    Natoms       = (double *) malloc(dust_nbins*sizeof(double));  
    dadt         = (double *) malloc(dust_nbins*sizeof(double));  
    dadt_fixed   = (double *) malloc(dust_nbins*sizeof(double));  
    number       = (double *) malloc(dust_nbins*sizeof(double)); 
    slope        = (double *) malloc(dust_nbins*sizeof(double)); 
    Mnew         = (double *) malloc(dust_nbins*sizeof(double)); 
    Nnew         = (double *) malloc(dust_nbins*sizeof(double)); 
    Snew         = (double *) malloc(dust_nbins*sizeof(double)); 
    
    dust_vrel = (double *) malloc(dust_nbins*sizeof(double)); 
        
    // determine where the cutoff point is for silicates
    if(fSi < 1.0){
        isilicone = Nabins;
    } else {
        isilicone = 0;
    }

    for(ibin = 0; ibin < dust_nbins; ibin++){
        if(ibin < isilicone){
            norm = rho_c/aveMatom_c;
            iabin = ibin;
        } else {
            norm = rho_s/aveMatom_s;
            iabin = ibin - isilicone;
        }
        Natoms[ibin] = volgrains[iabin] * norm;
    }

    // Set initial distribution
    getintegerdustpar("dust_initDist", &initDist);
    printf("INITDIST %d", initDist);
    if(initDist == 0){
        ierr = initDustPowerLaw();
        if(ierr < 0) {
            return -1;
        }
    }

    // get parameters for dust growth
    getintegerdustpar("dust_dadt_mode", &dadt_mode);
    getrealdustpar("dust_dadt_c", &dadt_c);
    getrealdustpar("dust_dadt_tscale", &dadt_tscale);

    // timestep limiters
    getrealdustpar("dust_dadt_lim", &dadt_lim);
    getintegerdustpar("dust_maxSubSteps", &dust_maxSubSteps);
    getintegerdustpar("dust_lowerBound_pileUp", &dust_lowerBound_pileUp);
    getintegerdustpar("dust_upperBound_pileUp", &dust_upperBound_pileUp);
    
    getintegerdustpar("dust_useRadiation", &dust_useRadiation);
    getintegerdustpar("dust_useSublimation", &dust_useSublimation);
    getintegerdustpar("dust_useSputtering", &dust_useSublimation);
    if(dust_useRadiation){
        
        // allocate  
        ierr = loadDustRadiationTables();
        if(ierr < 0) {
            return -1;
        }
        
        ida_tabQabs  = (int *) malloc(dust_nbins*sizeof(int));
        ida_tabQem   = (int *) malloc(dust_nbins*sizeof(int));
        for(ibin = 0; ibin < dust_nbins; ibin++){
            ida_tabQabs[ibin] = 0;
            ida_tabQem[ibin]  = 0;
        }
        for(idx = 0; idx < NdustBins; idx++){
            ibin = globalToLocalIndex(idx);
            if(ibin < isilicone){
                graphite = 1;
                iabin = ibin;
            } else {
                graphite = 0;
                iabin = ibin - isilicone;
            }
            ida_tabQabs[ibin] = getQabs_ida(abin_c[iabin], graphite);
            ida_tabQem[ibin]  = getQemAve_ida(abin_c[iabin], graphite);
        }

        if(dust_useSublimation){
            getrealdustpar("dust_dadt_min", &dadt_min);
            getintegerdustpar("dust_NtempBins", &NtempBins);
            initTemperatureDist(NtempBins);
        }
    }
    
    if(dust_useSputtering){
        ierr = loadDustSputteringTables();
        if(ierr < 0) {
            return -1;
        }
        ida_tabSput   = (int *) malloc(dust_nbins*sizeof(int));
        for(idx = 0; idx < NdustBins; idx++){
            ibin = globalToLocalIndex(idx);
            if(ibin < isilicone){
                graphite = 1;
                iabin = ibin;
            } else {
                graphite = 0;
                iabin = ibin - isilicone;
            }
            ida_tabSput[ibin] = getSputYield_ida(abin_c[iabin], graphite); 
        }
    }

    return 1;
}


// Set initial cell dust data, normalised such that everything in iabin >= 1 has a total mass matching the dustMass
int setCellInit(int icell, double dustMass){
    int idx, ibin, iabin, graphite;
    double mass;
    double Mtot_s, Mtot_c;
    double norm_s, norm_c, norm, Sj;

    Mtot_c = 0;
    Mtot_s = 0;

    for(idx = 0; idx < NdustBins; idx++){
        ibin = globalToLocalIndex(idx);
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
        ustate[icell*nvar + IDUST_START + idx*NdustVar]   = mass;
    }
    
    norm_s = fSi*dustMass/Mtot_s;
    norm_c = (1-fSi)*dustMass/Mtot_c;

    for(idx = 0; idx < NdustBins; idx++){
        ibin = globalToLocalIndex(idx);
        if(ibin < isilicone){
            graphite = 1;
        } else {
            graphite = 0;
        }

        // Mass (density) in bin
        mass = ustate[icell*nvar + IDUST_START + idx*NdustVar];
        Sj   = slope[ibin];
        // Normalize to match mass density
        if(graphite){
            mass = mass * norm_c;
            Sj   = Sj * norm_c; 
        } else {
            mass = mass * norm_s;
            Sj   = Sj * norm_s;
        }
        ustate[icell*nvar + IDUST_START + idx*NdustVar] = mass;
        ustate[icell*nvar + IDUST_START + idx*NdustVar+1] = Sj;
    }

    return 1;
}
#endif
