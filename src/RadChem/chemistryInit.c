#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "hydro.h"
#include "radchem.h"
#include "cgeneral.h"
#ifdef useDust
    #include "dust.h"
#endif

int nrealPars = 26;
real_list_t *chemDPars = NULL;
int nintPars = 13;
int_list_t *chemIPars = NULL;
//int nstrPars = 10;
//str_list_t chemSPars[];


double abar;
double mf_scale;
double abundHe, abundC, abundO, abundSi;
double ch_mH   = 1.672621637e-24;
double ch_muHe = 4.002602;
double ch_muC  = 12.011;
double ch_muO  = 15.9994;
double ch_muSi = 28.0855;
double ch_kb   = 1.38065e-16;
// set default parameter lists 
int setChemistryPars(){
    // Allocate spaces
    chemDPars = (real_list_t *) malloc(nrealPars * sizeof(real_list_t));
    chemIPars = (int_list_t *) malloc(nintPars * sizeof(int_list_t));

    // Real/Double parameters
    strcpy(chemDPars[0].name, "ch_deff");                  chemDPars[0].value = 1.0; 
    strcpy(chemDPars[1].name, "ch_abundHe");               chemDPars[1].value = 0.1; 
    strcpy(chemDPars[2].name, "ch_abundC");                chemDPars[2].value = 1.4e-4; 
    strcpy(chemDPars[3].name, "ch_abundO");                chemDPars[3].value = 3.2e-4;
    strcpy(chemDPars[4].name, "ch_abundSi");               chemDPars[4].value = 1.5e-5;
    strcpy(chemDPars[5].name, "ch_abundD");                chemDPars[5].value = 2.6e-5;
    strcpy(chemDPars[6].name, "ch_abundM");                chemDPars[6].value = 1e-7;
    strcpy(chemDPars[7].name, "ch_abundN");                chemDPars[7].value = 1.5e-5;
    strcpy(chemDPars[8].name, "ch_G0");                    chemDPars[8].value = 1.7e0;
    strcpy(chemDPars[9].name, "ch_f_rsc");                 chemDPars[9].value = 1e0;
    strcpy(chemDPars[10].name, "ch_phi_pah");              chemDPars[10].value = 0.5;
    strcpy(chemDPars[11].name, "ch_tdust");                chemDPars[11].value = 1e1;
    strcpy(chemDPars[12].name, "ch_dust_to_gas_ratio");    chemDPars[12].value = 1e0;
    strcpy(chemDPars[13].name, "ch_AV_conversion_factor"); chemDPars[13].value = 5.348e-22;
    strcpy(chemDPars[14].name, "ch_cosmic_ray_ion_rate");  chemDPars[14].value = 3e-17;
    strcpy(chemDPars[15].name, "ch_redshift");             chemDPars[15].value = 0e0;
    strcpy(chemDPars[16].name, "ch_NH_ext");               chemDPars[16].value = 1e20;
    strcpy(chemDPars[17].name, "ch_h2_form_ex");           chemDPars[17].value = 0e0;
    strcpy(chemDPars[18].name, "ch_h2_form_kin");          chemDPars[18].value = 0e0;
    strcpy(chemDPars[19].name, "ch_xray_scaling");         chemDPars[19].value = 1e0;
    strcpy(chemDPars[20].name, "ch_Z_atom");               chemDPars[20].value = 1e0;
    strcpy(chemDPars[21].name, "ch_max_ion_frac_change");  chemDPars[21].value = 1e-1;
    strcpy(chemDPars[22].name, "Tstar");  chemDPars[22].value = 5e3;
    strcpy(chemDPars[23].name, "Lstar");  chemDPars[23].value = 1e5;
    strcpy(chemDPars[24].name, "radEmin");  chemDPars[24].value = 1.6e-17;
    strcpy(chemDPars[25].name, "radEmax");  chemDPars[25].value = 1.6e-10;
    
    // Integer parameters
    strcpy(chemIPars[0].name, "ch_iphoto");               chemIPars[0].value =  0;
    strcpy(chemIPars[1].name, "ch_iflag_mn");             chemIPars[1].value =  1;
    strcpy(chemIPars[2].name, "ch_iflag_ad");             chemIPars[2].value =  1;
    strcpy(chemIPars[3].name, "ch_iflag_atom");           chemIPars[3].value =  3;
    strcpy(chemIPars[4].name, "ch_iflag_3bh2a");          chemIPars[4].value =  1;
    strcpy(chemIPars[5].name, "ch_iflag_3bh2b");          chemIPars[5].value =  1;
    strcpy(chemIPars[6].name, "ch_iflag_h3pra");          chemIPars[6].value =  1;
    strcpy(chemIPars[7].name, "ch_isrf_option");          chemIPars[7].value =  1;
    strcpy(chemIPars[8].name, "ch_no_chem");              chemIPars[8].value =  0;
    strcpy(chemIPars[9].name, "ch_use_photo_eqb_approx"); chemIPars[9].value =  0;
    strcpy(chemIPars[10].name, "radiationPressure");      chemIPars[10].value =  1;
    strcpy(chemIPars[11].name, "numBinsSubIon");       chemIPars[11].value = 30;
    strcpy(chemIPars[12].name, "numBinsFullIon");      chemIPars[12].value = 10;
    return 1;
}

int checkChemistryPars(char *name, char *value){
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
void getrealchemistrypar(char *name, double *value){
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
void getintegerchemistrypar(char *name, int *value){
    int ipar;
    for(ipar = 0; ipar < nintPars; ipar ++){
        if(compStr(name, chemIPars[ipar].name, 80) > 0){
            *value = chemIPars[ipar].value;
            return;
        }
    }
    printf("%s NOT FOUND\n", name);
}

int initChemistry(){
    getrealchemistrypar("ch_abundHe", &abundHe);
    getrealchemistrypar("ch_abundC", &abundC);
    getrealchemistrypar("ch_abundO", &abundO);
    getrealchemistrypar("ch_abundSi",&abundSi);
    
    abar = 1.0 + abundHe*ch_muHe + abundC*ch_muC + abundO*ch_muO + abundSi*ch_muSi ;
    mf_scale = 1.0 + abundC*ch_muC;
    
    //Call intit for fortran functions
    Chemistry_FortranInit();
    
    //Init radiation
    initRadiation();
#ifdef useDust 
    //Init Dust
    initDust();
#endif    
    return 1;

}


