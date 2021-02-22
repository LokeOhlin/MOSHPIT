#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "hydro.h"
#include "cgeneral.h"
#include "radchem.h"
#ifdef useDust
    #include "dust.h"
#endif

int nrealHydroPars = 15;
real_list_t *hydroDPars = NULL;
int nintHydroPars = 7;
int_list_t *hydroIPars = NULL;


double numd_init;
double temp_init;
double xH_init;
double xH2_init;
double xHp_init;
double xCO_init;

double adi;
double courant_number;

int geometry; // 0 = 1D carteesian
                  // 1 = 1D Spherical 
int setup;

int left_bound ;
double bdensL  ; 
double bvelL   ; 
double benerL  ; 

int right_bound ;
double bdensR   ; 
double bvelR    ; 
double benerR   ; 

int NCELLS;
int NGHOST;
int NINTER;

double *ustate = NULL;
double *pstate = NULL;
double *rs = NULL;
double *vol = NULL;
double *dr = NULL ;

int useHydro; 

int setHydroPars(){
    // Allocate spaces
    hydroDPars = (real_list_t *) malloc(nrealHydroPars * sizeof(real_list_t));
    hydroIPars = (int_list_t *) malloc(nintHydroPars * sizeof(int_list_t));
    
    // Real/Double parameters 
    strcpy(hydroDPars[0].name, "gamma"); hydroDPars[0].value = 5./3.; 
    strcpy(hydroDPars[1].name, "courant_number"); hydroDPars[1].value = 0.5; 

    strcpy(hydroDPars[2].name, "numd_init"); hydroDPars[2].value = 100.; 
    strcpy(hydroDPars[3].name, "temp_init"); hydroDPars[3].value = 100.; 
    strcpy(hydroDPars[4].name, "xH_init"); hydroDPars[4].value = .1; 
    strcpy(hydroDPars[5].name, "xH2_init"); hydroDPars[5].value = .4; 
    strcpy(hydroDPars[6].name, "xCO_init"); hydroDPars[6].value = 1.e-7;

    strcpy(hydroDPars[7].name, "rL"); hydroDPars[7].value = 0.; 
    strcpy(hydroDPars[8].name, "rR"); hydroDPars[8].value = 100.; 
    
    strcpy(hydroDPars[9].name, "bdensR"); hydroDPars[9].value = 1.; 
    strcpy(hydroDPars[10].name, "bvelR"); hydroDPars[10].value = 1.; 
    strcpy(hydroDPars[11].name, "benerR"); hydroDPars[11].value = 1.; 
    
    strcpy(hydroDPars[12].name, "bdensL"); hydroDPars[12].value = 1.; 
    strcpy(hydroDPars[13].name, "bvelL"); hydroDPars[13].value = 1.; 
    strcpy(hydroDPars[14].name, "benerL"); hydroDPars[14].value = 1.; 
     
    // Ints
    strcpy(hydroIPars[0].name, "ncells");      hydroIPars[0].value = 512;
    strcpy(hydroIPars[1].name, "nghost");      hydroIPars[1].value = 2;
    strcpy(hydroIPars[2].name, "geometry");    hydroIPars[2].value = 1;
    strcpy(hydroIPars[3].name, "left_bound");  hydroIPars[3].value = 0;
    strcpy(hydroIPars[4].name, "right_bound"); hydroIPars[4].value = 1;
    strcpy(hydroIPars[5].name, "setup");       hydroIPars[5].value = 0;
    strcpy(hydroIPars[6].name, "useHydro");    hydroIPars[6].value = 1;
   
    return 1;
}

int checkHydroPars(char *name, char *value){
    // Method to see if parameter matches any of the defined chemistry parameters
    int ipar;
    // Check reals
    for(ipar = 0; ipar < nrealHydroPars; ipar ++){
        if(compStr(name, hydroDPars[ipar].name, 128) > 0){
            hydroDPars[ipar].value = atof(value);
return 1;
        }
    }

    for(ipar = 0; ipar < nintHydroPars; ipar ++){
        if(compStr(name, hydroIPars[ipar].name, 128) > 0){
            hydroIPars[ipar].value = atoi(value);
            return 1;
        }
    }
    return -1;
}

// Method to get value of a double(real) parameter
void getrealhydropar(char *name, double *value){
    int ipar;
    for(ipar = 0; ipar < nrealHydroPars; ipar ++){
        if(compStr(name, hydroDPars[ipar].name, 80) > 0){
            *value = hydroDPars[ipar].value;
        }
    }
}

// Method to get value of an integer parameter
void getintegerhydropar(char *name, int *value){
    int ipar;
    for(ipar = 0; ipar < nintHydroPars; ipar ++){
        if(compStr(name, hydroIPars[ipar].name, 80) > 0){
            *value = hydroIPars[ipar].value;
        }
    }
}


// gets global parameters
int initHydro(){
    // Grid parameters
    getintegerhydropar("ncells", &NCELLS);
    getintegerhydropar("nghost", &NGHOST);
    // add ghost cells to total
    NCELLS = NCELLS + 2*NGHOST;
    NINTER = NCELLS - 3;

    // General hydro stuff
    getrealhydropar("gamma", &adi);
    getrealhydropar("courant_number", &courant_number);
    // Simulation selection
    getintegerhydropar("setup", &setup);

    // Boundary conditions
    getintegerhydropar("left_bound", &left_bound);
    if(left_bound == 2){ //fixed BCs : load BCs!
        getrealhydropar("bdensL", &bdensL);
        getrealhydropar("bvelL", &bvelL);
        getrealhydropar("benerL", &benerL);
    }

    getintegerhydropar("right_bound", &right_bound);
    if(right_bound == 2){ //fixed BCs : load BCs!
        getrealhydropar("bdensR", &bdensR);
        getrealhydropar("bvelR", &bvelR);
        getrealhydropar("benerR", &benerR);
    }

    getintegerhydropar("useHydro", &useHydro);
    return 1;
}

int init_uniform(){
    int icell,idx;
    double dens, eint, xH, xH2, rR, velinit= 0;
    getrealhydropar("numd_init", &numd_init);
    getrealhydropar("temp_init", &temp_init);
    getrealhydropar("xH_init", &xH_init);
    getrealhydropar("xH2_init", &xH2_init);
    getrealhydropar("xCO_init", &xCO_init);
    
    getrealhydropar("rR", &rR);
    
    dens = numd_init * ch_mH * abar;
    eint = numd_init*ch_kb*temp_init/(adi-1);

    for(icell = NGHOST; icell < NCELLS-NGHOST; icell++){
        idx = icell*nvar;
        ustate[idx] = dens;
        ustate[idx+1] = velinit*dens;
        ustate[idx+2] = eint + 0.5*dens*velinit*velinit;
#ifdef useChemistry
        xH  = xH_init ;// *(1-rs[icell]/rR);
        xH2 = xH2_init;// *rs[icell]/rR;
        ustate[idx+ICHEM_START]   = xH/mf_scale;
        ustate[idx+ICHEM_START+1] = 2*xH2/mf_scale;
        ustate[idx+ICHEM_START+2] = (1-xH-2*xH2)/mf_scale;
        ustate[idx+ICHEM_START+3] = xCO_init/mf_scale;
        ustate[idx+ICHEM_START+4] = (abundC-xCO_init)/mf_scale;
        ustate[idx+ICHEM_END] = 10.0;
#endif
    }
    return 1;
}


int init_leftwave(){
    int icell;
    double p0 = 1;
    double r1 = 1.24e4;
    double r0 = 1e5;
    double u0 = 0.5*sqrt(adi*p0/r0);
    double p1 = 1;
    double u1 = 0.5*sqrt(adi*p0/r0);

    left_bound = 2;
    bdensL = r0;
    bvelL  = r0*u0;
    benerL = r0*(p0/(r0*(adi-1))+0.5*u0*u0);

    right_bound = 2;
    bdensR = r1;
    bvelR  = r1*u1;
    benerR = r1*(p1/(r1*(adi-1))+0.5*u1*u1);
    geometry = 0; // 0 = 1D carteesian
    
    for(icell = 0; icell < NCELLS; icell++ ){
        if( icell > 2+ (NCELLS-4)/2){
            ustate[icell*nvar] = r1;
            ustate[icell*nvar + 1] = r1*u1;
            ustate[icell*nvar + 2] = r1*(p1/(r1*(adi-1))+0.5*u1*u1);
        }
        else{
            ustate[icell*nvar] = r0;
            ustate[icell*nvar +1 ] = r0*u0;
            ustate[icell*nvar +2 ] = r0*(p0/(r0*(adi-1))+0.5*u0*u0);
        }
    }
    return -1;
}

int init_sedov(){
    int icell, idx;
    double mH = 1.675e-24, kb = 1.3806e-16;
    double ESN=1e51;
    double num=100, dens;
    double temp=100, ener; 
    double xH, xH2;    
    getrealhydropar("xH_init", &xH_init);
    getrealhydropar("xH2_init", &xH2_init);
    getrealhydropar("xCO_init", &xCO_init);
    
    left_bound = 0;
    right_bound = 1;
    geometry = 1; // spherical
    
    ESN  = ESN;   // unit Msun*pc**2/(Myr)**2
    //ESN  = ESN/(4*M_PI*pow(dr[2]+dr[3],3)/3.); // make into energy density
    ESN  = ESN/(4*M_PI*pow(dr[2],3)/3.); // make into energy density
    dens = num*(mH*1.2); // units Msun/pc**3
    ener = num*kb*temp/(adi-1);
    for(icell = 0; icell < NCELLS; icell++ ){
        ustate[icell*nvar]   = dens;
        ustate[icell*nvar+1] = 0;
        ustate[icell*nvar+2] = ener;
        if(icell == 2){
            ustate[icell*nvar+2] = ener + ESN;
        }

#ifdef useChemistry
        idx = icell*nvar;
        xH  = xH_init ;// *(1-rs[icell]/rR);
        xH2 = xH2_init;// *rs[icell]/rR;
        ustate[idx+ICHEM_START]   = xH/mf_scale;
        ustate[idx+ICHEM_START+1] = 2*xH2/mf_scale;
        ustate[idx+ICHEM_START+2] = (1-xH-2*xH2)/mf_scale;
        ustate[idx+ICHEM_START+3] = xCO_init/mf_scale;
        ustate[idx+ICHEM_START+4] = (abundC-xCO_init)/mf_scale;
        ustate[idx+ICHEM_END] = 10.0;
#endif
    }
    return -1;
}

int init_grid(){
    int icell;
    double rL, rR, L, dr_const, rp, rm;
    getintegerhydropar("geometry", &geometry);
    getrealhydropar("rL", &rL);
    getrealhydropar("rR", &rR);
    L = rR -rL;
    dr_const = L/(NCELLS - 2*NGHOST);
    
    for(icell = 0; icell < NCELLS; icell++){
        rs[icell]  = dr_const*(icell - NGHOST +1./2.);

        rp = rs[icell] + 0.5*dr_const;
        rm = rs[icell] - 0.5*dr_const;
        
        dr[icell]  = dr_const;
        if(geometry == 1){
            vol[icell] = 4.0*M_PI*(pow(rp,3.0)-pow(rm,3.0))/3.0;
        } else {
            vol[icell] = dr_const*dr_const*dr_const;
        }
    }
   return 1; 
}
int init_domain(){
    int ierr = -1;
    int icell;
    // Allocate arrays
    ustate = (double *) malloc(NCELLS*nvar*sizeof(double));
    pstate = (double *) malloc(NCELLS*nvar*sizeof(double));
    rs     = (double *) malloc(NCELLS*sizeof(double));
    dr     = (double *) malloc(NCELLS*sizeof(double));
    vol    = (double *) malloc(NCELLS*sizeof(double));

    ierr = init_grid();
    if(setup == 0){
        ierr = init_uniform();
    } else if(setup == 1){
        ierr = init_leftwave();
    } else if(setup == 2){
        ierr = init_sedov();
    } else {
        printf("setup %d not recognised", setup);
        return -1;
    }
    if(ierr < 0) {
        return -1;
    }
#ifdef useDust
    double dustMass;
    double dust_to_gas_ratio;
    getrealchemistrypar("ch_dust_to_gas_ratio", &dust_to_gas_ratio);
    for(icell = NGHOST; icell < NCELLS - NGHOST; icell++){
        dustMass = dust_to_gas_ratio * ustate[icell*nvar];
        setCellInit(icell, dustMass);
    } 
#endif 
    return 1;
}





