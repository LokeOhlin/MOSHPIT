#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <hydro.h>
#include <cgeneral.h>
#include <radchem.h>
#ifdef useDust
    #include <dust.h>
#endif

int nrealHydroPars = 31;
real_list_t *hydroDPars = NULL;
int nintHydroPars = 10;
int_list_t *hydroIPars = NULL;

int init_numd;
double dens_init;
double numd_init;
double temp_init;
double velo_init;
double xH_init;
double xH2_init;
double xHp_init;
double xCO_init;
double cs_init;

double adi;
double courant_number;
double roe_p;

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

// threshhold fraction to use approximate Eint solution rather than total energy when
// Eint < Etot*hy_ethresh
double hy_ethresh;

int NCELLS;
int NGHOST;
int NINTER;

double *ustate = NULL;
double *pstate = NULL;
double *rs = NULL;
double *vol = NULL;
double *dr = NULL ;
double *varBuff = NULL;

int useHydro; 
int isothermal;

int setHydroPars(){
    // Allocate spaces
    hydroDPars = (real_list_t *) malloc(nrealHydroPars * sizeof(real_list_t));
    hydroIPars = (int_list_t *) malloc(nintHydroPars * sizeof(int_list_t));
    
    // Real/Double parameters 
    strcpy(hydroDPars[0].name, "gamma"); hydroDPars[0].value = 5./3.; 
    strcpy(hydroDPars[1].name, "courant_number"); hydroDPars[1].value = 0.5; 

    strcpy(hydroDPars[2].name, "dens_init");    hydroDPars[2].value = 100.; 
    strcpy(hydroDPars[3].name, "temp_init");    hydroDPars[3].value = 100.; 
    strcpy(hydroDPars[4].name, "velo_init");    hydroDPars[4].value = 0.; 
    strcpy(hydroDPars[5].name, "xH_init");      hydroDPars[5].value = .1; 
    strcpy(hydroDPars[6].name, "xH2_init");     hydroDPars[6].value = .4; 
    strcpy(hydroDPars[7].name, "xCO_init");     hydroDPars[7].value = 1.e-7;
    strcpy(hydroDPars[8].name, "cs_init");      hydroDPars[8].value = -1;

    strcpy(hydroDPars[9].name, "rL");           hydroDPars[9].value = 0.; 
    strcpy(hydroDPars[10].name, "rR");          hydroDPars[10].value = 100.; 
    strcpy(hydroDPars[11].name, "orth_extent"); hydroDPars[11].value = -1.0; 
    
    strcpy(hydroDPars[12].name, "bdensR");      hydroDPars[12].value = 1.; 
    strcpy(hydroDPars[13].name, "bvelR");       hydroDPars[13].value = 1.; 
    strcpy(hydroDPars[14].name, "benerR");      hydroDPars[14].value = 1.; 
    
    strcpy(hydroDPars[15].name, "bdensL");      hydroDPars[15].value = 1.; 
    strcpy(hydroDPars[16].name, "bvelL");       hydroDPars[16].value = 1.; 
    strcpy(hydroDPars[17].name, "benerL");      hydroDPars[17].value = 1.; 
    
    strcpy(hydroDPars[18].name, "hy_ethresh");  hydroDPars[18].value = 0.01; 
    
    strcpy(hydroDPars[19].name, "wave_vel_amplitude");  hydroDPars[19].value = 1.0; 
    strcpy(hydroDPars[20].name, "wave_rho_amplitude");  hydroDPars[20].value = 1.0; 
    strcpy(hydroDPars[21].name, "wave_vel_wavelength"); hydroDPars[21].value = 1.0; 
    strcpy(hydroDPars[22].name, "wave_rho_wavelength"); hydroDPars[22].value = 1.0; 
    strcpy(hydroDPars[23].name, "wave_vel_phase");      hydroDPars[23].value = 0.0; 
    strcpy(hydroDPars[24].name, "wave_rho_phase");      hydroDPars[24].value = 0.0; 
    
    strcpy(hydroDPars[25].name, "roe_p");      hydroDPars[25].value = 1.5; 

    strcpy(hydroDPars[26].name, "dens_initL");    hydroDPars[26].value = 100.; 
    strcpy(hydroDPars[27].name, "temp_initL");    hydroDPars[27].value = 100.; 
    strcpy(hydroDPars[28].name, "velo_initL");    hydroDPars[28].value = 0.; 
    strcpy(hydroDPars[29].name, "cs_initL");      hydroDPars[29].value = -1;
    strcpy(hydroDPars[30].name, "ESN_init");      hydroDPars[30].value = 1e51;
    
    // Ints
    strcpy(hydroIPars[0].name, "ncells");      hydroIPars[0].value = 512;
    strcpy(hydroIPars[1].name, "nghost");      hydroIPars[1].value = 2;
    strcpy(hydroIPars[2].name, "geometry");    hydroIPars[2].value = 1;
    strcpy(hydroIPars[3].name, "left_bound");  hydroIPars[3].value = 0;
    strcpy(hydroIPars[4].name, "right_bound"); hydroIPars[4].value = 1;
    strcpy(hydroIPars[5].name, "logspace_cells");       hydroIPars[5].value = 0;
    strcpy(hydroIPars[6].name, "setup");       hydroIPars[6].value = 0;
    strcpy(hydroIPars[7].name, "useHydro");    hydroIPars[7].value = 1;
    strcpy(hydroIPars[8].name, "isothermal");    hydroIPars[8].value = 0;
    strcpy(hydroIPars[9].name, "init_numd");    hydroIPars[8].value = 1;
   
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
#ifdef useDustDynamics
    if(NGHOST < 3){
        printf("NOTE: the ROE solver used for dust advection requires NGHOST>=3\n");
        printf("    : setting NGHOST = 3\n");
        NGHOST = 3;
    }
    getrealhydropar("roe_p", &roe_p);
#endif 
    //buffer used in outputs
    varBuff = (double *) malloc(NCELLS*sizeof(double));
    
    // add ghost cells to total
    NCELLS = NCELLS + 2*NGHOST;
    NINTER = NCELLS - 2*NGHOST + 1;

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
    
    // If one is periodic both better be periodic....
    if( (left_bound == 3 || right_bound == 3) && !(left_bound == 3 && right_bound == 3)){
        printf("One boundary condition was chosen as periodic and the other not.\n");
        printf("Need either both or none to be periodic... \n");
        exit(0);
    }
     
    getrealhydropar("hy_ethresh", &hy_ethresh);

    getintegerhydropar("useHydro", &useHydro);
    return 1;
}

int init_uniform(){
    int icell,idx;
    double eint, xH, xH2, rR;
    getintegerhydropar("init_numd", &init_numd);
    getrealhydropar("dens_init", &dens_init);
    getrealhydropar("temp_init", &temp_init);
    getrealhydropar("velo_init", &velo_init);
    getrealhydropar("xH_init", &xH_init);
    getrealhydropar("xH2_init", &xH2_init);
    getrealhydropar("xCO_init", &xCO_init);
    getrealhydropar("cs_init", &cs_init);
    
    getrealhydropar("rR", &rR);
    


#ifdef useDust
    double dustMass;
    double dust_to_gas_ratio;
    getrealchemistrypar("ch_dust_to_gas_ratio", &dust_to_gas_ratio);
#endif 
    if(init_numd){
        numd_init = dens_init;
#ifdef TESTATOMONLY 
        dens_init = dens_init * ch_mH;
#else
        dens_init = dens_init * ch_mH * abar;
#endif
    } else {
#ifdef TESTATOMONLY 
        numd_init = dens_init / ch_mH;
#else
        numd_init = dens_init / (ch_mH * abar);
#endif
    }
    
    if(isothermal){
        //if isothermal, keep token eint
        eint = dens_init * cs_init * cs_init /(adi*(adi-1));
    } else {
        eint = numd_init*ch_kb*temp_init/(adi-1);
    }
    for(icell = NGHOST; icell < NCELLS-NGHOST; icell++){
        idx = icell*nvar;
        ustate[idx] = dens_init;
        ustate[idx+1] = velo_init*dens_init;
        ustate[idx+2] = eint + 0.5*dens_init*velo_init*velo_init;
#ifdef useChemistry
        xH  = xH_init ;// *(1-rs[icell]/rR);
        xH2 = xH2_init;// *rs[icell]/rR;
        ustate[idx+ICHEM_START]   = dens_init*xH/mf_scale;
        ustate[idx+ICHEM_START+1] = dens_init*2*xH2/mf_scale;
        ustate[idx+ICHEM_START+2] = dens_init*(1-xH-2*xH2)/mf_scale;
        ustate[idx+ICHEM_START+3] = dens_init*xCO_init/mf_scale;
        ustate[idx+ICHEM_START+4] = dens_init*(abundC-xCO_init)/mf_scale;
        ustate[idx+ICHEM_END] = 10.0;
#endif
#ifdef useDust
        dustMass = dust_to_gas_ratio * 1e-2 * ustate[icell*nvar];
        setCellInit(icell, dustMass);
#endif
    }
    return 1;
}


int init_leftwave(){
    int icell;
    double dens_initL, temp_initL, velo_initL, cs_initL, numd_initL, eint, eintL; 
    getintegerhydropar("init_numd", &init_numd);
    
    getrealhydropar("dens_init", &dens_init);
    getrealhydropar("temp_init", &temp_init);
    getrealhydropar("velo_init", &velo_init);
    getrealhydropar("xH_init", &xH_init);
    getrealhydropar("xH2_init", &xH2_init);
    getrealhydropar("xCO_init", &xCO_init);
    getrealhydropar("cs_init", &cs_init);
    
    getrealhydropar("dens_initL", &dens_initL);
    getrealhydropar("temp_initL", &temp_initL);
    getrealhydropar("velo_initL", &velo_initL);
    getrealhydropar("cs_initL", &cs_initL);
    
    if(init_numd){
        numd_init  = dens_init;
        numd_initL = dens_initL;
#ifdef TESTATOMONLY 
        dens_init  = dens_init * ch_mH;
        dens_initL = dens_initL * ch_mH;
#else
        dens_init  = dens_init * ch_mH * abar;
        dens_initL = dens_initL * ch_mH * abar;
#endif
    } else {
#ifdef TESTATOMONLY 
        numd_init  = dens_init / ch_mH;
        numd_initL = dens_initL / ch_mH;
#else
        numd_init  = dens_init / (ch_mH * abar);
        numd_initL = dens_initL / (ch_mH * abar);
#endif
    }
    
    if(isothermal){
        //if isothermal, keep token eint
        eint  = dens_init * cs_init * cs_init /(adi*(adi-1));
        eintL = dens_initL * cs_initL * cs_initL /(adi*(adi-1));
    } else {
        eint  = numd_init *ch_kb * temp_init/(adi-1);
        eintL = numd_initL*ch_kb * temp_initL/(adi-1);
    }

    left_bound = 0;
    bdensL = dens_initL;
    bvelL  = dens_initL*velo_initL;
    benerL = eintL + 0.5*dens_initL*velo_initL*velo_initL;

    right_bound = 2;
    bdensR = dens_init;
    bvelR  = dens_init*velo_init;
    benerR = eint + 0.5*dens_init*velo_init*velo_init;
#ifdef useDust
    double dustMass, dust_velo_initL;
    double dust_to_gas_ratio;
    int idx, idust;
    getrealchemistrypar("ch_dust_to_gas_ratio", &dust_to_gas_ratio);
    getrealdustpar("dust_velo_initL", &dust_velo_initL);
#endif 
    
    for(icell = 0; icell < NCELLS; icell++ ){
        if( icell > NGHOST + (NCELLS-2*NGHOST)/2){
            ustate[icell*nvar] = dens_init;
            ustate[icell*nvar + 1] = dens_init*velo_init;
            ustate[icell*nvar + 2] = eint + 0.5*dens_init*velo_init*velo_init;
        }
        else{
            ustate[icell*nvar] = dens_initL;
            ustate[icell*nvar +1] = dens_initL*velo_initL;
            ustate[icell*nvar +2] = eintL + 0.5*dens_initL*velo_initL*velo_initL;
        }
#ifdef useChemistry
        idx = icell*nvar;
        ustate[idx+ICHEM_START]   = ustate[icell*nvar]*xH_init/mf_scale;
        ustate[idx+ICHEM_START+1] = ustate[icell*nvar]*2*xH2_init/mf_scale;
        ustate[idx+ICHEM_START+2] = ustate[icell*nvar]*(1-xH_init-2*xH2_init)/mf_scale;
        ustate[idx+ICHEM_START+3] = ustate[icell*nvar]*xCO_init/mf_scale;
        ustate[idx+ICHEM_START+4] = ustate[icell*nvar]*(abundC-xCO_init)/mf_scale;
        ustate[idx+ICHEM_END] = 10.0;
#endif
#ifdef useDust
        dustMass = dust_to_gas_ratio * 1e-2 * ustate[icell*nvar];
        setCellInit(icell, dustMass);
        // setCellInit already does the right cell correctly
        if(icell <= 2 + (NCELLS-2*NGHOST)/2){
            for(idust = 0; idust < NdustBins; idust++){
                idx = icell*nvar + IDUST_START + idust*NdustVar;
                ustate[idx + 2] = ustate[idx] * dust_velo_initL;
            }
        }
#endif
    }
    return 1;
}

int init_sedov(){
    int icell, idx;
    double ESN_init;
    double ener, eint; 
    getrealhydropar("dens_init", &dens_init);
    getrealhydropar("temp_init", &temp_init);
    getrealhydropar("velo_init", &velo_init);
    getrealhydropar("xH_init", &xH_init);
    getrealhydropar("xH2_init", &xH2_init);
    getrealhydropar("xCO_init", &xCO_init);
    getrealhydropar("cs_init", &cs_init);
    getrealhydropar("ESN_init", &ESN_init);

    if(init_numd){
        numd_init = dens_init;
#ifdef TESTATOMONLY 
        dens_init = dens_init * ch_mH;
#else
        dens_init = dens_init * ch_mH * abar;
#endif
    } else {
#ifdef TESTATOMONLY 
        numd_init = dens_init / ch_mH;
#else
        numd_init = dens_init / (ch_mH * abar);
#endif
    }
    
    
    if(isothermal){
        //if isothermal, keep token eint
        eint = dens_init * cs_init * cs_init /(adi*(adi-1));
    } else {
        eint = numd_init*ch_kb*temp_init/(adi-1);
    }
    ener = eint + 0.5*dens_init*velo_init*velo_init;
#ifdef useDust
    double dustMass;
    double dust_to_gas_ratio;
    getrealchemistrypar("ch_dust_to_gas_ratio", &dust_to_gas_ratio);
#endif 
    
    left_bound = 0;
    right_bound = 1;
    geometry = 1; // spherical
    
    
    ESN_init  = ESN_init/(4*M_PI*pow(dr[2],3)/3.); // make into energy density
    
    for(icell = 0; icell < NCELLS; icell++ ){
        ustate[icell*nvar]   = dens_init;
        ustate[icell*nvar+1] = velo_init;
        ustate[icell*nvar+2] = ener;
        if(icell == NGHOST){
            ustate[icell*nvar+2] = ener + ESN_init;
        }

#ifdef useChemistry
        idx = icell*nvar;
        ustate[idx+ICHEM_START]   = xH_init/mf_scale;
        ustate[idx+ICHEM_START+1] = 2*xH2_init/mf_scale;
        ustate[idx+ICHEM_START+2] = (1-xH_init-2*xH2_init)/mf_scale;
        ustate[idx+ICHEM_START+3] = xCO_init/mf_scale;
        ustate[idx+ICHEM_START+4] = (abundC-xCO_init)/mf_scale;
        ustate[idx+ICHEM_END] = 10.0;
#endif
#ifdef useDust
        dustMass = dust_to_gas_ratio * 1e-2 * ustate[icell*nvar];
        setCellInit(icell, dustMass);
#endif
    }
    return 1;
}


int init_wave(){
    int icell,idx;
    double eint, xH, xH2, rL, rR;
    double velo_cell, dens_cell;
    double wave_vel_amplitude, wave_rho_amplitude;
    double wave_vel_wavelength, wave_rho_wavelength;
    double wave_vel_phase, wave_rho_phase;
    
    getintegerhydropar("init_numd", &init_numd);
    getrealhydropar("dens_init", &dens_init);
    getrealhydropar("temp_init", &temp_init);
    getrealhydropar("velo_init", &velo_init);
    getrealhydropar("xH_init", &xH_init);
    getrealhydropar("xH2_init", &xH2_init);
    getrealhydropar("xCO_init", &xCO_init);

    getrealhydropar("wave_vel_amplitude", &wave_vel_amplitude);
    getrealhydropar("wave_rho_amplitude", &wave_rho_amplitude);
    getrealhydropar("wave_vel_wavelength", &wave_vel_wavelength);
    getrealhydropar("wave_rho_wavelength", &wave_rho_wavelength);
    getrealhydropar("wave_vel_phase", &wave_vel_phase);
    getrealhydropar("wave_rho_phase", &wave_rho_phase);
    
    getrealhydropar("rL", &rL);
    getrealhydropar("rR", &rR);
#ifdef useDust
    double dustMass;
    double dust_to_gas_ratio;
    int idust, idx2;
    getrealchemistrypar("ch_dust_to_gas_ratio", &dust_to_gas_ratio);
#endif 
    
    if(init_numd){
        numd_init = dens_init;
#ifdef TESTATOMONLY 
        dens_init = dens_init * ch_mH;
#else
        dens_init = dens_init * ch_mH * abar;
#endif
    } else {
#ifdef TESTATOMONLY 
        numd_init = dens_init / ch_mH;
#else
        numd_init = dens_init / (ch_mH * abar);
#endif
    }
    
    
    if(isothermal){
        //if isothermal, keep token eint
        eint = dens_init * cs_init * cs_init /(adi*(adi-1));
    } else {
        eint = numd_init*ch_kb*temp_init/(adi-1);
    }

    for(icell = NGHOST; icell < NCELLS-NGHOST; icell++){
        idx = icell*nvar;
        velo_cell = velo_init  + wave_vel_amplitude * sin( 2*M_PI*(rs[icell] - rL)/wave_vel_wavelength + wave_vel_phase);
        dens_cell = dens_init  + wave_rho_amplitude * sin( 2*M_PI*(rs[icell] - rL)/wave_rho_wavelength + wave_rho_phase);
        ustate[idx] = dens_cell;
        ustate[idx+1] = velo_cell*dens_cell;
        ustate[idx+2] = eint + 0.5*dens_cell*velo_cell*velo_cell;
#ifdef useChemistry
        xH  = xH_init ;// *(1-rs[icell]/rR);
        xH2 = xH2_init;// *rs[icell]/rR;
        ustate[idx+ICHEM_START]   = dens_cell*xH/mf_scale;
        ustate[idx+ICHEM_START+1] = dens_cell*2*xH2/mf_scale;
        ustate[idx+ICHEM_START+2] = dens_cell*(1-xH-2*xH2)/mf_scale;
        ustate[idx+ICHEM_START+3] = dens_cell*xCO_init/mf_scale;
        ustate[idx+ICHEM_START+4] = dens_cell*(abundC-xCO_init)/mf_scale;
        ustate[idx+ICHEM_END] = 10.0;
#endif
#ifdef useDust
        dustMass = dust_to_gas_ratio * 1e-2 * ustate[icell*nvar];
        setCellInit(icell, dustMass);
        for(idust = 0; idust < NdustBins; idust++){
            idx2 = idx + IDUST_START + idust*NdustVar;
            ustate[idx2 + 2] = ustate[idx2 + 2] + ustate[idx2] * wave_vel_amplitude * sin( 2*M_PI*(rs[icell] - rL)/wave_vel_wavelength + wave_vel_phase);
        }
#endif
    }
    return 1;
}

int init_grid(){
    int icell, logspace;
    double rL, rR, L, dr_const, rp, rm, orth_extent;
    getintegerhydropar("geometry", &geometry);
    getintegerhydropar("logspace_cells", &logspace);
    getrealhydropar("rL", &rL);
    getrealhydropar("rR", &rR);
    getrealhydropar("orth_extent", &orth_extent);
    
    getintegerhydropar("isothermal", &isothermal);
    if(isothermal == 1){
        getrealhydropar("cs_init", &cs_init);
    }

    if(logspace == 0){
        L = rR -rL;
        
        dr_const = L/(NCELLS - 2*NGHOST);
        // by default, assume cubical cells
        if(orth_extent < 0){
            orth_extent = dr_const;
        }

        for(icell = 0; icell < NCELLS; icell++){
            rs[icell]  = rL + dr_const*(icell - NGHOST +1./2.);
            rp = rs[icell] + 0.5*dr_const;
            rm = rs[icell] - 0.5*dr_const;
            
            dr[icell]  = dr_const;
            if(geometry == 1){
                vol[icell] = 4.0*M_PI*(pow(rp,3.0)-pow(rm,3.0))/3.0;
            } else {
                vol[icell] = dr_const*orth_extent*orth_extent;
            }
        }
    } else {
        // logspace. Distribute cell centers, with the caveat that the max and
        //  min cells will have their edges outside the bounds

        if(rL <= 0){
            printf("Left edge must be > 0 for logspace distributed cells \n");
            return -1;
        }
        dr_const = exp((log(rR) - log(rL))/(NCELLS - 2 * NGHOST));
        // distribute cell centers
        for(icell = NGHOST; icell < NCELLS; icell++){
            rs[icell] = rL*pow(dr_const, (icell - NGHOST));
        }
        // set lower ghost cells, reflect around 0
        for(icell = 0; icell < NGHOST; icell ++){
                rs[icell] = - rs[2*NGHOST - icell - 1]; 
        }

        // Calculate cell sizes and volumes
        // If we are using 1D cartesian, then we need to make sure that
        // the surface between cells are identical
        
        // by default, assume cells have the size in the orhogonal direction
        // same as for the first cell
        if(orth_extent < 0){
            orth_extent = dr[NGHOST];
        }
        
        for(icell = 1; icell < NCELLS - 1; icell++){
            rp = (rs[icell+1] + rs[icell])*0.5;
            rm = (rs[icell-1] + rs[icell])*0.5;
            dr[icell]  = rp - rm;
            if(geometry == 1){
                vol[icell] = 4.0*M_PI*(pow(rp,3.0)-pow(rm,3.0))/3.0;
            } else {
                vol[icell] = dr[icell]*orth_extent*orth_extent;
            }
        }
        rp = (rs[1] + rs[0])*0.5;
        rm = rs[0] + (rs[0]-rp);
        dr[0] = rp-rm;
        if(geometry == 1){
            vol[icell] = 4.0*M_PI*(pow(rp,3.0)-pow(rm,3.0))/3.0;
        } else {
            vol[icell] = dr[icell]*orth_extent*orth_extent;
        }

        rm = (rs[NCELLS-1] + rs[NCELLS-2])*0.5;
        rp = rs[NCELLS-1] + (rs[NCELLS-1]-rm);
        dr[NCELLS-1] = rp-rm;
        if(geometry == 1){
            vol[icell] = 4.0*M_PI*(pow(rp,3.0)-pow(rm,3.0))/3.0;
        } else {
            vol[icell] = dr[icell]*orth_extent*orth_extent;
        }


    }
   return 1; 
}
int init_domain(){
    int ierr = -1;
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
    } else if(setup == 3){
        ierr = init_wave();
    } else {
        printf("setup %d not recognised", setup);
        return -1;
    }
    if(ierr < 0) {
        return -1;
    }
    return 1;
}





