#ifndef HYDRO
#define HYDRO
#include "hydro.h"

// How large does our arrays have to be?
// eg how many variables and species
#ifdef useChemistry
    #define ICHEM_START  3
    #define ICHEM_END  8
    // No chemistry then no dust for now
    #ifdef useDust
        #ifndef NdustBins // if not user specified default 20
            # define NdustBins 20
        #endif
    
        #define IDUST_START (ICHEM_END+1)
        #define NdustVar 2
        // (5 + 2*NdustBins +1) new variables : H, H2, Hp, CO, Cp, dust mass and slope, Tdust from chemistry (not related to dust module)  
        #define nvar    (IDUST_START + NdustVar * NdustBins)  
    #else
        // 6 new variables: H, H2, Hp, CO, Cp, Tdust
        #define nvar   9
    #endif
#else
    #define ICHEM_END 3
    #define nvar 3
#endif
int toPrimitive();
int doHydroStep(double dt);
int getCFL(double *cfl);

int setHydroPars();
int checkHydroPars(char *name, char *value);
int initHydro();
int init_domain();

extern int useHydro;

// GRID GEOMETRY
extern int NGHOST; 
extern int NCELLS; 
extern int NINTER;
extern int geometry;

// BOUNDARY CONDITIONS
extern int left_bound;
extern int right_bound;
extern double bdensL, bvelL, benerL;
extern double bdensR, bvelR, benerR;

// HYDRO CONSTANTS
extern double adi;
extern double courant_number;

//Pointers to state arrays
extern double *ustate;
extern double *pstate;
extern double *rs;
extern double *vol;
extern double *dr;

#endif
