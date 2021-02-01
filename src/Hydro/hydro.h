#ifndef HYDRO
#define HYDRO
#include "hydro.h"
#ifdef useChemistry
    #define ICHEM_START  3
    #define ICHEM_END  8
    #define nvar   9
// 6 new variables: H, H2, Hp, CO, Cp, Tdus
// // 6 new variables: H, H2, Hp, CO, Cp, Tdust
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
