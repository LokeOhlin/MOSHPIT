#include <cgeneral.h>
#ifndef HYDRO
#define HYDRO
#include <hydro.h>

// Number of ghosts per side
#define NGHOST 4

// How large does our arrays have to be?
// eg how many variables and species
#define IHYDRO_START 0
#define IHYDRO_END 3
#define IADVECT_START 3
#ifdef useChemistry
    #define ICHEM_START  3
    #define ICHEM_END  8
    // No chemistry then no dust for now
    #ifdef useDust
        #ifndef NdustBins // if not user specified default 20
            # define NdustBins 20
        #endif
    
        #define IDUST_START (ICHEM_END+1)

        #if defined growthUpdateVelocities || defined useDustDynamics
            #define NdustVar 3
        #else
            #define NdustVar 2
        #endif
        // (5 + 2*NdustBins +1) new variables : H, H2, Hp, CO, Cp, dust mass and slope, Tdust from chemistry (not related to dust module)  
        #define nvar    (IDUST_START + NdustVar * NdustBins)  
    #else
        // 6 new variables: H, H2, Hp, CO, Cp, Tdust
        #define nvar   9
        #define IADVECT_END ICHEM_END
    #endif
    #if defined useDust 
        #if defined useDustDynamics
            #define IADVECT_END IDUST_START
        #else
            #define IADVECT_END nvar
        #endif
    #endif
#else
    #define ICHEM_START  3
    #define ICHEM_END 3
    #define nvar 3
    #define IADVECT_END   3
#endif
#define nFluxVar (nvar + 2)
int toPrimitive();
int doHydroStep(double dt);
int getCFL(double *cfl);

int setHydroPars();
int checkHydroPars(char *name, char *value);
int initHydro();
int init_domain();
double vanLeer(double a, double b);

extern int useHydro;
// ISOTHEMRMAL
extern int isothermal;
extern double cs_init;
// GRID GEOMETRY
extern int NCELLS; 
extern int NINTER;
extern int geometry;
// INTERPOLATION
extern int interpolator;

// BOUNDARY CONDITIONS
extern int left_bound;
extern int right_bound;
extern double bdensL, bvelL, benerL;
extern double bdensR, bvelR, benerR;

// HYDRO CONSTANTS
extern double adi;
extern double courant_number;
extern double hy_ethresh;
extern double roe_p;

//Pointers to state arrays
extern double *ustate;
extern double *pstate;
extern double *rs;
extern double *right_edge;
extern double *vol;
extern double *dr;

extern double *varBuff;
int Hydro_initIO();
int Hydro_output();

extern int nrealHydroPars, nintHydroPars;
extern int_list_t *hydroIPars;
extern real_list_t *hydroDPars;
#endif

