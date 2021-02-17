#ifndef radchem
#define radchem
#include "radchem.h"
void setRadiationData(double *radData, double dt);
void cellAbsorption(double *radData, double *specData, double numd, double Temp, double dr, double vol, double dt, double *absData);
int doChemistryStep(double dt, double *dt_chem);

int setChemistryPars();
int initChemistry();
int checkChemistryPars(char *name, char *value);
void getrealchemistrypar(char *name, double *value);
void getintegerchemistrypar(char *name, int *value);
#endif
extern double abar, mf_scale, abundC, abundO, abundSi, ch_mH, ch_kb, ch_muC;
extern int useRadiationPressure;
