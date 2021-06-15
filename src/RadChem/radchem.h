#ifndef radchem
#define radchem
#include <radchem.h>
void setRadiationData(double *radData, double dt);
void cellAbsorption(double *radData, double *specData, double numd, double Temp, double dr, double vol, double dt, double *absData);
int doChemistryStep(double dt, double *dt_chem);

int setChemistryPars();
int initChemistry();
void initRadiation();
int checkChemistryPars(char *name, char *value);
void getrealchemistrypar(char *name, double *value);
void getintegerchemistrypar(char *name, int *value);

void setFromStellarModel();
void setFromFile();
#endif
extern double abar, mf_scale, abundHe, abundC, abundO, abundSi, ch_mH, ch_kb, ch_muC;
extern int useRadiationPressure;
extern int numBinsSubIon;  // all bins below 11.2 eV
extern int numBinsFullIon; // all bins above 15.2 eV
extern int numRadiationBins;
extern int iE112; 
extern int iE136;
extern int iE152;
extern double radEmin;
extern double radEmax;
extern double dissE, ionE, ionEH2, sigmaLW, sion, sigma0, sigmaH0;
extern double *sionH, *sionH2, *EionH, *EionH2, *EbinEdges, *Nphots, *aveEphots;
