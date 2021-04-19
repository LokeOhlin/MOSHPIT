#ifndef dust
#define dust
#include "dust.h"


#define dust_nghost 1
// Initialisation methods
int initDust();
int setCellInit(int icell, double dustMass);

// User input paramterers
int setDustPars();
int checkDustPars(char *name, char *value);
void getrealdustpar(char *name, double *value);
void getintdustpar (char *name, double *value);

// Transfering data between cell and dust module
int setBinsCell(int icell, double *Mtot);
int getBinsCell(int icell, double *Mtot);

// Methods to calculate mass/number/slopes given the other two
// assuming piecewise linear distribution
double getSlope (double Nj, double Mj,int iabin);
double getNumber(double Sj, double Mj,int iabin);
double getMass  (double Nj, double Sj,int iabin);

// Method to limit slope do avoid dnda < 0 anywhere
int limitSlope(double *Njnew, double *Sjnew, double Nj, double Sj, double Mj, int iabin);

// change from global index (eg index in cell data array) 
// to what is used internally where we have ghost cells
int globalToLocalIndex(int idx);
int localToGlobalIndex(int ibin);

// Main call to update dust distribution
int dustCell(double *rpars, int *ipars, double dt_step);
#endif
// Scratch arrays
extern double *dadt, *number, *slope, *Mnew, *Nnew, *Snew;
// Values that are only dependent on bin size and are constant
extern double *abin_e, *abin_c;
extern double *NfactM, *NfactA;
extern double *SfactM, *SfactA; 

extern int dust_nbins, Nabins, isilicone;
// mass fraction of silicates
extern double fSi;
// density of silicates
extern double rho_s;
// density of graphites
extern double rho_c;
// Options for dust growth
extern int dadt_mode;
extern double dadt_c, dadt_tscale;
// Timestep limiter
extern double dadt_lim;
