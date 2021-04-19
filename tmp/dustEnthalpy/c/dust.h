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
extern int dust_nbins, Nabins, isilicone;

extern double *agrains;
extern double *Natoms;
extern double *volgrains;
extern double *pi_asquare;
extern int *ida_tabQabs;
extern int *ida_tabQem;


// mass fraction of silicates
extern double fSi;
// density of silicates
extern double rho_s;
// density of graphites
extern double rho_c;
