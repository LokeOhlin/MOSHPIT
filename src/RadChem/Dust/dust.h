#include <cgeneral.h>
#include <hdf5.h>
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

// methods used in sputtering 
// loading tables
int loadDustSputteringTables();
// getting indexes
int getSputYield_ida(double agrain, int graphite);
int getSputYield_idt(double tgas, int graphite);
// getting and calcualtint yields and sputtering rates
double getSputYield_t(double tgas, int graphite, int ida, int *idt);
double getSputYield_vel_t(double tgas, int graphite, int ida, int idv, int *idt);
double getSputYield(double agrain, double tgas, int graphite, int *ida, int *idt);
double getSputYield_vel(double agrain, double vdust, double tgas, int graphite, int *ida, int *idt);
double get_dadt_sputtering(double agrain, int ibin, int iabin, int graphite, double numd, double tgas, int *idt);
double get_dadt_sputtering_vel(double agrain, int ibin, int iabin, int graphite, double numd, double vdust, double tgas, int *idt);

// general method to calculate growth for each dustbin
// independent of substep
int set_dadt_fixed(double *rpars);
// dependent on substep
int set_dadt(double *rpars, double dt);


int Dust_initIO();
int Dust_outputCell(int icell, double dr);
int Dust_output();

// Methods for the dust velocity
double get_average_vdust(int ibin, int iabin, int graphite, double *rpars, double dt);
int update_vdust(int ibin, int iabin, int graphite, double *rpars, double dt);
int set_vdusts(int icell);
int get_vdust(int icell);


#endif
// Scratch arrays
extern double *dadt, *dadt_fixed, *number, *slope, *Mnew, *Nnew, *Snew, *vnew, *dust_vrel;
// Values that are only dependent on bin size and are constant
extern double *abin_e, *abin_c;
extern double *NfactM, *NfactA;
extern double *SfactM, *SfactA; 
extern double *pi_asquare;
extern double *volgrains;

extern int dust_nbins, Nabins, isilicone;
// mass fraction of silicates
extern double fSi;
// density and average atom mass of silicates 
extern double rho_s;
extern double aveMatom_s; //(~20 mH)
// density and average atom mass of graphites
extern double rho_c;
extern double aveMatom_c; //(~12 mH)
// Number of atoms in a grain
extern double *Natoms;

// Options for dust growth
extern int dadt_mode;
extern double dadt_c, dadt_tscale;
// Timestep limiter
extern double dadt_lim;

// flag for use of sputtering
extern int dust_useSputtering;
// index of dustgrain in sputtering table
extern int *ida_tabSput;

extern int dust_maxSubSteps;

// options for behaviour at upper and lower bound
extern int dust_lowerBound_pileUp;
extern int dust_upperBound_pileUp;

extern hid_t dustOutput;
extern int outputDust;
extern int outputNum;

extern int nrealDustPars, nintDustPars;
extern real_list_t *dustDPars;
extern int nintDustPars;
extern int_list_t *dustIPars;
