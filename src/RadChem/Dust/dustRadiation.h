#include <dust.h>
#ifndef DUSTRAD
#define DUSTRAD
#include <dustRadiation.h>

extern int nphotBins;
extern double *dust_Ephots;
extern double *dust_ELphots;
extern double *dust_ERphots;
extern double *dust_Nphots;
extern double *dust_dEphots;

extern double *tauDust;

extern int *ida_tabQabs;
extern int *ida_tabQem;
extern int dust_useRadiation;
extern int dust_useSublimation;
extern int NtempBins;

//mininmum sublimation rate, below which its assumed to be zero
extern double dadt_min;

int loadDustRadiationTables();
// get index in absorption and emission tables
int getQabs_ida(double agrain, int graphite);
int getQabs_idf(double freq, int graphite);
int getQemAve_ida(double agrain, int graphite);
int getQemAve_idt(double temp, int graphite);
// Methods to get interpolated values
double getQabs(double agrain, double freq, int graphite, int ida, int idf);
double getQsca(double agrain, double freq, int graphite, int ida, int idf);
double getgsca(double agrain, double freq, int graphite, int ida, int idf);
double getQemAve(double agrain, double temp, int graphite, int ida, int idt);
int initDustRadiation(int Nbins, double *Ebins);
int setRadiationBins(double *radData, double dt, double geoFact);

// method to get the total optical depth in one cell
double getDustOpticalDepth(int iEbin, double dr);

// Method to get the absorption of a dust grain from all radiation bins
double getAbsorption(int ibin, int graphite);
// Same but for all bins within range
double getAbsorption_range(int ibin, int graphite, double Emin, double Emax);
double getAbsorptionNum_range(int ibin, int graphite, double nuMin, double nuMax);
// Method to get energy density of photons of frequency freq (u_\nu)
double getRadEdens(double freq);

int initTemperatureDist(int Nbins);
int getTemperatureDist(int ibin, int iabin, int graphite, double Tmax);
double StableState_temperature(double agrain, double piaa, int ida, int graphite, double continous_absorption, double Tmin, double Tmax);

double get_Em(double agrain, double piaa, int ida, double temp, int graphite);
double getUpwardTransition(int ibin, int graphite, double Ui, double UiL, double UiR, double dUi,     double Uf, double UfL, double UfR, double dUf, double continousBound);
double intrabinUpwardTransition(int ibin, int graphite, double dUi, double Ui, double Uf);


double get_dadt_sublimation(int ibin, int iabin, int graphite);
int setDustOpticalDepthArray(double dr);
#endif
