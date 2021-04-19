#include "dust.h"
#ifndef DUSTRAD
#define DUSTRAD
#include "dustRadiation.h"
#define NdustBins 1

extern int nphotBins;
extern double *Ephots;
extern double *ELphots;
extern double *ERphots;
extern double *Nphots;
extern double *dEphots;

extern double parsec;
extern double clght;
extern double planck;
extern double planckInv;
extern double boltzmann;
extern double electronVolt;

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
// Method to get the absorption of a dust grain from all radiation bins
double getAbsorption(int ibin, int graphite);
// Same but for all bins within range
double getAbsorption_range(int ibin, int graphite, double nuMin, double nuMax);
double getAbsorptionNum_range(int ibin, int graphite, double nuMin, double nuMax);
// Method to get energy density of photons of frequency freq (u_\nu)
double getRadEdens(double freq);

int initTemperatureDist(int Nbins);
int getTemperatureDist(int ibin, int graphite, double Tmax, double *Ts, long double *Ps);
double StableState_temperature(double agrain, double piaa, int ida, int graphite, double continous_absorption, double Tmin, double Tmax);

double get_Em(double agrain, double piaa, int ida, double temp, int graphite);
double getUpwardTransition(int ibin, int graphite, double Ui, double UiL, double UiR, double dUi,     double Uf, double UfL, double UfR, double dUf, double continousBound);
#endif
