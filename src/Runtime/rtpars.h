#ifndef RTPARS
#define RTPARS
#include "rtpars.h"
int checkRuntimePars(char *name, char *value);
int initRuntimePars();
int setRuntimePars();
#endif
extern double t0;
extern double tend;
extern double dt_init, dt_max;
extern int imax;
