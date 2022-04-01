#include <cgeneral.h>
#ifndef RTPARS
#define RTPARS
#include <rtpars.h>
int checkRuntimePars(char *name, char *value);
int initRuntimePars();
int setRuntimePars();
int RuntimePars_initIO();
int RuntimePars_output();
#endif
extern double t0;
extern double tend;
extern double dt_init, dt_max;
extern double time, dt, timeOut;
extern int imax;
extern int nintRTPars, nrealRTPars;
extern int_list_t *RTIPars;
extern real_list_t *RTDPars;
