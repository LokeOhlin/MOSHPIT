#ifndef IO
#define IO
#include <IO.h>
int makeOutput(double t, char *outputName);

int setIOPars();
int checkIOPars(char *name, char *value);
int initIO();
#endif
extern int nstepOut;
