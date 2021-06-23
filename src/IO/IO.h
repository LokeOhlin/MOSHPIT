#include <hdf5.h>
#ifndef IO
#define IO
#include <IO.h>
int makeOutput(double t, char *outputName);

int setIOPars();
int checkIOPars(char *name, char *value);
int initIO();

int createOutputFile(int output_id);
int finalizeOutputFile();

int my_createDataset(char *datasetName, int rank, hsize_t *dims);
int my_writeToDataset(char *datasetName, double *varArr, int rank, hsize_t *start, hsize_t *stride, hsize_t *count);

int my_writeAttribute(const char *groupName, const char *name, const void *buff, hid_t type);
#endif
extern int nstepOut;
extern double dtOut;
extern hid_t h5Output;
