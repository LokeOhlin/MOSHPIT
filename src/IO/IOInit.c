#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <IO.h>
#include <cgeneral.h>

int nrealIOPars = 0;
real_list_t *IODPars = NULL;
int nintIOPars = 1;
int_list_t *IOIPars = NULL;

int nstepOut;

// set default parameters
int setIOPars(){
    IOIPars = (int_list_t * )malloc(nintIOPars * sizeof(int_list_t));
    strcpy(IOIPars[0].name, "nstepOut"); IOIPars[0].value = 10;
    return 1;
}

int checkIOPars(char *name, char *value){
    // Method to see if parameter matches any of the defined chemistry parameters
    int ipar;
    // Check reals
    for(ipar = 0; ipar < nrealIOPars; ipar ++){
        if(compStr(name, IODPars[ipar].name, 80) > 0){
            IODPars[ipar].value = atof(value);
            return 1;
        }
    }

    for(ipar = 0; ipar < nintIOPars; ipar ++){
        if(compStr(name, IOIPars[ipar].name, 80) > 0){
            IOIPars[ipar].value = atoi(value);
            return 1;
        }
    }
    return -1;
}

// Method to get value of a double(real) parameter
void getrealiopar(char *name, double *value){
    int ipar;
    for(ipar = 0; ipar < nrealIOPars; ipar ++){
        if(compStr(name, IODPars[ipar].name, 80) > 0){
            *value = IODPars[ipar].value;
        }
    }
}

// Method to get value of an integer parameter
void getintegeriopar(char *name, int *value){
    int ipar;
    for(ipar = 0; ipar < nintIOPars; ipar ++){
        if(compStr(name, IOIPars[ipar].name, 80) > 0){
            *value = IOIPars[ipar].value;
        }
    }
}


// gets global parameters
int initIO(){
    getintegeriopar("nstepOut", &nstepOut);
    return 1;
}

