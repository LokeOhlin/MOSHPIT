#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <rtpars.h>
#include <cgeneral.h>

int nrealRTPars = 4;
real_list_t *RTDPars = NULL;
int nintRTPars = 1;
int_list_t *RTIPars = NULL;

double t0, tend, dt_init, dt_max;
int imax;

// set default parameters
int setRuntimePars(){
    RTDPars = (real_list_t *) malloc(nrealRTPars * sizeof(real_list_t));
    RTIPars = (int_list_t *) malloc(nintRTPars * sizeof(int_list_t));
    strcpy(RTDPars[0].name, "t0"); RTDPars[0].value = 0.0;
    strcpy(RTDPars[1].name, "tend"); RTDPars[1].value = 1e3;
    strcpy(RTDPars[2].name, "dt_init"); RTDPars[2].value = 1.0;
    strcpy(RTDPars[3].name, "dt_max"); RTDPars[2].value = 1e10;
    
    strcpy(RTIPars[0].name, "imax"); RTIPars[0].value = 1000;
    return 1;
}

int checkRuntimePars(char *name, char *value){
    // Method to see if parameter matches any of the defined chemistry parameters
    int ipar;
    // Check reals
    for(ipar = 0; ipar < nrealRTPars; ipar ++){
        if(compStr(name, RTDPars[ipar].name, 80) > 0){
            RTDPars[ipar].value = atof(value);
            return 1;
        }
    }

    for(ipar = 0; ipar < nintRTPars; ipar ++){
        if(compStr(name, RTIPars[ipar].name, 80) > 0){
            RTIPars[ipar].value = atoi(value);
            return 1;
        }
    }
    return -1;
}

// Method to get value of a double(real) parameter
void getrealrtpar(char *name, double *value){
    int ipar;
    for(ipar = 0; ipar < nrealRTPars; ipar ++){
        if(compStr(name, RTDPars[ipar].name, 80) > 0){
            *value = RTDPars[ipar].value;
        }
    }
}

// Method to get value of an integer parameter
void getintegerrtpar(char *name, int *value){
    int ipar;
    for(ipar = 0; ipar < nintRTPars; ipar ++){
        if(compStr(name, RTIPars[ipar].name, 80) > 0){
            *value = RTIPars[ipar].value;
        }
    }
}


// gets global parameters
int initRuntimePars(){
    getrealrtpar("t0", &t0);
    getrealrtpar("tend", &tend);
    getrealrtpar("dt_init", &dt_init);
    getrealrtpar("dt_max", &dt_max);
    getintegerrtpar("imax", &imax);
    return 1;
}

