#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <mechanicalFeedback.h>
#include <cgeneral.h>

int nrealFBPars = 3;
real_list_t *FBDPars = NULL;
int nintFBPars = 1;
int_list_t *FBIPars = NULL;


int fb_useWind;
double fb_MdotWind, fb_vWind, fb_dEmaxWind;
// set default parameters
int setFeedbackPars(){
    FBDPars = (real_list_t * )malloc(nrealFBPars * sizeof(real_list_t));
    strcpy(FBDPars[0].name, "fb_MdotWind"); FBDPars[0].value = 5.89876e+15;
    strcpy(FBDPars[1].name, "fb_vWind"); FBDPars[1].value = 2.40608e+08;
    strcpy(FBDPars[2].name, "fb_dEmaxWind"); FBDPars[2].value = 0.01;

    FBIPars = (int_list_t * )malloc(nintFBPars * sizeof(int_list_t));
    strcpy(FBIPars[0].name, "fb_useWind"); FBIPars[0].value = 0;
    return 1;
}

int checkFeedbackPars(char *name, char *value){
    // Method to see if parameter matches any of the defined chemistry parameters
    int ipar;
    // Check reals
    for(ipar = 0; ipar < nrealFBPars; ipar ++){
        if(compStr(name, FBDPars[ipar].name, 80) > 0){
            FBDPars[ipar].value = atof(value);
            return 1;
        }
    }

    for(ipar = 0; ipar < nintFBPars; ipar ++){
        if(compStr(name, FBIPars[ipar].name, 80) > 0){
            FBIPars[ipar].value = atoi(value);
            return 1;
        }
    }
    return -1;
}

// Method to get value of a double(real) parameter
void getrealfbpar(char *name, double *value){
    int ipar;
    for(ipar = 0; ipar < nrealFBPars; ipar ++){
        if(compStr(name, FBDPars[ipar].name, 80) > 0){
            *value = FBDPars[ipar].value;
        }
    }
}

// Method to get value of an integer parameter
void getintegerfbpar(char *name, int *value){
    int ipar;
    for(ipar = 0; ipar < nintFBPars; ipar ++){
        if(compStr(name, FBIPars[ipar].name, 80) > 0){
            *value = FBIPars[ipar].value;
        }
    }
}


// gets global parameters
int initFeedback(){
    getintegerfbpar("fb_useWind", &fb_useWind);

    getrealfbpar("fb_MdotWind", &fb_MdotWind);
    getrealfbpar("fb_vWind", &fb_vWind);
    getrealfbpar("fb_dEmaxWind", &fb_dEmaxWind);
    return 1;
}

