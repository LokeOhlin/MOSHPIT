#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <mechanicalFeedback.h>
#include <cgeneral.h>
#include <hdf5.h>
#include <IO.h>
int Feedback_initIO(){
    int ipar, rpar, ierr;

    //  write module parameters    
    for(ipar = 0; ipar < nintFBPars; ipar ++){
        ierr = my_writeAttribute("/Parameters", FBIPars[ipar].name, &FBIPars[ipar].value, H5T_NATIVE_INT);
        if(ierr < 0){
            return -1;
        }
    }
    
    for(rpar = 0; rpar < nrealFBPars; rpar ++){
        ierr = my_writeAttribute("/Parameters", FBDPars[rpar].name, &FBDPars[rpar].value, H5T_NATIVE_DOUBLE);
        if(ierr < 0){
            return -1;
        }
    }

    return 1;
}

//Nothing is done here
int Feedback_output(){
    return 1;
}
