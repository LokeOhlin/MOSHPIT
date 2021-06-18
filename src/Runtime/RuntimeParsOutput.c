#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <rtpars.h>
#include <cgeneral.h>
#include <hdf5.h>
#include <IO.h>
int RuntimePars_initIO(){
    int ipar, rpar, ierr;

    //  write module parameters    
    for(ipar = 0; ipar < nintRTPars; ipar ++){
        ierr = my_writeAttribute("/Parameters", RTIPars[ipar].name, &RTIPars[ipar].value, H5T_NATIVE_INT);
        if(ierr < 0){
            return -1;
        }
    }
    
    for(rpar = 0; rpar < nrealRTPars; rpar ++){
        ierr = my_writeAttribute("/Parameters", RTDPars[rpar].name, &RTDPars[rpar].value, H5T_NATIVE_DOUBLE);
        if(ierr < 0){
            return -1;
        }
    }

    // write specific headers
    ierr = my_writeAttribute("/Headers", "time", &time, H5T_NATIVE_DOUBLE); 
    if(ierr < 0){
        return -1;
    }

    ierr = my_writeAttribute("/Headers", "dt", &dt, H5T_NATIVE_DOUBLE);
    if(ierr < 0){
        return -1;
    }
    
    return 1;
}

//Nothing done here
int RuntimePars_output(){
    return 1;
}
