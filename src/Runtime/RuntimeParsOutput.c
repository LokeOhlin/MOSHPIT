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

    
    return 1;
}

int RuntimePars_output(){
    int ierr;
    // write specific headers
    // This time has not yet been advanced. Add dt to make in sync
    ierr = my_writeAttribute("/Headers", "time", &timeOut, H5T_NATIVE_DOUBLE); 
    if(ierr < 0){
        return -1;
    }

    ierr = my_writeAttribute("/Headers", "dt", &dt, H5T_NATIVE_DOUBLE);
    if(ierr < 0){
        return -1;
    }
    return 1;
}
