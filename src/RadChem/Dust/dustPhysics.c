#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "dust.h"
#include "hydro.h"
#include "radchem.h"
#include "moshpit.h"


int setTdust_bin(double *rpars, int *ipars, int ibin, int iabin, int graphite){
    tdust[ibin] = 100;
    return 1;
}

int setTdust(double *rpars, int *ipars){
    int ibin, iabin, graphite, ierr;

    for(ibin = 0; ibin < dust_nbins; ibin++){
    
        if(ibin < isilicone){
            graphite = 1;
        } else {
            graphite = 0;
        }
        iabin = ibin - isilicone*(1-graphite);
        ierr  = setTdust_bin(rpars, ipars, ibin, iabin, graphite);
        if(ierr < 0){
            return -1;
        }
    }
    return 1;
}





