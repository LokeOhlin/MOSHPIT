#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <hydro.h>
#include <IO.h>
#ifdef useDust 
    #include <dust.h>
#endif
int makeOutput(double t, char *outputName){
    int icell,ivar,idx, ierr;
    FILE *fptr = fopen(outputName, "w");
    //write output info
    fprintf(fptr,"# time = %.4e\n", t);
    fprintf(fptr,"# nvar = %d\n", nvar);
#ifdef useDust
    double tmp;
    fprintf(fptr,"# IDUST_START = %d\n", IDUST_START);
    fprintf(fptr,"# NdustVar = %d\n", NdustVar);
    fprintf(fptr,"# fSi = %.4e\n", fSi);
    fprintf(fptr,"# Nabins = %d\n", Nabins-2*dust_nghost);
    if(fSi > 0.0 && fSi < 1.0){
        fprintf(fptr,"# isilicone = %d\n", isilicone-2*dust_nghost);
    } else {
        if(fSi > 0.0){
            fprintf(fptr,"# isilicone = %d\n", isilicone);
        } else {
            fprintf(fptr,"# isilicone = %d\n", isilicone- 2*dust_nghost);
        }
    }
    getrealdustpar("dust_amin", &tmp);
    fprintf(fptr,"# amin = %.4e\n", tmp); 
    getrealdustpar("dust_amax", &tmp);
    fprintf(fptr,"# amax = %.4e\n", tmp);
#endif

    ierr = toPrimitive();
    if(ierr < 0) {
        return -1;
    }
    for(icell = NGHOST; icell < NCELLS-NGHOST; icell++){
        fprintf(fptr, "%.8e\t", rs[icell]);
        idx = icell*nvar;
        for(ivar = 0; ivar < nvar; ivar++){
            fprintf(fptr, "%.8e\t", pstate[idx+ivar]);
        }
        fprintf(fptr,"\n");
    }
    fclose(fptr);
    return 1;
}
