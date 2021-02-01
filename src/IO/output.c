#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "hydro.h"
#include "IO.h"

int makeOutput(double t, char *outputName){
    int icell,ivar,idx, ierr;
    FILE *fptr = fopen(outputName, "w");

    //write output info
    fprintf(fptr,"# time = %.4e\n", t);

    ierr = toPrimitive();
    for(icell = 2; icell < NCELLS-2; icell++){
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
