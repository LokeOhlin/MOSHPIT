#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "hydro.h"
#include "IO.h"
#include "cgeneral.h"
#include "init.h"
#include "rtpars.h"
#ifdef useChemistry
#include "radchem.h"
#ifdef useDust
#include "dust.h"
#endif
#endif


int readParameterFile(){
    int ichar, endchar;
    int found;
    int containsPar;
    char line[128];
    char *name;
    char *value;
    FILE *fptr = fopen("simulation.par", "r");
    
    while(fgets(line, 128, fptr) !=NULL){
        containsPar = 0;
        for(ichar = 0; ichar < 128; ichar++){
            if(line[ichar]=='\n' || line[ichar]=='#'){
                endchar = ichar;
                break;
            }
            if(line[ichar]!= ' '){
                containsPar = 1;
            } 
        }

        if (containsPar <= 0){
            continue;
        }
        // endchar is either end of file or a comment. ether case: break line
        line[endchar] = '\0';
        sscanf(line, "%m[^=]=%ms", &name, &value); 
        if(value == NULL){
            printf("\nERROR : %s HAS NO VALUE ASSIGNED\n\n",name);
            continue;
        }
        trim(name);
        trim(value);
        // Checks modules for matching parameter name.
        // Continues when found.
        found = checkHydroPars(name, value);
        if(found > 0){
            continue;
        }
        found = checkIOPars(name, value);
        if(found > 0){
            continue;
        }
        found = checkChemistryPars(name, value);
        if(found > 0){
            continue;
        }
        found = checkDustPars(name, value);
        if(found > 0){
            continue;
        }
        found = checkRuntimePars(name, value);
        if(found > 0){
            continue;
        }
        printf("Parameter not found : %s \n",name);
    }
    return 1;
}
// Wrapper to init all other modules & initialise the domain
int initSimulation(){
    int ierr;
    
    printf("\tSetting all default parameters\n");
    ierr = setHydroPars();
    ierr = setRuntimePars();
    ierr = setIOPars();
    ierr = setChemistryPars();
    ierr = setDustPars();

    printf("\treading parameter file\n");
    ierr = readParameterFile();
    if(ierr < 0){
        return ierr;
    }
    printf("\tcalling initHydro\n");
    ierr = initHydro();
    if(ierr < 0){
        return ierr;
    }
    printf("\tcalling initRuntimePars\n");
    ierr = initRuntimePars();
    if(ierr < 0){
        return ierr;
    }
    printf("\tcalling initIO\n");
    ierr = initIO();
    if(ierr < 0){
        return ierr;
    }
    printf("\tcalling initChemistry\n");
    ierr = initChemistry();
    if(ierr < 0){
        return ierr;
    }
    return 1;
}



