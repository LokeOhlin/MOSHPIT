#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <hydro.h>
#include <cgeneral.h>
#include <hdf5.h>
#include <IO.h>
#include "constantsAndUnits.h"
#ifdef useChemistry
#include <radchem.h>
#endif 

int Hydro_initIO(){
    int ipar, rpar, rank, ierr;
    hsize_t dims_1d; 
    
    //  write module parameters    
    for(ipar = 0; ipar < nintHydroPars; ipar ++){
        ierr = my_writeAttribute("/Parameters", hydroIPars[ipar].name, &hydroIPars[ipar].value, H5T_NATIVE_INT);
        if(ierr < 0){
            return -1;
        }
    }
    
    for(rpar = 0; rpar < nrealHydroPars; rpar ++){
        ierr = my_writeAttribute("/Parameters", hydroDPars[rpar].name, &hydroDPars[rpar].value, H5T_NATIVE_DOUBLE);
        if(ierr < 0){
            return -1;
        }
    }

    rank = 1;
    dims_1d = NCELLS - 2*NGHOST;
    // create datasets
    
    // cell spatial data
    ierr = my_createDataset("coordinates", rank, &dims_1d);
    if(ierr < 0){
        return -1;
    }

    ierr = my_createDataset("deltax", rank, &dims_1d);
    if(ierr < 0){
        return -1;
    }
    
    ierr = my_createDataset("volume", rank, &dims_1d);
    if(ierr < 0){
        return -1;
    }
  
    // hydro data
    ierr = my_createDataset("density"       , rank, &dims_1d);
    if(ierr < 0){
        return -1;
    }

    ierr = my_createDataset("velocity"      , rank, &dims_1d);
    if(ierr < 0){
        return -1;
    }
    
    ierr = my_createDataset("pressure"      , rank, &dims_1d);
    if(ierr < 0){
        return -1;
    }
    
    ierr = my_createDataset("temperature"   , rank, &dims_1d);
    if(ierr < 0){
        return -1;
    }
    return 1;
}


int Hydro_output(){
    int icell, idx, ivar, rank, ierr;
    double dens, numd, velx, xH2, xHp;
    hsize_t start_1d, stride_1d, count_1d;
    
    rank = 1;
    start_1d = 0;
    stride_1d = 1;
    count_1d = NCELLS - 2*NGHOST;
    
    // cell spatial data

    // coordinate
    ierr = my_writeToDataset("coordinates", rs+NGHOST, rank, &start_1d, &stride_1d, &count_1d);
    if(ierr < 0){
        return -1;
    }
    
    // cell size
    ierr = my_writeToDataset("deltax", dr+NGHOST, rank, &start_1d, &stride_1d, &count_1d);
    if(ierr < 0){
        return -1;
    }

    // cell volume
    ierr = my_writeToDataset("volume", vol+NGHOST, rank, &start_1d, &stride_1d, &count_1d);
    if(ierr < 0){
        return -1;
    }
    


    // Density
    ivar = 0;
    for(icell = NGHOST; icell < NCELLS - NGHOST; icell++){
        idx = icell*nvar;
        varBuff[icell - NGHOST] = ustate[idx + ivar];
    }
    ierr = my_writeToDataset("density", varBuff, rank, &start_1d, &stride_1d, &count_1d);
    if(ierr < 0){
        return -1;
    }
    
    // velocity
    ivar = 1;
    for(icell = NGHOST; icell < NCELLS - NGHOST; icell++){
        idx = icell*nvar;
        dens = ustate[idx];
        varBuff[icell - NGHOST] = ustate[idx + ivar]/dens;
    }
    ierr = my_writeToDataset("velocity", varBuff, rank, &start_1d, &stride_1d, &count_1d);
    if(ierr < 0){
        return -1;
    }

    // pressure
    ivar = 2;
    for(icell = NGHOST; icell < NCELLS - NGHOST; icell++){
        idx = icell*nvar;
        dens = ustate[idx];
        velx = ustate[idx + 1]/dens;
        
        varBuff[icell - NGHOST] = (ustate[idx + ivar] - 0.5*dens*velx*velx) * (adi - 1);
    }

    ierr = my_writeToDataset("pressure", varBuff, rank, &start_1d, &stride_1d, &count_1d);
    if(ierr < 0){
        return -1;
    }

    // temperature 
    // we can reuse calculation in varBuff
    for(icell = NGHOST; icell < NCELLS - NGHOST; icell++){
        idx = icell*nvar;
        dens = ustate[idx];

#ifdef useChemistry
        // Get abundances
        xH2 = (ustate[idx+ICHEM_START+1]/dens) * mf_scale/2.0;
        xHp = (ustate[idx+ICHEM_START+2]/dens) * mf_scale;
 
        numd = dens*(1-xH2+xHp+abundHe)/mH/abar;
        
#else
        // just assume atomic hydogen + helium
        numd = dens*1.1/(mH*1.4);
#endif 
        varBuff[icell - NGHOST] = varBuff[icell - NGHOST]/(numd*boltzmann);
    }
    ierr = my_writeToDataset("temperature", varBuff, rank, &start_1d, &stride_1d, &count_1d);
    if(ierr < 0){
        return -1;
    }

    return 1;
}
