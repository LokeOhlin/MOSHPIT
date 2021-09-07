#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <hydro.h>
#include <radchem.h>
#include <cgeneral.h>
#include <hdf5.h>
#include <IO.h>
#ifdef useChemistry
int RadChem_initIO(){
    int ipar, rpar, ierr, rank;
    hsize_t dims_1d, start_1d, stride_1d, count_1d;
    hsize_t dims_2d[2];
    //  write module parameters    
    for(ipar = 0; ipar < nintChemPars; ipar ++){
        ierr = my_writeAttribute("/Parameters", chemIPars[ipar].name, &chemIPars[ipar].value, H5T_NATIVE_INT);
        if(ierr < 0){
            return -1;
        }
    }
    
    for(rpar = 0; rpar < nrealChemPars; rpar ++){
        ierr = my_writeAttribute("/Parameters", chemDPars[rpar].name, &chemDPars[rpar].value, H5T_NATIVE_DOUBLE);
        if(ierr < 0){
            return -1;
        }
    }

    // write specific headers
    ierr = my_writeAttribute("/Headers", "NphotBins", &numRadiationBins, H5T_NATIVE_INT);
    if(ierr < 0){
        return -1;
    }
     
    // create specific datasets, and fill those that are already known 
    rank = 1;
    dims_1d = numRadiationBins;
    ierr = my_createDataset("photonBins", rank, &dims_1d);
    if(ierr < 0){
        return -1;
    }
    
    rank = 2;
    dims_2d[0] = NCELLS - 2*NGHOST;
    dims_2d[1] = 5;

    ierr = my_createDataset("chemicalAbundances", rank, dims_2d);
    if(ierr < 0){
        return -1;
    }

#ifdef savePhotonFluxes
    rank = 2;
    dims_2d[0] = NCELLS - 2*NGHOST;
    dims_2d[1] = numRadiationBins;
    ierr = my_createDataset("photonFluxes", rank, dims_2d);
#endif
    
    //radiation bins
    rank = 1;
    start_1d  = 0;
    stride_1d = 1;
    count_1d  = numRadiationBins;
    ierr = my_writeToDataset("photonBins", EbinEdges, rank, &start_1d, &stride_1d, &count_1d);
    if(ierr < 0){
        return -1;
    }
    
    return 1;
}

//write chemical species at end of step
int RadChem_output(){
    int icell, idx, rank, ierr;
    double dens, xHI, xH2, xHp, xCO, xCp;
    hsize_t start_2d[2], stride_2d[2], count_2d[2];

    for(icell = NGHOST; icell < NCELLS - NGHOST; icell++){
        idx = icell * nvar;
        
        dens = ustate[idx];

        xHI = ustate[idx+ICHEM_START  ]/dens * mf_scale;
        xH2 = ustate[idx+ICHEM_START+1]/dens * mf_scale/2.0;
        xHp = ustate[idx+ICHEM_START+2]/dens * mf_scale;

        xCO = ustate[idx+ICHEM_START+3]/dens * mf_scale/ch_muC;
        xCp = ustate[idx+ICHEM_START+4]/dens * mf_scale/ch_muC;
        
        chemBuff[(icell - NGHOST) * 5    ] = xHI;
        chemBuff[(icell - NGHOST) * 5 + 1] = xH2;
        chemBuff[(icell - NGHOST) * 5 + 2] = xHp;
        chemBuff[(icell - NGHOST) * 5 + 3] = xCO;
        chemBuff[(icell - NGHOST) * 5 + 4] = xCp;
    }
    
    rank = 2;

    start_2d[0] = 0;
    start_2d[1] = 0;
    
    stride_2d[0] = 1;
    stride_2d[1] = 1;

    count_2d[0] = NCELLS - 2*NGHOST ;
    count_2d[1] = 5;
    
    ierr = my_writeToDataset("chemicalAbundances", chemBuff, rank, start_2d, stride_2d, count_2d);
    if(ierr < 0){
        return -1;
    }

#ifdef savePhotonFluxes
    rank = 2;

    start_2d[0] = 0;
    start_2d[1] = 0;
    
    stride_2d[0] = 1;
    stride_2d[1] = 1;

    count_2d[0] = NCELLS - 2*NGHOST;
    count_2d[1] = numRadiationBins;
    ierr = my_writeToDataset("photonFluxes", photonFluxes, rank, start_2d, stride_2d, count_2d);
    return 1;
#endif
}
#endif
