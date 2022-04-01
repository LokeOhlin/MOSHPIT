#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <hydro.h>
#include <radchem.h>
#include <cgeneral.h>
#include <dust.h>
#include <dustRadiation.h>
#include <constantsAndUnits.h>
#include <rtpars.h>
#include <hdf5.h>
#include <IO.h>
#ifdef useDust

int Dust_initIO(){
    int ipar, rpar, ierr, rank, ibuff;
    hsize_t dims_1d, dims_2d[2], dims_3d[3]; 
    hsize_t start_1d, stride_1d, count_1d;

    //  write module parameters    
    for(ipar = 0; ipar < nintDustPars; ipar ++){
        ierr = my_writeAttribute("/Parameters", dustIPars[ipar].name, &dustIPars[ipar].value, H5T_NATIVE_INT);
        if(ierr < 0){
            return -1;
        }
    }
    
    for(rpar = 0; rpar < nrealDustPars; rpar ++){
        ierr = my_writeAttribute("/Parameters", dustDPars[rpar].name, &dustDPars[rpar].value, H5T_NATIVE_DOUBLE);
        if(ierr < 0){
            return -1;
        }
    }
    
    // write specific headers
    ibuff = NdustBins;
    ierr = my_writeAttribute("/Headers", "NdustBins", &ibuff, H5T_NATIVE_INT);
    if(ierr < 0){
        return -1;
    }
    
    ibuff = Nabins - 2;
    ierr = my_writeAttribute("/Headers", "NdustSizeBins", &ibuff, H5T_NATIVE_INT);
    if(ierr < 0){
        return -1;
    }
    
    if(fSi > 0.0 && fSi < 1.0) { 
        ibuff = isilicone - 2;
    } else {
        if(fSi > 0.0){
            ibuff = 0; 
        } else {
            ibuff = isilicone - 2;
        }
    }
    ierr = my_writeAttribute("/Headers", "isilicone", &ibuff, H5T_NATIVE_INT);
    if(ierr < 0){
        return -1;
    }
     
    // create specific datasets, and fill those that are already known 
    rank = 1;
    dims_1d = Nabins - 2;
    ierr = my_createDataset("agrain" , rank, &dims_1d);
    if(ierr < 0){
        return -1;
    }
    
    dims_1d = Nabins - 1;
    ierr = my_createDataset("agrain_binEdges", rank, &dims_1d);
    if(ierr < 0){
        return -1;
    }
    
    // for dust we dont store all data, so datasets needs to be allocated now,
    // which is then filled later during the dust calculations
    
    // Data specific to each cell
    
    // Dust distribution and growth
    rank = 2;
    dims_2d[0] = NCELLS - 2*NGHOST;
    dims_2d[1] = NdustBins;
    ierr = my_createDataset("number", rank, dims_2d);
    if(ierr < 0){
        return -1;
    }
    
    ierr = my_createDataset("slope", rank, dims_2d);
    if(ierr < 0){
        return -1;
    }
    
    ierr = my_createDataset("dadt", rank, dims_2d);
    if(ierr < 0){
        return -1;
    }
#if defined useDustDynamics || defined  growthUpdateVelocities
    ierr = my_createDataset("dustVelocities", rank, dims_2d);
    if(ierr < 0){
        return -1;
    }
#endif
    
    // optical depth for each cell, for each dust bin, for each radiation bin
    if(dust_useRadiation){
        rank = 3;
        dims_3d[0] = NCELLS - 2*NGHOST;
        dims_3d[1] = NdustBins;
        dims_3d[2] = nphotBins;
        ierr = my_createDataset("opticalDepth", rank, dims_3d);
        if(ierr < 0){
            return -1;
        }
    }

    // agrains
    rank = 1;
    start_1d  = 0;
    stride_1d = 1;
    count_1d  = Nabins - 2;
    // +1 to avoid ghost bins
    ierr = my_writeToDataset("agrain", abin_c + 1, rank, &start_1d, &stride_1d, &count_1d);
    if(ierr < 0){
        return -1;
    }
   
    // bin edges
    rank = 1;
    start_1d  = 0;
    stride_1d = 1;
    count_1d  = Nabins - 1;
    ierr = my_writeToDataset("agrain_binEdges", abin_e + 1, rank, &start_1d, &stride_1d, &count_1d);
    if(ierr < 0){
        return -1;
    }
    
    return 1;
}

int Dust_outputCell_dadt(int icell, double dr){
    // We dont want to recalculate the growth rate so this is saved during the chemistry step
    hsize_t start_2d[2], stride_2d[2], count_2d[2];
    int ierr, rank;

    rank = 2;
    start_2d[0] = icell - NGHOST;
    start_2d[1] = 0;
    
    stride_2d[0] = 1;
    stride_2d[1] = 1;

    count_2d[0] = 1;
    count_2d[1] = Nabins-2;
   
    // + 1 to avoid ghost cells 
    ierr = my_writeToDataset("dadt"  , dadt   + 1, rank, start_2d, stride_2d, count_2d);
    if(ierr < 0){
        printf("ERROR cant write to dataset\n");
        return -1;
    }
    
    // if we have two species, we need to do these separate to avoid ghost cells
    if(fSi > 0.0 && fSi < 1.0){
        ierr = my_writeToDataset("dadt"  , dadt   + isilicone + 1, rank, start_2d, stride_2d, count_2d);
        if(ierr < 0){
            printf("ERROR cant write to dataset\n");
            return -1;
        }
    }
    return 1;
}
int Dust_outputCell(int icell, double dr){
    // start with grain distribution and evaporation rates
    hsize_t start_2d[2], stride_2d[2], count_2d[2];
    hsize_t start_3d[3], stride_3d[3], count_3d[3];
    int ierr, rank;

    rank = 2;
    start_2d[0] = icell - NGHOST;
    start_2d[1] = 0;
    
    stride_2d[0] = 1;
    stride_2d[1] = 1;

    count_2d[0] = 1;
    count_2d[1] = Nabins-2;
   
    // + 1 to avoid ghost cells 
    ierr = my_writeToDataset("number", number + 1, rank, start_2d, stride_2d, count_2d);
    if(ierr < 0){
        exit(0);
        printf("ERROR cant write to dataset\n");
        return -1;
    }
    ierr = my_writeToDataset("slope" , slope  + 1, rank, start_2d, stride_2d, count_2d);
    if(ierr < 0){
        printf("ERROR cant write to dataset\n");
        return -1;
    }
    
#if defined useDustDynamics || defined growthUpdateVelocities
    ierr = my_writeToDataset("dustVelocities", velocity + 1, rank, start_2d, stride_2d, count_2d);
    if(ierr < 0){
        printf("ERROR cant write to dataset\n");
        return -1;
    }
#endif
    // if we have two species, we need to do these separate to avoid ghost cells
    if(fSi > 0.0 && fSi < 1.0){
        start_2d[0] = icell - NGHOST;
        start_2d[1] = Nabins - 2;
        
        stride_2d[0] = 1;
        stride_2d[1] = 1;

        count_2d[0] = 1;
        count_2d[1] = Nabins - 2;
   
        // + 1 to avoid ghost cells 
        ierr = my_writeToDataset("number", number + isilicone + 1, rank, start_2d, stride_2d, count_2d);
        if(ierr < 0){
            printf("ERROR cant write to dataset\n");
            return -1;
        }
        
        ierr = my_writeToDataset("slope" , slope  + isilicone + 1, rank, start_2d, stride_2d, count_2d);
        if(ierr < 0){
            printf("ERROR cant write to dataset\n");
            return -1;
        }
        
        ierr = my_writeToDataset("dadt"  , dadt   + isilicone + 1, rank, start_2d, stride_2d, count_2d);
        if(ierr < 0){
            printf("ERROR cant write to dataset\n");
            return -1;
        }
#if defined useDustDynamics || defined growthUpdateVelocities
        ierr = my_writeToDataset("dustVelocities", velocity + isilicone + 1, rank, start_2d, stride_2d, count_2d);
        if(ierr < 0){
            printf("ERROR cant write to dataset\n");
            return -1;
        }
#endif
    }
    
    if(dust_useRadiation){
        ierr = setDustOpticalDepthArray(dr);

        rank = 3;
        start_3d[0] = icell - NGHOST;
        start_3d[1] = 0;
        start_3d[2] = 0;
        
        stride_3d[0] = 1;
        stride_3d[1] = 1;
        stride_3d[2] = 1;

        count_3d[0] = 1;
        count_3d[1] = NdustBins;
        count_3d[2] = nphotBins;
        ierr = my_writeToDataset("opticalDepth", tauDust, rank, start_3d, stride_3d, count_3d);
        if(ierr < 0){
            printf("ERROR cant write to dataset\n");
            return -1;
        }
    
    }
    return 1;

}


int Dust_output(){
    int icell, ierr;
    double rhoDust;
    //Loop over cells. set internal arrays and save variables to output
    for(icell = NGHOST; icell < NCELLS-NGHOST; icell++){
        ierr = setBinsCell(icell, &rhoDust);
        if(ierr < 0) {
            return -1; 
        }
        ierr = Dust_outputCell(icell, dr[icell]);
        if(ierr < 0) {
            return -1; 
        }
    }    
    return 1;
}

#endif
