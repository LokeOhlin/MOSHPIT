#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <hydro.h>
#include <radchem.h>
#include <cgeneral.h>
#include <dust.h>
#include <dustRadiation.h>
#include <constantsAndUnits.h>
#include <rtpars.h>
#include <hdf5.h>
#ifdef useDust
hid_t dustOutput;

int my_createDataset(char *datasetName, int rank, hsize_t *dims){
    hid_t dataset, datatype, dataspace;
    dataspace = H5Screate_simple(rank, dims, NULL);
    datatype  = H5Tcopy(H5T_NATIVE_DOUBLE);
    
    dataset   = H5Dcreate(dustOutput, datasetName, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    return 1;
}

int my_writeToDataset(char *datasetName, double *varArr, int rank, hsize_t *start, hsize_t *stride, hsize_t *count){
    hid_t dataset, dataspace, memspace;
    herr_t hstat;
    
    // open dataspace and dataset
    dataset = H5Dopen(dustOutput, datasetName, H5P_DEFAULT);
    dataspace = H5Dget_space(dataset);


    // open memory space and select slab
    memspace = H5Screate_simple(rank, count, NULL);
    hstat = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, stride, count, NULL);
    if(hstat < 0) {
        printf("Error in selecting hyperslab from %s \n", datasetName);
        return -1;
    } 
    
    hstat = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, varArr);
    if(hstat < 0) {
        printf("Error in writing data to %s \n", datasetName);
        return -1;
    } 
       
    // close all pointers
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Sclose(memspace);
    return 1;
}

int my_writeAttribute(hid_t handle, const char *name, const void *buff, hid_t type){
    hid_t dataspace, attribute;
    herr_t status;
    dataspace = H5Screate(H5S_SCALAR);
    attribute = H5Acreate(handle, name, type, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    status    = H5Awrite(attribute, type, buff);
    H5Aclose(attribute);
    H5Sclose(dataspace);   
    return 1;
}


int initDustOutput(double dt){
    char outname[13] = "";
    char fname[9] = "dustOut_";
    char fnum[5] = "";
    int ierr, rank;
    int ibuff;
    double dbuff;
    hid_t headerGroup;
    hsize_t dims_1d, dims_2d[2], dims_3d[3]; 
    hsize_t start_1d, stride_1d, count_1d;

    outputNum++;
    sprintf(fnum, "%04d", outputNum);
    strcat(outname, fname);
    strcat(outname, fnum);
    dustOutput = H5Fcreate(outname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    
    // create header data
    headerGroup = H5Gcreate(dustOutput, "/Header", 0, H5P_DEFAULT, H5P_DEFAULT);
    // store header attributes
    // data is saved at end of timestep
    dbuff = time + dt;
    ierr = my_writeAttribute(headerGroup, "time", &dbuff, H5T_NATIVE_DOUBLE);
    ierr = my_writeAttribute(headerGroup, "dt", &dt, H5T_NATIVE_DOUBLE);
    
    ibuff = NdustBins;
    ierr = my_writeAttribute(headerGroup, "NdustBins", &ibuff, H5T_NATIVE_INT);
    ibuff = Nabins - 1;
    ierr = my_writeAttribute(headerGroup, "NdustSizeBins", &ibuff, H5T_NATIVE_INT);
    ierr = my_writeAttribute(headerGroup, "NphotBins", &nphotBins, H5T_NATIVE_INT);
    
    
    ierr = my_writeAttribute(headerGroup, "fSilicone", &fSi, H5T_NATIVE_DOUBLE);
    if(fSi > 0.0 && fSi < 1.0) { 
        ibuff = isilicone - 2;
    } else {
        if(fSi > 0.0){
            ibuff = 0; 
        } else {
            ibuff = isilicone - 2;
        }
    }
    ierr = my_writeAttribute(headerGroup, "isilicone", &ibuff, H5T_NATIVE_INT);
     
    H5Gclose(headerGroup);
    
    
    //create all needed datasets
    // cell data
    rank = 1;
    dims_1d = NCELLS - 2*NGHOST;
    ierr = my_createDataset("radius" , rank, &dims_1d);
    
    // Data on dust bins & radiation bins independent of hydro cell
    rank = 1;
    dims_1d = Nabins - 2;
    ierr = my_createDataset("agrain" , rank, &dims_1d);
    dims_1d = Nabins - 1;
    ierr = my_createDataset("agrain_binEdges", rank, &dims_1d);
    
    if(dust_useRadiation){
        dims_1d = nphotBins;
        ierr = my_createDataset("photonBins", rank, &dims_1d);
    }

    // Data specific to each cell
    rank = 2;
    dims_2d[0] = NCELLS - 2*NGHOST;
    dims_2d[1] = NdustBins;
    ierr = my_createDataset("number", rank, dims_2d);
    ierr = my_createDataset("slope" , rank, dims_2d);
    ierr = my_createDataset("dadt"  , rank, dims_2d);
    
    // for each radiation bin
    if(dust_useRadiation){
        rank = 3;
        dims_3d[0] = NCELLS - 2*NGHOST;
        dims_3d[1] = NdustBins;
        dims_3d[2] = nphotBins;
        ierr = my_createDataset("opticalDepth", rank, dims_3d);
    }

    // fill up constant or cell independent dataspaces
    // cell radius
    rank = 1;
    start_1d = 0;
    stride_1d = 1;
    count_1d = NCELLS - 2*NGHOST;
    ierr = my_writeToDataset("radius", rs + NGHOST, rank, &start_1d, &stride_1d, &count_1d);


    // agrains
    rank = 1;
    start_1d  = 0;
    stride_1d = 1;
    count_1d  = Nabins - 2;
    // +1 to avoid ghost bins
    ierr = my_writeToDataset("agrain", abin_c + 1, rank, &start_1d, &stride_1d, &count_1d);
   
    // bin edges
    rank = 1;
    start_1d  = 0;
    stride_1d = 1;
    count_1d  = Nabins - 1;
    ierr = my_writeToDataset("agrain_binEdges", abin_e + 1, rank, &start_1d, &stride_1d, &count_1d);
    
    if(dust_useRadiation){
        //radiation bins
        rank = 1;
        start_1d  = 0;
        stride_1d = 1;
        count_1d  = nphotBins;
        ierr = my_writeToDataset("photonBins", dust_ELphots, rank, &start_1d, &stride_1d, &count_1d);
    }
    return 1;
}

int outputDustCell(int icell, double dr){
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
        printf("ERROR cant write to dataset\n");
        return -1;
    }
    
    ierr = my_writeToDataset("slope" , slope  + 1, rank, start_2d, stride_2d, count_2d);
    if(ierr < 0){
        printf("ERROR cant write to dataset\n");
        return -1;
    }
    
    ierr = my_writeToDataset("dadt"  , dadt   + 1, rank, start_2d, stride_2d, count_2d);
    if(ierr < 0){
        printf("ERROR cant write to dataset\n");
        return -1;
    }
    // if we hhave two species, we need to do these separate to avoid ghost cells
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


int finalizeDustOutput(){
    H5Fclose(dustOutput);
    return 1;
}
#endif
