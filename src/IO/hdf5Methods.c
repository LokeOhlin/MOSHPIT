#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <IO.h>
#include <hydro.h>
#include <rtpars.h>
#include <mechanicalFeedback.h>
#include <radchem.h>
#include <dust.h>
/*
 *  Create the input file and call all relevant modules for parameter and headers
 *  We crea
 */
hid_t h5Output;
int createOutputFile(int output_id){
    char outname [12] = "";
    char basename[8]  = "output_";
    char outnum[5]    = "";
    int ierr;

    sprintf(outnum, "%04d", output_id);
    strcat(outname, basename);
    strcat(outname, outnum);

    printf("Creating %s \n", outname);
    h5Output = H5Fcreate(outname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Runtime parameters
    ierr = RuntimePars_initIO();
    if(ierr < 0){
        printf("Error: cannot write runtime parameters\n");
        exit(0);
    }
    // Hydro
    ierr = Hydro_initIO();
    if(ierr < 0){
        printf("Error: cannot write hydro parameters\n");
        exit(0);
    }

    // Mechanical feedback
    ierr = Feedback_initIO();
    if(ierr < 0){
        printf("Error: cannot write mechanical feedback parameters\n");
        exit(0);
    }
#ifdef useChemistry 
    // Chemistry
    ierr = RadChem_initIO();
    if(ierr < 0){
        printf("Error: cannot write Radiation and chemistry parameters\n");
        exit(0);
    }
#ifdef useDust
    // Mechanical feedback
    ierr = Dust_initIO();
    if(ierr < 0){
        printf("Error: cannot write dust parameters\n");
        exit(0);
    }
#endif
#endif
    return 1;
}

int finalizeOutputFile(){
    int ierr;
    // Tell all modules to write their end of timestep output values
    //
    
    // Runtime parameters
    ierr = RuntimePars_output();
    if(ierr < 0){
        printf("Error: cannot write runtime\n");
        exit(0);
    }

    // Hydro
    ierr = Hydro_output();
    if(ierr < 0){
        printf("Error: cannot write hydro\n");
        exit(0);
    }

    // Mechanical feedback
    ierr = Feedback_output();
    if(ierr < 0){
        printf("Error: cannot write feedback\n");
        exit(0);
    }
#ifdef useChemistry
    // Chemistry
    ierr = RadChem_output();
    if(ierr < 0){
        printf("Error: cannot write radchem\n");
        exit(0);
    }
#ifdef useDust
    // Dust
    ierr = Dust_output();
    if(ierr < 0){
        printf("Error: cannot write dust\n");
        exit(0);
    }
#endif
#endif

    // close output
    H5Fclose(h5Output);
    printf("Output closed\n");
    return 1;
}

int my_createDataset(char *datasetName, int rank, hsize_t *dims){
    hid_t dataset, datatype, dataspace;
    dataspace = H5Screate_simple(rank, dims, NULL);
    datatype  = H5Tcopy(H5T_NATIVE_DOUBLE);
    
    dataset   = H5Dcreate(h5Output, datasetName, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    return 1;
}

int my_writeToDataset(char *datasetName, double *varArr, int rank, hsize_t *start, hsize_t *stride, hsize_t *count){
    hid_t dataset, dataspace, memspace;
    herr_t hstat;
    
    // open dataspace and dataset
    dataset = H5Dopen(h5Output, datasetName, H5P_DEFAULT);
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

int my_writeAttribute(const char *groupName, const char *name, const void *buff, hid_t type){
    hid_t group, dataspace, attribute;
    herr_t status;

    // see if group exists, otherwise create it 
    status = H5Eset_auto1(NULL,NULL);
    status = H5Gget_objinfo(h5Output, groupName, 0, NULL);
    if(status == 0){
        group = H5Gopen(h5Output, groupName, H5P_DEFAULT);
    } else {
        group = H5Gcreate(h5Output, groupName, 0, H5P_DEFAULT, H5P_DEFAULT);
    }
    status = H5Eset_auto1(NULL,NULL);

    dataspace = H5Screate(H5S_SCALAR);
    attribute = H5Acreate(group, name, type, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    status    = H5Awrite(attribute, type, buff);
    H5Aclose(attribute);
    H5Sclose(dataspace);
    H5Gclose(group);   
    return 1;
}
