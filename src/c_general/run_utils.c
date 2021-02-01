#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cgeneral.h>
#include <sys/types.h>
#include <sys/stat.h>
/*
 * Methods that are usefull in general usage, not exclusive to FLASH/HDF5 etc.
 */

// Search methods 
int binarySearch(double val, double *array, int size){
    int idx;
    int found = 0;
    int min = 0, max = size-1;     
    while(found == 0){
        idx = (max + min) / 2;

        // Check if at end points
        if(idx == size -1){
            if(val>=array[idx]){
                found = 1;
                break;
            }
        } else if (idx == 0){
            if(val < array[idx+1]){
                found = 1;
                break;
            }
        }

        if((val>=array[idx]) && (val<array[idx+1])){
            found = 1;
            break;
        }
        if((val>array[idx])){
            min = idx + 1;
        } else if(val<array[idx]) {
            max = idx - 1;
        }
    }
    return idx;
}   

int sequentialSearch(double val, double *array, int size){
    int i, idx;
    
    for(i=0 ; i < size; i++){
        if((val>=array[i]) && (val<array[i+1])){
            idx = i;
            break;
        }
    } 
    return idx;
}


// prints a progress bar
// Use inside loops
void DoProgress( int step, int total ){
    //progress width
    const int pwidth = 33;
    int i;
    //minus label len
    int width = pwidth;
    int pos = ( step * width ) / total ;

    
    int percent = ( step * 100 ) / total;

    printf( "[");

    //fill progress bar with =
    for ( i = 0; i < pos; i++ )  printf( "%c", '=' );

    //fill progress bar with spaces
    for ( i = 0; i<  (width - pos) ; i++) printf("%c", ' '); 
    
    printf( "]" );
    
    printf( " %3d%%\r", percent );

}

// Enforces periodic bountry conditions by checking 
// that  xmin < x < xmax, and otherwise
// moves x by L=(xmax-xmin) until that is true
double checkBnds(double x, double xmin, double xmax){
    double L = xmax - xmin;
    while(x<xmin){
        x = x+L;
    }
    while(x>xmax){
        x = x-L;
    }
    return x;
}
// 
// Recursive call to create dirs
//
int recmkdir(const char *dir) {
    char tmp[256];
    char *p = NULL;
    size_t len;
    int err;
    struct stat st = {0};
    snprintf(tmp, sizeof(tmp),"%s",dir);
    len = strlen(tmp);
    if(tmp[len - 1] == '/'){
        tmp[len - 1] = 0;
    }
    for(p = tmp + 1; *p; p++){
        if(*p == '/') {
            *p = 0;            
            if (stat(tmp, &st)==-1){
                err = mkdir(tmp, S_IRWXU);
                if(err < 0){
                    return -1;
                }
            }
            *p = '/';
        }
    }
    if (stat(tmp, &st)==-1){
        mkdir(tmp, S_IRWXU);
        if(err < 0){
            return -1;
        }
    }
}
