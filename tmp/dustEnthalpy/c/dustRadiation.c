#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "dustRadiation.h"
#include "cgeneral.h"
#include "dust.h"

int nphotBins;
double *Ephots = NULL;
double *Nphots = NULL;
double *dEphots = NULL;
double *ELphots = NULL;
double *ERphots = NULL;


int initDustRadiation(int Nbins, double *Ebins){
    // Initialise internal arrays for radiation bins in the dust model
    // In:
    //      Nbins   -   Number of radiation bins
    //      Ebins   -   edges of the radiation bins, size Nbins + 1
    int iphotBin;
    nphotBins = Nbins;
    Ephots  = (double *) malloc(nphotBins * sizeof(double));
    Nphots  = (double *) malloc(nphotBins * sizeof(double));
    dEphots = (double *) malloc(nphotBins * sizeof(double));
    ELphots = (double *) malloc(nphotBins * sizeof(double));
    ERphots = (double *) malloc(nphotBins * sizeof(double));
    
    // store binsizes
    for( iphotBin = 0; iphotBin < nphotBins; iphotBin++){
        ELphots[iphotBin] = Ebins[iphotBin];
        ERphots[iphotBin] = Ebins[iphotBin + 1];
        dEphots[iphotBin] = ERphots[iphotBin] - ELphots[iphotBin];
    }
    
    return 1;
}

int setRadiationBins(double *radData, double dt, double geoFact){
    // sets the energy bins internal to the dust module
    // In :
    //      radData  -  input data for the local radiation field 
    //                  order: 1)nphotBins of total number of photons [unitless]
    //                         2)nphotBins of average energy of photons in bin [erg]
    //                         3)Rest is unused in this module
    //      dt       -  size of current time step. Needed to convert number to rate
    //      geofact  -  inverted size of surface. Needed to convert rate to flux
    
    int iphotBin;
    FILE *fptr = fopen("raddata.dat", "w");
    for(iphotBin = 0; iphotBin < nphotBins; iphotBin++){
        Nphots[iphotBin] = radData[iphotBin] * geoFact / dt;  // make into flux
        Ephots[iphotBin] = radData[iphotBin + nphotBins];
        fprintf(fptr, "%.4e %.4e %.4e\n", (ELphots[iphotBin] + ERphots[iphotBin])*0.5,Nphots[iphotBin], Ephots[iphotBin]);
    }
    fclose(fptr);
    return 0;
}


//loop over all photon bins and sum up absorbed energy
double getAbsorption(int ibin, int graphite){
     int iphotBin;
     double aveNu, aveEphot, Qabs, EabsNu, EabsTot;
     EabsTot = 0;
     for (iphotBin = 0; iphotBin < nphotBins; iphotBin++){
          // average photon energy in bin
          aveEphot = Ephots[iphotBin];
          aveNu = aveEphot*planckInv;
          
          Qabs = getQabs(agrains[ibin], aveNu, graphite, ida_tabQabs[ibin], -1);
          EabsNu = Nphots[iphotBin] * aveEphot * Qabs;
          EabsTot = EabsTot + EabsNu;
    
     }
     return EabsTot * pi_asquare[ibin];
}

//loop over all photon bins and sum up absorbed energy from all within specified range
double getAbsorption_range(int ibin, int graphite, double nuMin, double nuMax){
     int iphotBin;
     double aveNu, aveEphot, Qabs, EabsNu, EabsTot;
     double Emax = nuMax*planck;
     double Emin = nuMin*planck;
     double leftEdge, rightEdge;
     EabsTot = 0;
     for (iphotBin = 0; iphotBin < nphotBins; iphotBin++){
          if(ELphots[iphotBin] > Emax){
              break;
          } 
          if(ERphots[iphotBin] < Emin){
              continue; // Ebins from lowest to highest
          }
          //Figure out how much of range is within the bin.
          leftEdge  = fmax(Emin, ELphots[iphotBin]);
          rightEdge = fmin(Emax, ERphots[iphotBin]);
          aveEphot = Ephots[iphotBin];
          aveNu = aveEphot*planckInv;
          Qabs = getQabs(agrains[ibin], aveNu, graphite, ida_tabQabs[ibin], -1);
          EabsNu = Nphots[iphotBin] * aveEphot * Qabs * (rightEdge - leftEdge)/dEphots[iphotBin];
          EabsTot = EabsTot + EabsNu;
     }
     EabsTot = EabsTot * pi_asquare[ibin];
     return EabsTot;
}

//loop over all photon bins and calculate transition from one dust energy bin to another
double getUpwardTransition(int ibin, int graphite, double Ui, double UiL, double UiR, double dUi, double Uf, double UfL, double UfR, double dUf, double continousBound){
    int iphotBin;
    double aveNu, aveEphot, Qabs, EabsNu, EabsTot;
    double leftEdge, rightEdge;
    double finiteWidthFactor;
    
    double Emin = fmax(UfL - UiR, continousBound);
    double Emax = UfR - UiL;
    //if(Emax < continousBound){
    //    return 0;
    //}
    
    double MaxEdgeDistance = fmax(UfR - UiR, UfL - UiL);
    double MinEdgeDistance = fmin(UfR - UiR, UfL - UiL);
   
    double uppRate = 0;
    for (iphotBin = 0; iphotBin < nphotBins; iphotBin++){
        if(ELphots[iphotBin] > Emax){
            break;
        } 
        if(ERphots[iphotBin] < Emin){
            continue; // Ebins from lowest to highest
        }
        
        aveEphot = Ephots[iphotBin];
        aveNu = aveEphot*planckInv;
        Qabs = getQabs(agrains[ibin], aveNu, graphite, ida_tabQabs[ibin], -1);
        // 1) 
        leftEdge  = fmax(Emin, ELphots[iphotBin]);
        if(leftEdge < MinEdgeDistance){
            rightEdge = fmin(MinEdgeDistance, ERphots[iphotBin]);
            if(aveEphot > Emin){
                finiteWidthFactor = aveEphot - Emin;
            } else {
                finiteWidthFactor = rightEdge - Emin;
            }
            EabsNu = Nphots[iphotBin] * aveEphot * Qabs * (rightEdge - leftEdge)/dEphots[iphotBin];
            uppRate = uppRate + finiteWidthFactor * EabsNu; 
        }
        if(ERphots[iphotBin] < MinEdgeDistance){
              continue;
        }
        // 2)
        leftEdge  = fmax(MinEdgeDistance, ELphots[iphotBin]);
        if(leftEdge < MaxEdgeDistance){
            rightEdge = fmin(MaxEdgeDistance, ERphots[iphotBin]);
            finiteWidthFactor = fmin(dUf, dUi);
            EabsNu = Nphots[iphotBin] * aveEphot * Qabs * (rightEdge - leftEdge)/dEphots[iphotBin];
            uppRate = uppRate + finiteWidthFactor * EabsNu; 
        }
        
        if(ERphots[iphotBin] < MaxEdgeDistance){
              continue;
        }
        // 3)

        leftEdge  = fmax(MaxEdgeDistance, ELphots[iphotBin]);
        rightEdge = fmin(Emax, ERphots[iphotBin]);
        
        if(aveEphot > Emax){ 
            finiteWidthFactor = Emax - MaxEdgeDistance;
        } else {
            finiteWidthFactor = Emax - aveEphot;
        }
        EabsNu = Nphots[iphotBin] * aveEphot * Qabs * (rightEdge - leftEdge)/dEphots[iphotBin];
        uppRate = uppRate + finiteWidthFactor * EabsNu; 
        
    }
    uppRate = uppRate * pi_asquare[ibin]  / dUi / (Uf - Ui);
    return uppRate;
}



//same but only number of photons
double getAbsorptionNum_range(int ibin, int graphite, double nuMin, double nuMax){
     int iphotBin;
     double aveNu, aveEphot, Qabs, EabsNu, EabsTot;
     double Emax = nuMax*planck;
     double Emin = nuMin*planck;
     double leftEdge, rightEdge;
     EabsTot = 0;
     for (iphotBin = 0; iphotBin < nphotBins; iphotBin++){
          if(ELphots[iphotBin] > Emax){
              break;
          } 
          if(ERphots[iphotBin] < Emin){
              continue; // Ebins from lowest to highest
          }
          //Figure out how much of range is within the bin.
          leftEdge  = fmax(Emin, ELphots[iphotBin]);
          rightEdge = fmin(Emax, ERphots[iphotBin]);
          aveEphot = Ephots[iphotBin];
          aveNu = aveEphot*planckInv;
          Qabs = getQabs(agrains[ibin], aveNu, graphite, ida_tabQabs[ibin], -1);
          EabsNu = Nphots[iphotBin] * Qabs * (rightEdge - leftEdge)/dEphots[iphotBin];
          EabsTot = EabsTot + EabsNu;
     }
     return EabsTot*pi_asquare[ibin];
}


// Method to get energy density of photons of frequency freq (u_\nu)
double getRadEdens(double freq){
     int iphotBin;
     double aveNu, aveEphot, Qabs, EabsNu, EabsTot;
     // Energy of desired photons
     double Ener = freq*planck;
     EabsTot = 0;
     // if(Ener > ERphots[nphotBins - 1]){
     //    return 0;
     // }
     // if(Ener < ELphots[0]){
     //    return 0;
     // }

     iphotBin = binarySearch(Ener, ELphots, nphotBins);
     return Nphots[iphotBin]*Ephots[iphotBin]*planck/clght/dEphots[iphotBin];
}



