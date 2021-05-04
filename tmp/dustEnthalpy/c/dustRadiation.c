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
double EbinMax;
double *Ephots = NULL;
double *Nphots = NULL;
double *NphotsE = NULL;
double *dEphots = NULL;
double *ELphots = NULL;
double *ERphots = NULL;
double *QabsBin = NULL;
double *NabsBinE = NULL;
double *EabsBinE = NULL;

int initDustRadiation(int Nbins, double *Ebins){
    // Initialise internal arrays for radiation bins in the dust model
    // In:
    //      Nbins   -   Number of radiation bins
    //      Ebins   -   edges of the radiation bins, size Nbins + 1
    int iphotBin;
    nphotBins = Nbins;
    Ephots  = (double *) malloc(nphotBins * sizeof(double));
    Nphots  = (double *) malloc(nphotBins * sizeof(double));
    NphotsE = (double *) malloc(nphotBins * sizeof(double));
    dEphots = (double *) malloc(nphotBins * sizeof(double));
    ELphots = (double *) malloc(nphotBins * sizeof(double));
    ERphots = (double *) malloc(nphotBins * sizeof(double));
    QabsBin = (double *) malloc(nphotBins * dust_nbins * sizeof(double));
    NabsBinE = (double *) malloc(nphotBins * dust_nbins * sizeof(double));
    EabsBinE = (double *) malloc(nphotBins * dust_nbins * sizeof(double));
    // store binsizes
    for( iphotBin = 0; iphotBin < nphotBins; iphotBin++){
        ELphots[iphotBin] = Ebins[iphotBin];
        ERphots[iphotBin] = Ebins[iphotBin + 1];
        dEphots[iphotBin] = ERphots[iphotBin] - ELphots[iphotBin];
    }
    EbinMax = ERphots[nphotBins -1];
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
    
    int iphotBin, idustBin, idx, graphite;
    double aveNu;
    FILE *fptr = fopen("raddata.dat", "w");
    for(iphotBin = 0; iphotBin < nphotBins; iphotBin++){
        Nphots[iphotBin] = radData[iphotBin] * geoFact / dt;  // make into flux
        NphotsE[iphotBin] = radData[iphotBin]/dEphots[iphotBin];  // make into flux
        Ephots[iphotBin] = radData[iphotBin + nphotBins];
        fprintf(fptr, "%.4e %.4e %.4e\n", (ELphots[iphotBin] + ERphots[iphotBin])*0.5,Nphots[iphotBin], Ephots[iphotBin]);
        aveNu = Ephots[iphotBin]*planckInv;
        for(idustBin = 0; idustBin < dust_nbins; idustBin++){
            graphite = (idustBin < isilicone);
            idx = idustBin * nphotBins + iphotBin;
            QabsBin[idx] = getQabs(agrains[idustBin], aveNu, graphite, ida_tabQabs[idustBin], -1);
            NabsBinE[idx] = QabsBin[idx] * NphotsE[iphotBin];
            EabsBinE[idx] = QabsBin[idx] * NphotsE[iphotBin] * Ephots[iphotBin];
        }
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
          
          Qabs = QabsBin[ibin*nphotBins + iphotBin]; //getQabs(agrains[ibin], aveNu, graphite, ida_tabQabs[ibin], -1);
          EabsNu = Nphots[iphotBin] * aveEphot * Qabs;
          EabsTot = EabsTot + EabsNu;
    
     }
     return EabsTot * pi_asquare[ibin];
}

//loop over all photon bins and sum up absorbed energy from all within specified range
double getAbsorption_range(int ibin, int graphite, double nuMin, double nuMax){
     int iphotBin, firstBin;
     double aveNu, aveEphot, Qabs, EabsNu, EabsTot;
     double Emax = nuMax*planck;
     double Emin = nuMin*planck;
     double leftEdge, rightEdge;
     EabsTot = 0;
     firstBin = binarySearch(Emin, ELphots, nphotBins);
     for (iphotBin = 0; iphotBin < nphotBins; iphotBin++){
          if(ERphots[iphotBin] < Emin){
              continue; // Ebins from lowest to highest
          }
          if(ELphots[iphotBin] > Emax){
              break;
          } 
          //Figure out how much of range is within the bin.
          leftEdge  = fmax(Emin, ELphots[iphotBin]);
          rightEdge = fmin(Emax, ERphots[iphotBin]);
          aveEphot = Ephots[iphotBin];
          Qabs = QabsBin[ibin*nphotBins + iphotBin]; //getQabs(agrains[ibin], aveNu, graphite, ida_tabQabs[ibin], -1);
          EabsNu = NphotsE[iphotBin] * aveEphot * Qabs * (rightEdge - leftEdge);
          EabsTot = EabsTot + EabsNu;
     }
     EabsTot = EabsTot * pi_asquare[ibin];
     return EabsTot;
}


//same but only number of photons
double getAbsorptionNum_range(int ibin, int graphite, double Emin, double Emax){
     int iphotBin, firstBin;
     double aveNu, aveEphot, Qabs, NabsNu, NabsTot;
     double leftEdge, rightEdge;
     if(EbinMax < Emin){
         return 0;
     }
     NabsTot = 0;
     int idx = ibin * nphotBins;
     if(Emin > 0){ 
         firstBin = binarySearch(Emin, ELphots, nphotBins);
     } else {
         firstBin = 0;
     }
     for (iphotBin = firstBin; iphotBin < nphotBins; iphotBin++){
          if(ERphots[iphotBin] < Emin){
              continue; // Ebins from lowest to highest
          }
          //Figure out how much of range is within the bin.
          leftEdge  = fmax(Emin, ELphots[iphotBin]);
          rightEdge = fmin(Emax, ERphots[iphotBin]);
          NabsNu = NabsBinE[idx + iphotBin] * (rightEdge - leftEdge);
          NabsTot = NabsTot + NabsNu;
     }
     return NabsTot*pi_asquare[ibin];
}


// Method to get energy density of photons of frequency freq (u_\nu)
double getRadEdens(double freq){
    int iphotBin;
    double aveNu, aveEphot, EabsNu, EabsTot;
    // Energy of desired photons
    double Ener = freq*planck;
    EabsTot = 0;
    if(Ener > ERphots[nphotBins - 1]){
       return 0;
    }
    if(Ener < ELphots[0]){
       return 0;
    }

    iphotBin = binarySearch(Ener, ELphots, nphotBins);
    return Nphots[iphotBin]*Ephots[iphotBin]*planck/clght/dEphots[iphotBin];
}


//loop over all photon bins and calculate transition from one dust energy bin to another
double getUpwardTransition(int ibin, int graphite, double Ui, double UiL, double UiR, double dUi, double Uf, double UfL, double UfR, double dUf, double continousBound){
    int iphotBin, firstBin, idx;
    double aveNu, aveEphot, Qabs, NabsE, EabsTot;
    double leftEdge, rightEdge, ER, EL;
    double finiteWidthFactor;
    double Emin = UfL - UiR;
    if(EbinMax < Emin){
        return 0;
    }
    
    double Emax = UfR - UiL;
    
    double MaxEdgeDistance = fmax(UfR - UiR, UfL - UiL);
    double MinEdgeDistance = fmin(UfR - UiR, UfL - UiL);
    double mindU = fmin(dUi, dUf); 
    double uppRate = 0;
    idx = ibin * nphotBins;
    if(Emin > 0){ 
        firstBin = binarySearch(Emin, ELphots, nphotBins);
    } else {
        firstBin = 0;
    }
    for (iphotBin = firstBin; iphotBin < nphotBins; iphotBin++){
        //if(ERphots[iphotBin] < Emin){ //SHOULD NOT BE NEEDED DUE TO FIRSTBIN
        //    continue; // Ebins from lowest to highest
        //} 
        EL = ELphots[iphotBin];
        if(EL > Emax){
            break;
        } 
        ER = ERphots[iphotBin];
        // We assume that the energy density of photons uE is constant in each photon bin.
        // We also assume that Qabs is constant within the bin which then means 
        //
        // c int G(E) uE * Cabs(E) dE = E_Eave Cabs(Eave) * int G(E) dE
        //
        // where G is the correction for dust energy bins being discreete 
        // And E_Eave is the number flux of photons with energy Eave (assumed to be constant)
        // (see eqs. 15 - 25 in Draine & Li 2001)
        
        NabsE = EabsBinE[idx + iphotBin]; // getQabs(agrains[ibin], aveNu, graphite, ida_tabQabs[ibin], -1);
        //Qabs     = QabsBin[ibin*nphotBins + iphotBin]; // getQabs(agrains[ibin], aveNu, graphite, ida_tabQabs[ibin], -1);
        //EabsNu   = NphotsE[iphotBin] * Qabs;
        //EabsNu   = Nphots[iphotBin] * aveEphot * Qabs / dEphots[iphotBin];
        
        // 1) 
        if(EL < MinEdgeDistance){
            leftEdge  = fmax(Emin, EL);
            rightEdge = fmin(MinEdgeDistance, ER);
            finiteWidthFactor = 0.5*(rightEdge*rightEdge  - leftEdge*leftEdge) - Emin * (rightEdge - leftEdge);
            uppRate = uppRate + finiteWidthFactor * NabsE; 
        }

        // 2)
        if(ER < MinEdgeDistance){
              continue;
        }
        if(EL < MaxEdgeDistance){
            leftEdge  = fmax(MinEdgeDistance, EL);
            rightEdge = fmin(MaxEdgeDistance, ER);
            finiteWidthFactor = mindU * (rightEdge - leftEdge);
            uppRate = uppRate + finiteWidthFactor * NabsE; 
        }

        // 3)
        if(ER < MaxEdgeDistance){
              continue;
        }
        leftEdge  = fmax(MaxEdgeDistance, EL);
        rightEdge = fmin(Emax, ER);
        finiteWidthFactor = Emax * (rightEdge - leftEdge) - 0.5 * (rightEdge*rightEdge - leftEdge*leftEdge);
        uppRate = uppRate + finiteWidthFactor * NabsE; 

    }
    uppRate = uppRate * pi_asquare[ibin]  / dUi / (Uf - Ui);
    return uppRate;
}

double intrabinUpwardTransition(int ibin, int graphite, double dUi, double Ui, double Uf){
    int iphotBin;
    double aveNu, aveEphot, Qabs, EabsNu, EabsTot;
    double leftEdge, rightEdge;
    double dUiInv = 1/dUi; 
    double uppRate = 0;
    int idx = ibin * nphotBins;
    for (iphotBin = 0; iphotBin < nphotBins; iphotBin++){
        leftEdge  = ELphots[iphotBin];
        if(leftEdge > dUi){
            break;
        } 
        EabsNu   = EabsBinE[idx + iphotBin]; //getQabs(agrains[ibin], aveNu, graphite, ida_tabQabs[ibin], -1);
        rightEdge = fmin(dUi, ERphots[iphotBin]); 
        uppRate = uppRate + EabsNu * ( (rightEdge - leftEdge) - 0.5 * (rightEdge*rightEdge - leftEdge*leftEdge)*dUiInv); 

    }
    uppRate = uppRate * pi_asquare[ibin]  / (Uf - Ui);
    return uppRate;
}


