#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dustRadiation.h>
#include <cgeneral.h>
#include <dust.h>
#include <constantsAndUnits.h>
#ifdef useDust
int nphotBins;

double EbinMax;

double *dust_Ephots = NULL;
double *dust_Nphots = NULL;
double *dust_NphotsE = NULL;
double *dust_dEphots = NULL;
double *dust_ELphots = NULL;
double *dust_ERphots = NULL;

double *tauDust = NULL;
double *QabsBin = NULL;
double *NabsBinE = NULL;
double *EabsBinE = NULL;

int initDustRadiation(int Nbins, double *Ebins){
    // Initialise internal arrays for radiation bins in the dust model
    // In:
    //      Nbins   -   Number of radiation bins
    //      Ebins   -   edges of the radiation bins, size Nbins + 1
    
    if(dust_useRadiation == 0) return 1;
    int iphotBin,idx;
    nphotBins = Nbins;
    dust_Ephots  = (double *) malloc(nphotBins * sizeof(double));
    dust_Nphots  = (double *) malloc(nphotBins * sizeof(double));
    dust_NphotsE = (double *) malloc(nphotBins * sizeof(double));
    dust_dEphots = (double *) malloc(nphotBins * sizeof(double));
    dust_ELphots = (double *) malloc(nphotBins * sizeof(double));
    dust_ERphots = (double *) malloc(nphotBins * sizeof(double));
    QabsBin      = (double *) malloc(nphotBins * dust_nbins * sizeof(double));
    tauDust      = (double *) malloc(nphotBins * NdustBins  * sizeof(double));
    NabsBinE     = (double *) malloc(nphotBins * dust_nbins * sizeof(double));
    EabsBinE     = (double *) malloc(nphotBins * dust_nbins * sizeof(double));
    // store binsizes
    for( iphotBin = 0; iphotBin < nphotBins; iphotBin++){
        dust_ELphots[iphotBin] = Ebins[iphotBin];
        dust_ERphots[iphotBin] = Ebins[iphotBin + 1];
        dust_dEphots[iphotBin] = dust_ERphots[iphotBin] - dust_ELphots[iphotBin];
        for(int idustBin = 0; idustBin < dust_nbins; idustBin ++){
            idx = idustBin * nphotBins + iphotBin;
            QabsBin[idx] = 0;
        }
    }
    EbinMax = dust_ERphots[nphotBins -1];
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
    
    if(dust_useRadiation == 0) return 1;
    int iphotBin, idustBin, idx, graphite, iabin;
    double aveNu;
    for(iphotBin = 0; iphotBin < nphotBins; iphotBin++){
        dust_Nphots[iphotBin] = radData[iphotBin] * geoFact / dt;  // make into flux
        dust_NphotsE[iphotBin] = radData[iphotBin]/dust_dEphots[iphotBin] * geoFact / dt;  // make into flux
        dust_Ephots[iphotBin] = radData[iphotBin + nphotBins];
        aveNu = dust_Ephots[iphotBin]*planckInv;
        for(idustBin = 0; idustBin < dust_nbins; idustBin++){
            if(idustBin < isilicone){
                graphite = 1;
                iabin = idustBin;
            } else {
                graphite = 0;
                iabin = idustBin - isilicone;
            }
            idx = idustBin * nphotBins + iphotBin;
            if(idx < 0 || idx > nphotBins*dust_nbins){
                printf("%d %d %d\n", iphotBin, idx, nphotBins*dust_nbins);
            }
            QabsBin[idx] = getQabs(abin_c[iabin], aveNu, graphite, ida_tabQabs[idustBin], -1);
            NabsBinE[idx] = QabsBin[idx] * dust_NphotsE[iphotBin];
            EabsBinE[idx] = QabsBin[idx] * dust_NphotsE[iphotBin] * dust_Ephots[iphotBin];
        }
    }
    return 1;
}

double getDustOpticalDepth(int iEbin, double dr){
    double tauDustTot=0;
    if(dust_useRadiation == 0) return 0.0;
    int ibin, idx, graphite, iabin;
    for(idx = 0; idx < NdustBins; idx++){
        ibin = globalToLocalIndex(idx);
        if(ibin < isilicone){
            graphite = 1;
        } else {
            graphite = 0;
        }
        iabin = ibin - isilicone * (1-graphite); 
        tauDustTot += number[ibin] * QabsBin[ibin*nphotBins + iEbin] * pi_asquare[iabin];
    }
    return tauDustTot * dr; 
}

int setDustOpticalDepthArray(double dr){
    // only used for outputs
    int iphotBin;
    int ibin, idx, graphite, iabin, idq, idq2;
    for(idx = 0; idx < NdustBins; idx++){
        ibin = globalToLocalIndex(idx);
        if(ibin < isilicone){
            graphite = 1;
        } else {
            graphite = 0;
        }
        iabin = ibin - isilicone * (1-graphite); 
        
        for(iphotBin = 0; iphotBin < nphotBins; iphotBin++){
            idq  = idx * nphotBins + iphotBin;
            idq2 = ibin * nphotBins + iphotBin;
            tauDust[idq] = number[ibin] * QabsBin[idq2] * pi_asquare[iabin];
        }
    }
    return 1;
}
//loop over all photon bins and sum up absorbed energy
double getAbsorption(int ibin, int graphite){
     int iphotBin, iabin;
     double aveEphot, Qabs, EabsNu, EabsTot;
     EabsTot = 0;
     for (iphotBin = 0; iphotBin < nphotBins; iphotBin++){
          // average photon energy in bin
          aveEphot = dust_Ephots[iphotBin];
          
          Qabs = QabsBin[ibin*nphotBins + iphotBin]; //getQabs(abin_c[ibin], aveNu, graphite, ida_tabQabs[ibin], -1);
          EabsNu = dust_Nphots[iphotBin] * aveEphot * Qabs;
          EabsTot = EabsTot + EabsNu;
    
     }
     if(ibin < isilicone){
        iabin = ibin;
     } else {
        iabin = ibin - isilicone;
     }
     return EabsTot * pi_asquare[iabin];
}

//loop over all photon bins and sum up absorbed energy from all within specified range
double getAbsorption_range(int ibin, int graphite, double nuMin, double nuMax){
     int iphotBin, firstBin, iabin;
     double aveEphot, Qabs, EabsNu, EabsTot;
     double Emax = nuMax*planck;
     double Emin = nuMin*planck;
     double leftEdge, rightEdge;
     EabsTot = 0;
     firstBin = binarySearch(Emin, dust_ELphots, nphotBins);
     for (iphotBin = firstBin; iphotBin < nphotBins; iphotBin++){
         if(dust_ERphots[iphotBin] < Emin){
              continue; // Ebins from lowest to highest
          }
          if(dust_ELphots[iphotBin] > Emax){
              break;
          } 
          //Figure out how much of range is within the bin.
          leftEdge  = fmax(Emin, dust_ELphots[iphotBin]);
          rightEdge = fmin(Emax, dust_ERphots[iphotBin]);
          aveEphot = dust_Ephots[iphotBin];
          Qabs = QabsBin[ibin*nphotBins + iphotBin]; //getQabs(abin_c[ibin], aveNu, graphite, ida_tabQabs[ibin], -1);
          EabsNu = dust_NphotsE[iphotBin] * aveEphot * Qabs * (rightEdge - leftEdge);
          EabsTot = EabsTot + EabsNu;
     }

     if(ibin < isilicone){
        iabin = ibin;
     } else {
        iabin = ibin - isilicone;
     }
     EabsTot = EabsTot * pi_asquare[iabin];
     return EabsTot;
}


//same but only number of photons
double getAbsorptionNum_range(int ibin, int graphite, double Emin, double Emax){
     int iphotBin, firstBin, iabin;
     double NabsNu, NabsTot;
     double leftEdge, rightEdge;
     if(EbinMax < Emin){
         return 0;
     }
     NabsTot = 0;
     int idx = ibin * nphotBins;
     if(Emin > 0){ 
         firstBin = binarySearch(Emin, dust_ELphots, nphotBins);
     } else {
         firstBin = 0;
     }
     for (iphotBin = firstBin; iphotBin < nphotBins; iphotBin++){
          if(dust_ERphots[iphotBin] < Emin){
              continue; // Ebins from lowest to highest
          }
          //Figure out how much of range is within the bin.
          leftEdge  = fmax(Emin, dust_ELphots[iphotBin]);
          rightEdge = fmin(Emax, dust_ERphots[iphotBin]);
          NabsNu = NabsBinE[idx + iphotBin] * (rightEdge - leftEdge);
          NabsTot = NabsTot + NabsNu;
     }
     if(ibin < isilicone){
        iabin = ibin;
     } else {
        iabin = ibin - isilicone;
     }
     return NabsTot*pi_asquare[iabin];
}


// Method to get energy density of photons of frequency freq (u_\nu)
double getRadEdens(double freq){
    int iphotBin;
    // Energy of desired photons
    double Ener = freq*planck;
    if(Ener > dust_ERphots[nphotBins - 1]){
       return 0;
    }
    if(Ener < dust_ELphots[0]){
       return 0;
    }

    iphotBin = binarySearch(Ener, dust_ELphots, nphotBins);
    return dust_Nphots[iphotBin]*dust_Ephots[iphotBin]*planck/clght/dust_dEphots[iphotBin];
}


//loop over all photon bins and calculate transition from one dust energy bin to another
double getUpwardTransition(int ibin, int graphite, double Ui, double UiL, double UiR, double dUi, double Uf, double UfL, double UfR, double dUf, double continousBound){
    int iphotBin, firstBin, idx, iabin;
    double NabsE;
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
        firstBin = binarySearch(Emin, dust_ELphots, nphotBins);
    } else {
        firstBin = 0;
    }
    for (iphotBin = firstBin; iphotBin < nphotBins; iphotBin++){
        //if(dust_ERphots[iphotBin] < Emin){ //SHOULD NOT BE NEEDED DUE TO FIRSTBIN
        //    continue; // Ebins from lowest to highest
        //} 
        EL = dust_ELphots[iphotBin];
        if(EL > Emax){
            break;
        } 
        ER = dust_ERphots[iphotBin];
        // We assume that the energy density of photons uE is constant in each photon bin.
        // We also assume that Qabs is constant within the bin which then means 
        //
        // c int G(E) uE * Cabs(E) dE = E_Eave Cabs(Eave) * int G(E) dE
        //
        // where G is the correction for dust energy bins being discreete 
        // And E_Eave is the number flux of photons with energy Eave (assumed to be constant)
        // (see eqs. 15 - 25 in Draine & Li 2001)
        
        NabsE = EabsBinE[idx + iphotBin]; // getQabs(abin_c[ibin], aveNu, graphite, ida_tabQabs[ibin], -1);
        //Qabs     = QabsBin[ibin*nphotBins + iphotBin]; // getQabs(abin_c[ibin], aveNu, graphite, ida_tabQabs[ibin], -1);
        //EabsNu   = dust_NphotsE[iphotBin] * Qabs;
        //EabsNu   = dust_Nphots[iphotBin] * aveEphot * Qabs / dust_dEphots[iphotBin];
        
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
    if(ibin < isilicone){
       iabin = ibin;
    } else {
       iabin = ibin - isilicone;
    }
    uppRate = uppRate * pi_asquare[iabin]  / dUi / (Uf - Ui);
    return uppRate;
}

double intrabinUpwardTransition(int ibin, int graphite, double dUi, double Ui, double Uf){
    int iphotBin, iabin;
    double EabsNu;
    double leftEdge, rightEdge;
    double dUiInv = 1/dUi; 
    double uppRate = 0;
    int idx = ibin * nphotBins;
    for (iphotBin = 0; iphotBin < nphotBins; iphotBin++){
        leftEdge  = dust_ELphots[iphotBin];
        if(leftEdge > dUi){
            break;
        } 
        EabsNu   = EabsBinE[idx + iphotBin]; //getQabs(abin_c[ibin], aveNu, graphite, ida_tabQabs[ibin], -1);
        rightEdge = fmin(dUi, dust_ERphots[iphotBin]); 
        uppRate = uppRate + EabsNu * ( (rightEdge - leftEdge) - 0.5 * (rightEdge*rightEdge - leftEdge*leftEdge)*dUiInv); 

    }
    if(ibin < isilicone){
       iabin = ibin;
    } else {
       iabin = ibin - isilicone;
    }
    uppRate = uppRate * pi_asquare[iabin]  / (Uf - Ui);
    return uppRate;
}

#endif
