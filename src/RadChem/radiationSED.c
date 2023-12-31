#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <radchem.h>
#include <constantsAndUnits.h>

#ifndef useDust
#include <cgeneral.h>
#include <constantDust.h>

double getTauDustConstant(double ELbin, double ERbin){
    int ifbin, ifbinp;
    double freqmin = ELbin*planckInv;
    double freqmax = ERbin*planckInv;
    double freq, freqp, leftEdge, rightEdge;

    ifbin = binarySearch(freqmin, MNR_abs_freqs, MNR_abs_Nfbins);
    ifbinp = ifbin + 1;
    double intTau = 0;
    double totFreq = 0;
    double delf, delTau, tau, taup;
    while(1) {
        if(ifbinp >= MNR_abs_Nfbins){
            printf("Energy bin exceeds the tabulated dust absorption values %.4e %.4e\n", freqmax, MNR_abs_freqs[MNR_abs_Nfbins]);
        }

        freq  = MNR_abs_freqs [ifbin];
        freqp = MNR_abs_freqs [ifbinp];
        tau   = MNR_abs_tauTot[ifbin]; 
        taup  = MNR_abs_tauTot[ifbinp];

        leftEdge  = fmin(freq, freqmin);
        rightEdge = fmax(freqp, freqmax);
        delf = rightEdge - leftEdge;
        totFreq += delf;

        delTau = taup - tau;
        intTau += tau*(rightEdge - leftEdge) + delTau/delf * (rightEdge*(0.5*rightEdge - freq) - leftEdge*(0.5*leftEdge - freq));
        if(freqp >= freqmax){
            break;
        }
        ifbin++;
        ifbinp++;
    }
    return intTau/totFreq;
}
#endif


//Piecewise H2 cross section
//{0, 0 <= nu <= 15.2}, {0.09*Mb, 15.2 < nu <= 15.45}, 
//{1.15*Mb, 15.45 < nu <= 15.7}, {3.00*Mb, 15.7 < nu <= 15.95}, 
//{5.00*Mb, 15.95 < nu <= 16.2}, {6.75*Mb, 16.2 < nu <= 16.4}, 
//{8.00*Mb, 16.4 < nu <= 16.65}, {9.0*Mb, 16.65 < nu <= 16.85}, 
//{9.5*Mb, 16.85 < nu <= 17.00}, {9.8*Mb, 17. < nu <= 17.2}, 
//{10.10*Mb, 17.2 < nu <= 17.65}, {9.85*Mb, 17.65 < nu < 18.1}, 
//{9.75*Mb*(18.1/nu)^3, 18.10 < nu}

double tab_es[12] = {2.43530838e-11, 2.47509e-11, 2.51514e-11,
                     2.55519e-11, 2.59524e-11, 2.62728e-11,
                     2.66733e-11, 2.69937e-11, 2.72340e-11,
                     2.75544e-11, 2.82753e-11, 2.89962e-11};


double tab_ss[12] = {0.09e-18, 1.15e-18, 3.0e-18, 
                     5.00e-18, 6.75e-18, 8.0e-18,
                     9.00e-18, 9.50e-18, 9.8e-18,
                     10.1e-18, 9.85e-18, 9.75e-18};

double get_sigmaH2(double ephot){
    int iEbin;
    if(ephot < ionEH2){
        return 0; //if below direct dissociation, return 0
    }
    if(ephot > tab_es[11]){
        return tab_ss[11] * pow(ephot/tab_es[11], -3);
    }
    for(iEbin = 0; iEbin < 11; iEbin++){
        if(ephot < tab_es[iEbin + 1]){
            return tab_ss[iEbin];
        }
    }
    return 0;
}

double get_sigmaH(double ephot){
    int iEbin;
    if(ephot < ionE){
        return 0;
    }
#ifdef FIXEDSIGMAHI
    printf("%.4e\n", FIXEDSIGMAHI);
    return FIXEDSIGMAHI;
#else
    return sigmaH0 * pow(ephot/ionE,-3.);
#endif
}

void setFromStellarModel(){
    double Edis, Eion, sion, sigma0_e;
    double Npeh0, Ndis0, Nion0, NionH20;
    double Teff, Lstar, RstarSQ;
    double tmp = 1.0;
    int iEbin;
    
    getrealchemistrypar("Tstar", &Teff);
    getrealchemistrypar("Lstar", &Lstar);

    RstarSQ = Lstar * Lsun/(4*M_PI*Teff*Teff*Teff*Teff*sigma_stefanBoltz);

    // treat the special bins
    // 13.6 to 15.2
    getsourceemission_(&Teff, &ionE, &ionEH2, &sigma0, &Eion, &Nion0, &sion);
    Nion0 = Nion0*4.0*M_PI*RstarSQ;
    sion  = sion*4.0*M_PI*RstarSQ/Nion0;

    Nphots[iE136]    = Nion0; 
    aveEphots[iE136] = Eion;
#ifdef FIXEDSIGMAHI
    sionH [0] = FIXEDSIGMAHI;
#else
    sionH [0] = sion;
#endif
    sionH2[0] = sigmaLW;
    
    Eion  = Eion - ionE;
    EionH [0] = Eion;
    EionH2[0] = Eion;
    
#ifndef useDust
    dustTau_perH[iE136] = getTauDustConstant(ionE, ionEH2);
#endif

    //11.2 - 13.6
    getemissionsigma_(&Teff, &dissE, &ionE, &tmp, &Edis, &Ndis0);
    Ndis0 = Ndis0*4.0*M_PI*RstarSQ;

    Nphots[iE112]    = Ndis0;
    aveEphots[iE112] = Edis; 
    
#ifndef useDust
    dustTau_perH[iE112] = getTauDustConstant(dissE, ionE);
#endif

    // 15.2 to Emax
    for(iEbin = 0; iEbin < numBinsFullIon; iEbin++){
        sigma0_e = sigmaH0 * pow(EbinEdges[iEbin + iE152] / ionEH2, -3);
        getsourceemission_(&Teff, &EbinEdges[iEbin + iE152], &EbinEdges[iEbin + 1 + iE152], &sigma0_e, &Eion, &NionH20, &sion);
        NionH20 = NionH20*4.0*M_PI*RstarSQ;
        Nphots[iEbin + iE152] = NionH20;
#ifdef FIXEDSIGMAHI
        sionH[iEbin + 1] = FIXEDSIGMAHI;
#else
        sionH[iEbin + 1] = sion*4.0*M_PI*RstarSQ/NionH20;
#endif
        aveEphots[iEbin + iE152] = Eion; 
        EionH[iEbin + 1] = Eion - ionE;


        directh2_(&Teff, &EbinEdges[iEbin + iE152], &EbinEdges[iEbin + 1 + iE152], &Eion, &sion);
        sionH2[iEbin+1] = sion*4.0*M_PI*RstarSQ/NionH20;
        EionH2[iEbin+1] = Eion - ionEH2;
        
#ifndef useDust
        dustTau_perH[iEbin+iE152] = getTauDustConstant(EbinEdges[iEbin + iE152], EbinEdges[iEbin + iE152 + 1]);
#endif
    } 

    //5.6 -11.2
    for(iEbin = 0; iEbin < numBinsSubIon; iEbin++){
        getemissionsigma_(&Teff, &EbinEdges[iEbin], &EbinEdges[iEbin+1], &tmp, &Edis, &Npeh0);
        Nphots[iEbin] = Npeh0*4.0*M_PI*RstarSQ;
        aveEphots[iEbin] = Edis;
#ifndef useDust
        dustTau_perH[iEbin] = getTauDustConstant(EbinEdges[iEbin], EbinEdges[iEbin + 1]);
#endif
    }
}


void setFromFile(){
    FILE *fptr = fopen("sed.dat", "r");
    char line[128];
    char delimiter[2] = ",";
    int nFileBins, iEbin;
    char *numphots_s, *eave_s;
    double numphots, eave;
    iEbin = 0;
    nFileBins = 0;
    while(fgets(line, 128, fptr) != NULL){
        numphots_s = strtok(line, delimiter);
        if(numphots_s == NULL){
            continue;
        } 
        eave_s = strtok(NULL, delimiter);
        if(eave_s == NULL){
            continue;
        }
        trim(numphots_s);
        trim(eave_s);
        numphots = atof(numphots_s);
        eave     = atof(eave_s);
        printf("%d %.4e,%.4e\n",iEbin,numphots,eave);
        nFileBins += 1;
        if(iEbin >= numRadiationBins){
            printf("ERROR: number of radiation bins in sed.dat (>%d) greater than defined number of bins (%d)\n", nFileBins, numRadiationBins);
            exit(-1);
        }

        Nphots[iEbin] = numphots;
        aveEphots[iEbin] = eave;
        if(eave < EbinEdges[iEbin] || eave > EbinEdges[iEbin + 1]){
            printf("WARNING: Missmatch between bins as defined and in file\n");
            printf("         Bin %d with Emin = %.4e, Emax = %.4e was given an average photon energy of %.4e \n", iEbin, EbinEdges[iEbin], EbinEdges[iEbin + 1], eave);
        }
        if(iEbin >= iE136){
            sionH [iEbin - iE136] = get_sigmaH(eave);
            sionH2[iEbin - iE136] = get_sigmaH2(eave);
        
            EionH [iEbin - iE136] = aveEphots[iEbin] - ionE;
            if(iEbin == iE136){
                EionH2[iEbin - iE136] = aveEphots[iEbin] - ionE;
            } else {
                EionH2[iEbin - iE136] = aveEphots[iEbin] - ionEH2;
            }

        }
#ifndef useDust
        dustTau_perH[iEbin] = getTauDustConstant(EbinEdges[iEbin], EbinEdges[iEbin + 1]);
#endif
        iEbin++;
    }
    if(nFileBins != numRadiationBins){
            printf("ERROR: number of radiation bins in sed.dat (%d) not equal to defined number of bins (%d)\n", nFileBins, numRadiationBins);
            exit(-1);
    }
}
