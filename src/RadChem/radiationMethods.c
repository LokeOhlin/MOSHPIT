#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
//#include "hydro.h"
#include <radchem.h>
#include <constantsAndUnits.h>
#ifdef useDust
#include <dust.h>
#include <dustRadiation.h>
#endif
int useRadiationPressure;
double radEmin, radEmax;
int numBinsSubIon;
int numBinsFullIon;
int numRadiationBins;
int iE112, iE136, iE152;

//double clght  = 29979245800.0;
double ionE   = 2.17863955e-11;
double ionEH2 = 2.43530838e-11;
double dissE  = 1.79443775e-11;
double heatE  = 9.02930155e-12;

double sigma0   = 6.30e-18;
double sigmaH0  = 4.5125966e-18;

double dustAtt5  = 2.5;
double dustAtt11 = 4.0;
double dustAtt13 = 5.0;
double dustAtt15 = 6.0;
double dust_to_gas_ratio;
double AV_conversion_factor;
double Temp_lim = 1e3;
double fsh_par1 = 82.3529411765;
double fpump = 6.94;


double *EbinEdges;
double *aveEphots;
double *Nphots;

double *EionH;
double *EionH2;

double sigmaLW = 2.46e-18;
double *sionH;
double *sionH2;

#ifndef useDust
double *dustTau_perH;
#endif


void initRadiation(){
    int iEbin, readSEDFromFile;
    double delE;
    
    getrealchemistrypar("ch_dust_to_gas_ratio", &dust_to_gas_ratio);
    getrealchemistrypar("ch_AV_conversion_factor", &AV_conversion_factor);
   
    getrealchemistrypar("ch_dust_to_gas_ratio", &dust_to_gas_ratio);
    getrealchemistrypar("ch_AV_conversion_factor", &AV_conversion_factor);
    
    getintegerchemistrypar("radiationPressure", &useRadiationPressure);
    getintegerchemistrypar("numBinsSubIon", &numBinsSubIon);
    getintegerchemistrypar("numBinsFullIon", &numBinsFullIon);
    getrealchemistrypar("radEmin", &radEmin);
    getrealchemistrypar("radEmax", &radEmax);

    getintegerchemistrypar("readSEDFromFile", &readSEDFromFile);
    
    numRadiationBins = numBinsSubIon + numBinsFullIon + 2;
    iE112 = numBinsSubIon;
    iE136 = iE112  + 1;
    iE152 = iE136  + 1;

    EbinEdges = (double * ) malloc((numRadiationBins + 1) * sizeof(double));
    // populate bins below 11.2 eV
    // Do it in logscale TODO: make user switchable to linear
    delE = exp((log(dissE) - log(radEmin))/(numBinsSubIon + 1));
    for(iEbin = 0; iEbin < iE112; iEbin++){
        EbinEdges[iEbin] = radEmin * pow(delE, iEbin);
    }

    // 11.2 eV
    EbinEdges[iE112] = dissE;
    
    // 13.6 eV
    EbinEdges[iE136] = ionE;

    // 15.2 eV and above
    delE = exp((log(radEmax) - log(ionEH2))/(numBinsFullIon + 1)); 
    for(iEbin = iE152; iEbin < numRadiationBins + 1; iEbin++){
        EbinEdges[iEbin] = ionEH2 * pow(delE, iEbin - iE152);
    }

    aveEphots = (double * ) malloc((numRadiationBins) * sizeof(double));
    Nphots    = (double * ) malloc((numRadiationBins) * sizeof(double));
    
    
    EionH  = (double * ) malloc((numBinsFullIon + 1) * sizeof(double));
    EionH2 = (double * ) malloc((numBinsFullIon + 1) * sizeof(double));
    sionH  = (double * ) malloc((numBinsFullIon + 1) * sizeof(double));
    sionH2 = (double * ) malloc((numBinsFullIon + 1) * sizeof(double));
#ifdef useDust
    // init local dust arrays for radiation SED and dust absorption
    int ierr = initDustRadiation(numRadiationBins, EbinEdges);
    if(ierr < 0){
        printf("Error when initialising dust radiation bins\n");
    }
#else
    dustTau_perH = (double * ) malloc(numRadiationBins * sizeof(double));
#endif 

    if(readSEDFromFile){
        setFromFile();
    } else {
        setFromStellarModel();
    }    
}

void setRadiationData(double *radData, double dt){
        // set intial package data
        // Number of photons 
        int iEbin;
        for(iEbin = 0; iEbin < numRadiationBins; iEbin++){
            radData[iEbin] = Nphots[iEbin] * dt;     
            radData[iEbin + numRadiationBins] = aveEphots[iEbin]; 
        }
        radData[2*numRadiationBins]     = 0;
        radData[2*numRadiationBins + 1] = 0;
        radData[2*numRadiationBins + 2] = 0;
}


void cellAbsorption(double *radData, double *specData, double numd, double Temp, double dr, double volcell, double dt, double *absData){
    int iEbin;
    // Radiation specific variables
    double Ndis, Edis, nphots, ephots;
    // Species data
    double numdr, xH, xH2, xHp;
    // Column densities
    double Hcolm, HcolmOld, H2colm, H2colmOld, Av, AvOld, normH2, normH2Old;
#ifndef useDust
    double dAv;
#endif
    // Optical depths
    double taud11In, taud11Out, DtauH, DtauH2, DtauD, tauTot;
    // Attenuation from dust and shielding from H2 in 5.2-13.6 bins
    double eAttavg, eAttd11_avg, eDAttd11_avg, fshieldIn, fshield_avg, fabs_in, Ndis_in, dNH211, dNdust11;
    double DNionH, DNionH2, DNdust, NabsTot, dNdust;
    // Rates and associated energies
    double kUV, kion, phih, hvphih, kdis, phih2, hvphih2, EtotPe;
    // Momentum injection
    double DustMom, H2Mom, HMom, DustEabs, Habs_est, H2abs_est;
    // Ugly
    double tmp, tmp2, norm;

    DustEabs  = 0;
    DustMom   = 0;
    EtotPe    = 0;
    H2Mom     = 0;
    HMom      = 0;
    H2abs_est = 0;
    Habs_est  = 0;

    kUV       = 0;
    phih      = 0;
    phih2     = 0;
    hvphih    = 0;
    hvphih2   = 0;
    // Calculate number and local column densite of relevant species
    xH  = fmax(specData[0],1e-10); // Lower limit to avoid division by zero
    xH2 = fmax(specData[1],1e-10); 
    xHp = fmax(specData[2],1e-10); 
    

    numdr = numd*dr;
    // Total up untill this point
    HcolmOld  = radData[2*numRadiationBins];
    H2colmOld = radData[2*numRadiationBins + 1];
    AvOld     = radData[2*numRadiationBins + 2];

    Hcolm  = HcolmOld  + numdr*(2*xH2+xH+xHp);
    H2colm = H2colmOld + numdr*xH2;

#ifndef useDust
    dAv = numdr*(2.0*xH2+xH+xHp)*exp(-Temp/Temp_lim);
    Av = AvOld + dAv;
#endif

    /////////////////////////////////////////////////////
    //
    //  Emin - 11.2 eV
    //
    for(iEbin = 0; iEbin < iE112; iEbin++){    
        nphots = radData[iEbin];
        if(nphots <= 0){
            continue;
        }
        ephots = radData[iEbin + numRadiationBins];
#ifdef useDust
    // Calculate the total optical depth from dust of these photons and how much they absorb
        DtauD = getDustOpticalDepth(iEbin, dr); 
#else
        //TODO: change dustAtt to array with dust size averaged cross sections.
        DtauD    = dustTau_perH[iEbin] * dAv   * dust_to_gas_ratio; 
#endif 
        if(DtauD > 1e-11){
            dNdust = nphots * (1 - exp(-DtauD));
        } else {
            dNdust = nphots * ( DtauD - 0.5 * DtauD * DtauD);
        }
        DustEabs = DustEabs + dNdust * ephots/dt;
        DustMom  = DustMom  + dNdust * ephots/clght;
        if(DtauD > 1e-11) { 
            eAttavg = 1.0/DtauD * ( 1 - exp(-DtauD));
        } else {
            eAttavg = (1.0 - 0.5*DtauD);
        }

        // Add to photo electric heating effect
        EtotPe += nphots*eAttavg*ephots/dt;

        // remove photons from bin
        nphots = nphots - dNdust;
        if(nphots < 0){
            nphots = 0;
        }
        radData[iEbin] = nphots;
    }
    /////////////////////////////////////////////////////
    //
    //      11.2 - 13.6 eV
    //
    //      To properly deal with H2 self shielding,
    //      we store the extinction of this bin
    //      
    //      Ndis = number of photons originally in this bin
    //      nphots = current number of photons.

    nphots = radData[iE112];
    if(nphots > 0){
        Ndis      = Nphots[iE112];
        Edis      = radData[iE112+numRadiationBins];
     
#ifdef useDust
        // Calculate the total optical depth from dust of these photons and how much they absorb
        taud11In  = AvOld; 
        DtauD = getDustOpticalDepth(iE112, dr);
        taud11Out = AvOld + DtauD; 
        Av = taud11Out;
#else
        taud11In  = dustTau_perH[iE112] * AvOld * dust_to_gas_ratio; 
        taud11Out = dustTau_perH[iE112] * Av    * dust_to_gas_ratio; 
        DtauD     = dustTau_perH[iE112] * dAv   * dust_to_gas_ratio; 
#endif 
    
        // Photons absorbed by H2 and dust, calculate in and out attenuation for both
        if(DtauD > 1e-11) { 
            eAttd11_avg  = 1.0/DtauD * (exp(-taud11In) - exp(-taud11Out));
            eDAttd11_avg = 1.0/DtauD * (1.0 - exp(-DtauD));
        } else {
            eAttd11_avg  = exp(-taud11In)*(1.0 - 0.5*DtauD);
            eDAttd11_avg = (1.0 - 0.5*DtauD);
        }

        normH2Old = H2colmOld/5e14;
        normH2    = H2colm/5e14;

        // H2 self shielding at beginning of cell:
        fshieldIn  = 0.965/pow(1.0+normH2Old,2) + 0.035/sqrt(1.0+normH2Old)*exp(-8.5e-4*sqrt(1.0+normH2Old));
    
    
        if ((H2colm-H2colmOld) > 1e11) {
        //this formula comes from averaging the shielding function over the H2 column density
            tmp  = (0.965/(normH2Old +1.0) + fsh_par1 * exp(-8.5e-4*sqrt(1.0 + normH2Old)));
            tmp2 = (0.965/(normH2    +1.0) + fsh_par1 * exp(-8.5e-4*sqrt(1.0 + normH2)));
            fshield_avg = 1./(H2colm - H2colmOld) * 5e14*(tmp-tmp2);     
        } else {
            fshield_avg = fshieldIn;
        }

   
        fabs_in = (1.0 + 0.0117 * normH2Old / (1.0 + normH2Old) - exp(-8.5e-4 * (sqrt(normH2Old + 1.0) - 1.0))) / 1.0117;
        Ndis_in = Ndis * exp(-taud11In) * (1.0 - fabs_in)*dt;
        // H2 ABSORPITON 
        //

        // Dissosiation rate of H2
        kUV = Ndis * eAttd11_avg * fshield_avg*sigmaLW*dr/volcell;
        //printf("%.4e %4e %.4e %.4e %.4e %.4e\n", Ndis, fshield_avg, eAttd11_avg, sigmaLW, dr/volcell, kUV);
        // Number of absorbed photons by H2
        dNH211 = (1 + fpump)* kUV * xH2*numd*volcell*dt;
    
        // radiation pressure on H2  
        H2Mom     += dNH211*Edis/clght;  
        H2abs_est += dNH211;
        // DUST ABSORPTION
    
        // Assume H2 takes president.. should probably be done at the same time, but thats not very easy
        if(DtauD > 1e-11){
            dNdust11 = fmax(Ndis_in - dNH211,0) * (1.0 - exp(-DtauD));
        } else {
            dNdust11 = fmax(Ndis_in - dNH211,0) * (DtauD-0.5*DtauD*DtauD);
        }

        DustEabs  = DustEabs  + dNdust11*Edis/dt; 
        DustMom   = DustMom   + dNdust11*Edis/clght; 
        EtotPe = EtotPe + fmax(Ndis_in-dNH211,0)/dt * eDAttd11_avg * Edis;
        // Update photons in bin
        nphots = nphots - dNdust11 - dNH211;
        if(nphots < 0){
            nphots = 0;
        }
        radData[iE112] = nphots;
        
    }
    
    /////////////////////////////////////////////////////
    //
    //      13.6 eV - Emax
    //
    for(iEbin = iE136; iEbin < numRadiationBins; iEbin++){
        nphots = radData[iEbin];
        if(nphots <= 0){
            continue;
        }
        ephots = radData[iEbin + numRadiationBins];
        // Calculate optical depths 
        DtauH  = numdr * sionH [iEbin - iE136] * xH;
        DtauH2 = numdr * sionH2[iEbin - iE136] * xH2;
#ifdef useDust
        DtauD = getDustOpticalDepth(iEbin, dr); 
#else
        DtauD = dustTau_perH[iEbin] * dAv * dust_to_gas_ratio; 
#endif
        tauTot = DtauH + DtauH2 + DtauD;
        if(tauTot == 0.0){
            continue;
        }
        norm = 1/tauTot;

        // total number of photons that are absorbed
        if(tauTot > 1e-11) {
            NabsTot = nphots * (1.0-exp(-tauTot));
        }else{
            NabsTot = nphots * (tauTot-0.5*tauTot*tauTot);
        }

        // Atomic hydrogen
        DNionH = NabsTot * DtauH * norm;
        kion   = DNionH/(numd * xH * volcell * dt);
        phih     += kion;
        hvphih   += kion*(ephots-ionE);
        HMom     += DNionH*ephots/clght;
        Habs_est += DNionH; 
        
        // Molecular hydrogen
        DNionH2 = NabsTot * DtauH2 * norm;
        kdis      = DNionH2/(numd * xH2 * volcell *dt);
        if(iEbin == iE136){
            kUV += kdis;
        } else {
            phih2     += kdis;
            hvphih2   += kdis*(ephots-ionEH2); 
        }
        H2Mom     += DNionH2*ephots/clght;
        H2abs_est += DNionH2; 
        
        // Dust
        DNdust = NabsTot * DtauD * norm;
        DustEabs += DNdust * ephots/dt;
        DustMom  += DNdust * DtauD * radData[iEbin + numRadiationBins]/clght;

        // Remove photons
        nphots = nphots - NabsTot;
        if(nphots < 0 ) {
            nphots = 0;
        }
        radData[iEbin] = nphots;
    }

    // Update arrays
    // Number of remaining photons 
    radData[2*numRadiationBins] = Hcolm;
    radData[2*numRadiationBins + 1] = H2colm;
    radData[2*numRadiationBins + 2] = Av;

    // Rates and energies
    absData[0] = EtotPe;
    absData[1] = kUV;
    absData[2] = phih;
    absData[3] = phih2;
    absData[4] = hvphih;
    absData[5] = hvphih2;
    absData[6] = DustEabs;

    // Momentum
    absData[7] = DustMom;
    absData[8] = HMom;
    absData[9] = H2Mom;

    // Estimates of number of absorbed photons
    absData[10] = Habs_est;
    absData[11] = H2abs_est;
}



