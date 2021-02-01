#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
//#include "hydro.h"
#include "radchem.h"

int useRadiationPressure;

double clght  = 29979245800.0;
double ionE   = 2.17863955e-11;
double ionEH2 = 2.43530838e-11;
double dissE  = 1.79443775e-11;
double heatE  = 9.02930155e-12;


double dustAtt5  = 2.5;
double dustAtt11 = 4.0;
double dustAtt13 = 5.0;
double dustAtt15 = 6.0;
double dust_to_gas_ratio;
double AV_conversion_factor;
double Temp_lim = 1e3;
double fsh_par1 = 82.3529411765;
double fpump = 6.94;

double Npeh0 = 0;
double Ndis0 = 0;
double Nion0 = 0;
double NionH20 = 0;

double Edis;
double Eion;
double EionH;
double EionH2;

double sigmaLW = 2.46e-18;
double sion;
double sionH;
double sionH2;

void initRadiation(){
        double Teff;
        double Lstar; 
        double RstarSQ;
        double inf = -1.0;
        double sigma0   = 6.30e-18;
        double sigma0H  = 4.5125966e-18;
        double boltzy = 5.6704e-5;
        double Lsun   = 3.838e+33;
        double tmp = 1.0, tmp2;
        
        getrealchemistrypar("ch_dust_to_gas_ratio", &dust_to_gas_ratio);
        getrealchemistrypar("ch_AV_conversion_factor", &AV_conversion_factor);
   
        getintegerchemistrypar("radiationPressure", &useRadiationPressure);
   
        getrealchemistrypar("Tstar", &Teff);
        getrealchemistrypar("Lstar", &Lstar);

        RstarSQ = Lstar * Lsun/(4*M_PI*Teff*Teff*Teff*Teff*boltzy);

        // 13.6 to 15.2
        getsourceemission_(&Teff, &ionE, &ionEH2, &sigma0, &Eion, &Nion0, &sion);
        Nion0 = Nion0*4.0*M_PI*RstarSQ;
        sion  = sion*4.0*M_PI*RstarSQ/Nion0;
        Eion  = Eion - ionE;

        // 15.2 to inf
        getsourceemission_(&Teff, &ionEH2, &inf, &sigma0H, &EionH, &NionH20, &sionH);
        NionH20 = NionH20*4.0*M_PI*RstarSQ;
        sionH = sionH*4.0*M_PI*RstarSQ/NionH20;
        EionH = EionH - ionE;

        directh2_(&Teff, &ionEH2, &inf, &EionH2, &sionH2);
        sionH2 = sionH2*4.0*M_PI*RstarSQ/NionH20;
        EionH2 = EionH2 - ionEH2;
        printf("Nphot > 13.6  = %.4e \n", NionH20 + Nion0);
        //11.2 - 13.6
        getemissionsigma_(&Teff, &dissE, &ionE, &tmp, &Edis, &Ndis0);
        Ndis0 = Ndis0*4.0*M_PI*RstarSQ;

        //5.6 -11.2
        getemissionsigma_(&Teff, &heatE, &dissE, &tmp, &tmp2, &Npeh0);
        Npeh0 = Npeh0*4.0*M_PI*RstarSQ*tmp2;
}
void setRadiationData(double *radData, double dt){
        // set intial package data
        // Number of photons 
        radData[0] = Npeh0;     // 5.6 - 11.2  eV no dt here...
        radData[1] = Ndis0*dt;     // 11.2 - 13.6 eV
        radData[2] = Nion0*dt;     // 13.6 - 15.2 eV
        radData[3] = NionH20*dt;   // 15.2 - inf  eV

        // Energies per photon  
        radData[4] = Edis;      // 11.2 - 13.6 eV
        radData[5] = Eion;      // 13.6 - 15.2 eV here its the ionisation energy eg Eion = Ephot - 13.6
        radData[6] = EionH;     // 15.2 - inf  eV ionisation of H  (EionH = Ephot - 13.6)
        radData[7] = EionH2;    // 15.2 - inf  eV ionisation of H2 (EionH2 = Ephot - 15.2)

        // Photon average cross sections
        radData[8]  = sion;       // 13.6 - 15.2 eV atomic hydrogen
        radData[9]  = sionH;      // 15.2 - inf  eV atomic hydrogen
        radData[10] = sionH2;     // 15.2 - inf  eV molecular hydrogen
    
        // Columdensities and cumulative attenuations
        radData[11] = 0;
        radData[12] = 0;
        radData[13] = 0;
}


void cellAbsorption(double *radData, double *specData, double numd, double Temp, double dr, double volcell, double dt, double *absData){
    // Radiation specific variables
    double Npeh, Ndis, Nion, NionH2, sion, sionH, sionH2, Edis, Eion, EionH, EionH2;
    // Species data
    double numdr, xH, xH2, xHp;
    // Column densities
    double Hcolm, HcolmOld, H2colm, H2colmOld, Av, AvOld, dAv, normH2, normH2Old;
    // Optical depths
    double taud5In, taud5Out, Dtaud5, taud11In, taud11Out, Dtaud11, DtauH, DtauH2, DtauD, tauTot;
    // Attenuation from dust and shielding from H2 in 5.2-13.6 bins
    double eAttd5_avg, eAttd11_avg, eDAttd11_avg, fshieldIn, fshield_avg, fabs_in, Ndis_in, dNH211, dNdust11;
    double DNionH, DNionH2, DNdust, NabsTot;
    // Rates and associated energies
    double kUV, kion, phih, hvphih, kdis, phih2, hvphih2, dEdust, EtotPe;
    // Momentum injection
    double DustMom, H2Mom, HMom, DustEabs, Habs_est, H2abs_est;
    // Ugly
    double tmp, tmp2, norm;

    // Get radiation data
    Npeh   = radData[0];
    Ndis   = radData[1];
    Nion   = radData[2];
    NionH2 = radData[3];
    Edis   = radData[4];
    Eion   = radData[5];
    EionH  = radData[6];
    EionH2 = radData[7];
    sion   = radData[8];
    sionH  = radData[9];
    sionH2 = 0;// radData[10];
    // initialise variables
    DustEabs  = 0;
    DustMom   = 0;
    EtotPe    = 0;
    H2Mom     = 0;
    HMom      = 0;
    H2abs_est = 0;
    Habs_est  = 0;

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
    HcolmOld  = radData[11];
    H2colmOld = radData[12];
    AvOld     = radData[13];

    Hcolm  = HcolmOld  + numdr*(2*xH2+xH+xHp);
    H2colm = H2colmOld + numdr*xH2;

#ifdef useDust
    // Calculate local absorption based on the size distribution of the local dust
    something_funny();
#else
    dAv = AV_conversion_factor*numdr*(2.0*xH2+xH2+xHp);//*exp(-Temp/Temp_lim);
#endif
    Av = AvOld + dAv;


    /////////////////////////////////////////////////////
    //
    //      5.6 - 11.2 eV
    //
        
#ifdef useDust
    // Calculate the total optical depth from dust of these photons and how much they absorb
    something_funnier();
#else
    taud5In  = dustAtt5 * AvOld * dust_to_gas_ratio; 
    taud5Out = dustAtt5 * Av    * dust_to_gas_ratio; 
    Dtaud5   = dustAtt5 * dAv   * dust_to_gas_ratio; 
    if(Dtaud5 > 1e-11){
        dEdust = Npeh * dt * (exp(-taud5In)- exp(-taud5Out));
    } else {
        dEdust = Npeh * dt * exp(-taud5In) * ( Dtaud5 - 0.5 * Dtaud5 * Dtaud5);
    }
    DustEabs = DustEabs + dEdust/dt;
#endif 
    DustMom = DustMom + dEdust/clght;
    if(Dtaud5 > 1e-11) { 
        eAttd5_avg = 1.0/Dtaud5 * (exp(-taud5In) - exp(-taud5Out));
    } else {
        eAttd5_avg = exp(-taud5In)*(1.0 - 0.5*Dtaud5);
    }

    // Add to photo electric heating effect
    EtotPe += Npeh*eAttd5_avg;
    
        
    /////////////////////////////////////////////////////
    //
    //      11.2 - 13.6 eV
    //
    
#ifdef useDust
    // Calculate the total optical depth from dust of these photons and how much they absorb
    something_funnier();
#else
    taud11In  = dustAtt11 * AvOld * dust_to_gas_ratio; 
    taud11Out = dustAtt11 * Av    * dust_to_gas_ratio; 
    Dtaud11   = dustAtt11 * dAv   * dust_to_gas_ratio; 
#endif 
    
    // Photons absorbed by H2 and dust, calculate in and out attenuation for both

    if(Dtaud11 > 1e-11) { 
        eAttd11_avg  = 1.0/Dtaud11 * (exp(-taud11In) - exp(-taud11Out));
        eDAttd11_avg = 1.0/Dtaud11 * (1.0 - exp(-Dtaud11));
    } else {
        eAttd11_avg  = exp(-taud11In)*(1.0 - 0.5*Dtaud11);
        eDAttd11_avg = (1.0 - 0.5*Dtaud11);
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
    Ndis_in = Ndis * exp(-taud11In) * (1.0 - fabs_in);
    
    // H2 ABSORPITON 
    //

    // Dissosiation rate of H2
    kUV = Ndis * eAttd11_avg * fshield_avg*sigmaLW*dr/(volcell*dt);
    // Number of absorbed photons by H2
    dNH211 = (1 + fpump)* kUV * xH2*numd*volcell*dt;
    
    // radiation pressure on H2  
    H2Mom     += dNH211*Edis/clght;  
    H2abs_est += dNH211;
    // DUST ABSORPTION
    
#ifdef useDust 
    something_funnier();
#else
    // Assume H2 takes president.. should probably be done at the same time, but thats not very easy
    if(Dtaud11 > 1e-11){
      dNdust11 = fmax(Ndis_in - dNH211,0) * (1.0 - exp(-Dtaud11));
    } else {
      dNdust11 = fmax(Ndis_in - dNH211,0) * (Dtaud11-0.5*Dtaud11*Dtaud11);
    }

    DustEabs = DustEabs + dNdust11*Edis/dt;
#endif
    DustMom  = DustMom  + dNdust11*Edis/clght; 
    
    EtotPe = EtotPe * fmax(Ndis_in-dNH211,0)/dt * eDAttd11_avg;
    
    /////////////////////////////////////////////////////
    //
    //      13.6 - 15.2 eV
    //
    
    if(Nion > 0) {
        // Calculate optical depths 
        DtauH  = numdr * sion * xH;
        DtauH2 = numdr * sigmaLW * xH2;
#ifdef useDust
        DtauD  = IMPLEMENT;
#else
        DtauD  = 0;//dustAtt13 * dAv * dust_to_gas_ratio; 
#endif
        tauTot = DtauH + DtauH2 + DtauD;

        if(tauTot == 0.0){
            goto noTau13;
        }
        norm = 1/tauTot;

        // total number of photons that are absorbed
        if(tauTot > 1e-11) {
            NabsTot = Nion * (1.0-exp(-tauTot));
        }else{
            NabsTot = Nion * (tauTot-0.5*tauTot*tauTot);
        }

        // Atomic hydrogen

        DNionH = NabsTot * DtauH * norm;
        kion = DNionH/(numd * xH * volcell * dt);
        phih     += kion;
        hvphih   += kion*Eion;
        HMom     += DNionH*(Eion + ionE)/clght;
        Habs_est += DNionH; 
        // Molecular hydrogen
        DNionH2 = NabsTot * DtauH2 * norm;
        kUV       += DNionH2/(numd * xH2 * volcell *dt);
        H2Mom     += DNionH2*(Eion + ionE)/clght;
        H2abs_est += DNionH2; 
        
        // Dust
        DNdust = NabsTot * DtauD * norm;
        
#ifdef useDust
        something_hilarious();
#else
        DustEabs = DustEabs + DNdust*(Eion + ionE)/dt;
#endif
        DustMom += DNdust * DtauD * (Eion + ionE)/clght;


        // Remove photons
        Nion = Nion - NabsTot;
        if(Nion < 0 ) {
            Nion = 0;
        }
    }
noTau13:
    
    /////////////////////////////////////////////////////
    //
    //      15.2 - inf eV
    //
    
    if(NionH2 > 0) {
        // Calculate optical depths 
        DtauH  = numdr * sionH  * xH;
        DtauH2 = numdr * sionH2 * xH2;
#ifdef useDust
        DtauD  = IMPLEMENT;
#else
        DtauD  = 0;dustAtt15 * dAv   * dust_to_gas_ratio; 
#endif
        tauTot = DtauH + DtauH2 + DtauD;

        if(tauTot == 0.0){
            goto noTau15;
        }
        norm = 1/tauTot;

        // total number of photons that are absorbed
        if(tauTot > 1e-11) {
            NabsTot = NionH2 * (1.0-exp(-tauTot));
        }else{
            NabsTot = NionH2 * (tauTot-0.5*tauTot*tauTot);
        }

        // Atomic hydrogen

        DNionH = NabsTot * DtauH * norm;
        kion = DNionH/(numd * xH * volcell * dt);
        phih     += kion;
        hvphih   += kion*EionH;
        HMom     += DNionH*(EionH + ionE)/clght;
        Habs_est += DNionH; 

        // Molecular hydrogen
        DNionH2 = NabsTot * DtauH2 * norm;
        kdis      = DNionH2/(numd * xH2 * volcell *dt);
        phih2     += kdis;
        hvphih2   += kdis*EionH2;
        H2Mom     += DNionH2*(EionH2 + ionEH2)/clght;
        H2abs_est += DNionH2; 
        
        // Dust
        DNdust = NabsTot * DtauD * norm;
        
#ifdef useDust
        something_hilarious();
#else
        DustEabs = DustEabs + DNdust*(EionH + ionE)/dt;
#endif
        DustMom += DNdust * (EionH + ionE)/clght;
        
        
        // Remove photons
        NionH2 = NionH2 - NabsTot;
        if(NionH2 < 0 ) {
            NionH2 = 0;
        }
    }
noTau15:


    // Update arrays

    // Number of remaining photons 
    radData[2] = Nion;     // 13.6 - 15.2 eV
    radData[3] = NionH2;   // 15.2 - inf  eV
    
    // Updated attenuation
    radData[11] = Hcolm;
    radData[12] = H2colm;
    radData[13] = Av;

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



