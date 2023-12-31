#ifdef useChemistry
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <hydro.h>
#include <radchem.h>
#include <cgeneral.h>
#include <rtpars.h>
#ifdef useDust
    #include <dust.h>
    #include <dustRadiation.h>
#endif

int doChemistryStep(double dt, double *dt_chem){
    int ierr;
    int icell, idx;
    double radData[2*numRadiationBins + 3];
    double absData[12];
    double abhtot, abctot;
    double specData[5], non_eq_species[3];
    double numd, Temp, xH, xHnew, xH2, xHp, xCO, xCp, volcell;
    double rho, etot, eint, vel, pres;
    double NionH0, NdisH2;
    double kUV, phih, hvphih, phih2, hvphih2, EtotPe;
    double DustMom, H2Mom, HMom, DustEabs, Habs_est, H2abs_est;
    double maxIon_Change;
    double Tdust, fshield_CO, fshield_H2, Av_mean, chi_mean, divv, redshift, tphoto,energy; 
    double dt_temp, sum_abs_est, sum_abs;
#ifdef useDust
    double rhoDust, rhoDustNew, rpars[24];
    int ipars[2], idust,ibin;
#endif

    getrealchemistrypar("ch_max_ion_frac_change", &maxIon_Change);
    // If we want radiation 
    setRadiationData(radData, dt);
    
    dt_chem[0] = 1e99;
    sum_abs = 0;
    sum_abs_est = 0;
    for(icell = NGHOST; icell < NCELLS-NGHOST; icell++){
        idx  = icell*nvar;
        // Get cell state
        rho  = ustate[idx];
        vel  = ustate[idx+1]/rho;
        etot = ustate[idx+2];
        eint = (etot - 0.5*rho*vel*vel);
        pres = eint*(adi-1);
         
        // Get abundances
        xH  = (ustate[idx+ICHEM_START]  ) * mf_scale;
        xH2 = (ustate[idx+ICHEM_START+1]) * mf_scale/2.0;
        xHp = (ustate[idx+ICHEM_START+2]) * mf_scale;
 
        numd = rho/ch_mH/abar;
        Temp = pres/(numd*ch_kb*(1-xH2+xHp+abundHe));
        
        abhtot = 2 * xH2 + xHp + xH;
        if(abhtot != 1.0){
            xH = xH/abhtot;
            xH2 = xH2/abhtot;
            xHp = xHp/abhtot;
        }

        xCO = (ustate[idx+ICHEM_START+3]) * mf_scale/ch_muC;
        xCp = (ustate[idx+ICHEM_START+4]) * mf_scale/ch_muC;
        abctot = xCO + xCp;
        if(abctot != abundC){
            xCO = xCO * abundC / abctot;
            xCp = xCp * abundC / abctot;
        }
        specData[0] = xH;
        specData[1] = xH2;
        specData[2] = xHp;
        specData[3] = xCO;
        specData[4] = xCp;
        
        volcell = vol[icell];
#ifdef useDust        
        // Set cell data in dust arrays
        // Needs to be done before radiative transfer
        ierr = setBinsCell(icell, &rhoDust);
        if(ierr < 0) {
            return -1; 
        }

        ierr = setRadiationBins(radData, dt, 1/(4*M_PI*pow((rs[icell]),2)));
        if(ierr < 0) {
            return -1; 
        }
#endif
#ifdef savePhotonFluxes
        // save incomming fluxes
        ierr = setPhotonFluxes(icell, radData, 1/(4*M_PI*pow((rs[icell]),2)));
#endif
        // calculate cell absorption
        
        cellAbsorption(radData, specData, numd, Temp, dr[icell], volcell, dt, absData);

#ifdef useDust
        rpars[0] = rho;
        rpars[1] = numd;
        rpars[2] = Temp;
        
        // approximate  mean molecular weight
        rpars[3] = ch_mH*abar/(1-xH2+xHp+abundHe);
        
        // gas velocity and acceleration from radiation pressure
        rpars[4] = vel;
        rpars[5] = (absData[8]+absData[9])/(dt*rho*volcell);
        
        // sound speed in cell
        rpars[6] = sqrt(pres * adi);
#ifdef growthUpdateVelocities 
        //if(icell == 4){
        //    verbose = 1;
        //    printf("vel %.8e\n", vel);
        //}
#endif
        ierr = dustCell(rpars, ipars, dt);
        if(ierr < 0) {
            return -1; 
        }

#ifdef growthUpdateVelocities   
        // If we calculate the dust and gas velocity in the growth module, get the updated gas velocity velocity
        vel = rpars[4];
        //if(icell == 4){
        //    printf("vel %.8e\n", vel);
        //}
#else
        // Otherwise update dust velocity here 
        ierr = dustCalcRadPress(dt);
        if(ierr < 0) {
            return -1; 
        }
#endif 
        if(outputDust){
            // If we want to output this step, we need to write the dadt's of the dust here
            ierr = Dust_outputCell_dadt(icell, dr[icell]);
            if(ierr < 0) {
                return -1; 
            }
        }
        ierr = getBinsCell(icell, &rhoDustNew);
        if(ierr < 0) {
            return -1; 
        }
#endif
        // get absorption variables
        //make into flux. take value at centre of cell
        if(geometry == 1){
            EtotPe    = absData[0]/(4*M_PI*pow((rs[icell]),2)); 
        } else {
            EtotPe    = absData[0]*dr[icell]/volcell; //make into flux. take value at centre of cell
        } 
        kUV       = absData[1];

        phih      = absData[2];
        phih2     = absData[3];
        hvphih    = absData[4];
        hvphih2   = absData[5];
        DustEabs  = absData[6]/volcell;
        
        DustMom   = absData[7]; 
        
        HMom      = absData[8]; 
        H2Mom     = absData[9]; 
        
        Habs_est  = absData[10];
        H2abs_est = absData[11];  
        // Get dust temperature for chemistry
        // Is not advected, but its set to equilibrium anyways, just need initial guess

        Tdust = ustate[idx+ICHEM_START+5];    
        // set variables
        non_eq_species[0] = xH2;
        non_eq_species[1] = xHp;
        non_eq_species[2] = xCO;
        NionH0 = 0;
        NdisH2 = 0;
        energy = eint;
        // FOR NOW NO EXTERNAL RADIATION
        Av_mean    = 100;
        chi_mean   = 0.0;
        fshield_H2 = 0.0;
        fshield_CO = 0.0;

        // other random stuff in call
        divv     = 0;
        redshift = 0;
        tphoto   = 0;
        if(noChemistry == 0){
            // TODO  make this interface more general... this might be comiler specific
            evolve_abundances_(&dt, &dr[icell], &numd, &divv, &energy, &redshift, 
                               non_eq_species, &fshield_H2, &fshield_CO, &Av_mean, 
                               &chi_mean, &Tdust, &phih, &hvphih, &phih2, &hvphih2, &kUV,
                               &EtotPe, &tphoto, &NionH0, &NdisH2, &DustEabs);  
    
            //printf(" %d : %.4e, %.4e, %4e %.4e \n", icell, eint, energy, (energy - eint)/eint, Tdust);
            sum_abs += NionH0*volcell;
            sum_abs_est += Habs_est; 
        //if( fabs(NionH0*volcell - Habs_est)/(NionH0*volcell) > 0.1) {
        //    printf("%d NionH0 %.4e %.4e   %.4e\n", icell, NionH0*volcell/dt, Habs_est/dt, xH);
        //}
#ifndef growthUpdateVelocities
            if(useRadiationPressure > 0){
#ifndef useDustDynamics   
                vel = vel + DustMom/(volcell*rho);
#endif
                NdisH2 = NdisH2 *  volcell ;
                NionH0 = NionH0 *  volcell ;
                if(H2abs_est > 1){
                    vel = vel + NdisH2/H2abs_est * H2Mom/(rho*volcell);
                }
                if(Habs_est > 1){
                    vel = vel + NionH0/Habs_est * HMom/(rho*volcell);
                }
            }
#endif
        }
        etot = energy + 0.5*rho*vel*vel;
        ustate[idx+1] = rho*vel;
        ustate[idx+2] = etot;
        ustate[idx + ICHEM_END] = Tdust;

        xHnew = 1-2*non_eq_species[0]-non_eq_species[1];
        if(fabs(xH-xHnew)>maxIon_Change){
            dt_temp =  maxIon_Change*dt/fabs(xHnew-xH);
            if(dt_chem[0] > dt_temp){
                dt_chem[0] = dt_temp;
            }
        } 
        // Unpack results of chemistry integration
        ustate[idx+ICHEM_START]   = (1-2*non_eq_species[0]-non_eq_species[1])/mf_scale;
        ustate[idx+ICHEM_START+1] = non_eq_species[0]*2.0/mf_scale;
        ustate[idx+ICHEM_START+2] = non_eq_species[1]/mf_scale;
        ustate[idx+ICHEM_START+3] = non_eq_species[2]*ch_muC/mf_scale;
        ustate[idx+ICHEM_START+4] = (abundC - non_eq_species[2])*ch_muC/mf_scale;
    }
    return 1;
}

#endif
