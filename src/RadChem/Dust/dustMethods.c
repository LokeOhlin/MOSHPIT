#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dust.h>
#include <hydro.h>
#include <radchem.h>
#include <moshpit.h>
#include <dustRadiation.h>
#include <constantsAndUnits.h>
#ifdef useDust
/////////////////////////////////////////
//
//  General utility methods
//
/////////////////////////////////////////

int globalToLocalIndex(int idx){
    int ibin;
    // idx = index in global dust bins (eg whats stored in cells)
    // ibin = index in local dust bins (eg what is used internally)
    // need to convert between the two
    ibin = idx + dust_nghost;
    // if ibin is outside of the number of size bins, 
    // that means that we are on a different species and therefore
    // need to account for 2 additional ghost cells
    if(ibin >= Nabins - dust_nghost){
        ibin += 2*dust_nghost;
    }

    return ibin;
}
int localToGlobalIndex(int ibin){
    int idx;
    // two options: 
    // 1) multiple species and we're on the second one, 
    //    need to remove three ghost cell
    // 2) single species or at the first species
    //    remove one ghost cell
    if(ibin > Nabins) {
        idx = ibin - 3*dust_nghost;
    } else {
        idx = ibin - dust_nghost;
    }
    return idx;
}

//////////////////////////////////////////
//
//  Methods to estimate the velocity of the dust grains
//
//////////////////////////////////////////

// get average dust velocity over the timestep relative to gas
double get_average_vdust(int ibin, int iabin, int graphite, double *rpars, double dt){

    // take the current velocity
    double vinit = dust_vrel[ibin];

    double numd    = rpars[0];
    double Tgas    = rpars[1]; 
    double mean_mw = rpars[2];
    double rhoGas    = numd * ch_kb * abar;
    // acceleration of the gas due to radiation pressure
    double gas_accel = rpars[3];

    // get the energy absorption rate
    double Eabs_dust = getAbsorption(ibin, graphite);
    // acceleration of the dust grain due to radiation pressure
    double mass;
    if(graphite){
        mass = aveMatom_c * Natoms[ibin] * aveMatom_c;
    } else {
        mass = aveMatom_s * Natoms[ibin] * aveMatom_s;
    }
    double dust_accel = Eabs_dust/clght/mass;

    // relative acceleration (always positive)
    double rad_accel = fabs(dust_accel - gas_accel);
    
    // the friction force from gas is proportional to the relative velocity
    // calculate the coefficient
    double fric_coeff = 4*sqrt(8*boltzmann*Tgas/(M_PI*mean_mw))*rhoGas*pi_asquare[iabin]/3;

    // calculate the average velocity based of the solution to
    // dv/dt + fric_coeff*v = rad_accell
    double mean_vdust = rad_accel/fric_coeff - (rad_accel/pow(fric_coeff,2) - vinit/fric_coeff)*(1-exp(-fric_coeff*dt))/dt;
    return mean_vdust;
}   

int update_vdust(int ibin, int iabin, int graphite, double *rpars, double dt){
    // take the current velocity
    double vinit = dust_vrel[ibin];

    double numd    = rpars[0];
    double Tgas    = rpars[1]; 
    double mean_mw = rpars[2];
    double rhoGas    = numd * ch_kb * abar;
    // acceleration of the gas due to radiation pressure
    double gas_accel = rpars[3];

    // get the energy absorption rate
    double Eabs_dust = getAbsorption(ibin, graphite);
    // acceleration of the dust grain due to radiation pressure
    double mass;
    if(graphite){
        mass = aveMatom_c * Natoms[ibin] * aveMatom_c;
    } else {
        mass = aveMatom_s * Natoms[ibin] * aveMatom_s;
    }
    double dust_accel = Eabs_dust/clght/mass;

    // relative acceleration (always positive)
    double rad_accel = fabs(dust_accel - gas_accel);
    
    // the friction force from gas is proportional to the relative velocity
    // calculate the coefficient
    double fric_coeff = 4*sqrt(8*boltzmann*Tgas/(M_PI*mean_mw))*rhoGas*pi_asquare[iabin]/3/mass;

    // solution to :
    // dv/dt + fric_coeff*v = rad_accell
    dust_vrel[ibin] = rad_accel*(1-exp(-fric_coeff*dt))/fric_coeff + vinit*exp(-fric_coeff*dt);
    return 1;

}


// loops over all local dustbins and sets them to their initial value for the global timestep
int set_vdusts(int icell){
    int idx, ibin;
    for(idx = 0; idx < dust_nbins; idx++){
        ibin = globalToLocalIndex(idx);
#ifdef trackDustVelocities
        //if the simulation is keeping track of the dust velocities we take these
        dust_vrel[ibin] = ustate[icell*nvar + IDUST_START + NdustVar*idx + 2];
#else
        // otherwise set this to zero
        dust_vrel[ibin] = 0.0;
#endif

    }
    // for safety, reset the velocities of the ghost bins
    dust_vrel[0] = 0.0;
    if(fSi < 1.0){
        dust_vrel[isilicone - 1] = 0.0;
    }
    if(fSi > 0.0){
        dust_vrel[isilicone] = 0.0;
    } 
    dust_vrel[dust_nbins - 1] = 0.0;
    return 1;
}

#ifdef trackDustVelocities
int get_vdust(icell){
    int idx, ibin;
    for(idx = 0; idx < dust_nbins; idx++){
        ibin = globalToLocalIndex(idx);
        //if the simulation is keeping track of the dust velocities we take these
        ustate[icell*nvar + IDUST_START + NdustVar*idx + 2] = dust_vrel[ibin];
    }
    return 1
}
#endif

//////////////////////////////////////////
//
//  Methods for evolution of dust distribution
//
//////////////////////////////////////////

// GRAIN DISTRIBUTION EVOLUTION
// GENERAL METHOD FROM McKinnon et al. 2018

// integral for Nj
double intNj(int ibin, int iabin, double a){
    return  number[ibin]*a/(abin_e[iabin+1] - abin_e[iabin]) + slope[ibin]*(a*a/2 - abin_c[iabin]*a);
}

// integral for Mj
double intMj(int ibin, int iabin, double a, double dt){
    double da  = dadt[ibin]*dt;
    double ac  = abin_c[iabin];
    double aep = abin_e[iabin+1];
    double ae  = abin_e[iabin];

    // integral (a+da)**3*(a-ac)
    double slope_fact = pow(a,5)/5. + (3*da-ac)*pow(a,4)/4.;
    slope_fact += da*(da-ac)*pow(a,3) + da*da*(da-3*ac)*a*a/2. +da*da*da*ac*a;

    return number[ibin]*pow(a+da,4)/(aep-ae)/4. + slope[ibin]*slope_fact;
}

//
// Assuming piecewise constant, solve for mass, slope or number based of the other two
// NOTE: WE DEAL WITH 4*pi*rho/3 LATER

// use Mj and Nj to solve for slope
double getSlope(double Nj, double Mj,int iabin){
    return (Mj - Nj*NfactM[iabin])/SfactM[iabin];

}
// use Mj and Sj to solve for Number
double getNumber(double Sj, double Mj,int iabin){
    if(Mj <= 0){
        return 0;
    }
    return (Mj - Sj*SfactM[iabin])/NfactM[iabin];

}

// use Nj and Sj to solve for Mass
double getMass(double Nj, double Sj,int iabin){
    if(Nj <= 0) {
        return 0;
    }
    return Nj*NfactM[iabin] + Sj*SfactM[iabin];
}

// Limit slope in cell if dnda goes negatie on either side.
// If both negative, something is wrong. Report and set all to zero
int limitSlope(double *Njnew, double *Sjnew, double Nj, double Sj, double Mj, int iabin){
    double ac  = abin_c[iabin];
    double aep = abin_e[iabin+1];
    double ae  = abin_e[iabin];
    double dac = aep - ae;
    double dap = aep - ac;
    double dam = ae  - ac;
    // Values at interfaces of the bins
    double dndap = Nj/dac + Sj*dap;
    double dnda  = Nj/dac + Sj*dam;
    // No grains, do nothing, set to zero
    if(Nj <= 0 || Mj <=0){
        *Njnew = 0;
        *Sjnew = 0;
    }
    // Both positive or both negative
    if( dndap*dnda > 0){
        //both negative set to zero
        if(dndap < 0){
            *Njnew = 0;
            *Sjnew = 0;
            return 1;
        }
        // Do nothing 
        *Njnew = Nj;
        *Sjnew = Sj;
        return 1;
    }

    double Nfact_tot;
    if(dndap < 0){
        Nfact_tot = (NfactM[iabin]-SfactM[iabin]/dap/dac);
        *Njnew = Mj/Nfact_tot;
        *Sjnew = -*Njnew/dac/dap;
    } else if(dnda < 0){
        Nfact_tot = (NfactM[iabin]-SfactM[iabin]/dam/dac);
        *Njnew = Mj/Nfact_tot;
        *Sjnew = -*Njnew/dac/dam;
    }
    return 1;    
}

// Rebinn such that the largest bin (amax,inf) contain zero stuff
// we keep the lowest bin size for the moment as a "resevoir" with zero opacity
// Assume all grains of size > amax have the size of amax. Update the average size of the grains in the bin
// And determine numbers and slope to satisfy both the new average size and new mass
int rebinn_upper(double Nn, double Nnp, double Sn, double Mn, double Mnp, double *Nntilde, double *Sntilde){
    double average_a, average_a_tilde;
    // Initial average
    average_a = NfactA[Nabins-2] + SfactA[Nabins-2]*Sn/Nn; 
    // Updated average
    average_a_tilde = (Nn*average_a + Nnp*abin_e[Nabins-1])/(Nn+Nnp); 

    // calculate Ntilde based on new mass
    *Nntilde = (Mn + Mnp) / (NfactM[Nabins-2] + (average_a_tilde - NfactA[Nabins-2])*SfactM[Nabins-2]/SfactA[Nabins-2]);
    // finally solve for Stilde
    *Sntilde = *Nntilde * (average_a_tilde - NfactA[Nabins-2])/SfactA[Nabins-2];  
    return 1; 
}
int rebinn_lower(double Nn, double Nnm, double Sn, double Mn, double Mnm, double *Nntilde, double *Sntilde){
    double average_a, average_a_tilde;
    // Initial average
    average_a = NfactA[1] + SfactA[1]*Sn/Nn; 
    // Updated average
    average_a_tilde = (Nn*average_a + Nnm*abin_e[1])/(Nn+Nnm); 

    // calculate Ntilde based on new mass
    *Nntilde = (Mn + Mnm) / (NfactM[1] + (average_a_tilde - NfactA[1])*SfactM[1]/SfactA[1]);
    // finally solve for Stilde
    *Sntilde = *Nntilde * (average_a_tilde - NfactA[1])/SfactA[1];  
    return 1; 
}
int rebinn(double dt){
    int ierr;
    double Mn, Mnp, Mnm, Nn, Nnp, Nnm, Nntilde, Sn, Sntilde;
    // If we have silicates
    if(fSi > 0.0){

        // upper boundary
        Mnp = Mnew[dust_nbins-1];
        if(Mnp > 0){ 
            // Does this boundary pile up? move dust grains to closest non-ghost cell 
            if(dust_upperBound_pileUp){
                Nnp = Mnp/pow(abin_e[Nabins-1], 3);
                Mn  = Mnew[dust_nbins-2];
                Nn  = Nnew[dust_nbins-2];
                Sn  = Snew[dust_nbins-2];
                // get updated bin values
                ierr = rebinn_upper(Nn, Nnp, Sn, Mn, Mnp, &Nntilde, &Sntilde);
                if(ierr < 0) {
                    return -1;
                }
                // Update the bins (and limit slope)
                Mnew[dust_nbins-2] = Mn + Mnp;        
                ierr = limitSlope(&Nnew[dust_nbins-2], &Snew[dust_nbins-2], Nntilde, Sntilde, Mn+Mnp, Nabins-2); 
                if(ierr < 0) {
                    return -1;
                }
            }
            Mnew[dust_nbins-1] = 0;      
            Nnew[dust_nbins-1] = 0;      
            Snew[dust_nbins-1] = 0;     
        } 
        
        
        // Lower boundary
        Mnm = Mnew[isilicone];
        
        if(Mnm > 0){ 
            if(dust_lowerBound_pileUp){
                Nnm = Mnm/pow(abin_e[1], 3);
                Mn  = Mnew[isilicone+1];
                Nn  = Nnew[isilicone+1];
                Sn  = Snew[isilicone+1];
                // get updated bin values
                ierr = rebinn_lower(Nn, Nnm, Sn, Mn, Mnm, &Nntilde, &Sntilde);
                if(ierr < 0) {
                    return -1;
                }

                // Update the bins (and limit slope)
                Mnew[isilicone+1] = Mn + Mnm;        
                ierr = limitSlope(&Nnew[isilicone+1], &Snew[isilicone+1], Nntilde, Sntilde, Mn+Mnm, 1); 
                if(ierr < 0) {
                    return -1;
                }
            } 
            
            Mnew[isilicone] = 0;      
            Nnew[isilicone] = 0;      
            Snew[isilicone] = 0;     
        } 

    }
    // if we have graphite grains
    if( isilicone > 0){
        // not same for carbonious grains
        Mnp = Mnew[isilicone-1];
        // Do we have grains in ghost bin? if so fix
        if(Mnp > 0){
            if(dust_upperBound_pileUp){
                Nnp = Mnp/pow(abin_e[Nabins-1],3);
                // get updated bin values
                Mn  = Mnew[isilicone-2];
                Nn  = Nnew[isilicone-2];
                Sn  = Snew[isilicone-2];
                
                ierr = rebinn_upper(Nn, Nnp, Sn, Mn, Mnp, &Nntilde, &Sntilde);
                if(ierr < 0) {
                    return -1;
                }
            
                 // Update the bins
                Mnew[isilicone-2] = Mn + Mnp;
                ierr = limitSlope(&Nnew[isilicone-2], &Snew[isilicone-2], Nntilde, Sntilde, Mn+Mnp, Nabins-2); 
                if(ierr < 0) {
                    return -1;
                }
            }
            Mnew[isilicone-1] = 0;      
            Nnew[isilicone-1] = 0;      
            Snew[isilicone-1] = 0;  
        }    
        
        Mnm = Mnew[0];

        if(Mnm > 0){ 
            if(dust_lowerBound_pileUp){
                Nnm = Mnm/pow(abin_e[1], 3);
                Mn  = Mnew[1];
                Nn  = Nnew[1];
                Sn  = Snew[1];
                // get updated bin values
                ierr = rebinn_lower(Nn, Nnm, Sn, Mn, Mnm, &Nntilde, &Sntilde);
                if(ierr < 0) {
                    return -1;
                }
                // Update the bins (and limit slope)
                Mnew[1] = Mn + Mnm;        
                ierr = limitSlope(&Nnew[1], &Snew[1], Nntilde, Sntilde, Mn+Mnm, 1); 
                if(ierr < 0) {
                    return -1;
                }
            }
            Mnew[0] = 0;      
            Nnew[0] = 0;      
            Snew[0] = 0;     
        } 
    }
    return 1; 
}

int dustCell(double *rpars, int *ipars, double dt_step){
    int idx, ibin, ibin2, iabin, iabin2, ierr;
    int rangeStart, rangeEnd, graphite;
    double x1, x2, intM, intN;
    double dt, time_dust, dadt_max_cell;
    double Nsum, Msum, Sest, NumTot_cell_g, NumTot_cell_s;
    double Mtot_old, Ntot_old, Mtot_new, Ntot_new, Ntot_mid;
    // Start by calculating change in size due to sublimation/sputtering/accretion
    //initialize by trying to take one full step
    time_dust = 0;
    // try to take all in one step
    dt = dt_step;
    
    // set the growth only dependent on current gas and radiation conditions.
    // eg does not change within loop
    ierr = set_dadt_fixed(rpars);
   

    // calculate total number of grains at start for normalisation
    NumTot_cell_g = 0;
    NumTot_cell_s = 0;
    for(ibin = 0; ibin < dust_nbins; ibin++){
        if(ibin < isilicone){
            NumTot_cell_g += number[ibin];
        } else {
            NumTot_cell_s += number[ibin];
        }
    }
    
    while(time_dust < dt_step){
        Mtot_old = 0;
        Ntot_old = 0;
        // update the growth rate
        ierr = set_dadt(rpars, dt);
        // Calculate change in each bin and limit timestep in needed
        for(ibin = 0; ibin < dust_nbins; ibin++){
            if(ibin < isilicone){
                graphite = 1;
            } else {
                graphite = 0;
            }
            // graphite and silicone have the same size bins. Use half array 
            iabin = ibin - isilicone*(1-graphite);    


            // Temperature of dust grains in this bin
            //temp         = getDustTemp(abin_c[iabin], rpars);
                
            //see if we have to limit the step
            if(number[ibin] > 0 && fabs(dadt[ibin]) > 0) {
                dadt_max_cell = dadt_lim * fabs((abin_e[iabin+1]-abin_e[iabin])/dadt[ibin]);
                if(dadt_max_cell < dt){
                    if(dt_step/dadt_max_cell > dust_maxSubSteps){
                        if(dadt[ibin]*dt > abin_c[iabin]){
                            number[ibin] = 0.0;
                            slope[ibin]  = 0.0;
                        } else {
                            dt = dt_step/dust_maxSubSteps; 
                        }
                    } else {
                        dt = dadt_max_cell;
                    }
                }
            }   
            // Debug total mass and number before step
            Ntot_old += number[ibin];
            Mtot_old += getMass(number[ibin], slope[ibin],iabin); 
        }
        
        // update the growth rate for the new stepsize
        ierr = set_dadt(rpars, dt);
        
        if(dt  < 0){
            return -1;
        }
        Ntot_new = 0;
        Ntot_mid = 0;
        Mtot_new = 0;
        // now loop over bins again and determine change between bins
        for(ibin = 0; ibin < dust_nbins; ibin++){
            // Graphite can only go into graphite bins
            // TODO: could be optimized with the dadt_lim limiter, since this restricts
            //       contributing bins to only a couple neighbours.
            if(ibin < isilicone){
                rangeStart = 0;
                rangeEnd   = isilicone;
                graphite = 1;
            } else {
                rangeStart = isilicone;
                rangeEnd   = dust_nbins;
                graphite = 0;
            }
            iabin = ibin - isilicone*(1-graphite);
            Nsum = 0;
            Msum = 0;
            for(ibin2 = rangeStart; ibin2 < rangeEnd; ibin2++){
                iabin2 = ibin2 - isilicone*(1-graphite);
                // no grains here, lets not risk grabbing any non-existent ones 
                if(number[ibin2] <= 0){
                    continue;
                }
                // See which range within ibin2 (if any) ends up in ibin1
                // left
                x1 = fmax(abin_e[iabin2]  ,abin_e[iabin]   - dadt[ibin2]*dt);
                // right
                x2 = fmax(fmin(abin_e[iabin2+1],abin_e[iabin+1] - dadt[ibin2]*dt),0.0);
                // If x2 > x1 then part of ibin2 will end up in ibin1
                if( x1 >= x2){
                    continue;
                }
                // Total number of grains and mass from ibin2 that ends up in ibin 
                intN = intNj(ibin2, iabin2, x2    ) - intNj(ibin2, iabin2, x1);
                // Multiply by 4pi/3 rho later. Currently unneccesary
                intM = intMj(ibin2, iabin2, x2, dt) - intMj(ibin2, iabin2, x1, dt);
                if(intN > 0 && intM >0){ 
                    Nsum = Nsum + intN; 
                    Msum = Msum + intM;
                }
            }
            
            if(Nsum < 0 || Msum < 0){
                Nnew[ibin] = 0;
                Snew[ibin] = 0;
                Mnew[ibin] = 0;
            } else {

                Sest = getSlope(Nsum, Msum, iabin);
                ierr = limitSlope(&Nnew[ibin], &Snew[ibin], Nsum, Sest, Msum, iabin); 
                Mnew[ibin] = Msum;

                Ntot_new += Nnew[ibin];
                Ntot_mid += Nsum;
                Mtot_new += Msum;
            }

        }
        ierr = rebinn(dt);
        if(ierr < 0){
            return -1;
        }
        // set arrays and recalculate number of grains in cell TODO: figure out reasonable "fail" limit and redo step at smaller stepsize if met
        NumTot_cell_s = 0;
        NumTot_cell_g = 0;
        for(ibin = 0; ibin < dust_nbins; ibin++){
            if(Nnew[ibin] < 0){
                printf("%d %.4e %.4e %.4e %.4e \n", ibin, Nnew[ibin], number[ibin], Snew[ibin], slope[ibin]);
            }
            number[ibin] = Nnew[ibin];
            slope[ibin]  = Snew[ibin];
            if(ibin < isilicone){
                NumTot_cell_g += number[ibin];
            } else {
                NumTot_cell_s += number[ibin];
            }
        }

#ifdef sputteringUseVelocities
        // update the velocities for the next timestep
        for(idx = 0; idx < NdustBins; idx++){
            ibin = globalToLocalIndex(idx);
            if(ibin < isilicone){
                graphite = 1;
            } else {
                graphite = 0;
            }
            iabin = ibin - isilicone * (1-graphite);
            ierr = update_vdust(ibin, iabin, graphite, rpars, dt);
        }
#endif

        // update time
        time_dust = time_dust + dt;
        // set timestep to attempt to do rest in one go
        dt = dt_step - time_dust;
    }
    return 1;
}



///////////////////////////////////////
//
//   Methods to communicate data between dust arrays and cells
//
///////////////////////////////////////
// Method to transfer cell data to smaller cells
int setBinsCell(int icell, double *Mtot){
    int idx, ibin, iabin, graphite;
    double norm, mass;
    
    *Mtot = 0;
    for( idx = 0; idx < NdustBins; idx++){
        ibin = globalToLocalIndex(idx);
        if(ibin < isilicone){
            graphite = 1;
        } else {
            graphite = 0;
        }

        // graphite and silicone have the same size bins. Use half array 
        iabin = ibin - isilicone*(1-graphite);    


        // Total number in bin stored as mass density 
        mass = ustate[icell*nvar + IDUST_START + NdustVar*idx];
        // current slope assuming piecewise linear dnumda
        slope[ibin]  = ustate[icell*nvar + IDUST_START + NdustVar*idx + 1];


        // convert to number via total mass in cell
        
        norm = 4*M_PI/3.;
        if(graphite){
            norm = norm * rho_c;
        } else {
            norm = norm * rho_s;
        }
    
        *Mtot = *Mtot + mass;
        
        number[ibin] = getNumber(slope[ibin], mass/norm, iabin);
    }
    return 1;
}

int getBinsCell(int icell, double *Mtot){
    int idx, ibin, iabin, graphite;
    double norm, mass;
    
    *Mtot = 0;

    for( idx = 0; idx < NdustBins; idx++){
        ibin = globalToLocalIndex(idx);
        if(ibin < isilicone){
            graphite = 1;
        } else {
            graphite = 0;
        }

        // graphite and silicone have the same size bins. Use half array 
        iabin = ibin - isilicone*(1-graphite);    

        // Conversion to mass density
        norm = 4*M_PI/3.;
        if(graphite){
            norm = norm * rho_c;
        } else {
            norm = norm * rho_s;
        }
        // Mass (density) in bin
        mass = getMass(number[ibin], slope[ibin], iabin);
        *Mtot = *Mtot + mass*norm;
        ustate[icell*nvar + IDUST_START + NdustVar*idx    ] = mass * norm;
        ustate[icell*nvar + IDUST_START + NdustVar*idx + 1] = slope[ibin];
    }
    return 1;
}

#endif
