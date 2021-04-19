#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "dust.h"
#include "hydro.h"
#include "radchem.h"
#include "moshpit.h"
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
//  Methods to calculate dust absorption and optical depth
//
//////////////////////////////////////////

int getDustOpticalDepth(){
    return 1;
}


//////////////////////////////////////////
//
//  Methods for evolution of dust distribution
//
//////////////////////////////////////////



//  Gas accretion and sputtering
//  sputtering : Tsai & Mathews 1995
double dadt_gas(double *pars){
//    double gas_density     = pars[0];
//    double H_density       = pars[1];
//    double gas_temperature = pars[2];
//    double Tsput           = pars[3];
//
//    double dadt = 0;
//
//    // Growth from accretion. 
//    dadt += (H_density/(1e-3))*sqrt(gas_temperature/10)*1e-4/(1e9*years); //micrometer per Gyr
//    // Destruction from sputtering
//    dadt -= (3.2e-18)*(gas_density/mass_p)/(pow(Tsput/gas_temperature,2.5)+1);
    return 0.;
}
// Sublimation rate from GUHATHAKURTA AND DRAINE 1989
// Assumes of the form
//  
//  dm/dt = -4*pi*a^2*m_part*prefact * exp (-(B-2e4*N^(-1/3)/kT_dust)
//  da/dt = -(m_part/rho)*prefact* exp
double dadt_sublimation(double ai, double temp_i, int graphite){
//    double Ni, prefact, B;
//    if(graphite == 1){
//        prefact = 2e14*m_carbon/rho_c;
//        Ni      = rho_c*4*M_PI*pow(ai,3)/3/m_carbon;
//        B       = 81200-20000*pow(Ni,-1./3.);
//    } else {
//        prefact = 2e15*m_silicate/rho_s;
//        Ni      = rho_s*4*M_PI*pow(ai,3)/3/m_silicon;
//        B       = 68100-20000*pow(Ni,-1./3.);
//    }
//
//    double dadt = -prefact*pow(ai,2)*exp(-B/temp_i);
    return 0.;
}



double get_dadt(double ai, double *rpars){
    if(dadt_mode == 1){
        return dadt_c;
    } else if(dadt_mode == 2) {
        return dadt_c * ai;
    } else if(dadt_mode == 3) {
        return dadt_c * sin(2*M_PI * time/dadt_tscale);
    }
    double dadt = 0;
    return dadt;
}

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
    return (Mj - Sj*SfactM[iabin])/NfactM[iabin];

}

// use Nj and Sj to solve for Mass
double getMass(double Nj, double Sj,int iabin){
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
    if(Nj <= 0 ){
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
        Mn  = Mnew[dust_nbins-2];
        Mnp = Mnew[dust_nbins-1];

        Nn  = Nnew[dust_nbins-2];
        Nnp = Mnp/pow(abin_e[Nabins-1], 3);
        
        Sn  = Snew[dust_nbins-2];
        
        if(Nnp > 0){ 
            // get updated bin values
            ierr = rebinn_upper(Nn, Nnp, Sn, Mn, Mnp, &Nntilde, &Sntilde);

            // Update the bins (and limit slope)
            Mnew[dust_nbins-2] = Mn + Mnp;        
            ierr = limitSlope(&Nnew[dust_nbins-2], &Snew[dust_nbins-2], Nntilde, Sntilde, Mn+Mnp, Nabins-2); 
            Mnew[dust_nbins-1] = 0;      
            Nnew[dust_nbins-1] = 0;      
            Snew[dust_nbins-1] = 0;     
        } 
        
        Mn  = Mnew[isilicone+1];
        Mnm = Mnew[isilicone];

        Nn  = Nnew[isilicone+1];
        Nnm = Mnm/pow(abin_e[1], 3);
        
        Sn  = Snew[isilicone+1];
        
        if(Nnm > 0){ 
            // get updated bin values
            ierr = rebinn_lower(Nn, Nnm, Sn, Mn, Mnm, &Nntilde, &Sntilde);

            // Update the bins (and limit slope)
            Mnew[isilicone+1] = Mn + Mnm;        
            ierr = limitSlope(&Nnew[isilicone+1], &Snew[isilicone+1], Nntilde, Sntilde, Mn+Mnm, 1); 
            Mnew[isilicone] = 0;      
            Nnew[isilicone] = 0;      
            Snew[isilicone] = 0;     
        } 

    }
    // if we have graphite grains
    if( isilicone > 0){
        // not same for carbonious grains
        Mn  = Mnew[isilicone-2];
        Mnp = Mnew[isilicone-1];

        Nn  = Nnew[isilicone-2];
        Nnp = Mnp/pow(abin_e[Nabins-1],3);
        Sn  = Snew[isilicone-2];
        // Do we have grains in ghost bin? if so fix
        if(Nnp > 0){
            // get updated bin values
            ierr = rebinn_upper(Nn, Nnp, Sn, Mn, Mnp, &Nntilde, &Sntilde);
            
            // Update the bins
            Mnew[isilicone-2] = Mn + Mnp;
            ierr = limitSlope(&Nnew[isilicone-2], &Snew[isilicone-2], Nntilde, Sntilde, Mn+Mnp, Nabins-2); 
            Mnew[isilicone-1] = 0;      
            Nnew[isilicone-1] = 0;      
            Snew[isilicone-1] = 0;  
        }    
        
        Mn  = Mnew[1];
        Mnm = Mnew[0];

        Nn  = Nnew[1];
        Nnm = Mnm/pow(abin_e[1], 3);
        
        Sn  = Snew[1];
        
        if(Nnm > 0){ 
            // get updated bin values
            ierr = rebinn_lower(Nn, Nnm, Sn, Mn, Mnm, &Nntilde, &Sntilde);
            // Update the bins (and limit slope)
            Mnew[1] = Mn + Mnm;        
            ierr = limitSlope(&Nnew[1], &Snew[1], Nntilde, Sntilde, Mn+Mnm, 1); 
            Mnew[0] = 0;      
            Nnew[0] = 0;      
            Snew[0] = 0;     
        } 
    }
    return 1; 
}

int dustCell(double *rpars, int *ipars, double dt_step){
    int ibin, ibin2, iabin, iabin2, ierr;
    int rangeStart, rangeEnd, graphite;
    double x1, x2;
    double dt, time_dust;
    double Nsum, Msum, Sest;
    double Mtot_old, Ntot_old, Mtot_new, Ntot_new, Ntot_mid;
    // Start by calculating change in size due to sublimation/sputtering/accretion
    //initialize by trying to take one full step
    time_dust = 0;
    // try to take all in one step
    dt = dt_step;
    
    while(time_dust < dt_step){
        // Calculate change in each bin and limit timestep in needed
        Mtot_old = 0;
        Ntot_old = 0;
        for(ibin = 0; ibin < dust_nbins; ibin++){
            if(ibin < isilicone){
                graphite = 1;
            } else {
                graphite = 0;
            }
            // graphite and silicone have the same size bins. Use half array 
            iabin = ibin - isilicone*(1-graphite);    

            // Debug total mass and number before step
            Ntot_old += number[ibin];
            Mtot_old += getMass(number[ibin], slope[ibin],iabin); 

            // Temperature of dust grains in this bin
            //temp         = getDustTemp(abin_c[iabin], rpars);
                
            dadt[ibin]   = get_dadt(abin_c[iabin], rpars);
            //see if we have to limit the step
            if(number[ibin] > 0 && dadt[ibin] > 0) {
                dt = fmin(dt, dadt_lim * fabs((abin_e[iabin+1]-abin_e[iabin])/dadt[ibin]));
            }   
        }
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
                x2 = fmin(abin_e[iabin2+1],abin_e[iabin+1] - dadt[ibin2]*dt);
                // If x2 > x1 then part of ibin2 will end up in ibin1
                if( x1 > x2){
                    continue;
                }
                // Total number of grains and mass from ibin2 that ends up in ibin 
                Nsum = Nsum + intNj(ibin2, iabin2, x2    ) - intNj(ibin2, iabin2, x1);
                // Multiply by 4pi/3 rho later. Currently unneccesary
                Msum = Msum + intMj(ibin2, iabin2, x2, dt) - intMj(ibin2, iabin2, x1, dt);
            } 
            Sest = getSlope(Nsum, Msum, iabin);
            ierr = limitSlope(&Nnew[ibin], &Snew[ibin], Nsum, Sest, Msum, iabin); 
            Mnew[ibin] = Msum;

            Ntot_new += Nnew[ibin];
            Ntot_mid += Nsum;
            Mtot_new += Msum;

        }
        ierr = rebinn(dt);
        if(ierr < 0){
            return -1;
        }
        // set arrays TODO: figure out reasonable "fail" limit and redo step at smaller stepsize if met
        for(ibin = 0; ibin < dust_nbins; ibin++){
            number[ibin] = Nnew[ibin];
            slope[ibin]  = Snew[ibin];
        }
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
    printf("initial dust density =  %.4e\n", *Mtot);
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
    printf("final dust density =  %.4e\n", *Mtot);
    return 1;

}


