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
    double ac  = abin_c[iabin];
    double aep = abin_e[iabin+1];
    double ae  = abin_e[iabin];
   
    
    double NfactP = pow(aep,4)/(4*(aep-ae));
    double Nfact  = pow(ae ,4)/(4*(aep-ae));

    double SfactP  = pow(aep,4)*(aep/5. - ac/4.);
    double Sfact   = pow(ae ,4)*(ae /5. - ac/4.);

    return (Mj - Nj*(NfactP - Nfact))/(SfactP -Sfact);

}
// use Mj and Sj to solve for Number
double getNumber(double Sj, double Mj,int iabin){
    double ac  = abin_c[iabin];
    double aep = abin_e[iabin+1];
    double ae  = abin_e[iabin];
   
    
    double NfactP = pow(aep,4)/(4*(aep-ae));
    double Nfact  = pow(ae ,4)/(4*(aep-ae));

    double SfactP  = pow(aep,4)*(aep/5. - ac/4.);
    double Sfact   = pow(ae ,4)*(ae /5. - ac/4.);

    return (Mj - Sj*(SfactP - Sfact))/(NfactP - Nfact);

}

// use Nj and Sj to solve for Mass
double getMass(double Nj, double Sj,int iabin){
    double ac  = abin_c[iabin];
    double aep = abin_e[iabin+1];
    double ae  = abin_e[iabin];
   
    
    double NfactP = pow(aep,4)/(4*(aep-ae));
    double Nfact  = pow(ae ,4)/(4*(aep-ae));

    double SfactP  = pow(aep,4)*(aep/5. - ac/4.);
    double Sfact   = pow(ae ,4)*(ae /5. - ac/4.);

    return Nj*(NfactP-Nfact) + Sj*(SfactP - Sfact);
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
//            printf("Warning : both left and right dust bin edges have dnda < 0. iabin =  %d \n", iabin);
            return 1;
        }
        // Do nothing 
        *Njnew = Nj;
        *Sjnew = Sj;
        return 1;
    }

    double Nfact_tot;
        
    double Nfact  = (pow(aep,4) - pow(ae ,4))/(4*(aep-ae));
    
    //double SfactP  = pow(aep,5)/5. - ac*pow(aep,4)/4.;
    double SfactP  = pow(aep,4)*(aep/5. - ac/4.);
    double Sfact   = pow(ae,4) *(ae /5. - ac/4.);
    //double Sfact   = pow(ae ,5)/5. - ac*pow(ae ,4)/4.;
    if(dndap < 0){
        Nfact_tot = (Nfact-(SfactP-Sfact)/dap/dac);
        *Njnew = Mj/Nfact_tot;
        *Sjnew = -*Njnew/dac/dap;
    } else if(dnda < 0){
        Nfact_tot = (Nfact-(SfactP-Sfact)/dam/dac);
        *Njnew = Mj/Nfact_tot;
        *Sjnew = -*Njnew/dac/dam;
    }
    
//    Nj[0] = Mj/();
//    Sj[0] = -Nj*div*div;
    return 1;    
}

// Rebinn such that the largest bin (amax,inf) contain zero stuff
// we keep the lowest bin size for the moment as a "resevoir" with zero opacity
// Assume all grains of size > amax have the size of amax. Update the average size of the grains in the bin
// And determine numbers and slope to satisfy both the new average size and new mass
int rebinn(double dt){
    double Mn, Mnp, Nn, Nnp, Nntilde, Sn, Sntilde;
    double ac  = abin_c[Nabins-1];
    double aep = abin_e[Nabins];
    double ae  = abin_e[Nabins-1];
    double average_a, average_a_tilde;

    double dac = aep - ae;
    // Factor in fron of N and S in mass equation
    double NfactM   = (pow(aep,4) - pow(ae ,4))/(4*(aep-ae));
    
    double SfactMP  = pow(aep,4)*(aep/5. - ac/4.);
    double SfactM   = pow(ae ,4)*(ae /5. - ac/4.);

    // Factors in front of N and S in average size equation
    double SfactAP   = pow(aep,2)*(aep/3. - ac/2.);
    double SfactA    = pow(ae ,2)*(ae /3. - ac/2.);
    double NfactA    = 0.5*(aep*aep-ae*ae)/dac; 


    // If we have silicates
    if(isilicone < NdustBins){
        Mn  = Mnew[NdustBins-2];
        Mnp = Mnew[NdustBins-1];
        Nn  = Nnew[NdustBins-2];
        Nnp = Nnew[NdustBins-1]; 
        Sn  = Nnew[NdustBins-2];
        if(Nnp > 0){ 
            // Initial average
            average_a = NfactA + (SfactAP - SfactA)*Sn/Nn; 
            // Updated average
            average_a_tilde = (Nn*average_a + Nnp*aep)/(Nn+Nnp); 

            // calculate Ntilde based on new mass
            Nntilde = (Mn+Mnp)/(NfactM + (average_a_tilde - NfactA)*(SfactMP-SfactM)/(SfactAP - SfactA));
            // finally solve for Stilde
            Sntilde = Nntilde * (average_a_tilde - NfactA)/(SfactAP-SfactA);  
            
            // Update the bins
            Mnew[NdustBins-2] = Mn + Mnp;        
            Nnew[NdustBins-2] = Nntilde;        
            Snew[NdustBins-2] = Sntilde; 
            
            Mnew[NdustBins-1] = 0;      
            Nnew[NdustBins-1] = 0;      
            Snew[NdustBins-1] = 0;     
        } 
    }
    // if we have graphite grains
    if( isilicone > 0){
        // not same for carbonious grains
        Mn  = Mnew[isilicone-2];
        Mnp = Mnew[isilicone-1];
        Nn  = Nnew[isilicone-2];
        Nnp = Nnew[isilicone-1]; 
        Sn  = Snew[isilicone-2];
        // Do we have grains in ghost bin? if so fix
        if(Nnp > 0){
            // Initial average
            average_a = NfactA + (SfactAP - SfactA)*Sn/Nn; 
            // Updated average
            average_a_tilde = (Nn*average_a + Nnp*aep)/(Nn+Nnp); 

            // calculate Ntilde based on new mass
            Nntilde = (Mn+Mnp)/(NfactM + (average_a_tilde - NfactA)*(SfactMP-SfactM)/(SfactAP - SfactA));
            // finally solve for Stilde
            Sntilde = Nntilde * (average_a_tilde - NfactA)/(SfactAP-SfactA);  
            // Update the bins
            Mnew[isilicone-2] = Mn + Mnp;        
            Nnew[isilicone-2] = Nntilde;        
            Snew[isilicone-2] = Sntilde; 
            
            Mnew[isilicone-1] = 0;      
            Nnew[isilicone-1] = 0;      
            Snew[isilicone-1] = 0;  
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
        for(ibin = 0; ibin < NdustBins; ibin++){
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
            if(ibin == 0) {
                printf("dadt = %.4e\n", dadt[ibin]);
            }    
            //see if we have to limit the step
            if(number[ibin] > 0) {
                dt = fmin(dt, dadt_lim * fabs((abin_e[iabin+1]-abin_e[iabin])/dadt[iabin]));
            }   
        }
        if(dt  < 0){
            return -1;
        }
        Ntot_new = 0;
        Ntot_mid = 0;
        Mtot_new = 0;
        // now loop over bins again and determine change between bins
        for(ibin = 0; ibin < NdustBins; ibin++){
            // Graphite can only go into graphite bins
            // TODO: could be optimized with the dadt_lim limiter, since this restricts
            //       contributing bins to only a couple neighbours.
            if(ibin < isilicone){
                rangeStart = 0;
                rangeEnd   = isilicone;
                graphite = 1;
            } else {
                rangeStart = isilicone;
                rangeEnd   = NdustBins;
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
        //printf("N: %.4e -> %.4e %.4e %.4e| M: %.4e -> %.4e %.4e \n", Ntot_old, Ntot_mid,Ntot_new,(Ntot_new-Ntot_old)/Ntot_old, Mtot_old, Mtot_new, (Mtot_new-Mtot_old)/Mtot_old);
        ierr = rebinn(dt);
        if(ierr < 0){
            return -1;
        }
        // set arrays TODO: figure out reasonable "fail" limit and redo step at smaller stepsize if met

        for(ibin = 0; ibin < NdustBins; ibin++){
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
    int ibin, iabin, graphite;
    double norm, mass;
    
    *Mtot = 0;
    for( ibin = 0; ibin < NdustBins; ibin++){
        if(ibin < isilicone){
            graphite = 1;
        } else {
            graphite = 0;
        }

        // graphite and silicone have the same size bins. Use half array 
        iabin = ibin - isilicone*(1-graphite);    


        // Total number in bin stored as mass density 
        mass = ustate[icell*nvar + IDUST_START + NdustVar*ibin];
        // current slope assuming piecewise linear dnumda
        slope[ibin]  = ustate[icell*nvar + IDUST_START + NdustVar*ibin + 1];


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
    int ibin, iabin, graphite;
    double norm, mass;
    
    *Mtot = 0;

    for( ibin = 0; ibin < NdustBins; ibin++){
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
        ustate[icell*nvar + IDUST_START + NdustVar*ibin    ] = mass * norm;
        ustate[icell*nvar + IDUST_START + NdustVar*ibin + 1] = slope[ibin];
    }
    printf("final dust density =  %.4e\n", *Mtot);
    return 1;

}


