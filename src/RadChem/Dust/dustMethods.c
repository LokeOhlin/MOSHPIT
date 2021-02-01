#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>


//////////////////////////////////////////
//
//  Methods to calculate dust absorption and optical depth
//
//////////////////////////////////////////

int getDustOpticalDepth(){
}


//////////////////////////////////////////
//
//  Methods for evolution of dust distribution
//
//////////////////////////////////////////



//  Gas accretion and sputtering
//  sputtering : Tsai & Mathews 1995
double dadt_gas(double *pars){
    double gas_density     = pars[0];
    double H_density       = pars[1];
    double gas_temperature = pars[2];
    double Tsput           = pars[3];

    double dadt = 0;

    // Growth from accretion. 
    dadt += (H_density/(1e-3))*sqrt(gas_temperature/10)*1e-4/(1e9*years); //micrometer per Gyr
    // Destruction from sputtering
    dadt -= (3.2e-18)*(gas_density/mass_p)/(pow(Tsput/gas_temperature,2.5)+1);
    return dadt
}
// Sublimation rate from GUHATHAKURTA AND DRAINE 1989
// Assumes of the form
//  
//  dm/dt = -4*pi*a^2*m_part*prefact * exp (-(B-2e4*N^(-1/3)/kT_dust)
//  da/dt = -(m_part/rho)*prefact* exp
double dadt_sublimation(double ai, double temp_i, int graphite){
    double Ni, prefact, B;
    if(graphite == 1){
        prefact = 2e14*m_carbon/rho_c;
        Ni      = rho_c*4*M_PI*pow(ai,3)/3/m_carbon;
        B       = 81200-20000*pow(Ni,-1./3.);
    } else {
        prefact = 2e15*m_silicate/rho_s;
        Ni      = rho_s*4*M_PI*pow(ai,3)/3/m_silicon;
        B       = 68100-20000*pow(Ni,-1./3.);
    }

    double dadt = -prefact*pow(ai,2)*exp(-B/temp_i);
    return dadt;
}



// GRAIN DISTRIBUTION EVOLUTION
// GENERAL METHOD FROM McKinnon et al. 2018

// integral for Nj
double intNj(int ibin, int iabin, double a){
    return  number[ibin]*a/(abin_e[iabin+1] - abin_e[iabin]) + slope[ibin]*(a*a/2 - abin_c[iabin]*a);
}

// integral for Mj
double intMj(int ibin, int iabin, double a, dt){
    double da  = dadt_i[iabin]*dt;
    double ac  = abin_c[iabin];
    double aep = abin_e[iabin+1];
    double ae  = abin_e[iabin];

    // integral (a+da)**3*(a-ac)
    double slope_fact = pow(a,5)/5. + (3*da-ac)*pow(a,4)/4.
    slope_fact += da(da-ac)*pow(a,3) + da*da*(da-3*ac)*a*a/2. +da*da*da*dc*a

    return number[ibin]*pow(a+da.4)/(aep-ae) + slope[ibin]*slope_fact;
}

// use Mj and Nj to solve for slope
double getSlope(double Nj, double Mj,int iabin){
    double ac  = abin_c[iabin];
    double aep = abin_e[iabin+1];
    double ae  = abin_e[iabin];
   
    
    double NfactP = pow(aep,4)/(4*(aep-ae));
    double Nfact  = pow(ae ,4)/(4*(aep-ae));

    double SfactP  = pow(aep,5)/5. - ac*pow(aep,4)/4.;
    double Sfact   = pow(ae ,5)/5. - ac*pow(ae ,4)/4.;

    return (Mj - Nj*(NfactP - Nfact))/(SfactP -Sfact);

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

    // Both positive or both negative
    if( dndap*dnda > 0){
        //both negative throw warning and set to zero
        if(dndap < 0){
            Njnew[0] = 0;
            Sjnew[0] = 0;
            printf("Warning : both left and right dust bin edges have dnda < 0. iabin =  %d \n", iabin);
            return 1;
        }
        // Do nothing 
        Njnew[0] = Nj;
        Sjnew[0] = Sj;
        return 1;
    }

    double Nfact_tot;
        
    double Nfact  = (pow(aep,4) - pow(ae ,4))/(4*(aep-ae));
    
    double SfactP  = pow(aep,5)/5. - ac*pow(aep,4)/4.;
    double Sfact   = pow(ae ,5)/5. - ac*pow(ae ,4)/4.;
    if(dndap < 0){
        Nfact_tot = (0.25*Nfact-(SfactP-Sfact)/dap)/dac;
        Njnew[0] = Mj/Nfact_tot;
        Sjnew[0] = -Njnew[0]/dac/dap;
    } else if(dnda < 0){
        Nfact_tot = (0.25*Nfact-(SfactP-Sfact)/dam)/dac;
        Njnew[0] = Mj/Nfact_tot;
        Sjnew[0] = -Njnew[0]/dac/dam;
    }
    
    Nj[0] = Mj/();
    Sj[0] = -Nj*div*div;
    return 1;    
}

// Rebinn such that the largest bin (amax,inf) contain zero stuff
// we keep the lowest bin size for the moment as a "resevoir" with zero opacity
// Assume all grains of size > amax have the size of amax. Update the average size of the grains in the bin
// And determine numbers and slope to satisfy both the new average size and new mass
int rebinn(dt){
    double Mn, Mnp, Nn, Nnp, Nntilde, Sn, Sntilde;
    double ac  = abin_c[NdustBins-2];
    double aep = abin_e[NdustBins-1];
    double ae  = abin_e[NdustBins-2];
    double average_a, average_a_tilde;

    double dac = aep - ae;
    double dap = aep - ac;
    double dam = ae  - ac;
    
    double Nfact_tot;
    // Factor in fron of N and S in mass equation
    double NfactM   = (pow(aep,4) - pow(ae ,4))/(4*(aep-ae));
    
    double SfactMP  = pow(aep,5)/5. - ac*pow(aep,4)/4.;
    double SfactM   = pow(ae ,5)/5. - ac*pow(ae ,4)/4.;

    // Factors in front of N and S in average size equation
    double SfactAP   = pow(aep,3)/3. - ac*pow(aep,2)/2.;
    double SfactA    = pow(ae ,3)/3. - ac*pow(ae ,2)/2.;
    double NfactA    = 0.5*(aep*aep-ae*ae)/dac; 


    // If we have silicates
    if(isilicone < NdustBins){
        Mn  = Mnew[NdustBins-2];
        Mnp = Mnew[NdustBins-1];
        Nn  = Nnew[NdustBins-2];
        Nnp = Nnew[NdustBins-1]; 
        Sn  = Nnew[NdustBins-2];
        
        // Initial average
        average_a = NfactA + (SfactAP - SfactA)*Sn/Nn; 
        // Updated average
        average_a_tilde = (Nn*averagea + Nnp*aep)/(Nn+Nnp); 

        // calculate Ntilde based on new mass
        Nntilde = (Mn+Mnp)/(NfactM + (average_a_tilde - NfactA)*(SfactMP-SfactM)/(SfactAP - SfactA));
        // finally solve for Stilde
        Sntilde = Nntilde * (average_a_tilde - NfactA)/(SfactAp-SfactA);  
        
        // Update the bins
        Mnew[NdustBins-2] = Mn + Mnp;        
        Nnew[NdustBins-2] = Nntilde;        
        Snew[NdustBins-2] = Sntilde; 
        
        Mnew[NdustBins-1] = 0;      
        Nnew[NdustBins-1] = 0;      
        Snew[NdustBins-1] = 0;      
    }

    // not same for carbonious grains
    Mn  = Mnew[isilicone-2];
    Mnp = Mnew[isilicone-1];
    Nn  = Nnew[isilicone-2];
    Nnp = Nnew[isilicone-1]; 
    Sn  = Nnew[isilicone-2];
    
    // Initial average
    average_a = NfactA + (SfactAP - SfactA)*Sn/Nn; 
    // Updated average
    average_a_tilde = (Nn*averagea + Nnp*aep)/(Nn+Nnp); 

    // calculate Ntilde based on new mass
    Nntilde = (Mn+Mnp)/(NfactM + (average_a_tilde - NfactA)*(SfactMP-SfactM)/(SfactAP - SfactA));
    // finally solve for Stilde
    Sntilde = Nntilde * (average_a_tilde - NfactA)/(SfactAp-SfactA);  
    
    // Update the bins
    Mnew[isilicone-2] = Mn + Mnp;        
    Nnew[isilicone-2] = Nntilde;        
    Snew[isilicone-2] = Sntilde; 
    
    Mnew[isilicone-1] = 0;      
    Nnew[isilicone-1] = 0;      
    Snew[isilicone-1] = 0;      
    return 1; 
}


int dustCell(int icell, double rpars, int ipars){
    int ibin, ibin2, iabin;
    // Start by unpacking dust data and calculate change in size due to sublimation/sputtering/accretion
    
    for( ibin = 0; ibin < NdustBins; ibin++){
        if(ibin < isilicone){
            graphite = 1;
        } else {
            graphite = 0;
        }

        // graphite and silicone have the same size bins. Use half array 
        iabin = ibin -isilicone*(1-graphite);    


        // Total number in bin
        number[ibin] = ustate[icell + IDUST_BEGIN + ndustVar*ibin];
        // current slope assuming piecewise linear dnumda
        slope[ibin]  = ustate[icell + IDUST_BEGIN + ndustVar*ibin + 1];
        // Temperature of dust grains in this bin
        temp         = getDustTemp(abin_c[iabin],rpar);
            
        dadt_i[ibin]   = dadt_gas(rpars);
        dadt_i[ibin]  += dadt_sublimation(abin_c[iabin], graphite)
    } 
    
    // now loop over bins again and determine change between bins
    for(ibin = 0; ibin < NdustBins; ibin++){
        // Graphite can only go into graphite bins
        if(ibin < isilicone){
            rangeStart = 0;
            rangeEnd   = isilicione;
            graphite = 1;
        } else {
            rangeStart = isilicone;
            rangeENd   = NdustBins;
            graphite = 0;
        }
        iabin = ibin - isilicone(1-graphite);
        Nsum = 0;
        Msum = 0;
        for(ibin2 = rangeStart; ibin2 < rangeEnd; ibin2++){
            iabin2 = ibin2 - isilicone*(1-graphite);
            // See which range within ibin2 (if any) ends up in ibin1
            // left
            x1 = fmax(abin_e[iabin2],abin_e[iabin1]-dadt_i[ibin2]);
            // right
            x2 = fmin(abin_e[iabin2+1],abin_e[iabin1+1]-dadt_i[ibin2]);
            // If x2 > x1 then part of ibin2 will end up in ibin1
            if( x1 > x2){
                continue;
            }
            
            Nsum = Nsum + intNj(ibin2, iabin2, x2) - intNj(ibin2, iabin2, x1);
            // Multiply by 4pi/3 rho later. Currently unneccesary
            Msum = Msum + intMj(ibin2, iabin2, x2) - intMj(ibin2, iabin2, x1);

        } 

        Sest = getSlope(Nsum, Msum, iabin);
        ierr = limitSlope(&Nnew[ibin], &Snew[ibin], Nsum, Sest, Msun, iabin); 
        Mnew[ibin] = Msum;

    }

    ierr = rebinn(dt);
    return 1;
}



