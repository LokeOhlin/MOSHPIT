#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dust.h>
#include <dustRadiation.h>
#include <hydro.h>
#include <radchem.h>
#include <moshpit.h>
#ifdef useDust

int set_dadt(double *rpars){
    int idx, ibin, iabin, graphite;
    int idt; 

    // start by setting everything to constant (default zero)
    for(ibin = 0; ibin < dust_nbins; ibin++){
        if(ibin < isilicone){
            graphite = 1;
            iabin = ibin;
        } else {
            graphite = 0;
            iabin = ibin - isilicone;
        }

        if(dadt_mode == 1){
            dadt[ibin] = dadt_c;   
        } else if(dadt_mode == 2) {
            dadt[ibin] = dadt_c * abin_c[iabin];   
        } else if(dadt_mode == 3) {
            dadt[ibin] = dadt_c * sin(2*M_PI * time/dadt_tscale);
        } else {
            dadt[ibin] = 0;
        }
    }

    // Start with sputtering
    if(dust_useSputtering){
        double numd = rpars[0];
        double Tgas = rpars[1];
        // only over physical range
        for(idx = 0; idx < NdustBins; idx++){
            ibin = globalToLocalIndex(idx);
            if(ibin < isilicone){
                graphite = 1;
                iabin = ibin;
            } else {
                graphite = 0;
                iabin = ibin - isilicone;
            }

            dadt[ibin] += get_dadt_sputtering(abin_c[iabin], ibin, iabin, graphite, numd, Tgas, &idt);
        }
    }

    // Next do sublimation (currently only possible with radiation)
    if(dust_useSublimation && dust_useRadiation){
        double dadt_sub;
        for(idx = 0; idx < NdustBins; idx ++){
            ibin = globalToLocalIndex(idx);
            if(ibin < isilicone){
                graphite = 1;
                iabin = ibin;
            } else {
                graphite = 0;
                iabin = ibin - isilicone;
            }
            dadt_sub = get_dadt_sublimation(ibin, iabin, graphite);
            dadt[ibin] += dadt_sub;
            // if sublimation for this bin is lower than 10% of total, or less than a minimum, 
            //printf("%d %.4e %.4e %.4e \n", ibin, abin_c[iabin], dadt[ibin], dadt_sub);
            // we no longer have to calculate for the larger grains.
            if(fabs(dadt_sub) < fabs(0.1*dadt[ibin]) || fabs(dadt_sub) < abin_c[iabin]*dadt_min){
                // if we are looping over graphite grains, jump to silicates
                // otherwise break.
                if(graphite){
                    idx = localToGlobalIndex(isilicone)-2;
                } else {
                    break;
                }
            }
        }
    }
    return 1;
}





#endif
