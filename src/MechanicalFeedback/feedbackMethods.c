#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <mechanicalFeedback.h>
#include <cgeneral.h>
#include <hydro.h>
#include <constantsAndUnits.h>

int doFeedbackStep(double dt, double *dt_feedback){
    // start with wind
    //
    *dt_feedback = 1e99;
    if(fb_useWind) {
        double dens0     = ustate[NGHOST*nvar];
        double mom0      = ustate[NGHOST*nvar + 1];
        double etot0     = ustate[NGHOST*nvar + 2];

        double volcell = vol[NGHOST];
        double densWind = fb_MdotWind*dt/volcell;
        double dens = dens0 + densWind;
        
        double momWind  = densWind * fb_vWind;
        double mom  = mom0 + momWind;

        double Twind  = 1e4;
        double eintWind = 1.5 * densWind * boltzmann * Twind / (0.5 * mH);
        double etotWind  = 0.5*densWind * fb_vWind*fb_vWind + eintWind;
        double etot = etot0 + etotWind;

    
        ustate[NGHOST*nvar]     = dens;
        ustate[NGHOST*nvar + 1] = mom ;
        ustate[NGHOST*nvar + 2] = etot; 

        if(etotWind > etot0*fb_dEmaxWind){
            *dt_feedback = etot0*fb_dEmaxWind * dt / etotWind;
        }
    }
    // possibly other stuff 
    return 1;
}
