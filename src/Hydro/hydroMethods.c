#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "hydro.h"


int getCFL(double *cfl){
    int icell;
    double rho, rhoinv, pre, cs, vel, wspeed;
    
    cfl[0]=0; // initially inverted 
    for(icell = 2; icell < NCELLS-2; icell++){
        rho = ustate[icell*nvar];
        rhoinv = 1/rho;

        vel = ustate[icell*nvar+1]*rhoinv;
        pre = (ustate[icell*nvar+2]*rhoinv-0.5*vel*vel)*(adi-1);
        
        cs = sqrt(adi*pre);
        wspeed = cs + fabs(vel);
        cfl[0] = fmax(wspeed/dr[icell], cfl[0]);
    }
    cfl[0] = courant_number/cfl[0];
    return 1;
}

int setBoundary(){
    int ivar;
    // Left boundary
    for( ivar = 0; ivar < nvar; ivar++){
        if(left_bound == 0){ // reflective
            if(ivar == 1){
                ustate[ivar] = -ustate[3*nvar+ivar];
                ustate[nvar+ivar] = -ustate[2*nvar+ivar];
            }else{
                ustate[ivar] = ustate[3*nvar+ivar];
                ustate[nvar+ivar] = ustate[2*nvar+ivar];
            }
        } else if(left_bound == 1){ // outflow
            if(ivar == 1){
                ustate[ivar] = fmin(ustate[2*nvar+ivar],0);
                ustate[nvar+ivar] = fmin(ustate[2*nvar+ivar],0);
            }else{
                ustate[ivar] = ustate[2*nvar+ivar];
                ustate[nvar+ivar] = ustate[2*nvar+ivar];
            }
        } else if(left_bound == 2){ // fixed
            if(ivar == 0){
                ustate[ivar] = bdensL;
                ustate[nvar + ivar] = bdensL;
            } else if(ivar == 1){
                ustate[ivar] = bvelL;
                ustate[nvar + ivar] = bvelL;
            } else if(ivar == 2){
                ustate[ivar] = benerL;
                ustate[nvar + ivar] = benerL;
            } else { // all others variables set as in outflow
                ustate[ivar] = ustate[2*nvar+ivar];
                ustate[nvar + ivar] = ustate[2*nvar+ivar];
            }
        }
    }

    // Right boundary
    for (ivar = 0; ivar < nvar;ivar++){
        if(right_bound == 0){ // reflective
            if(ivar == 1){
                ustate[(NCELLS-1)*nvar + ivar] = -ustate[(NCELLS-4)*nvar +ivar];
                ustate[(NCELLS-2)*nvar + ivar] = -ustate[(NCELLS-3)*nvar +ivar];
            } else {
                ustate[(NCELLS-1)*nvar + ivar] = ustate[(NCELLS-4)*nvar +ivar];
                ustate[(NCELLS-2)*nvar + ivar] = ustate[(NCELLS-3)*nvar +ivar];
            }
        } else if(right_bound == 1) { // outflow
            if(ivar == 1){
                ustate[(NCELLS-1)*nvar + ivar] = fmax(ustate[(NCELLS-3)*nvar +ivar],0);
                ustate[(NCELLS-2)*nvar + ivar] = fmax(ustate[(NCELLS-3)*nvar +ivar],0);
            } else {
                ustate[(NCELLS-1)*nvar + ivar] = ustate[(NCELLS-3)*nvar +ivar];
                ustate[(NCELLS-2)*nvar + ivar] = ustate[(NCELLS-3)*nvar +ivar];
            }
        } else if(right_bound == 2) { //fixed
            if(ivar == 0){
                ustate[(NCELLS-1)*nvar+ivar] = bdensR;
                ustate[(NCELLS-2)*nvar+ivar] = bdensR;
            } else if(ivar == 1){
                ustate[(NCELLS-1)*nvar+ivar] = bvelR;
                ustate[(NCELLS-2)*nvar+ivar] = bvelR;
            } else if(ivar == 2){
                ustate[(NCELLS-1)*nvar+ivar] = benerR;
                ustate[(NCELLS-2)*nvar+ivar] = benerR;
            } else { // all others variables set as in outflow
                ustate[(NCELLS-1)*nvar + ivar] = ustate[(NCELLS-3)*nvar +ivar];
                ustate[(NCELLS-2)*nvar + ivar] = ustate[(NCELLS-3)*nvar +ivar];
            }
        }
    }
    return 1;
}
int toPrimitive(){
    int icell, idx, ivar;
    double rho, rhoinv, velx, ekin,eint;
    for(icell = 0; icell < NCELLS; icell++){
        //Index of cell in array
        idx=icell*nvar;
        // Density -> density
        rho = ustate[idx];
        rhoinv = 1/rho;
        pstate[idx] = rho;
        // Momentum -> velocity
        velx = ustate[idx + 1]*rhoinv;
        pstate[idx + 1] = velx;
        // Total energy -> internal energy -> Pressure
        // specific ekin
        ekin = 0.5*velx*velx;
        eint = ustate[idx + 2]*rhoinv - ekin;
        pstate[idx + 2]= rho * eint * (adi-1);
#ifdef useChemistry
        // copy over mass fractions
        for(ivar = ICHEM_START; ivar < ICHEM_END; ivar++){
            pstate[idx + ivar] = ustate[idx + ivar];
        }
#endif 
    }
    return 1;
}

int fixState(){
    int icell, idx, ivar;
    double val;
    for(icell = 0; icell < NCELLS; icell++){
        idx = icell *nvar;
#ifdef useChemistry
        for(ivar = ICHEM_START; ivar < ICHEM_END; ivar++){
            val = fmax(ustate[idx + ivar], 1e-10);
            ustate[idx + ivar] = val;
        }
#endif 
    }
    return 1;
}

double minmod(double a, double b){
    double sign;
    if( a*b > 0) {
        return 0;
    }
    if( fabs(a) > 0 ){
        sign = a/fabs(a);
    } else {
        sign = 1;
    }

    return sign * fmin(fabs(a),fabs(b));
}
double maxmod(double a, double b){
    double sign;
    if( a*b > 0) {
        return 0;
    }
    if( fabs(a) > 0 ){
        sign = a/fabs(a);
    } else {
        sign = 1;
    }

    return sign * fmax(fabs(a),fabs(b));
}
int getTVDSlopes(double *dq, double *dr, double dt){
    int icell, idx, idxm, idxp, ivar;
    double dlft, drgt, slope, sigm1, sigm2;

    // only need slopes on the first ghost cells 
    for(icell = 1; icell < NCELLS-1; icell++){
        idx = nvar*icell;
        idxp = nvar*(icell+1);
        idxm = nvar*(icell-1);
        for(ivar = 0; ivar < nvar; ivar++){
            // Superbee slope limiter
            dlft = (pstate[idx  + ivar] - pstate[idxm +ivar])/dr[icell];
            drgt = (pstate[idxp + ivar] - pstate[idx  +ivar])/dr[icell];
            
            sigm1 = minmod(  drgt,2*dlft);
            sigm2 = minmod(2*drgt,  dlft);
            slope = maxmod(sigm1, sigm2);
            dq[idx + ivar] = slope; 
        }
    }
    // set dq to zero on outmost boundary
    for(ivar = 0; ivar < nvar; ivar++){
        dq[ivar] = 0;
        dq[(NCELLS-1)*nvar+ivar]=0;
    }
    return 1;
}

int getRiemannStates(double *dq, double *qP, double *qM, double *rs, double *dr, double dt){
    int icell, idx, ivar;
    double rho, drhox, Srho;
    double vel, dvelx, Svel;
    double pre, dprex, Spre;
    double qi , dqix , Sqi;
    double dtdx;
    // Geometric terms
    double twooR, drp, drm;
    for(icell = 1; icell < NCELLS-1; icell++){
        dtdx = dt/dr[icell];
        idx = nvar*icell;
        // Get hydro variables
        rho = pstate[idx];
        vel = pstate[idx+1];
        pre = pstate[idx+2];
        
        // Slopes
        drhox = dq[idx];
        dvelx = dq[idx+1];
        dprex = dq[idx+2];

        //Source terms 
        Srho = -vel*drhox - dvelx*rho;
        Svel = -vel*dvelx - dprex/rho;
        Spre = -vel*dprex - dvelx*adi*pre;

        // Hydro left & right
        qM[idx]   = rho - 0.5*drhox + 0.5*Srho*dtdx;
        qM[idx+1] = vel - 0.5*dvelx + 0.5*Svel*dtdx;
        qM[idx+2] = pre - 0.5*dprex + 0.5*Spre*dtdx;
                                                   
        qP[idx]   = rho + 0.5*drhox + 0.5*Srho*dtdx;
        qP[idx+1] = vel + 0.5*dvelx + 0.5*Svel*dtdx;
        qP[idx+2] = pre + 0.5*dprex + 0.5*Spre*dtdx;

        //Geometric source terms if needed
        if(geometry == 1){
            drp= rs[icell] + dr[icell]/2;
            drm= rs[icell] - dr[icell]/2;
            twooR = 2/rs[icell];

            // add geometric source terms to left and right states
            qM[idx]     = qM[idx]     - rho*vel*twooR * 0.5*dt;
            //qM[idx + 1] = qM[idx + 1] - pre/rho  * 0.5*dt;
            qM[idx + 2] = qM[idx + 2] - pre*vel*twooR  * 0.5*dt;
            
            qP[idx]     = qP[idx]     - rho*vel*twooR * 0.5*dt;
            //qP[idx + 1] = qP[idx + 1] - pre/rho  * 0.5*dt;
            qP[idx + 2] = qP[idx + 2] - pre*vel*twooR* 0.5*dt;

            if(drm < dr[icell/2]){ // No velocity at zero
                qM[idx + 1 ] = 0;
            }
            if(drp <= dr[icell]/2){ // No velocity at zero
                qP[idx + 1 ] = 0;
            }
        }
        if(qP[idx] <= 0){
            printf("\n idx = %d \n", icell);
        }
        if(qM[idx] <= 0){
            printf("\n idx = %d \n", icell);
        }
        // other variables left and right
#ifdef useChemistry
        for(ivar = ICHEM_START; ivar < ICHEM_END; ivar++){
            // cell centered value
            qi = pstate[idx + ivar];
            // slope
            dqix = dq[idx + ivar];
            //source terms
            Sqi  = -vel*dqix;

            qM[idx+ivar] = qi - 0.5*dqix + 0.5*Sqi*dtdx;
            qP[idx+ivar] = qi + 0.5*dqix + 0.5*Sqi*dtdx;
        }
#endif
    }
    return 1;  
}
int getHLLCFlux(double *qL, double *qR, double*flux) {
    double rhoL, velL, preL, etotL, eL; 
    double rhoR, velR, preR, etotR, eR;
    double rho0, vel0, pre0, var0, etot0;
    double velS, preS;
    double rhoSL, etotSL;
    double rhoSR, etotSR;
    double SL, SR;
    double entho;
    double cfastL, rcL;
    double cfastR, rcR;
    int ivar;

    entho = 1/(adi-1);
    
    // Get left and right state variables
    rhoL = qL[0]; 
    velL = qL[1];
    preL = qL[2]; 
    eL   = preL*entho;
    etotL = eL + 0.5*rhoL*velL*velL;
    
    rhoR = qR[0]; 
    velR = qR[1];
    preR = qR[2]; 
    eR   = preR*entho;
    if(rhoR == 0){
        return -1;
    }
    etotR = eR + 0.5*rhoR*velR*velR;
    
    // Largest eigenvalues 
    cfastL = sqrt(adi* preL/rhoL);
    cfastR = sqrt(adi* preR/rhoR);
    
    SL=fmin(velL,velR)-fmax(cfastL,cfastR);
    SR=fmax(velL,velR)+fmax(cfastL,cfastR);
    
    rcL = rhoL*(velL-SL);
    rcR = rhoR*(SR-velR);
    
    // compute star states
    velS = (rcR*velR + rcL*velL +         (preL-preR))/(rcR+rcL);
    preS = (rcR*preL + rcL*preR + rcL*rcR*(velL-velR))/(rcR+rcL);
    
    // left star states
    rhoSL  = rhoL*(SL-velL)                           /(SL-velS);
    etotSL =     ((SL-velL)*etotL-preL*velL+preS*velS)/(SL-velS);
    //eSL    =   eL*(SL-velL)                           /(SL-velS); 
    
    // right star states
    rhoSR  = rhoR*(SR-velR)                           /(SR-velS);
    etotSR =     ((SR-velR)*etotR-preR*velR+preS*velS)/(SR-velS);
    //eSR    =   eR*(SR-velR)                           /(SR-velS); 

    // determine solution
    if(SL > 0){
        rho0  = rhoL;
        vel0  = velL;
        pre0  = preL;
        etot0 = etotL;
        //e0    = eL; 
    } else if(velS  > 0){
        rho0  = rhoSL;
        vel0  = velS;
        pre0  = preS;
        etot0 = etotSL;
        //e0    = eSL;
    } else if(SR > 0){
        rho0  = rhoSR;
        vel0  = velS;
        pre0  = preS;
        etot0 = etotSR;
        //e0    = eSR; 
    } else {
        rho0  = rhoR;
        vel0  = velR;
        pre0  = preR;
        etot0 = etotR;
        //e0    = eR; 
    }

    // Calculate HLLC fluxes 
    flux[0] = rho0*vel0;
    flux[1] = rho0*vel0*vel0+pre0;
    flux[2] = (etot0+pre0)*vel0;
    // Advected variables
#ifdef useChemistry
    for(ivar = ICHEM_START; ivar < ICHEM_END; ivar++){
        if(velS > 0){
            var0 = qL[ivar];
        } else {
            var0 = qR[ivar];
        }   
        flux[ivar] = var0*vel0;
    }
#endif
    return 1;
}

int getFluxes(double *qM, double *qP, double *fluxes){
    int iinter, idx, ivar, ierr;
    double qL[nvar], qR[nvar], flux[nvar];

    // Loop over interfaces
    for(iinter=0; iinter < NINTER; iinter++){
        // left states
        idx = (iinter+1)*nvar;
        for(ivar = 0; ivar < nvar; ivar++){
            qL[ivar] = qP[idx+ivar];
        }
        // right states 
        idx = (iinter+2)*nvar;
        for(ivar = 0; ivar < nvar; ivar++){
            qR[ivar] = qM[idx+ivar];
        }

        // now get flux at interface
        ierr = getHLLCFlux(qL, qR, flux);
        if(ierr < 0){
            printf("\n HLLC FLUX ERROR AT INTERFACE %d involving cells %d & %d\n",iinter,iinter+1, iinter+2);
            return -1;
        }
        // save fluxes
        idx = iinter*nvar;
        for(ivar = 0; ivar < nvar; ivar++){
            fluxes[idx+ivar] = flux[ivar];
        }
    }
    return -1;
}

int doHydroStep(double dt){
    int ierr;
    int icell, idx,ifp,ifn, ivar;
    double dq[nvar*NCELLS], qM[nvar*NCELLS], qP[nvar*NCELLS];
    double fluxes[nvar*NINTER], unew;
    double rm,rp, surfm = 1, surfp = 1, vol;
    double smom;
    ierr = setBoundary();
    if(ierr < 0){
        return -1;
    }
    // update the primitive variables 
    ierr = toPrimitive();
    if(ierr < 0){
        return -1;
    }
    if(useHydro > 0) { 
        // calculate TVD Slope
        ierr = getTVDSlopes(dq, dr, dt);
        if(ierr < 0){
            return -1;
        }

        // calculate left and right riemann states
        ierr = getRiemannStates(dq, qP, qM, rs, dr, dt);
        if(ierr < 0){
            return -1;
        }
        // calculate Fluxes 
        ierr = getFluxes(qM, qP, fluxes);
        if(ierr < 0){
            return -1;
        }

        // Update the conservative variables
        for(icell = 2; icell < NCELLS-2; icell++){
            idx  = icell*nvar;
            ifp = (icell-1)*nvar;
            ifn = (icell-2)*nvar;

            if(geometry == 1){
                rm = rs[icell] - dr[icell]*0.5;
                rp = rs[icell] + dr[icell]*0.5;
                
                surfm = rm*rm;//*dx_sph;
                surfp = rp*rp;//*dx_sph;
                
                vol = (rp*rp*rp -rm*rm*rm)/3.;
            } else {
                vol = dr[icell];
            }

            for(ivar = 0; ivar<nvar; ivar++){
                unew = ustate[idx+ivar]+(fluxes[ifn+ivar]*surfm-fluxes[ifp+ivar]*surfp)*dt/vol;

                if((ivar==0 || ivar==2) && unew<0 ){
                    printf("\nNEGATIVE IVAR = %d\n",ivar);
                    printf("u-1 = %.4e, u=%.4e , u+1=%.4e \n",ustate[(icell-1)*nvar+ivar], ustate[(icell)*nvar+ivar], ustate[(icell+1)*nvar+ivar]);
                    printf("u-1 = %.4e, u=%.4e , u+1=%.4e \n",ustate[(icell-1)*nvar+ivar], ustate[idx+ivar], ustate[(icell+1)*nvar+ivar]);
                    printf("FL = %.4e FR= %.4e \n",surfm*fluxes[ifn+ivar]*dt/vol,surfp*fluxes[ifp+ivar]*dt/vol);
                    printf("icell = %d \n",icell); 
                    return -1;
                }
                ustate[idx + ivar] = unew;
            }
            // add geometric source terms
            if(geometry == 1){
                // Adjusted pressure such that P = constant & u = 0 => static
                // qstate still not updated
                smom = (surfp - surfm) * pstate[idx + 2]/vol;
                //smom = 2*qstate[idx + 2]/rs[icell];
                ustate[idx + 1] = ustate[idx + 1] + smom*dt;
            }

        }
    }
    fixState();
    return 1;
}

