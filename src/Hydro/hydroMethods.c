#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <hydro.h>
#ifdef useDust
#include <dust.h>
#endif



int getCFL(double *cfl){
    int icell, idust, idx;
    double rho, rhoinv, pre, cs, vel, wspeed;
    if(useHydro == 0){ // if no hydro no cfl
        cfl[0] = 1e99;
        return 1;
    }
    *cfl=0; // initially inverted 
    for(icell = NGHOST; icell < NCELLS - NGHOST; icell++){
        rho = ustate[icell*nvar];
        rhoinv = 1/rho;

        vel = ustate[icell*nvar+1]*rhoinv;
        if(isothermal == 1){
            cs = cs_init;
        } else {
            pre = (ustate[icell*nvar+2]*rhoinv-0.5*vel*vel)*(adi-1);
            cs = sqrt(adi*pre);
        }
        wspeed = cs + fabs(vel);
        *cfl = fmax(wspeed/dr[icell], *cfl);
#if defined useDustDynamics || defined growthUpdateVelocities
        for(idust = 0; idust < NdustBins; idust ++){
            idx =  icell*nvar + IDUST_START + idust*NdustVar;
            if(ustate[idx] > 0){
                vel = fabs(ustate[idx + 2] / ustate[idx]);
                if(vel != vel){
                    printf("ERROR: NAN VELOCITY \n");
                    exit(0);
                }
                if(vel/dr[icell] > *cfl){
                    printf("%.4e %.4e %.4e\n", vel, ustate[idx], ustate[idx+2]);
                }
                *cfl = fmax(vel/dr[icell], *cfl);
            } 
        }
#endif 
    }
    *cfl = courant_number/cfl[0];
    return 1;
}

int setBoundary(){
    int ivar, ighost;
    
    // Loop over variables
    for( ivar = 0; ivar < nvar; ivar++){
        // Loop over ghost cells
        for( ighost = 0; ighost < NGHOST; ighost++){
            // Left boundary
            if(left_bound == 0){ // reflective
                if(ivar == 1){
                    ustate[nvar*ighost+ivar] = -ustate[(2*NGHOST - 1 - ighost) *nvar+ivar];
                }else{
                    ustate[nvar*ighost+ivar] =  ustate[(2*NGHOST - 1 - ighost) *nvar+ivar];
                }
            } else if(left_bound == 1){ // outflow
                if(ivar == 1){
                    ustate[nvar*ighost+ivar] = fmin(ustate[NGHOST*nvar+ivar],0);
                }else{
                    ustate[nvar*ighost+ivar] = ustate[NGHOST*nvar+ivar];
                }
            } else if(left_bound == 2){ // fixed
                if(ivar == 0){
                    ustate[nvar*ighost + ivar] = bdensL;
                } else if(ivar == 1){
                    ustate[nvar*ighost + ivar] = bvelL;
                } else if(ivar == 2){
                    ustate[nvar*ighost + ivar] = benerL;
                } else { // all others variables set as in outflow
                    ustate[nvar*ighost + ivar] = ustate[NGHOST*nvar+ivar];
                }
            } else if(left_bound == 3){ // Periodic
                ustate[nvar*ighost + ivar] = ustate[(NCELLS - 2*NGHOST + ighost)*nvar  + ivar];
            }
            // Right boundary
            if(right_bound == 0){ // reflective
                if(ivar == 1){
                    ustate[(NCELLS - 1 - ighost)*nvar + ivar] = -ustate[(NCELLS - 2*NGHOST + ighost)*nvar +ivar];
                } else {
                    ustate[(NCELLS - 1 - ighost)*nvar + ivar] = ustate[(NCELLS - 2*NGHOST + ighost)*nvar +ivar];
                }
            } else if(right_bound == 1) { // outflow
                if(ivar == 1){
                    ustate[(NCELLS - 1 - ighost)*nvar + ivar] = fmax(ustate[(NCELLS - NGHOST - 1)*nvar +ivar],0);
                } else {
                    ustate[(NCELLS - 1 - ighost)*nvar + ivar] = ustate[(NCELLS - NGHOST - 1)*nvar +ivar];
                }
            } else if(right_bound == 2) { //fixed
                if(ivar == 0){
                    ustate[(NCELLS - 1 - ighost)*nvar + ivar] = bdensR;
                    ustate[(NCELLS-2)*nvar+ivar] = bdensR;
                } else if(ivar == 1){
                    ustate[(NCELLS - 1 - ighost)*nvar + ivar] = bvelR;
                    ustate[(NCELLS-2)*nvar+ivar] = bvelR;
                } else if(ivar == 2){
                    ustate[(NCELLS - 1 - ighost)*nvar + ivar] = benerR;
                } else { // all others variables set as in outflow
                    ustate[(NCELLS - 1 - ighost)*nvar + ivar] = ustate[(NCELLS - NGHOST - 1)*nvar +ivar];
                }
            } else if(right_bound == 3){ // Periodic
                ustate[(NCELLS - 1 - ighost)*nvar + ivar] = ustate[(2*NGHOST - 1 - ighost)*nvar + ivar];
            }
        } 
    }
    return 1;
}
int toPrimitive(){
    int icell, idx, ivar;
#ifdef useDust
    int idust;
#endif
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
        if(isothermal){
            pstate[idx + 2] = cs_init*cs_init * rho;
        } else {
            // Total energy -> internal energy -> Pressure
            // specific ekin
            ekin = 0.5*velx*velx;
            eint = ustate[idx + 2]*rhoinv - ekin;
            pstate[idx + 2]= rho * eint * (adi-1);
        }
#ifdef useChemistry
        // copy over mass fractions
        for(ivar = ICHEM_START; ivar < ICHEM_END; ivar++){
            pstate[idx + ivar] = ustate[idx + ivar];
        }
#endif 
#ifdef useDust
        for(idust = 0; idust < NdustBins; idust++){
            idx = icell*nvar + IDUST_START + idust*NdustVar;
            // denisty -> density
            pstate[idx] = ustate[idx];
            // slope   -> slope 
            pstate[idx + 1] = ustate[idx + 1];
#ifdef useDustDynamics
            // momentum -> velocity
            pstate[idx + 2] = ustate[idx + 2]/ustate[idx];
#elif defined(growthUpdateVelocitites)
            // momentum -> momentum (as scalar field)
            pstate[idx + 2] = ustate[idx + 2];
#endif
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
            val = fmax(ustate[idx + ivar], 1e-10*ustate[idx]);
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

    // Slopes only made for the inner cells  + ghost cells
    for(icell = 1; icell < NCELLS-1; icell++){
        idx = nvar*icell;
        idxp = nvar*(icell+1);
        idxm = nvar*(icell-1);
        for(ivar = 0; ivar < nvar; ivar++){
            // Superbee slope limiter
            dlft = (pstate[idx  + ivar] - pstate[idxm + ivar])/dr[icell];
            drgt = (pstate[idxp + ivar] - pstate[idx  + ivar])/dr[icell];
            
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
            twooR = 2/fabs(rs[icell]);

            // add geometric source terms to left and right states
            qM[idx]     = qM[idx]     - rho*vel*twooR * 0.5*dt;
            qM[idx + 2] = qM[idx + 2] - pre*vel*twooR * 0.5*dt;
            
            qP[idx]     = qP[idx]     - rho*vel*twooR * 0.5*dt;
            qP[idx + 2] = qP[idx + 2] - pre*vel*twooR * 0.5*dt;
        }
        
        // other variables left and right
#ifdef useChemistry
        for(ivar = ICHEM_START; ivar < ICHEM_END; ivar++){
            // cell centered value
            qi = pstate[idx + ivar];
            // slope
            dqix = dq[idx + ivar];
            //source terms
            Sqi  = -vel*dqix - dvelx*qi;

            qM[idx+ivar] = qi - 0.5*dqix + 0.5*Sqi*dtdx;
            qP[idx+ivar] = qi + 0.5*dqix + 0.5*Sqi*dtdx;
            // Geometric source term same as for density
            if(geometry == 1){
                qM[idx + ivar] = qM[idx + ivar] - qi*vel*twooR * 0.5 * dt;
                qP[idx + ivar] = qP[idx + ivar] - qi*vel*twooR * 0.5 * dt;
            }
        }
#ifdef useDust
#ifndef useDustDynamics
        for(ivar = IDUST_START; ivar < nvar; ivar++){
            // cell centered value
            qi = pstate[idx + ivar];
            // slope
            dqix = dq[idx + ivar];
            //source terms
            Sqi  = -vel*dqix - dvelx*qi;

            qM[idx+ivar] = qi - 0.5*dqix + 0.5*Sqi*dtdx;
            qP[idx+ivar] = qi + 0.5*dqix + 0.5*Sqi*dtdx;
            // Geometric source term same as for density
            if(geometry == 1){
                qM[idx + ivar] = qM[idx + ivar] - qi*vel*twooR * 0.5 * dt;
                qP[idx + ivar] = qP[idx + ivar] - qi*vel*twooR * 0.5 * dt;
            }
        }
#else 
        int ierr = getRiemannStates_dust(dq, qP, qM, dt, dtdx, icell);
        if(ierr < 0){
            return ierr;
        }     
#endif //useDynamicDust
#endif //useDust
#endif //useChemistry 
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
    rhoSL  = rhoL*(SL-velL)/(SL-velS);
    etotSL = ((SL-velL)*etotL-preL*velL+preS*velS)/(SL-velS);
    //eSL    =   eL*(SL-velL)/(SL-velS); 
    
    // right star states
    rhoSR  = rhoR*(SR-velR)/(SR-velS);
    etotSR = ((SR-velR)*etotR-preR*velR+preS*velS)/(SR-velS);
    //eSR    =   eR*(SR-velR)/(SR-velS); 

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
        if(SL > 0){
            var0  = qL[ivar];
        } else if(velS  > 0){
            var0  = qL[ivar]*(SL-velL)/(SL-velS); ;
        } else if(SR > 0){
            var0  = qR[ivar]*(SR-velR)/(SR-velS);
        } else {
            var0  = qR[ivar];
        }
        flux[ivar] = var0*vel0;
    }
#endif
#ifdef useDust
#ifndef useDustDynamics
    for(ivar = IDUST_START; ivar < nvar; ivar++){
        if(SL > 0){
            var0  = qL[ivar];
        } else if(velS  > 0){
            var0  = qL[ivar]*(SL-velL)/(SL-velS); ;
        } else if(SR > 0){
            var0  = qR[ivar]*(SR-velR)/(SR-velS);
        } else {
            var0  = qR[ivar];
        }

        flux[ivar] = var0*vel0;
    }
#endif
#endif
    // transport of internal energy. Needed for eintswitch when Ekin ~ Etot
    // based on li 2008 "A Simple Dual Implementation to Track Pressure Accurately"
    // average pressure
    if(hy_ethresh >= 0.0){
        if(flux[0] == 0.0){
            flux[nvar] = 0.0;
            flux[nvar + 1] = 0.0;
        } else { 
            if(flux[0] > 0){ // if flux flowing from left to right
                flux[nvar] = flux[0]/qL[0] * qL[2]/(adi-1); 
                flux[nvar + 1] = flux[0]/qL[0];
            } else {
                flux[nvar] = flux[0]/qR[0] * qR[2]/(adi-1);
                flux[nvar + 1] = flux[0]/qR[0];
            }
        }
    }
    return 1;
}

int getFluxes(double *qM, double *qP, double *fluxes, double dt){
    int iinter, idx, ivar, ierr;
    double qL[nvar], qR[nvar], flux[nFluxVar];

    // Loop over interfaces
    for(iinter=0; iinter < NINTER; iinter++){
        // left states
        idx = (iinter + NGHOST - 1)*nvar;
        for(ivar = 0; ivar < nvar; ivar++){
            qL[ivar] = qP[idx+ivar];
        }
        // right states 
        idx = (iinter + NGHOST)*nvar;
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
        idx = iinter*nFluxVar;
        for(ivar = 0; ivar < nFluxVar; ivar++){
            fluxes[idx+ivar] = flux[ivar];
        }
    }
    printf("%.4e, %.4e \n", fluxes[(NINTER -2)*nFluxVar], fluxes[(NINTER -1)*nFluxVar]); 
#ifdef useDustDynamics
    ierr = getFluxes_dust(qM, qP, fluxes, dt);
#endif
    return 1;
}

int doHydroStep(double dt){
    int ierr;
    int icell, idx,ifp,ifn, ivar;
    double dq[nvar*NCELLS], qM[nvar*NCELLS], qP[nvar*NCELLS];
    double fluxes[nFluxVar*NINTER], unew, uold[nvar];
    double rm,rp, surfm = 1, surfp = 1, vol;
    double smom;
    double Ethermal, Ekinetic, presStar;
    int idxm, idxp;
    ierr = setBoundary();
    if(ierr < 0){
        printf("Error in setBoundary\n");
        return -1;
    }
    // update the primitive variables 
    ierr = toPrimitive();
    if(ierr < 0){
        printf("Error in toPrimitive\n");
        return -1;
    }
    if(useHydro > 0) { 
        // calculate TVD Slope
        ierr = getTVDSlopes(dq, dr, dt);
        if(ierr < 0){
            printf("Error in getTVDslopes\n");
            return -1;
        }

        // calculate left and right riemann states
        ierr = getRiemannStates(dq, qP, qM, rs, dr, dt);
        if(ierr < 0){
            printf("Error in getRiemannStates\n");
            return -1;
        }
        // calculate Fluxes 
        ierr = getFluxes(qM, qP, fluxes, dt);
        if(ierr < 0){
            printf("Error in getFluxes\n");
            return -1;
        }

        // Update the conservative variables
        for(icell = NGHOST; icell < NCELLS-NGHOST; icell++){
            idx  = icell*nvar;
            ifp = (icell - NGHOST + 1)*nFluxVar;
            ifn = (icell - NGHOST)*nFluxVar;

            if(geometry == 1){
                rm = rs[icell] - dr[icell]*0.5;
                rp = rs[icell] + dr[icell]*0.5;
                
                surfm = rm*rm;
                surfp = rp*rp;
                
                vol = (rp*rp*rp -rm*rm*rm)/3.;
            } else {
                vol = dr[icell];
            }
            
            if(hy_ethresh >= 0.0) {
                // get initial internal temperature
                // Ekin = 0.5*rho v**2 = 0.5 * momx*momx/rho
                Ethermal = ustate[idx + 2] - 0.5*ustate[idx + 1] * ustate[idx + 1] /ustate[idx];
                // predicted half step pressure of at cell center
                presStar = (surfp*qP[idx + 2] + surfm*qM[idx + 2])/(surfm+surfp);
                
                Ethermal = Ethermal +            (fluxes[ifn + nvar]*surfm     - fluxes[ifp + nvar]*surfp)*dt/vol;
                Ethermal = Ethermal + presStar * (fluxes[ifn + nvar + 1]*surfm - fluxes[ifp + nvar + 1]*surfp)*dt/vol;
                Ethermal = fmax(Ethermal, 1e-40);
            }
            
            for(ivar = 0; ivar<nvar; ivar++){
                if(isothermal && ivar == 2){
                    continue;
                }
                unew = ustate[idx+ivar]+(fluxes[ifn+ivar]*surfm-fluxes[ifp+ivar]*surfp)*dt/vol;
#ifdef useDust
                if((ivar==0 || ivar==2 || (ivar >= IDUST_START && (ivar - IDUST_START)%NdustVar == 0)) && unew<0 ){
#else
                if((ivar==0 || ivar==2) && unew<0 ){
#endif
                    printf("\nNEGATIVE IVAR = %d\n",ivar);
                    printf("unew=%.4e\n", unew);
                    printf("u-1 = %.4e, u=%.4e , u+1=%.4e \n",ustate[(icell-1)*nvar+ivar], ustate[(icell)*nvar+ivar], ustate[(icell+1)*nvar+ivar]);
                    printf("u-1 = %.4e, u=%.4e , u+1=%.4e \n",ustate[(icell-1)*nvar+ivar], ustate[idx+ivar], ustate[(icell+1)*nvar+ivar]);
                    printf("FL = %.4e FR= %.4e \n",surfm*fluxes[ifn+ivar]*dt/vol,surfp*fluxes[ifp+ivar]*dt/vol);
                    printf("icell = %d / %d\n",icell, NCELLS - NGHOST - 1); 
                    return -1;
                }
                uold[ivar] = ustate[idx+ivar];
                ustate[idx + ivar] = unew;
            }

            // add geometric source terms
            if(geometry == 1){
                // Adjusted pressure such that P = constant & u = 0 => static
                // pstate still not updated so we can freely use
                smom = (surfp - surfm) * pstate[idx + 2]/vol;
                //smom = 2*qstate[idx + 2]/rs[icell];
                ustate[idx + 1] = ustate[idx + 1] + smom*dt;
            }
             
            if(hy_ethresh>=0.0 && isothermal == 0){
                if(Ethermal < ustate[idx + 2]*hy_ethresh){
                    Ekinetic = 0.5*ustate[idx + 1] * ustate[idx + 1] /ustate[idx];
                    //printf("%.4e %.4e \n", Ethermal, Ekinetic);
                    ustate[idx + 2] = Ekinetic + Ethermal;
                }
            }
            
            if(isothermal == 0) { 
                if(ustate[idx + 2] - 0.5 * ustate[idx + 1] * ustate[idx + 1]/ustate[idx] < 0){
                    idxm = (icell-1)*nvar;
                    idxp = (icell+1)*nvar;
                    printf("NEGATIVE INTERNAL ENERGY\n");
                    printf("cell = %d \n", icell);
                    printf("dens = %.4e %.4e %.4e \n", ustate[idxm], uold[0], ustate[idxp]); 
                    printf("velx = %.4e %.4e %.4e \n", ustate[idxm+1]/ustate[idxm], uold[1]/uold[0], ustate[idxp+1]/ustate[idxp]); 
                    printf("etot = %.4e %.4e %.4e \n", ustate[idxm+2], uold[2], ustate[idxp+2]); 
                    printf("enew = %.4e %.4e %.4e \n", ustate[idxm+2], ustate[idx + 2], ustate[idxp+2]); 
                    printf("ekin = %.4e\n",  0.5 * ustate[idx + 1] * ustate[idx + 1]/ustate[idx]); 
                    printf("eint = %.4e\n", ustate[idx+2] - 0.5 * ustate[idx + 1] * ustate[idx + 1]/ustate[idx]); 
                    printf("rhoFluxes = %.4e \t %.4e\n", fluxes[ifn + 0]*surfm, fluxes[ifp+0]*surfp);
                    printf("momFluxes = %.4e \t %.4e\n", fluxes[ifn + 1]*surfm, fluxes[ifp+1]*surfp);
                    printf("eneFluxes = %.4e \t %.4e\n", fluxes[ifn + 2]*surfm, fluxes[ifp+2]*surfp);

                    if(ustate[idx+2] - 0.5 * ustate[idx + 1] * ustate[idx + 1]/ustate[idx] < 0){
                        return -1;
                    }
                }
            }
        }
    }
    fixState();
    return 1;
}

