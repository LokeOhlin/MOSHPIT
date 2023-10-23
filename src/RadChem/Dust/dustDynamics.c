#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dust.h>
#include <hydro.h>
#include <radchem.h>
#include <moshpit.h>
#include <rtpars.h>
#include <dustRadiation.h>
#include <constantsAndUnits.h>
#ifdef useDust


/*
    METHOD TO UPDATE DUST VELOCITY FROM RADIATION PRESSURE 
 */

int dustCalcRadPress(double dt){
    int idx, ibin, graphite, imax;
    double mass, Eabs_dust, dust_accel, velmax;
    velmax = -1e99;
    for(idx = 0; idx < NdustBins; idx++){
        
        // get index in the local arrays
        ibin = globalToLocalIndex(idx);
        // If no grains here, dont accelerate them
        if(number[ibin] <=0){
            continue;
        } 

        if(ibin < isilicone){
            graphite = 1;
        } else {
            graphite = 0;
        }
        
        if(graphite){
            mass = aveMatom_c * Natoms[ibin];
        } else {
            mass = aveMatom_s * Natoms[ibin];
        }
        // get the energy absorption rate
        Eabs_dust = getAbsorption(ibin, graphite);
        // Calcualte acceleration
        dust_accel = Eabs_dust/clght/mass;
        // update velocity
        velocity[ibin] = velocity[ibin] + dust_accel*dt;
        if(velocity[ibin] > velmax){
            velmax = velocity[ibin];
            imax = idx;
        }
    }
    return 1;

}


/*
    DUST GAS DRAG

    METHOD FROM Benítez-LLambay & Krapp 2019 / Krapp &  Benítez-LLambay 2020
 */

double drag_coefficient(int iabin, int graphite, double rhoDust, double velDust, double rho, double vel, double cs){
    double rhoGrain;
    // If there are no dust grains left in this bin, set drag constant to zero
    if(rhoDust < 1e-40){
        return 0;
    }
    if(drag_mode == 0){
        if(graphite){
            rhoGrain = rho_c;
        } else {
            rhoGrain = rho_s;
        }
        return (2.*sqrt(2) * rho * cs)/(sqrt(M_PI*adi) * abin_c[iabin] * rhoGrain) * sqrt(1. + 9.*M_PI/128. * pow(abs(velDust - vel)/cs, 2)) ;
    } else if (drag_mode == 1){
        return drag_const / rhoDust;

    } else if (drag_mode == 2){
        return drag_const / rhoDust * pow(fabs(velDust - vel), drag_par);

    } else if (drag_mode == 3){
        return drag_const / rhoDust * (1 + drag_par * pow(fabs(velDust - vel),2));
    
    } else {
        return drag_const / rhoDust * sqrt(1 + drag_par * pow(fabs(velDust - vel),2)); 
    }
}

int set_drag_coefficient(double rho, double vel, double cs, double *dtmax){
    int idx, ibin, iabin, graphite;
    double velDust, rhoDust, norm;  
    // set to stupid high value
    *dtmax = 1e99;
    // loop over dust species
    for(idx = 0; idx < NdustBins; idx++){

        ibin = globalToLocalIndex(idx);
        if(ibin < isilicone){
            graphite = 1;
            norm = 4.*M_PI/3. * rho_c;
        } else {
            graphite = 0;
            norm = 4.*M_PI/3. * rho_s;
        }
        iabin = ibin - isilicone * (1-graphite);

        rhoDust = norm * getMass(number[ibin], slope[ibin], iabin);
        velDust = velocity[ibin];
        // get the drag coefficient for the current dust species from the gas
        dragCoef[ibin] = drag_coefficient(iabin, graphite, rhoDust, velDust, rho, vel, cs);
        *dtmax = fmin(*dtmax, drag_dtmax_fact/dragCoef[ibin]);
    }
    return 1;
}

int dustCalcGasDrag(double rho, double *vel, double cs, double dt_drag){
    int idx, ibin, iabin, graphite;
    double rhoDust, velDust, contrib_fact;
    double rhoTot, momTot;
    double norm;

    // Initialize the density and momentum quantities described in Krapp & Benítez-LLambay 2020
    rhoTot = rho;
    momTot = *vel * rho;
    // Loop over dust species
    for(idx = 0; idx < NdustBins; idx++){
    
        // get index in the local arrays
        ibin = globalToLocalIndex(idx);
        if(ibin < isilicone){
            graphite = 1;
            norm = 4.*M_PI/3. * rho_c;
        } else {
            graphite = 0;
            norm = 4.*M_PI/3. * rho_s;
        }
        iabin = ibin - isilicone * (1-graphite);
        

        // Density and Momentum of the dust species
        rhoDust = norm * getMass(number[ibin], slope[ibin], iabin) ;
        velDust = velocity[ibin];
    
        // Calculate the contribution to the momentum and density
        contrib_fact = dt_drag*dragCoef[ibin]/(1+dt_drag*dragCoef[ibin]);
        rhoTot = rhoTot + rhoDust*contrib_fact;
        momTot = momTot + rhoDust*velDust*contrib_fact;
    }
    
    //printf("%.4e %.4e %.4e\n", rhoTot, momTot, momTot/rhoTot);
    // Update dust momentum
    for(idx = 0; idx < NdustBins; idx++){
        // get index in the local arrays
        ibin = globalToLocalIndex(idx);
        if(ibin < isilicone){
            graphite = 1;
            norm = 4.*M_PI/3. * rho_c;
        } else {
            graphite = 0;
            norm = 4.*M_PI/3. * rho_s;
        }
        iabin = ibin - isilicone * (1-graphite);
        

        // Density and Momentum of the dust species
        rhoDust = norm * getMass(number[ibin], slope[ibin], iabin) ;
        velDust = velocity[ibin];
    
        contrib_fact = dt_drag*dragCoef[ibin]/(1+dt_drag*dragCoef[ibin]);
        velocity[ibin] = (velDust/(1+dt_drag*dragCoef[ibin]) + contrib_fact * momTot/rhoTot);
        if(velocity[ibin] != velocity[ibin]){
            printf("ibin %d : vold %.4e, vnew %.4e, dragCoef %.4e, contrib_fact %.4e, momTot %.4e rhoTot %.4e \n", ibin, velDust, velocity[ibin], dragCoef[ibin], contrib_fact, momTot, rhoTot);     
    
            exit(0);
        }
    }
    *vel = momTot/rhoTot;
    return 1;
}
int doDustGasDrag_cell(double rho, double *vel, double cs, double subcycle_step){
    double time_drag = 0;
    double dt_drag = 1e99;
    int ierr;
    int istep = 0;  
    // Subcycle in dt
    while (time_drag < subcycle_step){
        // set the drag coefficients and estimate a desired dt
        ierr = set_drag_coefficient(rho, *vel, cs, &dt_drag);
        if(ierr < 0){
            printf("ERROR in set_drag_coefficient");
            return -1;
        }
        dt_drag = fmin(dt_drag, subcycle_step - time_drag);
        // calulate the effect of the drag 
        ierr = dustCalcGasDrag(rho, vel, cs, dt_drag);
        if(ierr < 0){
            printf("ERROR in dustCalcGasDrag");
            return -1;
        }
        time_drag = time_drag + dt_drag;
        istep = istep + 1;
        //if(time_drag >= 4.5e-2){
        //    printf(" %.4e %.4e %.4e %.4e\n", time_drag, dt_drag, *vel, velocity[globalToLocalIndex(0)]);
        //    exit(0);
        //}
    }
    
    return istep;
}

int doDustGasDrag(double dt_step){
    int icell, ierr;
    double rho, vel, pre, cs, Mtot;
    
    for(icell = NGHOST; icell < NCELLS - NGHOST; icell++){
        rho = ustate[icell*nvar];
        vel = ustate[icell*nvar + 1]/rho;
        pre = (ustate[icell*nvar + 2]/rho - 0.5*vel*vel) * (adi-1);
        cs = sqrt(adi*pre);
        // set dust arrays
        ierr = setBinsCell(icell, &Mtot);     
        
        // Update dust velocities
        ierr = doDustGasDrag_cell(rho, &vel, cs, dt_step);
        
        // set momentum and update energy 
        ustate[icell*nvar + 1] = vel*rho;
        ustate[icell*nvar + 2] = rho*(pre/(adi-1) + 0.5*vel*vel);
    
        // get dust arrays
        ierr = getBinsCell(icell, &Mtot);     
        
    }
    return ierr;
}



/*
    DUST TRANSPORT
*/




/*
    We do the same calcualtion of the riemann states as one would do with gas.
    We treat both dust mass density and slope as conserved quantities
    
    Using these riemann states the fluxes are determined using the method of
    Paardekooper & Mellema 2006

 */



int addGeometricSourceTerms_interfaces_dust(int icell, double *qM, double *qP, double dt){
    double Sgeo, twoOverR;
    int idx, idust;
    if(geometry == 0){
        return 1;
    }
    twoOverR = 2./rs[icell];
    for(icell = 1; icell < NCELLS -1; icell++){
        for(idust = 0; idust < NdustBins; idust++){
            idx = icell*nvar + IDUST_START + idust*NdustVar;
            
            //dust density    = -2rho_dust * v_dust / r
            Sgeo = -pstate[idx]*pstate[idx+2]*twoOverR;
            qP[idx] += 0.5*dt*Sgeo;
            qM[idx] += 0.5*dt*Sgeo;
        }
    }
    return 1;
}


int getTVDslopes_dust(double *dpstate){
    int ivar, idust, idx, idxp, idxm, icell;
    double dminus, dplus;
    for(icell = 1; icell < NCELLS -1; icell++){
        for(idust = 0; idust < NdustBins; idust++){
            for(ivar=0; ivar < NdustVar; ivar++){
                idx  = icell*nvar + IDUST_START + idust*NdustVar + ivar;
                idxp = (icell+1)*nvar + IDUST_START + idust*NdustVar + ivar;
                idxm = (icell-1)*nvar + IDUST_START + idust*NdustVar + ivar;
                
                dplus   = (pstate[idxp] - pstate[idx] );
                dminus  = (pstate[idx]  - pstate[idxm]);

                dpstate[idx] = vanLeer(dminus, dplus); 
            }
        }
    }
    return 1;
}



int reconstructStates_dust(double *qM, double *qP, double *dpstate, double dt){
    int ivar, idust, idx, idxm, idxp, icell;
    double vel_dust;
    double dtdx, halfdtdx, twothirddtdx, fourthirddtdx_sq;
    double dPrim, prim6, primLeft, primRight, dminusPrim, dplusPrim, min_eigenVal, max_eigenVal;
    for(icell = 1; icell < NCELLS -1; icell++){
        dtdx = dt/dr[icell];
        halfdtdx = 0.5*dtdx;
        twothirddtdx = 2./3.*dtdx;
        fourthirddtdx_sq = 4./3.*pow(dtdx,2.);
        for(idust = 0; idust < NdustBins; idust++){
            idx  = icell*nvar + IDUST_START + idust*NdustVar;
            idxp = (icell+1)*nvar + IDUST_START + idust*NdustVar;
            idxm = (icell-1)*nvar + IDUST_START + idust*NdustVar;
                

            vel_dust = pstate[idx + 2]; 
            max_eigenVal = fmax(vel_dust,0);
            min_eigenVal = fmin(vel_dust,0);
            for(ivar=0; ivar < NdustVar; ivar++){
                if(interpolator == 0){ //Piecewise linear muscle
                    qP[idx + ivar] = pstate[idx + ivar] + (0.5 - max_eigenVal*halfdtdx)*dpstate[idx + ivar];
                    qM[idx + ivar] = pstate[idx + ivar] - (0.5 + min_eigenVal*halfdtdx)*dpstate[idx + ivar];
                } else {
                    primRight = 0.5*(pstate[idxp + ivar] + pstate[idx  + ivar]) - (dpstate[idxp + ivar] + dpstate[idx  + ivar])/6.;
                    primLeft  = 0.5*(pstate[idx  + ivar] + pstate[idxm + ivar]) - (dpstate[idx  + ivar] + dpstate[idxm + ivar])/6.;

                    primRight = fmax(fmin(pstate[idx + ivar], pstate[idxp + ivar]), fmin(fmax(pstate[idx + ivar], pstate[idxp + ivar]), primRight));
                    primLeft  = fmax(fmin(pstate[idx + ivar], pstate[idxm + ivar]), fmin(fmax(pstate[idx + ivar], pstate[idxm + ivar]), primLeft));                    
                    
                    // Ensure monoticity
                    if( (primRight - pstate[idx + ivar])*(pstate[idx + ivar] - primLeft) <= 0){
                        primRight = pstate[idx + ivar];
                        primLeft  = pstate[idx + ivar];
                    }
                    if( 6*(primRight - primLeft)*(pstate[idx + ivar] - 0.5*(primLeft + primRight)) > pow((primRight - primLeft), 2.)){
                        primLeft = 3*pstate[idx + ivar] - 2*primRight;
                    }
                    if( 6*(primRight - primLeft)*(pstate[idx + ivar] - 0.5*(primLeft + primRight)) < -pow((primRight - primLeft), 2.)){
                        primRight = 3*pstate[idx + ivar] - 2*primLeft;
                    }

                    dPrim = (primRight - primLeft);
                    prim6 = 6*(pstate[idx + ivar] - 0.5*(primRight + primLeft));


                    qP[idx + ivar] = primRight - halfdtdx*max_eigenVal*(dPrim - (1 - max_eigenVal*twothirddtdx)*prim6);
                    qM[idx + ivar] = primLeft  + halfdtdx*min_eigenVal*(dPrim + (1 + min_eigenVal*twothirddtdx)*prim6);
                }
            }  // loop over dust variable
        } // loop over dust species
    } // loop over cells
    return 1;
}

double limiter(double a, double b){
    return fmax(0, fmin(roe_p * a, fmax(b, fmin(a, roe_p * b)))) + \
            fmin(0, fmax(roe_p * a,   fmin(b,   fmax(a,   roe_p * b))));
}

double limit_2nd_order_flux(double rhoL, double rhoR, double F1rho, double F2rho, double dtdx){
    double frac_p, frac_m;
    double rho1, rhom, rhop;
    double epsd = fmin(0.0, fmin(rhoL,rhoR)); 
    frac_p = 1.0;
    rhop = 0.5*(rhoL - 2.0*dtdx*F2rho);
    if(rhop < epsd){
        rho1 = 0.5*(rhoL - 2.0*dtdx*F1rho);
        frac_p = fmax(0, (rho1 - epsd)/(rho1 - rhop));
        printf("frac_p %.4e %.4e\n", frac_p, rho1);
    }

    frac_m = 1.0;
    rhom = 0.5*(rhoR + 2.0*dtdx*F2rho);
    if(rhom < epsd){
        rho1 = 0.5*(rhoR + 2.0*dtdx*F1rho);
        frac_m = fmax(0, (rho1 - epsd)/(rho1-rhom));
        printf("frac_m %.4e %.4e\n", frac_m, rho1);
    }
    return fmin(frac_m, frac_p);
}

int getRoeFlux_dust(double *qL, double *qR, double *am, double *a0, double *ap, double *UL, double *UR, double dtdxl, double dtdxr, double *flux) {
    double UrhoL, rhoL, slopeL, velL;
    double UrhoR, rhoR, slopeR, velR;
    
    double FrhoL, FslopeL, FmomL;
    double FrhoR, FslopeR, FmomR;
    double F1rho, F1slope, F1mom;
    double F2rho, F2slope, F2mom;
    double bm, bp;

    double a1, a2, a3, a1u, a2u, a3u; 
    double b1, b2, b3;
    double limiter_1, limiter_2, limiter_3;
    double order2_frac;
    double velS, dtdx;
    int idust, idx;

    for(idust = 0; idust < NdustBins; idust++){ 

        // Left and right states to the interface 
        idx = idust*NdustVar;
        rhoL   = qL[idx]; 
        slopeL = qL[idx+1];
        velL   = qL[idx+2];
        
        rhoR   = qR[idx]; 
        slopeR = qR[idx+1];
        velR   = qR[idx+2];

        
        // If velocities are diverging, set flux to zero
        if(velL < 0 && velR > 0){
            flux[idx] = 0;
            flux[idx + 1] = 0;
            flux[idx + 2] = 0;
            continue;
        }
        // Right and left fluxes 
        FrhoL   =   rhoL*velL;
        FslopeL = slopeL*velL;
        FmomL   =   rhoL*velL*velL;
        
        FrhoR   =   rhoR*velR;
        FslopeR = slopeR*velR;
        FmomR   =   rhoR*velR*velR;

        // Compute singularity velocity
        velS = (sqrt(rhoL)*velL +sqrt(rhoR)*velR)/(sqrt(rhoL) + sqrt(rhoR));
       
        // Signal speeds
        bm = fmin(0.0, fmin(velS,velL));
        bp = fmax(0.0, fmax(velS,velR));

        // First order fluxes
        F1rho = (bp*FrhoL - bm*FrhoR + bm*mp*(rhoR - rhoL))/(bp-bm);
        F1slope = (bp*FslopeL - bm*FslopeR + bm*mp*(slopeR - slopeL))/(bp-bm);
        F1mom = (bp*FmomL - bm*FmomR + bm*mp*(rhoR*velR - rhoL*velL))/(bp-bm);

        // first order coefficients
        a1 = a0[idx];
        a2 = a0[idx + 1];
        a3 = a0[idx + 2];

        // upwind
        if(velS > 0){
            a1u = am[idx];
            a2u = am[idx + 1];
            a3u = am[idx + 2];
            dtdx = dtdxl;
        } else {
            a1u = ap[idx];
            a2u = ap[idx + 1];
            a3u = ap[idx + 2];
            dtdx = dtdxr;
        }

        // Limiters
        limiter_1 = limiter(a1, a1u);
        limiter_2 = limiter(a2, a2u);
        limiter_3 = limiter(a3, a3u);

        // second order coefficiencts
        b1  = -fabs(velS)*a1 + limiter_1*(fabs(velS) - dtdx*velS*velS);
        b2  = -fabs(velS)*a2 + limiter_2*(fabs(velS) - dtdx*velS*velS);
        b3  = -fabs(velS)*a3 + limiter_3*(fabs(velS) - dtdx*velS*velS);

        // second order fluxes
        F2rho   = 0.5*(FrhoL + FrhoR + b1); 
        F2slope = 0.5*(FslopeL + FslopeR + b2); 
        F2mom   = 0.5*(FmomL + FmomR + b3); 
    
        // To avoid rho <= 0, we limit the flux
        
        // density of left and right cell
        UrhoL = UL[idx];
        UrhoR = UR[idx];
        
        order2_frac = limit_2nd_order_flux(UrhoL, UrhoR, F1rho, F2rho, dtdx);

        flux[idx]     = order2_frac * F2rho   + (1 - order2_frac) * F1rho;
        flux[idx + 1] = order2_frac * F2slope + (1 - order2_frac) * F1slope;
        flux[idx + 2] = order2_frac * F2mom   + (1 - order2_frac) * F1mom;

    }
    return 1;
}

int getAdvectionFlux(double *qL, double *qR, double *flux){
    int idust, idx;
    double velS;
    double rhoL, slopeL, velL;
    double rhoR, slopeR, velR;
    for(idust = 0; idust < NdustBins; idust++){ 

        // Left and right states to the interface 
        idx = idust*NdustVar;
        rhoL   = qL[idx]; 
        slopeL = qL[idx+1];
        velL   = qL[idx+2];
        
        rhoR   = qR[idx]; 
        slopeR = qR[idx+1];
        velR   = qR[idx+2];
       
        if(velL < 0 && velR > 0){
            flux[idx] = 0;
            flux[idx  + 1] = 0;
            flux[idx  + 2] = 0;
            continue;
        } 
        // central velocity
        velS = (rhoR*velR + rhoL*velL)/(rhoR + rhoL);

        if( velS > 0){
            flux[idx]     = velL*rhoL;
            flux[idx + 1] = velL*slopeL;
            flux[idx + 2] = flux[idx]*velL; 
        } else {
            flux[idx]     = velR*rhoR;
            flux[idx + 1] = velR*slopeR;
            flux[idx + 2] = flux[idx]*velR; 
        } 


    }
    return 1;
}
int getCoefficients(double *qL, double *qR, double *as){
    double rhoL, slopeL, velL;
    double rhoR, slopeR, velR;
    
    double FrhoL, FslopeL, FmomL;
    double FrhoR, FslopeR, FmomR;
    double dFrho, dFslope, dFmom;
    double velS;
    int idust, idx;

    for(idust = 0; idust < NdustBins; idust++){ 
        // Left and right states 
        idx = idust*NdustVar;
        rhoL   = qL[idx]; 
        slopeL = qL[idx+1];
        velL   = qL[idx+2];
        
        rhoR   = qR[idx]; 
        slopeR = qR[idx+1];
        velR   = qR[idx+2];

        // Right and left fluxes 
        FrhoL   =   rhoL*velL;
        FslopeL = slopeL*velL;
        FmomL   =   rhoL*velL*velL;
        
        FrhoR   =   rhoR*velR;
        FslopeR = slopeR*velR;
        FmomR   =   rhoR*velR*velR;

        // Differences
        dFrho   = FrhoR - FrhoL;
        dFslope = FslopeR - FslopeL;
        dFmom   = FmomR - FmomL;

        velS = (sqrt(rhoL)*velL +sqrt(rhoR)*velR)/(sqrt(rhoL) + sqrt(rhoR));
        
        // first order coefficients
        as[idx]     = dFrho / velS;
        as[idx + 1] = dFslope / velS;
        as[idx + 2] = dFmom / velS;
    }
    return 1;
}

int getFluxes_dust(double *qM, double *qP, double *fluxes, double dt){
    int iinter, idxmm, idxm, idx, idxp, idxpp, ivar, ierr;
    double dtdxl, dtdxr;
    double am[NdustBins*NdustVar], a0[NdustBins*NdustVar], ap[NdustBins*NdustVar], flux[NdustBins*NdustVar];
   
    // Loop over interfaces
    for(iinter=0; iinter < NINTER; iinter++){
        // left, current and right cell
        idxmm = (iinter+NGHOST - 2)*nvar + IDUST_START;
        idxm  = (iinter+NGHOST - 1)*nvar + IDUST_START;
        idx   = (iinter+NGHOST    )*nvar + IDUST_START;
        idxp  = (iinter+NGHOST + 1)*nvar + IDUST_START;
        idxpp = (iinter+NGHOST + 2)*nvar + IDUST_START;
        
        //if first step populate coefficient matrices
        if(iinter == 0){
            ierr =  getCoefficients(qP + idxmm, qM + idxm, am);
            ierr =  getCoefficients(qP + idxm , qM + idx, a0);
        } else {
            // Otherwise copy the ones already calcuated
            for(ivar = 0; ivar <NdustBins*NdustVar; ivar ++){
                am[ivar] = a0[ivar];
                a0[ivar] = ap[ivar];
            }
        }
        // Calculate the coefficient on the next interface
        ierr = getCoefficients(qP + idx, qM + idxp, ap);
        
        dtdxl = dt/dr[iinter + NGHOST - 1];
        dtdxr = dt/dr[iinter + NGHOST];
        
        // now get flux at interface
        //ierr = getRoeFlux_dust(qP + idxm, qM + idx, am, a0, ap, ustate + idxm, ustate + idx, dtdxl, dtdxr, flux);
        ierr = getAdvectionFlux(qP+idxm, qM+idx, flux);
        if(ierr < 0){
            printf("\n ROE DUST FLUX ERROR AT INTERFACE %d involving cells %d & %d\n",iinter + NGHOST - 1,iinter + NGHOST, iinter + NGHOST +1);
            return -1;
        }
        // save fluxes
        idx = iinter*nFluxVar + IDUST_START;
        for(ivar = 0; ivar < NdustBins*NdustVar; ivar++){
            fluxes[idx+ivar] = flux[ivar];
        }
    }
    return 1;
}


int updateDustVars(int icell, double *fluxes, double *uold, double surfm, double surfp, double vol, double dt){
    int idust, ivar, ifp, ifn, idx;
    
    for(idust = 0; idust < NdustBins; idust++){
        idx = icell*nvar + IDUST_START + idust*NdustVar;
        ifp = (icell - NGHOST + 1)*nFluxVar + IDUST_START + idust*NdustVar;
        ifn = (icell - NGHOST    )*nFluxVar + IDUST_START + idust*NdustVar;
        ustate[idx] = ustate[idx] + (fluxes[ifn]*surfm - fluxes[ifp]*surfp)*dt/vol;
        ustate[idx+1] = ustate[idx+1] + (fluxes[ifn+1]*surfm - fluxes[ifp+1]*surfp)*dt/vol;
        ustate[idx+2] = ustate[idx+2] + (fluxes[ifn+2]*surfm - fluxes[ifp+2]*surfp)*dt/vol;
    }
    return 1;    
}
#endif
