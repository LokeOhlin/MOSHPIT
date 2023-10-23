#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <hydro.h>
#include <moshpit.h>
#ifdef useDust
#include <dust.h>
#endif



int getCFL(double *cfl){
    int icell, idx;
#if defined useDustDynamics
    int idust;
#endif
    double rho, rhoinv, pre, cs, vel, wspeed;
    if(useHydro == 0){ // if no hydro no cfl
        cfl[0] = 1e99;
        return 1;
    }
    *cfl=0; // initially inverted 
    for(icell = NGHOST; icell < NCELLS - NGHOST; icell++){
        idx = icell*nvar;
        rho = ustate[idx];
        rhoinv = 1/rho;
        
        vel = ustate[idx+1]*rhoinv;
        if(isothermal == 1){
            cs = cs_init;
        } else {
            pre = (ustate[idx+2]*rhoinv-0.5*vel*vel)*(adi-1);
            cs = sqrt(adi*pre);
        }
        wspeed = cs + fabs(vel);
        *cfl = fmax(wspeed/dr[icell - 1],fmax(wspeed/dr[icell+1], fmax(wspeed/dr[icell], *cfl)));
#if defined useDustDynamics || defined growthUpdateVelocities
        for(idust = 0; idust < NdustBins; idust ++){
            idx =  icell*nvar + IDUST_START + idust*NdustVar;
            if(ustate[idx] > 0){
                vel = fabs(ustate[idx + 2] / ustate[idx]);
                if(vel != vel){
                    printf("ERROR: NAN VELOCITY \n");
                    exit(0);
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
                ustate[nvar*ighost + ivar] = ustate[MAX(NCELLS - 2*NGHOST + ighost, NGHOST)*nvar  + ivar];
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
                ustate[(NCELLS - 1 - ighost)*nvar + ivar] = ustate[MIN(2*NGHOST - 1 - ighost, NCELLS-NGHOST-1)*nvar + ivar];
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
            if(pstate[idx + 2] <= 0){
                printf("pstate < 0\n");
                printf("icell %d  \n", icell);
                printf("eint %.4e \n", eint);
                printf("etot %.4e \n", ustate[idx + 2]*rhoinv);
                printf("ekin %.4e \n", ekin);
                printf("rho  %.4e \n", rho);
                printf("%.4e\t%.4e\t%.4e\t\n", ustate[idx], ustate[idx+1], ustate[idx+2]);
                return -1;
            }
        }
        // copy over advected quantities
        for(ivar = IADVECT_START; ivar < IADVECT_END; ivar++){
            pstate[idx + ivar] = ustate[idx + ivar];
        }
#ifdef useDustDynamics
        for(idust = 0; idust < NdustBins; idust++){
            idx = icell*nvar + IDUST_START + idust*NdustVar;
            // denisty -> density
            pstate[idx] = ustate[idx];
            // slope   -> slope 
            pstate[idx + 1] = ustate[idx + 1];
            // momentum -> velocity
            if(ustate[idx] > 0){
                pstate[idx + 2] = ustate[idx + 2]/ustate[idx];
            } else {
                pstate[idx + 2] = 0.0;
            }
        }
#endif
    }
    return 1;
}


int fixState(){
#ifdef useChemistry
    int icell, idx, ivar;
    double val;
    for(icell = 0; icell < NCELLS; icell++){
        idx = icell *nvar;
        for(ivar = ICHEM_START; ivar < ICHEM_END; ivar++){
            val = fmax(ustate[idx + ivar], 1e-10);
            ustate[idx + ivar] = val;
        }
    }
#endif 
    return 1;
}

double minmod(double a, double b){
    double sign;
    if( a*b <= 0) {
        return 0;
    }
    if( fabs(a) > 0 ){
        sign = a/fabs(a);
    } else if (fabs(b) > 0){
        sign = b/fabs(b);
    } else {
        sign = 0.0;
    }

    return sign * fmin(fabs(a),fabs(b));
}
double maxmod(double a, double b){
    double sign;
    if( a*b <= 0) {
        return 0.0;
    }
    printf("a*b = %.4e\n", a,b);
    if( fabs(a) > 0 ){
        sign = a/fabs(a);
    } else if (fabs(b) > 0){
        sign = b/fabs(b);
    } else {
        sign = 0;
    }

    return sign * fmax(fabs(a),fabs(b));
}

double vanLeer(double a, double b){
    if (a*b <= 0){
        return 0.0;
    }
    return 2.* a*b/(a+b);
}

int addGeometricSourceTerms_interfaces(int icell, double *qM, double *qP, double dt){
    double Sgeo, cssq, twoOverR;
    int idx, ivar, ierr;
    if(geometry == 0){
        return 1;
    }
    idx = icell*nvar;
    twoOverR = 2./fabs(rs[icell]);
    if(twoOverR != twoOverR){
        printf(" rs = %.4e\n",rs[icell]);
        return -1;
    }  
    
    // density = -2rho v/r
    Sgeo = -pstate[idx]*pstate[idx + 1]*twoOverR;
    qP[idx] += 0.5*dt*Sgeo;
    qM[idx] += 0.5*dt*Sgeo;

    // velocity = 0

    // Pressure = -2 cs^2 rho v /r
    cssq = adi*pstate[idx + 2]/pstate[idx];
    Sgeo = Sgeo * cssq;
    qP[idx + 2] += 0.5*dt*Sgeo;
    qM[idx + 2] += 0.5*dt*Sgeo;

    //// Advected quantities = - 2 var vel/r
    //for(ivar = IADVECT_START; ivar < IADVECT_END; ivar++){
    //    Sgeo = -pstate[idx + ivar] * pstate[idx + 1] * twoOverR;
    //    qP[idx + ivar] += 0.5*dt*Sgeo;
    //    qM[idx + ivar] += 0.5*dt*Sgeo;
    //}

#ifdef useDustDynamics
    ierr = addGeometricSourceTerms_interfaces_dust(icell, qM, qP, dt);    
#endif 

    return 1;
}   


int addHydroGeometricSourceTerms(int icell, double *qM, double *qP, double dt){
    double smom, surfp,surfm, rm, rp, vol;
    int idx;
    if(geometry == 0){
        return 1;
    }
    rm = right_edge[icell - 1];
    rp = right_edge[icell];
    
    surfm = rm*rm;
    surfp = rp*rp;
    
    vol = (rp*rp*rp -rm*rm*rm)/3.;
    idx = icell*nvar;
    // add geometric source terms
    // Adjusted pressure such that P = constant & u = 0 => static
    
    smom = (surfp*qP[idx+2] + surfm*qM[idx+2])/(surfm + surfp);
    // pstate still not updated so we can freely use
    //smom = (surfp - surfm) * pstate[idx + 2]/vol;
    
    //smom = 2*qstate[idx + 2]/rs[icell];
    //if(icell < 2*NGHOST){
    //    printf("%d  %.4e %.4e %.4e %.4e \n",icell, qM[idx+2], qP[idx+2],smom, rs[icell]);
    //}
    ustate[idx + 1] = ustate[idx + 1] + 2.0*smom*dt*dr[icell]*rs[icell]/vol;
    
    return 1;
}
int getEigenVects(int icell, double *leftEigenVects, double *rightEigenVects, double *eigenVals){
    int idx;
    double rho, vel, pre, cs;
    idx = icell*nvar;
    rho = pstate[idx];
    vel = pstate[idx + 1];
    pre = pstate[idx + 2];
    cs  = sqrt(adi*pre/rho);
    if(cs != cs || cs <= 0){
        printf("cs <= 0 \n");
        printf("cs = %.4e \n", cs);
        printf("pre = %.4e \n", pre);
        printf("rho = %.4e \n", rho);
        return -1; 
    } 
    // lambda = vel
    eigenVals[0] = vel;
    leftEigenVects[0*3 + 0] = 1.;
    leftEigenVects[0*3 + 1] = 0;
    leftEigenVects[0*3 + 2] = -1/(cs*cs);
    
    rightEigenVects[0*3 + 0] = 1.;
    rightEigenVects[1*3 + 0] = 0;
    rightEigenVects[2*3 + 0] = 0;

    // lambda = vel + cs
    eigenVals[1] = vel + cs;
    leftEigenVects[1*3 + 0] = 0;
    leftEigenVects[1*3 + 1] = 0.5*rho/cs;
    leftEigenVects[1*3 + 2] = 0.5/(cs*cs);
    
    rightEigenVects[0*3 + 1] = 1.;
    rightEigenVects[1*3 + 1] = cs/rho;
    rightEigenVects[2*3 + 1] = cs*cs;
    
    // lambda = vel - cs
    eigenVals[2] = vel - cs;
    leftEigenVects[2*3 + 0] = 0;
    leftEigenVects[2*3 + 1] = -0.5*rho/cs;
    leftEigenVects[2*3 + 2] = 0.5/(cs*cs);
    
    rightEigenVects[0*3 + 2] = 1.;
    rightEigenVects[1*3 + 2] = -cs/rho;
    rightEigenVects[2*3 + 2] = cs*cs;
    return 1;
}


int getTVDslope(double *dpstate){
#ifdef useDustDynamics
    int ierr;
#endif
    int icell, idx, idxm, idxp, ivar, icvar;
    double leftEigenVects[9], rightEigenVects[9], eigenVals[3];
    double dcenterPrim[3], dminusPrim[3], dplusPrim[3];
    double dcenterChar, dminusChar, dplusChar, dChar[3];
    //double sign; 
    // Slopes only made for the inner cells  + ghost cells
    for(icell = 1; icell < NCELLS-1; icell++){
        idx = nvar*icell;
        idxp = nvar*(icell+1);
        idxm = nvar*(icell-1);
        

        /*
         *      Hydrodynamical properties
         *  
         *
         */
        for(ivar = 0; ivar < 3; ivar++){
            dcenterPrim[ivar] = 0;
            dminusPrim[ivar] = 0;
            dplusPrim[ivar] = 0;
            dChar[ivar] = 0;
        }

        // Get Eigenvectors and eigenvalues
        getEigenVects(icell, leftEigenVects, rightEigenVects, eigenVals);

        // Set diferences in hydro vars
        for(ivar = 0; ivar < 3; ivar++){
            dcenterPrim  [ivar]   = 0.5*(pstate[idxp + ivar] - pstate[idxm + ivar]);
            dminusPrim[ivar] = pstate[idx  + ivar] - pstate[idxm + ivar];
            dplusPrim [ivar] = pstate[idxp + ivar] - pstate[idx  + ivar];
        }

        // Project into characteristic variables and apply limiter
        for(icvar = 0; icvar < 3; icvar++){
            dcenterChar = 0; 
            dminusChar  = 0;
            dplusChar   = 0;
            for(ivar = 0; ivar < 3; ivar ++){
                dcenterChar += leftEigenVects[icvar*3 + ivar]*dcenterPrim[ivar]; 
                dminusChar  += leftEigenVects[icvar*3 + ivar]*dminusPrim[ivar];
                dplusChar   += leftEigenVects[icvar*3 + ivar]*dplusPrim[ivar];
            }
            // slope limiting
            //if(dcenterChar >= 0){
            //    sign = 1;
            //} else {
            //    sign = -1;
            //}

            dChar[icvar] = vanLeer(dminusChar, dplusChar); //sign * fmin(2*fabs(dminusChar), fmin(2*fabs(dplusChar), fabs(dcenterChar)));
            if(dChar[icvar] != dChar[icvar]){
                printf("Error: Reconstruction gave NaN characteristic slopes\n");
                printf("icvar = %d\n", icvar);
                printf("dChar = %.4e\n", dChar[icvar]);
                printf("dminusChar = %.4e\n", dminusChar);
                printf("dplusChar = %.4e\n", dplusChar);
                printf("left eigenvect : ");
                for(ivar = 0; ivar < 3; ivar++){
                    printf("%.4e\t", leftEigenVects[icvar*3 + ivar]);
                }
                printf("\n");
                return -1;
            }
        }

        // Project back into primitive 
        for(ivar = 0; ivar < 3; ivar++){
            dpstate[idx + ivar] = 0;
            for(icvar = 0; icvar < 3; icvar++){
                dpstate[idx + ivar] += rightEigenVects[ivar*3 + icvar]*dChar[icvar];
            }
            if(dpstate[idx + ivar] != dpstate[idx+ ivar]){
                printf("Error: Reconstruction gave NaN slopes\n");
                printf("ivar = %d\n", ivar);
                printf("dpstate = %.4e\n", dpstate[idx + ivar]);
                printf("right eigenvect : ");
                for(icvar = 0; icvar < 3; icvar++){
                    printf("%.4e\t", rightEigenVects[ivar*3 + icvar]);
                }
                printf("\n");
                printf("dChar : ");
                for(icvar = 0; icvar < 3; icvar++){
                    printf("%.4e\t", dChar[icvar]);
                }
                printf("\n");
                return -1;
            }
        }
        
        
        /*
         *
         *      Advected mass scalars properties
         *      for these, characteristic quanteties are the same as the primitive and 
         *      so no eigenvectors are neccesary, and the eigen value is just the velocity
         *      
         */

        for(ivar = IADVECT_START; ivar < IADVECT_END; ivar++){
            // Differences 
            dminusPrim[0]  = pstate[idx  + ivar] - pstate[idxm + ivar];
            dplusPrim [0]  = pstate[idxp + ivar] - pstate[idx  + ivar];
            
            dpstate[idx+ivar] = vanLeer(dminusPrim[0], dplusPrim[0]);

        }
    }
#ifdef useDustDynamics
    ierr = getTVDslopes_dust(dpstate);
    if(ierr < 0){
        return ierr;
    }
#endif

    return 1;
}

int reconstructStates(double *qP, double *qM, double dt){
    int icell, idx, idxm, idxp, ivar, icvar, ierr;
    double leftEigenVects[9], rightEigenVects[9], eigenVals[3], max_eigenVal, min_eigenVal;
    double dPrim[nvar], prim6[nvar], primLeft, primRight;
    double dotProd, dotProdA, dotProdB; 
    double dtdx, halfdtdx, twothirddtdx, fourthirddtdx_sq;
    double dpstate[nvar*NCELLS];

    //
    // Get initial dW
    //
    ierr = getTVDslope(dpstate);
    if(ierr < 0){
        return ierr;
    }
            

    //
    // Characteristic tracing 
    //

    for(icell = 1; icell < NCELLS - 1; icell++){
        idx = nvar*icell;
        idxp = nvar*(icell+1);
        idxm = nvar*(icell-1);
        
        // Get Eigenvectors and eigenvalues
        getEigenVects(icell, leftEigenVects, rightEigenVects, eigenVals);
        

        dtdx = dt/dr[icell];
        halfdtdx = 0.5*dtdx;
        twothirddtdx = 2./3.*dtdx;
        fourthirddtdx_sq = 4./3.*pow(dtdx,2.);
        
        if(interpolator == 0) { // Piecewise linear MUSCL
            for(ivar = 0; ivar < IADVECT_END; ivar ++){
                qP[idx + ivar] = pstate[idx + ivar];
                qM[idx + ivar] = pstate[idx + ivar];
            }
            // Loop over eigen vectors 
            for(icvar = 0; icvar < 3; icvar ++){
                if(eigenVals[icvar] == 0){
                    continue;
                }
                // convert back to characteristic
                dotProd = 0;
                for(ivar = 0; ivar < 3; ivar++){
                    dotProd += leftEigenVects[icvar*3 + ivar]*dpstate[idx + ivar];
                }
                
                // add to interfaces
                for(ivar = 0; ivar <3; ivar ++){
                    qP[idx + ivar]  += ( 0.5 - halfdtdx*eigenVals[icvar])*rightEigenVects[ivar*3+icvar]*dotProd;
                    qM[idx + ivar]  += (-0.5 - halfdtdx*eigenVals[icvar])*rightEigenVects[ivar*3+icvar]*dotProd;
                    if((qP[idx + ivar] != qP[idx + ivar]) ||(qP[idx + ivar] != qP[idx + ivar]) ){
                        printf("Error: Reconstruction gave NaN states\n");
                        printf("ivar = %d\n", ivar);
                        printf("dpstate = %.4e\n", dpstate[idx + ivar]);
                        printf("qM = %.4e\n", qM[idx + ivar]);
                        printf("qP = %.4e\n", qP[idx + ivar]);
                        return -1;
                    }
                }
            }
            
            max_eigenVal = fmax(eigenVals[0],0);
            min_eigenVal = fmin(eigenVals[0],0);
            for(ivar = IADVECT_START; ivar < IADVECT_END; ivar ++){
                qP[idx + ivar] = pstate[idx + ivar] + (0.5 - max_eigenVal*halfdtdx)*dpstate[idx +ivar];
                qM[idx + ivar] = pstate[idx + ivar] - (0.5 + min_eigenVal*halfdtdx)*dpstate[idx +ivar];
            }

        } else {  // PPM

            for(ivar = 0; ivar< IADVECT_END; ivar++){
                primRight = 0.5*(pstate[idxp + ivar] + pstate[idx  + ivar]) - (dpstate[idxp + ivar] + dpstate[idx  + ivar])/6.;
                primLeft  = 0.5*(pstate[idx  + ivar] + pstate[idxm + ivar]) - (dpstate[idx  + ivar] + dpstate[idxm + ivar])/6.;
            
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

                dPrim[ivar] = (primRight - primLeft);
                prim6[ivar] = 6*(pstate[idx + ivar] - 0.5*(primRight + primLeft)); 
                
                if(ivar < 3){
                    max_eigenVal = fmax(eigenVals[1],0);
                    min_eigenVal = fmin(eigenVals[2],0);
                } else {
                    max_eigenVal = fmax(eigenVals[0],0);
                    min_eigenVal = fmin(eigenVals[0],0);
                }
                
                if(ivar < 3){
                    qP[idx + ivar] = pstate[idx + ivar] + prim6[ivar]/12.;
                    qM[idx + ivar] = pstate[idx + ivar] + prim6[ivar]/12.;
                } else { 
                    qP[idx + ivar] = primRight - halfdtdx*max_eigenVal*(dPrim[ivar] - (1 - max_eigenVal*twothirddtdx)*prim6[ivar]);
                    qM[idx + ivar] = primLeft  + halfdtdx*min_eigenVal*(dPrim[ivar] + (1 - min_eigenVal*twothirddtdx)*prim6[ivar]);
                }
            }

            /*
             *
             *  Characteristic traicing of Hydro variables (advected variables only have one eigenvalue, velocity)
             *
             */

            
            // Loop over eigen vectors 
            for(icvar = 0; icvar < 3; icvar ++){
                // Dot product between leftEigens and dPrim and prim6
                dotProdA = 0;
                dotProdB = 0;
                for(ivar = 0; ivar < 3; ivar++){
                    dotProdA += leftEigenVects[icvar*3 + ivar]*dPrim[ivar];
                    dotProdB -= leftEigenVects[icvar*3 + ivar]*prim6[ivar];
                }
               
                for(ivar = 0; ivar <3; ivar ++){
                    qP[idx + ivar] += 0.5*(1.  - dtdx*eigenVals[icvar])*rightEigenVects[ivar*3+icvar]*dotProdA;
                    qP[idx + ivar] += 0.25*(1. - 2.*dtdx*eigenVals[icvar] 
                                               + fourthirddtdx_sq*pow(eigenVals[icvar],2.))*rightEigenVects[ivar*3+icvar]*dotProdB;
                    
                    qM[idx + ivar] += 0.5*(-1. - dtdx*eigenVals[icvar])*rightEigenVects[ivar*3+icvar]*dotProdA;
                    qM[idx + ivar] += 0.25*(1. + 2.*dtdx*eigenVals[icvar] 
                                               + fourthirddtdx_sq*pow(eigenVals[icvar],2.))*rightEigenVects[ivar*3+icvar]*dotProdB;
                }
            }
        }
    }    
#ifdef useDustDynamics
    ierr = reconstructStates_dust(qM, qP, dpstate, dt);
    if(ierr < 0){
        return ierr;
    }
#endif 


    /*
     *  
     *  Geometric source terms
     *
     */
    for(icell = 1; icell< NCELLS -1; icell++){
        ierr = addGeometricSourceTerms_interfaces(icell,  qM, qP, dt);
        if(ierr < 0){
            return ierr;
        }
    }

    return 1;
}


int getHLLCFlux(double *qL, double *qR, double *flux) {
    double rhoL, velL, preL, etotL, eL; 
    double rhoR, velR, preR, etotR, eR;
    double rho0, vel0, pre0, etot0;
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
    /*
    if(flux[0] != flux[0]){
        printf("density flux is NaN \n");
        printf("rhoL %.4e \n", rhoL);
        printf("rhoSL %.4e \n", rhoSL);
        printf("rhoSR %.4e \n", rhoSR);
        printf("rhoR %.4e \n", rhoR);
        printf("vel0 %.4e \n", vel0);
        for(ivar = 0; ivar  < 3; ivar ++){
            printf("ivar %d qL %.4e qR %.4e \n", ivar, qL[ivar], qR[ivar]);
        }
    }
    */
    // Advected variables
    for(ivar = IADVECT_START; ivar < IADVECT_END; ivar++){
        if(flux[0]> 0){
            flux[ivar] = qL[ivar]*flux[0];
        } else {
            flux[ivar] = qR[ivar]*flux[0];
        }
    }
    
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
#ifdef useDustDynamics
    ierr = getFluxes_dust(qM, qP, fluxes, dt);
#endif
    return 1;
}


int printCellHydro(int icell, int allVars, double *fluxes, double *qP, double *qM, double dt){
    int idx = icell*nvar, idxm = (icell - 1)*nvar, idxp = (icell + 1)*nvar;
    int ifp = (icell - NGHOST + 1) * nFluxVar;
    int ifn = (icell - NGHOST    ) * nFluxVar;
    double rm, rp, surfm = 1, surfp = 1, vol;
    double Ethermal, presStar;

    printf("--------- icell %d ------ \n", icell);
    if(geometry == 1){
        rm = right_edge[icell - 1];
        rp = right_edge[icell];
        
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

    printf("Density: \n");
    printf("\t u-1\tu\tu+1\n\t %.4e\t%.4e\t%.4e\n", ustate[idxm], ustate[idx], ustate[idxp]);
    printf("\t f-1/2\tf+1/2 -> du\n\t %.4e\t%.4e\t%.4e\n", fluxes[ifn]*surfm, fluxes[ifp]*surfp, (fluxes[ifn]*surfm -  fluxes[ifp]*surfp)*dt/vol);
    printf("Momentum :\n");
    printf("\t u-1\tu\tu+1\n\t %.6e\t%.6e\t%.6e\n", ustate[idxm + 1], ustate[idx + 1], ustate[idxp + 1]);
    printf("\t p-1\tp\tp+1\n\t %.6e\t%.6e\t%.6e\n", pstate[idxm + 1], pstate[idx + 1], pstate[idxp + 1]);
    if(geometry == 1){
        presStar = (surfp*qP[idx + 2] + surfm*qM[idx + 2])/(surfm+surfp);
        double smom = 2.0*presStar*dr[icell]*rs[icell]/vol;
        printf("\t f-1/2\tf+1/2 -> du\n\t %.6e\t%.6e\t%.6e\n", fluxes[ifn + 1]*surfm, fluxes[ifp + 1]*surfp, (fluxes[ifn + 1]*surfm - fluxes[ifp + 1]*surfp)*dt/vol+ smom*dt);
        
    } else {
        printf("\t f-1/2\tf+1/2 -> du\n\t %.6e\t%.6e\t%.6e\n", fluxes[ifn + 1]*surfm, fluxes[ifp + 1]*surfp, (fluxes[ifn + 1]*surfm - fluxes[ifp + 1]*surfp)*dt/vol);
    }
    printf("Energy : \n");
    printf("\t u-1\tu\tu+1\n\t %.4e\t%.4e\t%.4e\n", ustate[idxm + 2], ustate[idx + 2], ustate[idxp + 2]);
    printf("\t p-1\tp\tp+1\n\t %.4e\t%.4e\t%.4e\n", pstate[idxm + 2], pstate[idx + 2], pstate[idxp + 2]);
    printf("\t f-1/2\tf+1/2 -> du\n\t %.4e\t%.4e\t%.4e\n", fluxes[ifn + 2]*surfm, fluxes[ifp + 2]*surfp, (fluxes[ifn + 2]*surfm - fluxes[ifp + 2]*surfp)*dt/vol);

    if(hy_ethresh >= 0.0) {
        // get initial internal temperature
        // Ekin = 0.5*rho v**2 = 0.5 * momx*momx/rho
        Ethermal = ustate[idx + 2] - 0.5*ustate[idx + 1] * ustate[idx + 1] /ustate[idx];
        // predicted half step pressure of at cell center
        presStar = (surfp*qP[idx + 2] + surfm*qM[idx + 2])/(surfm+surfp);
        
        printf("Ethermal : \n");
        printf("\t Ethermal, presStar\t %.4e\t%.4e \n", Ethermal, presStar);
        printf("\t f-1/2\tf+1/2 -> du\n\t %.4e\t%.4e\t%.4e \n", (fluxes[ifn + nvar] + presStar*fluxes[ifn + nvar + 1])*surfm, (fluxes[ifp+nvar] + presStar*fluxes[ifp + nvar])*surfp, fluxes[ifn + nvar]*surfm - fluxes[ifp + nvar]*surfp + presStar*(fluxes[ifn + nvar + 1]*surfm - fluxes[ifp + nvar + 1]*surfp)*dt/vol);

    }
    if(ustate[idx] != ustate[idx]){
        exit(0);
    }
    return 1;

}




int updateHydroVars(int icell, double *fluxes, double *uold, double *qM, double *qP, double surfm, double surfp, double vol, double dt){
    int ifp, ifn, idx, ierr;
    int ivar;
    double unew[3];
    double Ethermal, Ekinetic, presStar;
    
    ifp = (icell - NGHOST + 1)*nFluxVar;
    ifn = (icell - NGHOST)*nFluxVar;
    idx = icell*nvar;
    
    // Update Ethermal before anything changes
    if(hy_ethresh>=0.0 && isothermal == 0){
        // get initial internal temperature
        // Ekin = 0.5*rho v**2 = 0.5 * momx*momx/rho
        Ethermal = uold[2] - 0.5*uold[1] * uold[1] /uold[0];
        // predicted half step pressure of at cell center
        presStar = (surfp*qP[idx + 2] + surfm*qM[idx + 2])/(surfm+surfp);
        
        Ethermal = Ethermal +            (fluxes[ifn + nvar]*surfm     - fluxes[ifp + nvar]*surfp)*dt/vol;
        Ethermal = Ethermal + presStar * (fluxes[ifn + nvar + 1]*surfm - fluxes[ifp + nvar + 1]*surfp)*dt/vol;
        Ethermal = fmax(Ethermal, 1e-40);
    }

    // update density momentum and velocity
    for(ivar = 0; ivar < 3; ivar++){
        if(isothermal && ivar == 2){
            continue;
        }
        unew[ivar] = ustate[idx+ivar]+(fluxes[ifn+ivar]*surfm-fluxes[ifp+ivar]*surfp)*dt/vol;
        if(((ivar==0 || ivar==2) && unew[ivar]<0 ) || unew[ivar] != unew[ivar]){
            printf("\nNEGATIVE IVAR = %d\n",ivar);
            printf("dold=%.4e\n", uold[0]);
            printf("dnew=%.4e\n", unew[0]);
            printf("eold=%.4e\n", uold[2]);
            printf("enew=%.4e\n", unew[2]);
            printf("FL = %.4e FR= %.4e \n",fluxes[ifn+ivar],fluxes[ifp+ivar]);
            printf("FL = %.4e FR= %.4e \n",surfm*fluxes[ifn+ivar]*dt/vol,surfp*fluxes[ifp+ivar]*dt/vol);
            printf("icell = %d / %d\n",icell, NCELLS - NGHOST - 1); 
            printf("rsm = %.4e rs = %.4e rsp = %.4e \n", rs[icell-1], rs[icell], rs[icell+1]);
            printf("nvar = %d \n", nvar);
            return -1;
        }

        ustate[idx + ivar] = unew[ivar];
    }
   
    ierr = addHydroGeometricSourceTerms(icell, qM, qP, dt);
    if(ierr < 0){
        return ierr;
    } 
    // If kinetic energy dominates to the point where thermal energy is lost to machine precision, explicitly add that back in 
    if(isothermal == 0){
        Ekinetic = 0.5*ustate[idx + 1] * ustate[idx + 1] /ustate[idx];
        if(ustate[idx+2] - Ekinetic < Ekinetic*hy_ethresh){
            ustate[idx + 2] = Ekinetic + Ethermal;
        } 
    }
    

    return 1;
}

int updateAdvectedVars(int icell, double *fluxes, double *uold, double *qM, double *qP, double surfm, double surfp, double vol, double dt){
    int ifp, ifn, idx;
    int ivar;
    
    ifp = (icell - NGHOST + 1)*nFluxVar;
    ifn = (icell - NGHOST)*nFluxVar;
    idx = icell*nvar; 

    for(ivar = IADVECT_START; ivar < IADVECT_END; ivar++){
        // phi_new = (phi_old*rho_old + dmass_flux)/rho_new
        ustate[idx+ivar] = (ustate[idx+ivar]*uold[0] + (fluxes[ifn+ivar]*surfm-fluxes[ifp+ivar]*surfp)*dt/vol)/ustate[idx];
    }


    return 1;
}

int doHydroStep(double dt){
    int ierr;
    int icell, idx, ivar;
    double qM[nvar*NCELLS], qP[nvar*NCELLS];
    double fluxes[nFluxVar*NINTER], uold[nvar];
    double rm,rp, surfm = 1, surfp = 1, vol;
    
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
        // calculate left and right riemann states
        ierr = reconstructStates(qP, qM, dt);
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
#ifdef DEBUG_CELL
        ierr = printCellHydro(DEBUG_CELL, 1, fluxes, qP,qM, dt);
#endif
        // Update the conservative variables
        for(icell = NGHOST; icell < NCELLS-NGHOST; icell++){
            idx  = icell*nvar;

            if(geometry == 1){
                rm = right_edge[icell - 1];
                rp = right_edge[icell];
                
                surfm = rm*rm;
                surfp = rp*rp;
                
                vol = (rp*rp*rp -rm*rm*rm)/3.;
            } else {
                vol = dr[icell];
            }

            // Save old state 
            for(ivar = 0;  ivar < nvar; ivar++){
                uold[ivar] = ustate[idx + ivar];
            }

            ierr = updateHydroVars(icell, fluxes, uold, qM, qP, surfm, surfp, vol, dt);
            if(ierr < 0){
                return ierr;
            }

            ierr = updateAdvectedVars(icell, fluxes, uold, qM, qP, surfm, surfp, vol, dt);
            if(ierr < 0){
                return ierr;
            }
#ifdef useDustDynamics
            ierr = updateDustVars(icell, fluxes, uold, surfm, surfp, vol, dt);
            if(ierr < 0){
                return ierr;
            }
#endif

        }
    }
    // convert from their conserved forms to whatever was before
    fixState();
    return 1;
}

