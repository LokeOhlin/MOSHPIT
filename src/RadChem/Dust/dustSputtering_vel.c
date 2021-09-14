#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <cgeneral.h>
#include <dust.h>
#ifdef useDust
#ifdef sputteringUseVelocities 
double *sput_yield_tab_agrain   = NULL;
double *sput_yield_tab_tgas     = NULL;
double *sput_yield_tab_vdust_g  = NULL;
double *sput_yield_tab_vdust_s  = NULL;
double *sput_yield_tab_yield_g  = NULL;
double *sput_yield_tab_yield_s  = NULL;


int s_tab_na;
int s_tab_nt;
int s_tab_nv;

double sput_yield_tab_maxAgrain;
double sput_yield_tab_minAgrain;
double sput_yield_tab_maxTgas;
double sput_yield_tab_minTgas;



int loadYieldTable(char *yieldTable){
    double yield_s, yield_g, agrain, tgas, vel_g, vel_s;
    int iline;
    char line[128];
    char *name, *value;
    FILE *fptr = fopen(yieldTable, "r");
    int tab_na = -1;
    int tab_nt = -1;
    int tab_nv = -1;
    // read header data
    while(fgets(line, 128, fptr) != NULL){
        if(line[0] != '#'){
            // end if not header
            break;
        }
        sscanf(line, "#%m[^=]=%ms", &name, &value); 
        trim(name); 
        trim(value);

        // find table dimensions
        // number of grain sizes
        if(compStr("ngrain", name, 80) > 0){
            tab_na = atoi(value);
        }
        // number of frequencies
        if(compStr("ntemps", name, 80) > 0){
            tab_nt = atoi(value);
        }
        if(compStr("nvels", name, 80) > 0){
            tab_nv = atoi(value);
        }
    } 
    if(tab_na < 1 || tab_nt < 1){
        printf("Table dimensions not set!");
    }
    sput_yield_tab_agrain = (double *) malloc(tab_na*sizeof(double));  
    sput_yield_tab_tgas   = (double *) malloc(tab_nt*sizeof(double));
    
    sput_yield_tab_vdust_g   = (double *) malloc(tab_na*tab_nv*sizeof(double)); 
    sput_yield_tab_vdust_s   = (double *) malloc(tab_na*tab_nv*sizeof(double)); 
    
    sput_yield_tab_yield_g  = (double *) malloc(tab_na*tab_nt*tab_nv*sizeof(double)); 
    sput_yield_tab_yield_s  = (double *) malloc(tab_na*tab_nt*tab_nv*sizeof(double));    
    
    s_tab_na = tab_na;
    s_tab_nt = tab_nt;
    s_tab_nv = tab_nv;
    // rewind, should not be necesserary but looks cleaner
    rewind(fptr);
    int idx = 0;
    int ida, idt, idv;
    for(iline = 0; iline < tab_na*tab_nv*tab_nt; iline++){
        fgets(line, 128, fptr);
        //skip header data
        while(line[0] == '#'){
            fgets(line, 128, fptr);
        }
        idt = idx%tab_nt;
        idv = (idx/tab_nt)%tab_nv;
        ida = idx/tab_nt/tab_nv;
        sscanf(line, "%lf %lf %lf %lf %lf %lf ", &agrain, &vel_g, &vel_s, &tgas, &yield_g, &yield_s);
        sput_yield_tab_agrain[ida] = agrain;
        sput_yield_tab_tgas[idt]   = tgas;

        sput_yield_tab_vdust_g[ida*tab_nv + idv]   = vel_g;
        sput_yield_tab_vdust_s[ida*tab_nv + idv]   = vel_s;
        
        sput_yield_tab_yield_g[idx]  = yield_g;  
        sput_yield_tab_yield_s[idx]  = yield_s;
      idx = idx + 1;  
    }
    sput_yield_tab_minAgrain = sput_yield_tab_agrain[0];
    sput_yield_tab_maxAgrain = sput_yield_tab_agrain[tab_na - 1];

    sput_yield_tab_minTgas= sput_yield_tab_tgas[0];
    sput_yield_tab_maxTgas= sput_yield_tab_tgas[tab_nt - 1];
    fclose(fptr);
    return 1;
}

int loadDustSputteringTables(){
    int ierr = 1;

    printf("Loading sputtering data\n");
    char *fname_q = "sputtering_yield.dat";
    if(access(fname_q, F_OK) != 0){
        printf("ERROR: could not find sputtering yield table %s\n", fname_q);
        exit(-1);
    }
    ierr = loadYieldTable(fname_q);
    if(ierr < 0){
        return ierr;
    }
    return ierr;
}



// Methods to get index of grain size and temperature in arrays
int getSputYield_ida(double agrain, int graphite){
    int ida;
    if(agrain < sput_yield_tab_minAgrain){
        ida  = 0;
    } else if(agrain >= sput_yield_tab_maxAgrain){
        ida  = s_tab_na - 2;
    } else{
        ida = binarySearch(agrain, sput_yield_tab_agrain, s_tab_na);
    }
    return ida;
}

int getSputYield_idt(double tgas, int graphite){
    int idt;
    if(tgas >= sput_yield_tab_maxTgas){
        idt  = s_tab_nt-2;
    } else {
        idt = binarySearch(tgas, sput_yield_tab_tgas, s_tab_nt);
    }
    return idt;
}

int getSputYield_idv(double vdust, int ida, int graphite) {
    int idv;
    double *vdusts;
    if(graphite){
        vdusts = sput_yield_tab_vdust_g;
    } else {
        vdusts = sput_yield_tab_vdust_s;
    }

    if(vdust >= vdusts[ida*s_tab_nv + s_tab_nv - 1]){
        idv = s_tab_nv - 2;
    } else {
        idv = binarySearch(vdust, vdusts + ida*s_tab_nv, s_tab_nv);
    }
    return idv;

}

// Linear interpolation of sputtering yield
double getSputYield_vel_t(double tgas, int graphite, int ida, int idv, int *idt){
    int idtm, idtp;
    int idxm, idxp;
    // point to correct tables
    // if below tabulated values take as 0 (Should be valid for Tmin < 5000K)
    if (tgas < sput_yield_tab_minTgas){
        return 0;
    } 
    // if above, use simple linear extrapolation (should probably be
    if(tgas > sput_yield_tab_maxTgas){
        idxm = ida*s_tab_nt*s_tab_nv + idv*s_tab_nt + s_tab_nt - 1;
        if(graphite){
            return sput_yield_tab_yield_g[idxm];
        } else {
            return sput_yield_tab_yield_s[idxm];
        }
    }
    if(*idt > 0){
        idtm = *idt;
    } else {
        idtm = getSputYield_idt(tgas, graphite);
        *idt = idtm;
    }
    idtp = idtm + 1;

    idxm = ida*s_tab_nt*s_tab_nv + idv*s_tab_nt + idtm;
    idxp = ida*s_tab_nt*s_tab_nv + idv*s_tab_nt + idtp;

    double tm = sput_yield_tab_tgas[idtm];
    double dt = sput_yield_tab_tgas[idtp] - tm;
    double Ym, dY;
    if(graphite){ 
        Ym = sput_yield_tab_yield_g[idxm];
        dY = sput_yield_tab_yield_g[idxp] - Ym;
    } else {
        Ym = sput_yield_tab_yield_s[idxm];
        dY = sput_yield_tab_yield_s[idxp] - Ym;
    }
    return Ym + (tgas - tm) * dY/dt;
}

double getSputYield_vel_v(double vdust, double tgas, int graphite, int ida, int *idt){
    int idvm, idvp;
    double *vdusts;
    if(graphite){
        vdusts = sput_yield_tab_vdust_g;
    } else {
        vdusts = sput_yield_tab_vdust_s;
    }
    idvm = getSputYield_idv(vdust, ida, graphite);
    idvp = idvm + 1;
    
    // if we are greater or equal to the last velocity, assume 0
    if(idvp + 1 == s_tab_nv){
        printf(".... %.4e %.4e %d %d\n", vdust, vdusts[ida*s_tab_nv + s_tab_nv - 1], idvm, s_tab_nv);
        return 0;
    }
    double vm = vdusts[ida*s_tab_nv + idvm];
    double dv = vdusts[ida*s_tab_nv + idvp] - vm;
    
    double Ym = getSputYield_vel_t(tgas, graphite, ida, idvm, idt);
    double dY = getSputYield_vel_t(tgas, graphite, ida, idvp, idt);

    return Ym + (vdust - vm)* dY/dv;
}


double getSputYield_vel(double agrain, double vdust, double tgas, int graphite, int *ida, int *idt){
    /////////////////////////////////////////////////////
    //Method to get sputtering yields based of interpolation (and extrapolation) of tabulated values
    //
    //In :
    //      double agrain - size of dust grain
    //      double vdust  - relative velocity between dust and gas
    //      double tgas   - temperature of gas
    //      int graphite  - Switch between graphite (1) or siliciate (0) grains
    //      int ida       - Pointer to index of dust grain size in table. If <0 we recalculate this with a binary search
    //      int idt       - Pointer Index of dust temperature in table. If <0 we recalculate this with a binary search
    //Out:
    //      double yield  - <Yv> where Y is the sputtering yeild and the average is over the maxwellian distribution for a given temperature
    //////////////////////////////////////////////////
    int idam, idap;
    // if outside range simple linear extrapolation
    if(*ida > 0){
        idam = *ida;
    } else {
        idam = getSputYield_ida(agrain, graphite);
        *ida = idam;
    }
    idap = idam + 1;
    double am = sput_yield_tab_agrain[idam];
    double da = sput_yield_tab_agrain[idap] - am;

    double Ym = getSputYield_vel_v(vdust, tgas, graphite, idam, idt);
    double dY = getSputYield_vel_v(vdust, tgas, graphite, idap, idt) - Ym;
    return Ym + (agrain - am) * dY/da; 
}


double get_dadt_sputtering_vel(double agrain, int ibin, int iabin, int graphite, double numd, double vdust, double tgas, int *idt){
    double prefact;
    double yieldv = getSputYield_vel(agrain, vdust, tgas, graphite, &ida_tabSput[ibin], idt);
    
    if(graphite) {
        prefact = aveMatom_c/(2*rho_c);   // 12 mH
    } else {
        prefact = aveMatom_s/(2*rho_s);   // 20 mH  
    }
    return -yieldv*prefact*numd;
}
#endif
#endif
