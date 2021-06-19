#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <dustRadiation.h>
#include <cgeneral.h>
double *Qabs_g_tab_Qabs = NULL;
double *Qabs_g_tab_Qsca = NULL;
double *Qabs_g_tab_gsca = NULL;
double *Qabs_g_tab_agrain = NULL;
double *Qabs_g_tab_freqs = NULL;

double *Qabs_s_tab_Qabs = NULL;
double *Qabs_s_tab_Qsca = NULL;
double *Qabs_s_tab_gsca = NULL;
double *Qabs_s_tab_agrain = NULL;
double *Qabs_s_tab_freqs = NULL;

double *Qem_tab_agrain = NULL;
double *Qem_tab_temp   = NULL;
double *Qem_tab_Qem_g  = NULL;
double *Qem_tab_Qem_s  = NULL;


int g_tab_na;
int g_tab_nf;

int s_tab_na;
int s_tab_nf;

int q_tab_na;
int q_tab_nt;

double Qabs_g_tab_maxAgrain;
double Qabs_g_tab_minAgrain;
double Qabs_g_tab_maxFreq;
double Qabs_g_tab_minFreq;

double Qabs_s_tab_maxAgrain;
double Qabs_s_tab_minAgrain;
double Qabs_s_tab_maxFreq;
double Qabs_s_tab_minFreq;

double Qem_tab_maxAgrain;
double Qem_tab_minAgrain;
double Qem_tab_maxTemp;
double Qem_tab_minTemp;



int loadQabsTable(char *qabsTable, int graphite){
    double Qabs, Qsca, agrain, freq, gsca;
    int iline;
    char line[128];
    char *name, *value;
    FILE *fptr = fopen(qabsTable, "r");
    int tab_na = -1;
    int tab_nf = -1;
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
        if(compStr("nfreqs", name, 80) > 0){
            tab_nf = atoi(value);
        }
    } 
    if(tab_na < 1 || tab_nf < 1){
        printf("Table dimensions not set!");
    }
    if(graphite){
        Qabs_g_tab_Qabs = (double *) malloc(tab_na*tab_nf*sizeof(double));
        Qabs_g_tab_Qsca = (double *) malloc(tab_na*tab_nf*sizeof(double));
        Qabs_g_tab_gsca = (double *) malloc(tab_na*tab_nf*sizeof(double));
        Qabs_g_tab_agrain = (double *) malloc(tab_na*sizeof(double));
        Qabs_g_tab_freqs  = (double *) malloc(tab_nf*sizeof(double));
        g_tab_na = tab_na;
        g_tab_nf = tab_nf;
    
    } else {
        Qabs_s_tab_Qabs = (double *) malloc(tab_na*tab_nf*sizeof(double));
        Qabs_s_tab_Qsca = (double *) malloc(tab_na*tab_nf*sizeof(double));
        Qabs_s_tab_gsca = (double *) malloc(tab_na*tab_nf*sizeof(double));
        Qabs_s_tab_agrain = (double *) malloc(tab_na*sizeof(double));
        Qabs_s_tab_freqs  = (double *) malloc(tab_nf*sizeof(double));
        s_tab_na = tab_na;
        s_tab_nf = tab_nf;
    }
    // rewind, should not be necesserary but looks cleaner
    rewind(fptr);
    int idx = 0;
    int ida, idf;
    for(iline = 0; iline < tab_na*tab_nf; iline++){
        fgets(line, 128, fptr);
        //skip header data
        while(line[0] == '#'){
            fgets(line, 128, fptr);
        }
        idf = idx%tab_nf;
        ida = (idx - idf)/tab_nf;
        sscanf(line, "%lf %lf %lf %lf %lf", &agrain, &freq, &Qabs, &Qsca, &gsca);
        
        if(graphite){
            Qabs_g_tab_agrain[ida] = agrain;
            Qabs_g_tab_freqs[idf]  = freq;
            Qabs_g_tab_Qabs[idx] = Qabs;
            Qabs_g_tab_Qsca[idx] = Qsca;
            Qabs_g_tab_gsca[idx] = gsca;
        } else {
            Qabs_s_tab_agrain[ida] = agrain;
            Qabs_s_tab_freqs[idf]  = freq;
            Qabs_s_tab_Qabs[idx] = Qabs;
            Qabs_s_tab_Qsca[idx] = Qsca;
            Qabs_s_tab_gsca[idx] = gsca;
        }
        idx = idx + 1;
    }
    if(graphite){
        Qabs_g_tab_minAgrain = Qabs_g_tab_agrain[0];
        Qabs_g_tab_maxAgrain = Qabs_g_tab_agrain[tab_na - 1];
        
        Qabs_g_tab_minFreq = Qabs_g_tab_freqs[0];
        Qabs_g_tab_maxFreq = Qabs_g_tab_freqs[tab_nf - 1];
    } else {
        Qabs_s_tab_minAgrain = Qabs_s_tab_agrain[0];
        Qabs_s_tab_maxAgrain = Qabs_s_tab_agrain[tab_na - 1];
        
        Qabs_s_tab_minFreq = Qabs_s_tab_freqs[0];
        Qabs_s_tab_maxFreq = Qabs_s_tab_freqs[tab_nf - 1];
    }
    fclose(fptr);
    return 1;
}

int loadQemTable(char *qemTable){
    double Qem_s, Qem_g, agrain, temp;
    int iline;
    char line[128];
    char *name, *value;
    FILE *fptr = fopen(qemTable, "r");
    int tab_na = -1;
    int tab_nt = -1;
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
    } 
    if(tab_na < 1 || tab_nt < 1){
        printf("Table dimensions not set!");
    }
    Qem_tab_agrain = (double *) malloc(tab_na*sizeof(double));  
    Qem_tab_temp   = (double *) malloc(tab_nt*sizeof(double)); 
    Qem_tab_Qem_g  = (double *) malloc(tab_na*tab_nt*sizeof(double)); 
    Qem_tab_Qem_s  = (double *) malloc(tab_na*tab_nt*sizeof(double));    
    
    q_tab_na = tab_na;
    q_tab_nt = tab_nt;
    // rewind, should not be necesserary but looks cleaner
    rewind(fptr);
    int idx = 0;
    int ida, idt;
    for(iline = 0; iline < tab_na*tab_nt; iline++){
        fgets(line, 128, fptr);
        //skip header data
        while(line[0] == '#'){
            fgets(line, 128, fptr);
        }
        idt = idx%tab_nt;
        ida = (idx - idt)/tab_nt;
        sscanf(line, "%lf %lf %lf %lf ", &agrain, &temp, &Qem_g, &Qem_s);
        Qem_tab_agrain[ida] = agrain;
        Qem_tab_temp[idt]   = temp;
        Qem_tab_Qem_g[idx]  = Qem_g;  
        Qem_tab_Qem_s[idx]  = Qem_s;
      idx = idx + 1;  
    }
    Qem_tab_minAgrain = Qem_tab_agrain[0];
    Qem_tab_maxAgrain = Qem_tab_agrain[tab_na - 1];

    Qem_tab_minTemp = Qem_tab_temp[0];
    Qem_tab_maxTemp = Qem_tab_temp[tab_nt - 1];
    
    fclose(fptr);
    return 1;
}

int loadDustRadiationTables(){
    int ierr;
    printf("Loading graphite data\n");
    char *fname_g = "dustAbsorption_graphite.txt";
    if(access(fname_g, F_OK) != 0){
        printf("ERROR: could not find graphite dust absorption table %s\n", fname_g);
        exit(-1);
    }
    ierr = loadQabsTable(fname_g, 1);
    if(ierr < 0){
        return ierr;
    }
    printf("Loading silicone data\n");
    char *fname_s = "dustAbsorption_silicone.txt";
    if(access(fname_s, F_OK) != 0){
        printf("ERROR: could not find silicate dust absorption table %s\n", fname_s);
        exit(-1);
    }
    ierr = loadQabsTable(fname_s, 0);
    if(ierr < 0){
        return ierr;
    }
    
    printf("Loading <Qabs> data\n");
    char *fname_q = "QemAve.dat";
    if(access(fname_q, F_OK) != 0){
        printf("ERROR: could not find dust emission table %s\n", fname_q);
        exit(-1);
    }
    ierr = loadQemTable(fname_q);
    if(ierr < 0){
        return ierr;
    }
    
    return ierr;
}


// Get corresponding index of photon frequency in the Qabs tables
int getQabs_idf(double freq, int graphite){
    double *freqs;
    int idf;
    int tab_nf;
    double minFreq, maxFreq;
    // point to correct tables
    if(graphite){
        freqs = Qabs_g_tab_freqs;
        tab_nf = g_tab_nf;
        minFreq = Qabs_g_tab_minFreq;
        maxFreq = Qabs_g_tab_maxFreq;
    } else {
        freqs = Qabs_s_tab_freqs;
        tab_nf = s_tab_nf;
        minFreq = Qabs_s_tab_minFreq;
        maxFreq = Qabs_s_tab_maxFreq;
    }
    // if below tabulated values take functional continuation (no need to care for index)
    if (freq < minFreq){
        idf = -1;
    } 
    // if above, use simple linear extrapolation (should probably be done better...)
    if(freq >= maxFreq){
        idf  = tab_nf-2;
    } else {
        idf = binarySearch(freq, freqs, tab_nf);
        //if(freqs[idf] > freq){
        //    idf = idf - 1;
        //}
    }
    return idf;
}
// Get corresponding index of dust grain size in the Qabs tables
int getQabs_ida(double agrain, int graphite){
    double *agrains;
    int ida;
    int tab_na;
    double maxAgrain, minAgrain;
    // point to correct table depending on grain type
    if(graphite){
        agrains = Qabs_g_tab_agrain;
        minAgrain = Qabs_g_tab_minAgrain;
        maxAgrain = Qabs_g_tab_maxAgrain;
        tab_na = g_tab_na;
    } else {
        agrains = Qabs_s_tab_agrain;
        minAgrain = Qabs_s_tab_minAgrain;
        maxAgrain = Qabs_s_tab_maxAgrain;
        tab_na = s_tab_na;
    }

    // if outside range simple linear extrapolation
    if(agrain < minAgrain){
        ida  = 0;
    } else if(agrain >= maxAgrain){
        ida  = tab_na - 2;
    } else{
        ida = binarySearch(agrain, agrains, tab_na);
        //if(agrains[ida] > agrain){
        //    ida = ida - 1;
        //}
    }
    return ida;
}

// Linear interpolation of Qabs
double getQabs_nu(double freq, int graphite, int ida, int *idf){
    double *Qabs;
    double *freqs;
    int idfm, idfp;
    int idx, idxp;
    int tab_nf;
    double minFreq, maxFreq;
    // point to correct tables
    if(graphite){
        freqs = Qabs_g_tab_freqs;
        tab_nf = g_tab_nf;
        Qabs  = Qabs_g_tab_Qabs;
        minFreq = Qabs_g_tab_minFreq;
        maxFreq = Qabs_g_tab_maxFreq;
    } else {
        freqs = Qabs_s_tab_freqs;
        tab_nf = s_tab_nf;
        Qabs  = Qabs_s_tab_Qabs;
        minFreq = Qabs_s_tab_minFreq;
        maxFreq = Qabs_s_tab_maxFreq;
    }
    // if below tabulated values take functional continuation
    if (freq < minFreq){
        if(graphite) {
            return 1.1126500560536185e-23 * freq * freq * Qabs_g_tab_agrain[ida];
        } else {
            return 1.5577100784750660e-23 * freq * freq * Qabs_s_tab_agrain[ida];
        } 
    } 
    if(freq > maxFreq){
        idx = ida*tab_nf + tab_nf -1;
        return Qabs[idx];
    }
    if(*idf < 0){
        idfm  = getQabs_idf(freq, graphite);
        *idf = idfm;
    } else {
        idfm = *idf;
    }
    idfp = idfm + 1; 
    
    idx  = ida*tab_nf + idfm;
    idxp = ida*tab_nf + idfp;

    double fm = freqs[idfm];
    double df = freqs[idfp] - fm;
    
    double Qm = Qabs[idx];
    double dQ = Qabs[idxp] - Qm;

    return Qm + (freq - fm) * dQ/df;
}

double getQabs(double agrain, double freq, int graphite, int ida, int idf){
    /////////////////////////////////////////////////////
    //Method to get Qabs based of interpolation (and extrapolation) of tabulated values
    //
    //In :
    //      double agrain - size of dust grain
    //      double freq   - frequency of absorbed photons
    //      int graphite  - Switch between graphite (1) or siliciate (0) grains
    //      int ida       - Index of dust grain size in table. If <0 we recalculate this with a binary search
    //      int idf       - Index of photon frequency in table. If <0 we recalculate this with a binary search
    //Out:
    //      double Qabs   - Reduced cross section Qabs of dust grain for a given frequency
    //////////////////////////////////////////////////
    double *agrains;
    int idap, idam;
    // point to correct table depending on grain type
    if(graphite){
        agrains = Qabs_g_tab_agrain;
    } else {
        agrains = Qabs_s_tab_agrain;
    }
    if(ida < 0){
        idam  = getQabs_ida(agrain, graphite);
    } else {
        idam = ida;
    }
    idap = idam + 1;
    double am = agrains[idam];
    double da = agrains[idap] - am;

    double Qm = getQabs_nu(freq, graphite, idam, &idf);
    double dQ = getQabs_nu(freq, graphite, idap, &idf) - Qm;

    return Qm + (agrain - am) * dQ/da;
}

// Methods to get index of grain size and temperature in arrays
int getQemAve_ida(double agrain, int graphite){
    int ida;
    if(agrain < Qem_tab_minAgrain){
        ida  = 0;
    } else if(agrain >= Qem_tab_maxAgrain){
        ida  = q_tab_na - 2;
    } else{
        ida = binarySearch(agrain, Qem_tab_agrain, q_tab_na);
    }
    return ida;
}

int getQemAve_idt(double temp, int graphite){
    int idt;
    if(temp >= Qem_tab_maxTemp){
        idt  = q_tab_nt-2;
    } else {
        idt = binarySearch(temp, Qem_tab_temp, q_tab_nt);
    }
    return idt;
}

// Linear interpolation of Qem
double getQemAve_t(double temp, int graphite, int ida, int *idt){
    int idtm, idtp;
    int idxm, idxp;
    // point to correct tables
    // if below tabulated values take as 0 (temperature should be low)
    if (temp < Qem_tab_minTemp){
        return 0;
    } 
    // if above, use simple linear extrapolation (should probably be
    if(temp > Qem_tab_maxTemp){
        idxm = ida*q_tab_nt + q_tab_nt - 1;
        if(graphite){
            return Qem_tab_Qem_g[idxm];
        } else {
            return Qem_tab_Qem_s[idxm];
        }
    }
    if(*idt > 0){
        idtm = *idt;
    } else {
        idtm = getQemAve_idt(temp, graphite);
        *idt = idtm;
    }
    idtp = idtm + 1;

    idxm = ida*q_tab_nt + idtm;
    idxp = ida*q_tab_nt + idtp;

    double tm = Qem_tab_temp[idtm];
    double dt = Qem_tab_temp[idtp] - tm;
    double Qm, dQ;
    if(graphite){ 
        Qm = Qem_tab_Qem_g[idxm];
        dQ = Qem_tab_Qem_g[idxp] - Qm;
    } else {
        Qm = Qem_tab_Qem_s[idxm];
        dQ = Qem_tab_Qem_s[idxp] - Qm;
    }
    return Qm + (temp - tm) * dQ/dt;
}

double getQemAve(double agrain, double temp, int graphite, int ida, int idt){
    /////////////////////////////////////////////////////
    //Method to get Qabs based of interpolation (and extrapolation) of tabulated values
    //
    //In :
    //      double agrain - size of dust grain
    //      double temp   - temperature of dust grains
    //      int graphite  - Switch between graphite (1) or siliciate (0) grains
    //      int ida       - Index of dust grain size in table. If <0 we recalculate this with a binary search
    //      int idt       - Index of dust temperature in table. If <0 we recalculate this with a binary search
    //Out:
    //      double Qem   - Reduced frequency averaged cross section Qem of dust grain for a given temperature
    //////////////////////////////////////////////////
    int idam, idap;
    // if outside range simple linear extrapolation
    if(ida > 0){
        idam = ida;
    } else {
        idam = getQemAve_ida(agrain, graphite);
    }
    idap = idam + 1;
    double am = Qem_tab_agrain[idam];
    double da = Qem_tab_agrain[idap] - am;

    double Qm = getQemAve_t(temp, graphite, idam, &idt);
    double dQ = getQemAve_t(temp, graphite, idap, &idt) - Qm;
    return Qm + (agrain - am) * dQ/da; 
}




