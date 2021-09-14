#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <init.h>
#include <hydro.h>
#include <IO.h>
#include <rtpars.h>
#ifdef useChemistry
#include <radchem.h>
#endif
#ifdef useDust
#include <dust.h>
#endif
#include <mechanicalFeedback.h>
#include <moshpit.h>

double time, dt;

int mainLoop(){
    int istep = 0, io=0, iostep=0, ierr;
    int doOutput;
    double dt_new, dt_chem=1e99, dt_feedback, ioTime=0.0;
    time=t0;
    dt = dt_init;
    // First output
    ierr = createOutputFile(io);
    ierr = finalizeOutputFile();
    io = io + 1;
    iostep = iostep + nstepOut;
    ioTime = ioTime + dtOut;

    while( time < tend ){
        // check timestep
        ierr = getCFL(&dt_new);
        if(dt_new < dt){
            dt = dt_new;
        } else if (fmin(dt_new,dt_chem) > 2*dt){
            dt = 2*dt;
        } else{
            dt = fmin(dt_new, dt_chem);
        }
        if(dt > dt_max){
            dt = dt_max;
        }
        printf("step %d time = %.4e dt = %.4e || dt_CFL = %.4e  dt_chem= %.4e\n", istep, time ,dt, dt_new,dt_chem);
        // check if we are outputing on this step
        if((istep + 1 >= iostep) || (time + dt >= ioTime)){
            // we create file here since some modules (dust) does not store much of the data, so must therefore write during the step
            ierr = createOutputFile(io);
            doOutput = 1;
        } else {
            doOutput = 0;
        }
        
        //Hydro
        ierr = doHydroStep(dt);
        if(ierr < 0){
            printf("error in hydro\n");
            return -1;
        }
#ifdef useChemistry 
        // Check if we want to write dust data
#ifdef useDust
        if(doOutput){
            outputDust = 1;    
        } else {
            outputDust = 0;
        }

#endif       
        // Chemistry
        ierr = doChemistryStep(dt, &dt_chem);
        if(ierr < 0){
            return -1;
        }
#endif
        // Feedback
        ierr = doFeedbackStep(dt, &dt_feedback);
        if(ierr < 0){
            return -1;
        }
        
        // finalize the output
        if(doOutput){
            ierr = finalizeOutputFile();
            iostep = istep + nstepOut;
            ioTime = time + dtOut;
            printf("next output at time = %.4e or step = %d\n",ioTime, iostep);
            io = io + 1;
        }
        
        istep = istep + 1;
        time = time+dt;
        if(istep > imax){
            printf("\nMAXIMUM STEPS REACHED\n");
            break;
        }

        //update time steps
#ifdef useChemistry
        if(dt_chem<dt){
            dt = dt_chem;
        }

#endif
        if(dt_feedback < dt){
            dt = dt_feedback;
        }

    }
    // final output
    ierr = createOutputFile(io);
    ierr = finalizeOutputFile();
    io = io + 1;
    return 1;
}



int main(){
    int ierr;
    printf("\t\tMY HYDRO !\n"); 
    printf("\ncalling initSimulation\n");
    ierr = initSimulation(); 
    if(ierr < 0){
        return -1;
    }
    printf("Done\n"); 
   
     
    printf("\ncalling init_domain\n");
    ierr = init_domain();
    if(ierr < 0){
        return -1;
    }
    printf("done\n"); 

    printf("\nEntering main loop\n");
    ierr = mainLoop();
    if(ierr < 0){
        return -1;
    }
    printf("done\n"); 
    
}
    
