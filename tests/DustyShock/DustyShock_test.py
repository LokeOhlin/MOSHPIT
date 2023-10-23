from moshpit_utils.test_utils.SodShockTube import solve_sodshock, soundspeed_dust
import moshpit_utils.units.cgs as cgs
import moshpit_utils.dust_utils as dust_utils
from moshpit_utils import set_parameter, get_parameter, current_time
from moshpit_utils.python_utils import color_scheme, label_scheme, get_figure_parameters

import h5py as hp
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import subprocess
import argparse
import glob 



def check_makefile(path):
    with open(path, "r") as f:
        if "moshpit" in f.read():
            return True
        else:
            return False

def get_exact(time, P0, P1, dtg, gamma, xplot, x0 = 0.0):
    rpar = (gamma, x0, dtg)

    dens, velo, pres = solve_sodshock(P0, P1, rpar, time, xplot, cs_method = soundspeed_dust)
    
    eint = pres/dens/(gamma - 1)
    return dens, velo, pres, eint

def get_simulated(h5file):
    gamma = h5file["Parameters"].attrs.get("gamma")
    pos = h5file["coordinates"][:]
    rhogas = h5file["density"][:]
    rhodust = dust_utils.dustDensity(h5file)
    vgas = h5file["velocity"][:]
    vdust = h5file["dustVelocities"][:,0]
    pres = h5file["pressure"][:]
    eint = pres/rhogas/(gamma - 1)
    return pos, rhogas, vgas, pres, eint, rhodust, vdust

def run_simulation(kamp, args, simname):
    # Make sure rundir exists
    if not os.path.exists("rundir/"):
        os.makedirs("rundir/")
    
    # copy default simulation.par
    subprocess.run(["cp", args.moshpit_dir + "tests/DustyShock/simulation.par", "rundir/"])
    
    # move into rundir
    cwd = os.getcwd()
    os.chdir("rundir")
    
    # Set parameters 
    P1 = [0.125, 0, 0.1]
    P0 = [1,0,1]
    gamma = args.gamma
    set_parameter("simulation.par", "tend", args.tmax)
    set_parameter("simulation.par", "dtOut", args.dtout)
    set_parameter("simulation.par", "gamma", args.gamma)
    set_parameter("simulation.par", "dust_drag_dtmax_fact", args.drag_dt)
    set_parameter("simulation.par", "dust_drag_const", kamp)
    
    # remove any previous outputs
    files = glob.glob("output_*")
    if len(files) > 0:
        subprocess.run(["rm", "-f"]+ files)
    
    # run
    subprocess.run([args.moshpit_dir + "/moshpit"])
    # process outputs
    files = np.atleast_1d(sorted(glob.glob("output_*"))[-1])
    nplot = 1024
    for ifile, fname in enumerate(files):
        hfile = hp.File(fname, "r")
        dust_to_gas = get_parameter(hfile, "ch_dust_to_gas_ratio")
        time = current_time(hfile)
        # Simulated values
        xplot, rhogas, vgas, pgas, egas, rhodust, vdust = get_simulated(hfile)
        #plot point frequency
        dplot = max(int(len(xplot)/nplot), 1)
        xplot = xplot[::dplot]
        rhogas = rhogas[::dplot]
        vgas = vgas[::dplot]
        pgas = pgas[::dplot]
        egas = egas[::dplot]
        rhodust = rhodust[::dplot]
        vdust = vdust[::dplot]
        # exact values
        dens, velo, pres, eint = get_exact(time, P0, P1, dust_to_gas/100, gamma, xplot)
       
        nrows = 2
        ncols = 2

        figsize, subplots_pars = get_figure_parameters(ncols,nrows, sharex = True)
        fig, axes = plt.subplots(nrows = nrows, ncols = ncols, figsize = figsize, sharex = True)
        plt.subplots_adjust(**subplots_pars)
        axes[0,0].plot(xplot,rhogas,  c = color_scheme["gas"], marker = "o")
        axes[0,0].plot(xplot,rhodust, c = color_scheme["dust"], marker = "x")
        axes[0,0].plot(xplot,dens, c = color_scheme["gas_analytic"])
        axes[0,0].set_ylim(-0.1, 1.1)
        axes[0,0].set_ylabel(r"$%s,\,%s$ [code units]"%(label_scheme["rhogas"],label_scheme["rhodust"]))

        axes[0,1].plot(xplot,vgas,  c = color_scheme["gas"], marker = "o")
        axes[0,1].plot(xplot,vdust, c = color_scheme["dust"], marker = "x")
        axes[0,1].plot(xplot,velo, c = color_scheme["gas_analytic"])
        axes[0,1].set_ylim(-0.1, 0.9)
        axes[0,1].set_ylabel(r"$%s,\,%s$ [code units]"%(label_scheme["vgas"],label_scheme["vdust"]))
        
        axes[1,0].plot(xplot,pgas,  c = color_scheme["gas"], marker = "o")
        axes[1,0].plot(xplot,pres, c = color_scheme["gas_analytic"])
        axes[1,0].set_ylabel(r"$%s$ [code units]"%label_scheme["pres"])
        axes[1,0].set_ylim(-0.1, 1.1)
        
        axes[1,1].plot(xplot,egas,  c = color_scheme["gas"], marker = "o")
        axes[1,1].plot(xplot,eint, c = color_scheme["gas_analytic"])
        axes[1,1].set_xlim(-0.55, 0.55)
        axes[1,1].set_ylim(0.85, 2.3)
        axes[1,1].set_ylabel(r"$%s$ [code units]"%label_scheme["eint"])
        

        axes[0,1].text(-0.4, 0.8, "$t = %.2f$"%time)
        axes[-1,0].set_xlabel(r"$x$ [code units]")
        axes[-1,1].set_xlabel(r"$x$ [code units]")

        print(simname+"_"+fname[-4:]+ ".png")
        plt.savefig(simname+"_"+fname[-4:]+ ".pdf")
        plt.show()
        plt.close(fig)
    os.chdir(cwd)
# set labels and legends

cwd = os.getcwd()

# Parameters
ap=argparse.ArgumentParser()
# Location of moshpit root directory
ap.add_argument('--moshpit_dir', default = cwd+"/../../")
ap.add_argument('--gamma', default = 5/3, type = float)
ap.add_argument('--drag_dt', default = 1e-2, type = float)
ap.add_argument('--tmax', default = 0.2, type = float)
ap.add_argument('--dtout', default = 0.2, type = float)
ap.add_argument('--no_recompile', action = "store_true")
args=ap.parse_args()

#####
# Compile the code
#####
if not args.no_recompile:
    # Check that there is a makefile
    if not os.path.exists(args.moshpit_dir + "/Makefile"):
        sys.exit("Supplied moshpit_dir \"%s\" does not contain a makefile.\n use --moshpit_dir path/to/moshpit_dir/ as an argument to point to the moshpit root directory"%(args.moshpit_dir))
    
    if not check_makefile(args.moshpit_dir + "/Makefile"):
        sys.exit("Makefile \"%s\" does not contain moshpit.\n use --moshpit_dir path/to/moshpit_dir/ as an argument to point to the moshpit root directory"%(args.moshpit_dir + "/Makefile"))
        
    # Move the config file
    subprocess.run(["cp",args.moshpit_dir+"/tests/DustyShock/Config", args.moshpit_dir])
    # Change to the moshpit root directory and compile
    os.chdir(args.moshpit_dir)
    # Make clean
    subprocess.run(["make", "clean"])
    # Make
    subprocess.run(["make"])
    # Change back to cwd
    os.chdir(cwd)


# Transient stage
run_simulation(1, args, "transient")
# Stationary regime 
run_simulation(1000, args, "stationary")
