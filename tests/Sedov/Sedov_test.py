from moshpit_utils.test_utils.SedovTaylor import get_solution, get_position, get_shockparams
import moshpit_utils.units.cgs as cgs
from moshpit_utils import set_parameter, current_time

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


def get_simulated(h5file):
    pos = h5file["coordinates"][:]
    rhogas = h5file["density"][:]
    vgas = h5file["velocity"][:]
    pres = h5file["pressure"][:]
    return pos, rhogas, vgas, pres

def run_simulation(kamp, args, simname):
    # Make sure rundir exists
    if not os.path.exists("rundir/"):
        os.makedirs("rundir/")
    
    # copy default simulation.par
    subprocess.run(["cp", args.moshpit_dir + "tests/Sedov/simulation.par", "rundir/"])
    
    # move into rundir
    cwd = os.getcwd()
    os.chdir("rundir")
    
    # Set parameters 
    gamma = args.gamma
    set_parameter("simulation.par", "tend", args.tmax)
    set_parameter("simulation.par", "dtOut", args.dtout)
    set_parameter("simulation.par", "gamma", args.gamma)
    set_parameter("simulation.par", "ncells", args.ncells)
    
    if args.logspace:
        set_parameter("simulation.par", "logspace_cells", 1)
        set_parameter("simulation.par", "rL", 1e-3)
    else:
        set_parameter("simulation.par","logspace_cells", 0)
        set_parameter("simulation.par", "rL", 0)

    
    # remove any previous outputs
    files = glob.glob("output_*")
    if len(files) > 0:
        subprocess.run(["rm", "-f"]+ files)
    
    # run
    subprocess.run([args.moshpit_dir + "/moshpit"])
    # process outputs
    files = sorted(glob.glob("output_*"))[1:]
    print(files)
    
    for ifile, fname in enumerate(files):
        hfile = hp.File(fname, "r")
        time = current_time(hfile)
        # Simulated values
        pos, rhogas, vgas, pgas = get_simulated(hfile)
        # semi analytic values
        pos_sa, rho_sa, vel_sa, pre_sa = get_solution(time, 1, 1, args.gamma)

        rad_shock, dennorm, vel_shock, prenorm = get_shockparams(time, 1, 1, args.gamma)
        pos_shock = get_position(time, 1, 1)

        fig, axes = plt.subplots(nrows = 3, ncols = 1, sharex = True, figsize =(3.46, 9) )
        axes[0].plot(pos/rad_shock,rhogas/dennorm,  c = "navy")
        axes[0].plot(pos_sa/rad_shock,rho_sa/dennorm, c = "k", ls = "--")
        axes[0].set_ylabel(r"$\rho / \rho_\mathrm{2}$")
        axes[0].set_ylim(0,None)

        axes[1].plot(pos/rad_shock,vgas/vel_shock,  c = "navy")
        axes[1].plot(pos_sa/rad_shock,vel_sa/vel_shock, c = "k", ls = "--")
        axes[1].set_ylabel(r"$u/u_\mathrm{sh}$")
        axes[1].set_ylim(0,None)
        
        axes[2].plot(pos/rad_shock,pgas/prenorm,  c = "navy")
        axes[2].plot(pos_sa/rad_shock,pre_sa/prenorm, c = "k", ls = "--")
        axes[2].set_ylabel(r"$P/P_\mathrm{sh}$")
        axes[2].set_ylim(0,None)
        
        

        axes[2].set_xlabel(r"$r/r_\mathrm{sh}$")
        axes[2].set_xlim(0,1.2) 
        

        print(simname+"_"+fname[-4:]+ ".png")
        plt.subplots_adjust(hspace = 0.025)
        for ax in axes:
            labels = ax.get_yticklabels()
            labels[0] = ""
            ax.set_yticklabels(labels)

        plt.savefig(simname+"_"+fname[-4:]+ ".png")
        plt.show()
        plt.close(fig)
    os.chdir(cwd)
# set labels and legends

cwd = os.getcwd()

# Parameters
ap=argparse.ArgumentParser()
# Location of moshpit root directory
ap.add_argument('--moshpit_dir', default = cwd+"/../../")
ap.add_argument('--logspace', action = "store_true")
ap.add_argument('--ncells', default = 512, type = int)
ap.add_argument('--gamma', default = 5/3, type = float)
ap.add_argument('--tmax', default = 0.1, type = float)
ap.add_argument('--dtout', default = 0.01, type = float)
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
    subprocess.run(["cp",args.moshpit_dir+"/tests/Sedov/Config", args.moshpit_dir])
    # Change to the moshpit root directory and compile
    os.chdir(args.moshpit_dir)
    # Make clean
    subprocess.run(["make", "clean"])
    # Make
    subprocess.run(["make"])
    # Change back to cwd
    os.chdir(cwd)


run_simulation(1, args, "sedov")
