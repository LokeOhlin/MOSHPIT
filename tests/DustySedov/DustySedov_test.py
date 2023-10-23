from moshpit_utils.test_utils.SedovTaylor import get_solution, get_shockparams
import moshpit_utils.units.cgs as cgs
import moshpit_utils.dust_utils as dust_utils
from moshpit_utils import set_parameter, current_time
from moshpit_utils.python_utils import fmt_10, color_scheme, label_scheme, get_figure_parameters

import h5py as hp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
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
    rhodust = dust_utils.dustDensity(h5file)
    vgas = h5file["velocity"][:]
    vdust = h5file["dustVelocities"][:,0]
    pres = h5file["pressure"][:]
    return pos, rhogas, vgas, pres, rhodust, vdust

def run_simulation(kamp, args, simname):
    # Make sure rundir exists
    if not os.path.exists("rundir/"):
        os.makedirs("rundir/")
    
    # copy default simulation.par
    subprocess.run(["cp", args.moshpit_dir + "tests/DustySedov/simulation.par", "rundir/"])
    
    # move into rundir
    cwd = os.getcwd()
    os.chdir("rundir")
    
    # Set parameters 
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
    files = sorted(glob.glob("output_*"))[1:-1]
    print(files)
    
    for ifile, fname in enumerate(files):
        hfile = hp.File(fname, "r")
        time = current_time(hfile)
        # Simulated values
        pos, rhogas, vgas, pgas, rhodust, vdust = get_simulated(hfile)
        # semi analytic values
        pos_sa, rho_sa, vel_sa, pre_sa = get_solution(time, 1, 1, args.gamma)  
        pos_shock, dennorm, vel_shock, prenorm = get_shockparams(time, 1, 1, args.gamma)

        nrows = 2
        ncols = 1

        figsize, subplots_pars = get_figure_parameters(ncols,nrows, sharex = True)
        fig, axes = plt.subplots(nrows = nrows, ncols = ncols, figsize = figsize, sharex = True)
        plt.subplots_adjust(**subplots_pars)
        axes[0].text(0.05,3.5, r"$t=%.1f$"%time)
        axes[0].plot(pos,rhogas,  c = color_scheme["gas"])
        axes[0].plot(pos_sa,rho_sa, ls = "--", c = "k")
        axes[0].plot(pos, rhodust*100, c = color_scheme["dust"])
        axes[0].set_ylabel(r"$%s,\,100%s $ [code units]"%(label_scheme["rhogas"], label_scheme["rhodust"]))
        axes[0].set_ylim(0.0,4.5)
        axes[0].set_yticks([0.0,1.5,3.0, 4.5])
        

        axes[1].plot(pos, vgas/vel_shock,  c = color_scheme["gas"], label = "Gas")
        axes[1].plot(pos, vdust/vel_shock, c = color_scheme["dust"], label = "Dust")
        axes[1].plot(pos_sa,vel_sa/vel_shock, ls = "--", c = "k", label = "SA")
        axes[1].set_ylabel(r"$%s/U_{sh},\,%s/U_{sh}$"%(label_scheme["vgas"], label_scheme["vdust"]))
         

        axes[1].set_xlabel(r"$r$ [code units]")
        axes[1].set_ylim(0.0, 0.8)
        axes[1].legend()
        axes[1].set_yticks([0,0.2,0.4, 0.6])
        axes[1].set_xlim(0, np.max(pos))
        print(simname+"_"+fname[-4:]+ ".pdf")
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
ap.add_argument('--tmax', default = 0.1, type = float)
ap.add_argument('--dtout', default = 0.1, type = float)
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
    subprocess.run(["cp",args.moshpit_dir+"/tests/DustySedov/Config", args.moshpit_dir])
    # Change to the moshpit root directory and compile
    os.chdir(args.moshpit_dir)
    # Make clean
    subprocess.run(["make", "clean"])
    # Make
    subprocess.run(["make"])
    # Change back to cwd
    os.chdir(cwd)


run_simulation(1, args, "sedov")
