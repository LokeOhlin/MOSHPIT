import moshpit_utils.test_utils.DustyWave.dustywaves as dustywaves
import moshpit_utils.units.cgs as cgs
import moshpit_utils.dust_utils as dust_utils
from moshpit_utils import set_parameter, current_time
from moshpit_utils.python_utils import color_scheme, label_scheme, get_figure_parameters

import h5py as hp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
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
def get_exact(time, ampl, cs, Kdragin, lam, x0, rhogeq, rhodeq, xplot):
    vgaso = np.zeros(xplot.shape)
    vdusto = np.zeros(xplot.shape)
    rhogaso = np.zeros(xplot.shape)
    rhodusto = np.zeros(xplot.shape)
    ierr = 0
    vgaso,vdusto,rhogaso,rhodusto,ierr = dustywaves.exact_dustywave(time,ampl,cs,Kdragin,lam,x0,rhogeq,rhodeq, xplot)
    return vgaso, vdusto, rhogaso, rhodusto

def get_simulated(h5file):
    pos = h5file["coordinates"][:]
    rhogas = h5file["density"][:]
    rhodust = dust_utils.dustDensity(h5file)
    vgas = h5file["velocity"][:]
    vdust = h5file["dustVelocities"][:,0]
    return pos, vgas, vdust, rhogas, rhodust

cwd = os.getcwd()

# Parameters
ap=argparse.ArgumentParser()
# Location of moshpit root directory
ap.add_argument('--moshpit_dir', default = cwd+"/../../")
ap.add_argument('--Kamp', default = 1, type = float)
ap.add_argument('--drag_dt', default = 1e-2, type = float)
ap.add_argument('--rhog', default = 1, type = float)
ap.add_argument('--rhod', default = 1, type = float)
ap.add_argument('--cs'  , default = 1, type = float)
ap.add_argument('--lam' , default = 1, type = float)
ap.add_argument('--amp' , default = 1e-4, type = float)
ap.add_argument('--xmax' , default = 1, type = float)
ap.add_argument('--tmax', default = 2, type = float)
ap.add_argument('--dtout', default = 1, type = float)
ap.add_argument('--no_recompile', action = "store_true")
ap.add_argument('--no_rerun', action = "store_true")
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
    subprocess.run(["cp",args.moshpit_dir+"/tests/DustyWave/Config", args.moshpit_dir])
    # Change to the moshpit root directory and compile
    os.chdir(args.moshpit_dir)
    # Make clean
    subprocess.run(["make", "clean"])
    # Make
    subprocess.run(["make"])
    # Change back to cwd
    os.chdir(cwd)

    
#unpack all parameters
Kamp = args.Kamp
drag_dt  = args.drag_dt
rhog = args.rhog
rhod = args.rhod
cs = args.cs
lam = args.lam
amp = args.amp

xmax = args.xmax
tmax = args.tmax
dtout = args.dtout

# Get dust to gas ratio (normalized to 1%)
dust_to_gas = rhod/rhog * 100
# calculate constant
drag_const = Kamp
# Make sure rundir exists
if not os.path.exists("rundir/"):
    os.makedirs("rundir/")


# move into rundir
os.chdir("rundir")
    
if not args.no_rerun:
    # copy default simulation.par
    subprocess.run(["cp", args.moshpit_dir + "tests/DustyWave/simulation.par", "./"])
    # Set parameters 
    set_parameter("simulation.par", "tend", tmax)
    set_parameter("simulation.par", "dtOut", dtout)
    set_parameter("simulation.par", "dens_init", rhog)
    set_parameter("simulation.par", "cs_init", cs)
    set_parameter("simulation.par", "ch_dust_to_gas_ratio", dust_to_gas)
    set_parameter("simulation.par", "dust_drag_const", drag_const)
    
    set_parameter("simulation.par", "wave_vel_amplitude", amp)
    set_parameter("simulation.par", "wave_rho_amplitude", amp)
    set_parameter("simulation.par", "wave_vel_wavelength", lam)
    set_parameter("simulation.par", "wave_rho_wavelength", lam)
    
    set_parameter("simulation.par", "rR", xmax)
    
    # remove any previous outputs
    files = glob.glob("output_*")
    if len(files) > 0:
        subprocess.run(["rm", "-f"]+ files)
    
    # run
    subprocess.run([args.moshpit_dir + "/moshpit"])

# process outputs
files = sorted(glob.glob("output_*"))
files = files[:-1]
nfiles = len(files)
if nfiles <= 3:
    nrows = nfiles
    ncols = 1
else:
    nrows = 3
    ncols = 1


figsize, subplots_pars = get_figure_parameters(ncols,nrows, sharex = True, sharey =True)
fig, axes = plt.subplots(nrows = nrows, ncols = ncols, figsize = figsize, sharex = True, sharey =True)
plt.subplots_adjust(**subplots_pars)
#if len(axes.shape) == 1:
#    axes = [axes]

# number of markers
nmarks = 25
for ifile, fname in enumerate(files):
    hfile = hp.File(fname, "r")
    time = current_time(hfile)
    # Simulated values
    xplot, vgas_s, vdust_s, rhogas_s, rhodust_s = get_simulated(hfile)
    # exact values
    vgas_e, vdust_e, rhogas_e, rhodust_e = get_exact(time, amp, cs, Kamp, lam, 0, rhog, rhod, xplot)
    
    dmark = len(xplot)//nmarks
    # Velocity
    axes[ifile].plot(xplot, vgas_e/1e-4, c = color_scheme["gas_analytic"])
    axes[ifile].plot(xplot[::dmark], vgas_s[::dmark]/1e-4, c = color_scheme["gas"],ls = "", marker = "x")
    axes[ifile].plot(xplot, vdust_e/1e-4, c = color_scheme["dust_analytic"])
    axes[ifile].plot(xplot[::dmark], vdust_s[::dmark]/1e-4, c = color_scheme["dust"], ls = "", marker = "x")
    axes[ifile].text(0.05, -0.8, "time = %d"%int(time))
   
    # Density
    #axes[ifile,1].plot(xplot[::dmark], rhogas_s[::dmark], c = color_scheme["gas"], marker = "x")
    #axes[ifile,1].plot(xplot, rhogas_e, c = color_scheme["gas_analytic"])
    #axes[ifile,1].plot(xplot[::dmark], rhodust_s[::dmark], c = color_scheme["dust"], marker = "x")
    #axes[ifile,1].plot(xplot, rhodust_e, c = color_scheme["dust_analytic"])


    #axes[ifile,1].set_ylim(rhog-amp*1.05, rhog+amp*1.05)

# set labels and legends
# xlabel
axes[-1].set_ylim(-1.15, 1.15)
axes[-1].set_xlabel(r"$x$ [code units]")

# ylabel
axes[1].set_ylabel(r"$%s/A,\,%s/A$ [code units]"%(label_scheme["vgas"], label_scheme["vdust"]))
#axes[1,0 ].set_ylabel("density [code units]")

# legend
axes[0].plot([],[], marker = "x", c = color_scheme["gas_analytic"], mfc = color_scheme["gas"], mec = color_scheme["gas"], label = "Gas")
axes[0].plot([],[], marker = "x", c = color_scheme["dust_analytic"], mfc = color_scheme["dust"], mec = color_scheme["dust"], label = "Dust")
axes[0].legend()
plt.savefig("DustyWave.pdf")
plt.show()

#for ifile, fname in enumerate(files):
#    hfile = hp.File(fname, "r")
#    time = current_time(hfile)
#    # Simulated values
#    xplot, vgas_s, vdust_s, rhogas_s, rhodust_s = get_simulated(hfile)
#    # exact values
#    vgas_e, vdust_e, rhogas_e, rhodust_e = get_exact(time, amp, cs, Kamp, lam, 0, rhog, rhod, xplot)
#
#    fig, axes = plt.subplots(figsize = (6,4), nrows = 1, ncols = 2, sharex = True)
#    print(ifile, iax, fname)
#    dmark = len(xplot)//nmarks
#    # Velocity
#    axes[0].plot(xplot[::dmark], vgas_s[::dmark], c = "navy", marker = "x")
#    axes[0].plot(xplot, vgas_e, c = "navy")
#    axes[0].plot(xplot[::dmark], vdust_s[::dmark], c = "orange", marker = "x")
#    axes[0].plot(xplot, vdust_e, c = "orange")
#    
#    # Density
#    axes[1].plot(xplot[::dmark], rhogas_s[::dmark], c = "navy", marker = "x")
#    axes[1].plot(xplot, rhogas_e, c = "navy")
#    axes[1].plot(xplot[::dmark], rhodust_s[::dmark], c = "orange", marker = "x")
#    axes[1].plot(xplot, rhodust_e, c = "orange")
#
#    axes[0].text(np.max(xplot)*0.05, max(np.max(vgas_e), np.max(vgas_s))* 0.9, "t = %.4e"%time)
#
#    axes[0].set_ylim(-amp*1.05, amp*1.05)
#    axes[1].set_ylim(rhog-amp*1.05, rhog+amp*1.05)
#
#    plt.savefig("plot_"+fname[-4:]+".png")
#    plt.close(fig)
