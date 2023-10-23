import numpy as np
import matplotlib.pyplot as plt
import argparse
import h5py
import subprocess
import glob
import os

from moshpit_utils.dust_utils import dustDensity
from moshpit_utils import set_parameter, current_time
from moshpit_utils.python_utils import color_scheme, label_scheme, get_figure_parameters

def check_makefile(path):
    with open(path, "r") as f:
        if "moshpit" in f.read():
            return True
        else:
            return False

def get_analytic(time, K, par, vbar, deltaV0, dens, dustDens, mode):
    arg = K*(1/dustDens + 1/dens)*time
    if mode == 1:
        delta = deltaV0 * np.exp(-arg)
    elif mode == 2:
        delta = deltaV0/(1+par*np.abs(deltaV0)**par*arg)**(1/par)
    elif mode == 3:
        delta = deltaV0 * np.exp(-arg) /np.sqrt(1+par*deltaV0**2*(1- np.exp(-2*arg)))
    else:
        delta = (np.sign(deltaV0)/np.sqrt(par))*np.sqrt(((np.sinh(arg) +np.sqrt(1+par*deltaV0**2)*np.cosh(arg)) / 
                                          (np.cosh(arg) +np.sqrt(1+par*deltaV0**2)*np.sinh(arg)))**2 - 1)

    return vbar + dustDens*delta/(dens + dustDens), vbar - dens*delta/(dens + dustDens)

def get_simulated(files):
    vdusts = np.zeros(len(files))
    vgass = np.zeros(len(files))
    times = np.zeros(len(files))
    for i, fname in enumerate(files):
        hfile = h5py.File(fname,"r")
    
        vdusts[i] = hfile["dustVelocities"][0,0]
        vgass[i] = hfile["velocity"][0]
        times[i] = current_time(hfile)

    return times, vgass, vdusts


cwd = os.getcwd()

# Parameters
ap=argparse.ArgumentParser()
# Location of moshpit root directory
ap.add_argument('--moshpit_dir', default = cwd+"/../../")
ap.add_argument('--Kamp', default = [1],nargs = "+", type = float)
ap.add_argument('--drag_par', default = [1e-1],nargs = "+", type = float)
ap.add_argument('--drag_mode', default = [1],nargs = "+", type = int)

ap.add_argument('--drag_dt', default = 1e-2, type = float)
ap.add_argument('--rhog', default = 1, type = float)
ap.add_argument('--rhod', default = 1, type = float)
ap.add_argument('--velg', default = 0, type = float)
ap.add_argument('--veld', default = 1, type = float)


ap.add_argument('--tmax', default =2.0, type = float)
ap.add_argument('--dtout', default = 2.5e-1, type = float)
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
    subprocess.run(["cp",args.moshpit_dir+"/tests/DustyBox/Config", args.moshpit_dir])
    # Change to the moshpit root directory and compile
    os.chdir(args.moshpit_dir)
    # Make clean
    subprocess.run(["make", "clean"])
    # Make
    subprocess.run(["make"])
    # Change back to cwd
    os.chdir(cwd)


# Make sure rundir exists
if not os.path.exists("rundir/"):
    os.makedirs("rundir/")

# Copy simulation file
subprocess.run(["cp", args.moshpit_dir + "tests/DustyBox/simulation.par", "rundir/"])

# move into rundir
os.chdir("rundir")

# Set fixed parameters
dust_to_gas = args.rhod/args.rhog*100
vbar = (args.rhog*args.velg + args.rhod*args.veld)/(args.rhog + args.rhod)
set_parameter("simulation.par", "tend", args.tmax)
set_parameter("simulation.par", "dtOut", args.dtout)
set_parameter("simulation.par", "dt_init", args.dtout)
set_parameter("simulation.par", "dt_max", args.dtout)

set_parameter("simulation.par", "dens_init", args.rhog)
set_parameter("simulation.par", "velo_init", args.velg)
set_parameter("simulation.par", "ch_dust_to_gas_ratio", dust_to_gas)
set_parameter("simulation.par", "dust_velo_init", args.veld)
set_parameter("simulation.par", "dust_drag_dtmax_fact", args.drag_dt)



Ks    = [1,1]
pars  = [1,10]
modes = [1,4] 

nruns = 1

if len(Ks) == 1:
    kamp = Ks[0]
else:
    nruns = len(Ks)
if len(pars) == 1:
    par = pars[0]
else:
    if nruns == 1:
        nruns = len(pars)
    if len(pars) != nruns:
        sys.exit("Kamp, drag_pars and drag_mode must all be of same length or length 1")

if len(modes) == 1:
    mode = modes[0]
else:
    if nruns == 1:
        nruns = len(modes)
    if len(modes) != nruns:
        sys.exit("Kamp, drag_pars and drag_mode must all be of same length or length 1")



line_styles = ["solid", "dotted", "dashed", "dashdot"]
markers = ["o", "x", ".", "+"]


figsize, subplots_pars = get_figure_parameters(1,1)
fig, axes = plt.subplots(nrows = 1, ncols = 1, figsize = figsize)
plt.subplots_adjust(**subplots_pars)
axes = [axes]
tana = np.linspace(0, args.tmax)
for irun in range(nruns):
    if len(Ks) > 1:
        kamp = Ks[irun]
    if len(pars) > 1:
        par = pars[irun]
    if len(modes) > 1:
        mode = modes[irun]

    files = glob.glob("output_*")
    if len(files) > 0:
        subprocess.run(["rm", "-f"]+ files)
    set_parameter("simulation.par", "dust_drag_const", kamp)
    set_parameter("simulation.par", "dust_drag_par", par)
    set_parameter("simulation.par", "dust_drag_mode", mode)
    
    subprocess.run([args.moshpit_dir + "/moshpit"])
    files = glob.glob("output_*")
    files.sort()
    times, vgas_simulated, vdust_simulated = get_simulated(files)
    vgas_analytic, vdust_analytic = get_analytic(tana, kamp, par, vbar, args.velg-args.veld, args.rhog, args.rhod, mode)
    
    axes[0].plot(tana, vgas_analytic, ls = line_styles[irun], c = color_scheme["gas_analytic"])
    axes[0].plot(times, vgas_simulated, ls = "", marker = markers[irun], c =  color_scheme["gas"])
    axes[0].plot(tana, vdust_analytic, ls = line_styles[irun], c = color_scheme["dust_analytic"])
    axes[0].plot(times, vdust_simulated, ls = "", marker = markers[irun], c = color_scheme["dust"])
    #err = np.max((vdust_simulated - vdust_analytic)/vdust_analytic)
    #axes[1].scatter(irun, err, c = P.get_color())    
axes[0].set_xlabel(r"$t$ [code units]")
axes[0].set_ylabel(r"$%s,\,%s$ [code units]"%(label_scheme["vgas"], label_scheme["vdust"]))
axes[0].plot([],[], marker = "o", label  = "Gas", c = color_scheme["gas_analytic"], mfc = color_scheme["gas"], mec = color_scheme["gas"])
axes[0].plot([],[], marker = "o", label  = "Dust", c = color_scheme["dust_analytic"], mfc = color_scheme["dust"], mec = color_scheme["dust"])
axes[0].legend()
plt.savefig("dustyBox.pdf")
plt.show()



# move back
os.chdir(cwd)

