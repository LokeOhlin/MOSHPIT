import numpy as np
import matplotlib.pyplot as plt
import argparse
import h5py as hp
import subprocess
import glob
import os

from moshpit_utils.dust_utils import dustDensity
from moshpit_utils import set_parameter, get_parameter, current_time
from moshpit_utils.python_utils import color_scheme, label_scheme, get_figure_parameters

def check_makefile(path):
    with open(path, "r") as f:
        if "moshpit" in f.read():
            return True
        else:
            return False

def numdist(a, Ni, Si, da, ac):
    return Ni/da + Si * (a - ac)

def get_analytic(a, time, amin, amax, plaw, dadt, tscale, rho_grain, mtot):
    norm = (plaw + 1)/(amax**(plaw+1) - amin**(plaw+1))
    mfirst = 4*np.pi *  rho_grain/3 * norm * (amax**(plaw+4)-amin**(plaw+4))/(plaw+4)
    amax = amax - dadt * tscale * (np.cos(2*np.pi*time/tscale)-1)/(2*np.pi)
    amin = amin - dadt * tscale * (np.cos(2*np.pi*time/tscale)-1)/(2*np.pi)
    norm = (plaw + 1)/(amax**(plaw+1) - amin**(plaw+1))
    norm = norm * mtot / mfirst
    
    dist = np.zeros(a.shape)
    dist[ (a<amin) ] = 0
    dist[ (a>amax) ] = 0
    mask = (a>amin)*(a<amax) 
    dist[mask] = norm*a[mask]**plaw
    return dist

def get_min_max(amin, amax, dadt, tscale):
    maximum = amax + dadt * tscale / np.pi
    minimum = amin 

    return minimum, maximum


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
ap.add_argument('--tmax', default =2.0, type = float)
ap.add_argument('--dtout', default = 0.5, type = float)
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
    subprocess.run(["cp",args.moshpit_dir+"/tests/DustSinGrowth/Config", args.moshpit_dir])
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
subprocess.run(["cp", args.moshpit_dir + "tests/DustSinGrowth/simulation.par", "rundir/"])

# move into rundir
os.chdir("rundir")

# Set fixed parameters
set_parameter("simulation.par", "tend", args.tmax)
set_parameter("simulation.par", "dtOut", args.dtout)

# remove any previous outputs
files = glob.glob("output_*")
if len(files) > 0:
    subprocess.run(["rm", "-f"]+ files)

# run
subprocess.run([args.moshpit_dir + "/moshpit"])

files = glob.glob("output_*")
files.sort()
files = files[:-1]
ncols = 1
nrows = len(files)
figsize, subplots_pars = get_figure_parameters(2.5,nrows, sharex = True, sharey =True, ax_height_ratio= 1.5/5)
fig, axes = plt.subplots(nrows = nrows, ncols = ncols, figsize = figsize, sharex = True, sharey =True)
plt.subplots_adjust(**subplots_pars)




for ifile, fname in enumerate(files):
    hfile = hp.File(fname,"r")
    amin   = get_parameter(hfile, "dust_amin")
    amax   = get_parameter(hfile, "dust_amax")
    cmin   = get_parameter(hfile, "dust_cmin")
    cmax   = get_parameter(hfile, "dust_cmax")
    dadt   = get_parameter(hfile, "dust_dadt_c")
    tscale = get_parameter(hfile, "dust_dadt_tscale")
    
    dist_min, dist_max = get_min_max(cmin, cmax, dadt, tscale)
    time = current_time(hfile)
    dust_numd = hfile["number"][:,:]
    dust_slope = hfile["slope"][:,:]
    dust_density = hfile["density"][0]*get_parameter(hfile, "ch_dust_to_gas_ratio")/100
    abin_c = hfile["agrain"][:]
    abin_e = hfile["agrain_binEdges"][:]
    nabins = dust_numd.shape[1]
    aplot = np.zeros(3)
    nplot = np.zeros(3)

    a_analytic = np.logspace(np.log10(amin), np.log10(amax), 300)
    analytic = get_analytic(a_analytic, time, cmin, cmax, 0, dadt, tscale, 2.26, dust_density)
    norm = np.max(analytic)
    #print(dust_density)
    for ia in range(nabins):
        aplot[0] = abin_e[ia]
        aplot[1] = abin_c[ia]
        aplot[2] = abin_e[ia+1]

        nplot[:] = numdist(aplot, dust_numd[0,ia], dust_slope[0,ia], abin_e[ia+1] - abin_e[ia], abin_c[ia])
        axes[ifile].plot(aplot/1e-4, nplot/norm, c = color_scheme["dust"])
        axes[ifile].plot(aplot[1]/1e-4, nplot[1]/norm, ls = "", marker = "o", c = color_scheme["dust"])
    

    axes[ifile].plot(a_analytic/1e-4, analytic/norm, ls = ":", c = "k")
    axes[ifile].text(1.5*dist_min/2.5/1e-4, 0.8, r"$t=%.1fT$"%time)

axes[-1].set_xlim(dist_min/2.5/1e-4, dist_max*2.5/1e-4)
axes[-1].set_xscale("log")
axes[-1].set_xlabel(r"$a$ [$\mu$m]")
axes[len(files)//2].set_ylabel(r"$\partial n/\partial a$ [normalized]")
plt.savefig("dustSinGrowth.pdf")
plt.show()
os.chdir(cwd)


