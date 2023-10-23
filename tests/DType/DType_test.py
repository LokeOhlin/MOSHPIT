import numpy as np
import matplotlib.pyplot as plt
import argparse
import h5py as hp
import subprocess
import glob
import os
from scipy.optimize import curve_fit

from moshpit_utils import current_time, set_parameter
import moshpit_utils.units.cgs as cgs

def check_makefile(path):
    with open(path, "r") as f:
        if "moshpit" in f.read():
            return True
        else:
            return False

adi = 5./3.
RST = 0 

meanMol = 0.5

def spitzer(time, Temp):
    cs = np.sqrt(5/3 * cgs.kb * Temp / (meanMol*cgs.mH))

    return RST * (1 + 7/4 * cs * time / RST)**(4/7)
def rfront(tnorm):
    return (1 - np.exp(-tnorm))**(1/3)

def Stromgren(Nphot,alpha,numdens):
    return (3*Nphot/(4*np.pi*alpha))**(1/3)*numdens**(-2/3)
cwd = os.getcwd()

# Parameters
ap=argparse.ArgumentParser()
# Location of moshpit root directory
ap.add_argument('--moshpit_dir', default = cwd+"/../../")
ap.add_argument('--Nion', default = 1e49, type = float)
ap.add_argument('--Eion', default = 2.49912e-11, type = float)
ap.add_argument('--ncells', default = 256, type = int)
ap.add_argument('--numd_init', default = 1e-3, type = float)
ap.add_argument('--rR', default = 5, type = float)
ap.add_argument('--tend', default = 4, type = float)

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
    subprocess.run(["cp",args.moshpit_dir+"/tests/DType/Config", args.moshpit_dir])
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
subprocess.run(["cp", args.moshpit_dir + "tests/DType/simulation.par", "rundir/"])

# move into rundir
os.chdir("rundir")

# Energies and fluxes
fluxes = [0,0,0,0]
energies = [1.121526e-11, 1.922616e-11, 2.243052e-11, 2.563488e-11]
if args.Eion > 2.4353136e-11 :
    energies[-1] = args.Eion
    fluxes[-1] = args.Nion
elif args.Eion > 2.1789648e-11:
    energies[2] = args.Eion
    fluxes[2] = args.Nion
else:
    print("ERROR: supplied Eion is not ionizing")
    sys.exit(0)


# Create SED file
with open("sed.dat", "w") as f:
    for ie in range(len(energies)):
        f.write("%.8e, %.8e\n"%(fluxes[ie], energies[ie]))


# Calculate Strömgren radius and convert max distance
recomb = 2.6e-13
RST = Stromgren(args.Nion, recomb, args.numd_init)
trec = 1/(args.numd_init*recomb)
rR = args.rR*RST
tend = args.tend * trec
dtOut = tend/200


# Set parameters
set_parameter("simulation.par", "rR", rR)
set_parameter("simulation.par", "tend", tend)
set_parameter("simulation.par", "dt_max", dtOut/10)
set_parameter("simulation.par", "dtOut", dtOut)
set_parameter("simulation.par", "dt_init", 1e-3)
set_parameter("simulation.par", "ncells", args.ncells)
set_parameter("simulation.par", "dens_init", args.numd_init)

files = glob.glob("output_*")
if len(files) > 0:
    subprocess.run(["rm", "-f"]+ files)
subprocess.run([args.moshpit_dir + "/moshpit"])
files = sorted(glob.glob("output_*"))

times = np.zeros(len(files))
dts   = np.zeros(len(files))
radIF = np.zeros(len(files))
Tmin  = np.zeros(len(files))
Tavg  = np.zeros(len(files))
Tmax  = np.zeros(len(files))
print(files)
for i, f in enumerate(files):
    
    h5out = hp.File(f, "r")
    abundHe = h5out["Parameters"].attrs.get("ch_abundHe")
    meanMol = (1 + abundHe*4)/(2+2*abundHe)
    times[i]  = current_time(h5out)

    rad  = h5out["coordinates"][:]
    dens = h5out["density"][:]
    temp = h5out["temperature"][:]
    pres = h5out["pressure"][:]
    xH0  = h5out["chemicalAbundances"][:,0]
    xH2  = h5out["chemicalAbundances"][:,1]
    xHp  = h5out["chemicalAbundances"][:,2]
   

    ionReg = xHp > 0.5
    if(np.sum(ionReg)==0):
        continue
     
    abar = 1 # We assume no other species here
    numd = dens/(abar*cgs.mH)
    Tmin   = np.min (temp[ionReg])
    Tmax   = np.max (temp[ionReg])
    Tmean  = np.mean(temp[ionReg])
    Tpres  = temp[pres.argmax()]
    print(Tmin, Tmax, Tmean) 
    dr = rad[0] # constant dr, r[0] = dr
    rm = rad - dr
    rp = rad + dr

    vol = 4*np.pi*(rp**3-rm**3)/3

    Tavg[i] = np.average(temp[ionReg], weights = vol[ionReg]*dens[ionReg])

    
    radIF[i] = np.max(rad[ionReg])
    
    print(np.argmin(np.abs(radIF[i] - rad)))
    h5out.close()

idx = np.abs(radIF - RST).argmin()
plt.plot(times/trec, radIF/RST)
plt.show()
while (radIF[idx] > RST): 
    idx = idx - 1
tm = times[idx]
tp = times[idx + 1]
rm = radIF[idx]
rp = radIF[idx + 1]

if rp == rm:
    tSP = tp
else:
    tSP = tm + (tp - tm) * (RST - rm)/(rp - rm)
print(idx, tm, tp, rm, rp, tSP)
print(rm, rp)

tfit = times[times > tSP] - tSP
Rfit = radIF[times > tSP]
pars = curve_fit(spitzer, tfit, Rfit, 5e3)
Tfit = pars[0][0]
print(Tfit, Tmin, Tmean, Tmax, Tpres)
plt.plot(times/trec,radIF/RST)
tplot = np.linspace(tSP, np.max(times))
#plt.plot(tplot/trec, spitzer(tplot-tSP, Tfit)/RST , ls = '--')
plt.plot(tplot/trec, spitzer(tplot-tSP, Tmin)/RST, ls = '--')
plt.plot(tplot/trec, spitzer(tplot-tSP, Tmean)/RST, ls = '--')
plt.plot(tplot/trec, spitzer(tplot-tSP, Tmax)/RST, ls = '--')
plt.plot(tplot/trec, spitzer(tplot-tSP, Tpres)/RST, ls = '--')
tplot = np.linspace(0, np.max(times)/trec, 1000)
plt.plot(tplot, rfront(tplot), ls = ':')
plt.xlim(0, None)
plt.ylim(0, None)
plt.show()
