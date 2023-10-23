import matplotlib.pyplot as plt
import matplotlib.colors as mcol
import numpy as np

from moshpit_utils.units import  cgs
import sputtering_utils as utils
# projectile properties
proj_mass  = np.array([1, 4  , 12  , 14  , 16  ])
proj_atom  = np.array([1, 2  , 6   , 7   , 8   ])
proj_abund = np.array([1, 0.1, 1e-4, 1e-5, 3e-4]) 

# carbon
Md_c = 12
Zd_c = 6
K_c  = 0.61
U0_c = 4.0

#silicate
Md_s = 20
Zd_s = 10
K_s  = 0.1
U0_s = 5.8



table_path = "sputtering_yield.dat"
ngrain = 0
nvels  = 0
ntemps = 0
with open(table_path, "r") as f:
    lines = f.readlines()
    for line in lines:
        if line[0] != "#":
            break
        strings = line.split(" ")
        if strings[1] == "ngrain":
            ngrain = int(strings[-1])
        elif strings[1] == "nvels":
            nvels = int(strings[-1])
        elif strings[1] == "ntemps":
            ntemps = int(strings[-1])

print(ngrain, nvels, ntemps)

agrains = np.zeros((ngrain))
temps   = np.zeros((ntemps))
vels_g  = np.zeros((ngrain, nvels))
vels_s  = np.zeros((ngrain, nvels))
sputtering_rate_g = np.zeros((ngrain, nvels, ntemps))
sputtering_rate_s = np.zeros((ngrain, nvels, ntemps))

with open(table_path, "r") as f:
    lines = f.readlines()
    i = 0
    for line in lines:
        if line[0] == "#":
            continue
        it = i%ntemps
        iv = (i//ntemps)%nvels
        ia = (i//ntemps)//nvels
        
        strings = line.split(" ")
        agrains[ia] = float(strings[0])
        vels_g[ia, iv] = float(strings[1])
        vels_s[ia, iv] = float(strings[2])
        temps[it] = float(strings[3])
        sputtering_rate_g[ia, iv, it] = float(strings[4])
        sputtering_rate_s[ia, iv, it] = float(strings[5])
        i = i+1


sputtering_rate_g = sputtering_rate_g * (Md_c*cgs.mH/(2*2.5)) * cgs.yr/1e-4
sputtering_rate_s = sputtering_rate_s * (Md_s*cgs.mH/(2*3.2)) * cgs.yr/1e-4

ncols = 3
color_schemes = ["Blues", "Oranges", "Greens"]


ias = [int(ngrain*0.25), int(ngrain*0.5), int(ngrain*0.75)]
ivs = [int(nvels*0.25), int(nvels*0.5), int(nvels*0.75)]
its = [int(ntemps*0.25), int(ntemps*0.5), int(ntemps*0.75)]


fig, axes = plt.subplots(ncols = 3, nrows = 2)

def get_color(value, vmin, vmax, cmap):
    norm = mcol.LogNorm(vmin = vmin, vmax = vmax)
    return plt.get_cmap(cmap)(norm(value))

# color v, x axis T
for ip in range(ncols):
    ia = ias[ip]
    for iv in range(nvels):
        axes[0,0].plot(temps, sputtering_rate_g[ia, iv, :], c = get_color(vels_g[ia,iv], vels_g[ia,1]/10, vels_g[ia, -1], color_schemes[ip]))
axes[0,0].set_xscale("log")
axes[0,0].set_xlim(temps[0], temps[-1])
axes[0,0].set_yscale("log")
axes[0,0].set_ylim(1e-9, 1e-4)


# color T, x axis v
maxv = 0
for ip in range(ncols):
    ia = ias[ip]
    minvel, maxvel = utils.getVlims(10*agrains[ia], agrains[ia], proj_mass, proj_atom, Md_c, U0_c, 1, minv_fact = 1)
    maxv = max(vels_g[ia,-1]/minvel, maxv)
    for it in range(ntemps):
        axes[1,0].plot(vels_g[ia,:]/minvel, sputtering_rate_g[ia, :, it], c = get_color(temps[it], temps[0]/10, temps[-1], color_schemes[ip]))
axes[1,0].set_xscale("log")
axes[1,0].set_xlim(0.1, maxv)
axes[1,0].set_yscale("log")
axes[1,0].set_ylim(1e-9, 1e-4)

# color v, x axis a
for ip in range(ncols):
    it = its[ip]
    for iv in range(nvels):
        axes[0,1].plot(agrains, sputtering_rate_g[:, iv, it], c = get_color(vels_g[ia,iv], vels_g[ia,1]/10, vels_g[ia, -1], color_schemes[ip]))
axes[0,1].set_xscale("log")
axes[0,1].set_xlim(np.min(agrains), np.max(agrains))
axes[0,1].set_yscale("log")
axes[0,1].set_ylim(1e-9, 1e-4)


# color a, x axis v
maxv = 0
for ip in range(ncols):
    it = its[ip]
    maxv = max(vels_g[ia,-1]/minvel, maxv)
    for ia in range(ngrain):
        minvel, maxvel = utils.getVlims(10*agrains[ia], agrains[ia], proj_mass, proj_atom, Md_c, U0_c, 1, minv_fact = 1)
        axes[1,1].plot(vels_g[ia,:]/minvel, sputtering_rate_g[ia, :, it], c = get_color(agrains[ia], agrains[0], agrains[-1], color_schemes[ip]))
axes[1,1].set_xscale("log")
axes[1,1].set_xlim(0.1, maxv)
axes[1,1].set_yscale("log")
axes[1,1].set_ylim(1e-9, 1e-4)

# color a, x axis T
for ip in range(ncols):
    iv = ivs[ip]
    for ia in range(ngrain):
        axes[0,2].plot(temps, sputtering_rate_g[ia, iv, :], c = get_color(agrains[ia], agrains[0], agrains[-1], color_schemes[ip]))
axes[0,2].set_xscale("log")
axes[0,2].set_xlim(temps[0], temps[-1])
axes[0,2].set_yscale("log")
axes[0,2].set_ylim(1e-9, 1e-4)

# color T, x axis a
for ip in range(ncols):
    iv = ivs[ip]
    for it in range(ntemps):
        axes[1,2].plot(agrains, sputtering_rate_g[:, iv, it], c = get_color(temps[it], temps[0]/10, temps[-1], color_schemes[ip]))
axes[1,2].set_xscale("log")
axes[1,2].set_xlim(np.min(agrains), np.max(agrains))
axes[1,2].set_yscale("log")
axes[1,2].set_ylim(1e-9, 1e-4)
plt.show()



#for ia in range(ngrain):
#    fig, axes = plt.subplots(ncols = 2, nrows = 1)
#    # color v, x axis T
#    for iv in range(nvels):
#        axes[0].plot(temps, sputtering_rate_g[ia, iv, :], c = get_color(vels_g[ia,iv], vels_g[ia,1]/10, vels_g[ia, -1], color_schemes[0]))
#    axes[0].set_xscale("log")
#    axes[0].set_xlim(temps[0], temps[-1])
#    axes[0].set_yscale("log")
#    axes[0].set_ylim(1e-9, 1e-4)
#    
#    
#    # color T, x axis v
#    maxv = 0
#    minvel, maxvel = utils.getVlims(10*agrains[ia], agrains[ia], proj_mass, proj_atom, Md_c, U0_c, 1, minv_fact = 1)
#    maxv = max(vels_g[ia,-1]/minvel, maxv)
#    for it in range(ntemps):
#        axes[1].plot(vels_g[ia,:]/minvel, sputtering_rate_g[ia, :, it], c = get_color(temps[it], temps[0]/10, temps[-1], color_schemes[0]))
#    axes[1].set_xscale("log")
#    axes[1].set_xlim(0.1, maxv)
#    axes[1].set_yscale("log")
#    axes[1].set_ylim(1e-9, 1e-4)
#    plt.show()


for ia in range(ngrain):
    fig, axes = plt.subplots(ncols = 2, nrows = 1)
    # color v, x axis T
    for iv in range(nvels):
        axes[0].plot(temps, sputtering_rate_s[ia, iv, :], c = get_color(vels_g[ia,iv], vels_s[ia,1]/10, vels_s[ia, -1], color_schemes[0]))
    axes[0].set_xscale("log")
    axes[0].set_xlim(temps[0], temps[-1])
    axes[0].set_yscale("log")
    axes[0].set_ylim(1e-13, 1e-8)
    
    
    # color T, x axis v
    maxv = 0
    minvel, maxvel = utils.getVlims(10*agrains[ia], agrains[ia], proj_mass, proj_atom, Md_c, U0_c, 1, minv_fact = 1)
    maxv = max(vels_g[ia,-1]/minvel, maxv)
    for it in range(ntemps):
        axes[1].plot(vels_s[ia,:]/minvel, sputtering_rate_s[ia, :, it], c = get_color(temps[it], temps[0]/10, temps[-1], color_schemes[0]))
    axes[1].set_xscale("log")
    axes[1].set_xlim(0.1, maxv)
    axes[1].set_yscale("log")
    axes[1].set_ylim(1e-13, 1e-8)
    plt.show()
