# MOSHPIT <br/> <font size='4'>MOdeling Silicates and Hydrodynamics in PhotoIonized Terrains </font> </br>

MOSHPIT is a dust dynamics and evolution code, meant to simulate the dynamical impact of dust grains inside of HII regions close to star formation. It couples a simple 1D hydrodynamics code to a set of modules that specify the dynamical evolution of dust grains, including the direct coupling to the gas via drag, acceleration due to radiation pressure, and the evolution of the size distribution of the grains via sublimation and sputtering.

## Developers

Loke Lönnblad Ohlin - Main developer (loke.lonnblad@gmail.com)

## Install

### Prerequisites 
C, Fortran77 and Fortran90 compilers (see eg. [GNU](https://gcc.gnu.org/))

[HDF5 Library](https://www.hdfgroup.org/solutions/hdf5/)

### Installation
In the root directory of MOSHPIT, edit the Makefile to fit your installation of compilers and HDF5. Then copy the template_Config into Config, and choose included physics (by commenting or de-commenting specified lines), then install via the command-line
```console
$ make
```
This will create an executable "moshpit".

## Running
To run moshpit, first create the desired executable, then create a parameter file called "simulation.par" in the directory where you wish to run MOSHPIT. In this parameter file write out all the desired parameters to configure the run. Then execute moshpit:
```console
$ /path/to/moshpit
```


## Python utility packages
Inside of the tools directory, there are three packages/scripts. Two of these are to generate look-up tables used when calculating dust absorption coefficients and sputtering yields. The last package, moshpit_utils, is a python package meant to be used to easily read the output hdf5 files, and run tests.

### moshpit-utils installation
the moshpit-utils package is pip installable via running the following command
```
$ pip install -e moshpit_utils
```
inside of the tools directory.