# Compilers
CCOMP = gcc
FCOMP = gfortran

# Libraries
LIBS = -lm -lgfortran -lquadmath -lstdc++
INCS = 

# hdf5
HDF5_DIR = /usr/local/hdf5

LIBS += -L$(HDF5_DIR)/lib -lhdf5
INCS += -I$(HDF5_DIR)/include

# compiler flags
FLAGS = -g 
CFLAGS = -fbounds-check -Wall $(FLAGS)
FFLAGS = -fdefault-real-8 -std=legacy $(FLAGS)

# precompilaton Flags and variables used in code
include ./Config

# add globa dependencies
CDEPS = $(INCS)
FDEPS =


VPATH = ./src/Main:./src/Init:./src/Constants:./src/c_general/:./src/IO:./src/Hydro:./src/RadChem:./src/RadChem/Chemistry:./src/RadChem/Dust:./src/Runtime:./src/MechanicalFeedback



OBJDIR=./bin
# Code stuff

# Each directory determines what object it thinks it needs
CGENOBJ = 
include ./src/c_general/Makefile
CDEPS += -I./src/c_general/

HYDROBJ = 
include ./src/Hydro/Makefile
CDEPS += -I./src/Hydro/

CHEMOBJ = 
include ./src/RadChem/Chemistry/Makefile
FDEPS += -I./src/RadChem/Chemistry/

include ./src/RadChem/Makefile
CDEPS += -I./src/RadChem/

include ./src/RadChem/Dust/Makefile
CDEPS += -I./src/RadChem/Dust/

IOOBJ   = 
include ./src/IO/Makefile
CDEPS += -I./src/IO/

FBOBJ   = 
include ./src/MechanicalFeedback/Makefile
CDEPS += -I./src/MechanicalFeedback/

INITOBJ   = 
include ./src/Init/Makefile
CDEPS += -I./src/Init/

RTOBJ   = 
include ./src/Runtime/Makefile
CDEPS += -I./src/Runtime/

CONSOBJ   = 
include ./src/Constants/Makefile
CDEPS += -I./src/Constants/

MAINOBJ = 
include ./src/Main/Makefile
CDEPS += -I./src/Main/

# Gather modules
OBJ  = 
OBJ += $(CHEMOBJ)
OBJ += $(CGENOBJ)
OBJ += $(MAINOBJ)
OBJ += $(HYDROBJ)
OBJ += $(IOOBJ)
OBJ += $(FBOBJ)
OBJ += $(RTOBJ)
OBJ += $(INITOBJ)
OBJ += $(CONSOBJ)

# Put them in object directory
OBJ := $(OBJ:%.o=$(OBJDIR)/%.o)
CFLAGS += $(CDEPS)
FFLAGS += $(FDEPS)
.SUFFIXES: .c .F90 .F .o 
$(OBJDIR)/%.o: %.F $(OBJDIR)
	$(FCOMP) -c -o $@ $< $(FLAGS) $(FFLAGS)
$(OBJDIR)/%.o: %.F90 $(OBJDIR)
	$(FCOMP) -c -o $@ $< $(FLAGS) $(FFLAGS)
$(OBJDIR)/%.o: %.c $(OBJDIR)
	$(CCOMP) -c -o $@ $< $(FLAGS) $(CFLAGS)
moshpit: $(OBJ)
	$(CCOMP) -o $@ $^ $(FLAGS) $(CFLAGS) $(LIBS) 
$(OBJDIR): 
	mkdir -p $(OBJDIR)
clean :
	rm -f $(OBJDIR)/*.o
