#-------------------------------------------------------------------------
#
#   ppm example client Makefile
#
#   Author: Birte Schrader and Philipp Boenhof
#
#-------------------------------------------------------------------------
#  Parallel Particle Mesh Library (PPM)
#  ETH Zurich
#  CH-8092 Zurich, Switzerland
#-------------------------------------------------------------------------
PROGRAM  = lennardjones
PROG_SRC = lj.f90
SRCS =	lj_init.f90 \
	lj_set_particles.f90 \
	lj_newtopo.f90 \
	lj_move_particles.f90 \
	lj_remap.f90 \
	lj_change_pos_vel.f90 \
	lj_finalize.f90 \
	$(PROG_SRC)
MODSRCS =	lj_module_global.f90 \
		lj_module_writeout.f90 

# Directory where to move the program when running 'make install'
DESTDIR=.

PPMDIR = /home/danny/workspace/deploy/PPM_Core
PPMNUMDIR = /home/danny/workspace/deploy/PPM_Numerics
METISDIR = /home/danny/workspace/deploy/Metis
# PPMDIR = /home/user/ppmcore/libppm/
# PPMNUMDIR = /home/user/ppmnumerics/libppmnumerics/
# METISDIR = /home/user/metis/intel/

# Compiler flags and settings for my machine
FC   = mpif90
CPP  = cpp
LIBS =  -L$(PPMDIR)/lib64 -L$(PPMNUMDIR)/lib64 -L$(METISDIR) -lppm -lppmnumerics -lmetis
OPT  = $(COPT)
INCL = -I$(PPMDIR)/include/ppm -I$(PPMDIR)/src/ -I$(PPMNUMDIR)/include/ppm

<<<<<<< HEAD
# Flags and settings for ALL machines
COPT= -O
=======
LIBS_METIS    = -L$(DIR_LINK)/Metis -lmetis
# LIBS_FFTW     = -L$(DIR_LINK)/FFTW/lib64 -lfftw3
LIBS_FFTW     = -L$(DIR_LINK)/FFTW/lib -lfftw3
# LIBS_PPM_CORE = -L$(DIR_LINK)/PPM_Core/lib64 -lppm
LIBS_PPM_CORE = -L$(DIR_LINK)/PPM_Core/lib -lppm
# LIBS_PPM_NUM  = -L$(DIR_LINK)/PPM_Numerics/lib64 -lppmnumerics
LIBS_PPM_NUM  = -L$(DIR_LINK)/PPM_Numerics/lib -lppmnumerics
LIBS          = $(LIBS_FFTW) $(LIBS_METIS) $(LIBS_PPM_CORE) $(LIBS_PPM_NUM)
>>>>>>> 17cdcfd95617e420fe3c2e425710c8d9fde1c530

### Usually you should not change these variables below.

# Collecting compilerflags in one variable.
FFLAGS= $(OPT) $(REP) $(COMPATFLAGS)

###### You shouldn't need to change anything below this line #############

# generic variables
OBJS   = $(SRCS:.f90=.o)
MODULES = $(MODSRCS:.f90=.o)
W_TIME = $(shell echo $$(($(TIME)/60)):$$(($(TIME)%60)))

####### Build and Install Targets ########################################
all: $(PROGRAM)

$(PROGRAM): $(MODULES) $(OBJS)
	$(FC) $(FFLAGS) $(MODULES) $(OBJS) -o $@ $(LIBS) $(INCL) 

install: $(PROGRAM) 
	mv $(PROGRAM) $(DESTDIR) 

%.o: %.f90
	$(CPP) -P $*.f90 __$*.f90
	$(FC) -c $(FFLAGS) __$*.f90 -o $*.o $(LIBS) $(INCL)

###### Cleaning up #######################################################
clean:
	$(RM) *.o *.a *.mod __*.f90 core $(OUTF) $(ERRF)

cleanall: clean
	$(RM) $(PROGRAM)
