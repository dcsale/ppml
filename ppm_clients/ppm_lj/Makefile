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


#----------------------------------------------------------------------
#  User parameters (set directories or true/false)
#----------------------------------------------------------------------
DIR_LINK = /home/danny/workspace/deploy
build_Parallel = true
# build_Debug    = true

#----------------------------------------------------------------------
#  Compiler settings
#----------------------------------------------------------------------
INC_FFTW     = -I$(DIR_LINK)/FFTW/include
# INC_PPM_CORE = -I$(DIR_LINK)/PPM_Core/include/ppm
INC_PPM_CORE = -I$(DIR_LINK)/PPM_Core/include/ppm -I$(DIR_LINK)/PPM_Core/src
INC_PPM_NUM  = -I$(DIR_LINK)/PPM_Numerics/include/ppm 
INC          = $(INC_FFTW) $(INC_PPM_CORE) $(INC_PPM_NUM)

LIBS_FFTW     = -L$(DIR_LINK)/FFTW/lib64 -lfftw3
LIBS_METIS    = -L$(DIR_LINK)/Metis -lmetis
LIBS_PPM_CORE = -L$(DIR_LINK)/PPM_Core/lib64 -lppm
LIBS_PPM_NUM  = -L$(DIR_LINK)/PPM_Numerics/lib64 -lppmnumerics
LIBS          = $(LIBS_FFTW) $(LIBS_METIS) $(LIBS_PPM_CORE) $(LIBS_PPM_NUM)

ifdef build_Parallel
	DEFINE = -D__MPI
	FC = mpif90
else
	FC = gfortran
endif
CPP = cpp

ifdef build_Debug
	FFLAGS = $(INC) -ffree-form -ffree-line-length-none -fimplicit-none -g -Wall -fbounds-check -fbacktrace	
else
	FFLAGS = $(INC) -ffree-form -ffree-line-length-none -fimplicit-none -O3
endif

###### You shouldn't need to change anything below this line #############

# generic variables
OBJS    = $(SRCS:.f90=.o)
MODULES = $(MODSRCS:.f90=.o)
W_TIME  = $(shell echo $$(($(TIME)/60)):$$(($(TIME)%60)))

####### Build and Install Targets ########################################
all: $(PROGRAM)

$(PROGRAM): $(MODULES) $(OBJS)
	$(FC) $(FFLAGS) $(MODULES) $(OBJS) -o $@ $(LIBS) $(INC) 

install: $(PROGRAM) 
	mv $(PROGRAM) $(DESTDIR) 

%.o: %.f90
	$(CPP) -P $*.f90 __$*.f90
	$(FC) -c $(FFLAGS) __$*.f90 -o $*.o $(LIBS) $(INC)

###### Cleaning up #######################################################
clean:
	$(RM) *.o *.a *.mod __*.f90 core $(OUTF) $(ERRF)

cleanall: clean
	$(RM) $(PROGRAM)
