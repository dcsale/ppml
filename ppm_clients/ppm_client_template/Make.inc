
#----------------------------------------------------------------------
# PATHS (need to be adjusted)
#----------------------------------------------------------------------
PPMDIR = /home/danny/workspace/deploy/PPM_Core
PPMNUMDIR = /home/danny/workspace/deploy/PPM_Numerics
METISDIR = /home/danny/workspace/deploy/Metis

# LIBS =  -L$(PPMDIR)/lib64 -L$(PPMNUMDIR)/lib64 -L$(METISDIR) -lppm -lppmnumerics -lmetis
LIBS =  -L$(PPMDIR)/lib -L$(PPMNUMDIR)/lib -L$(METISDIR) -lppm -lppmnumerics -lmetis
INC = -I$(PPMDIR)/include/ppm -I$(PPMDIR)/src/ -I$(PPMNUMDIR)/include/ppm

#----------------------------------------------------------------------
#  Compiler settings
#----------------------------------------------------------------------
# Build settings (comment/uncomment)
build_Debug = true
# build_Parallel = true

# Compiler flags and settings for my machine
CPP  = cpp
ifdef build_Parallel
	F77 = mpif90
	FC  = mpif90
else
	F77 = gfortran
	FC  = gfortran
endif
ifdef build_Debug
	FFLAGS = $(INC) -ffree-form -ffree-line-length-none -fimplicit-none -g -Wall -fbounds-check -fbacktrace
else
	FFLAGS = $(INC) -ffree-form -ffree-line-length-none -fimplicit-none -O3
endif

CFLAGS      = -O3
LD          = $(F77)
LDFLAGS     = $(FFLAGS) -static
CPP         = cpp
LDFLAGS     = $(LIBS)





