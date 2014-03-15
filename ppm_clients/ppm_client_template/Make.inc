
#----------------------------------------------------------------------
# PATHS (need to be adjusted)
#----------------------------------------------------------------------
# PPM Library
PPMDIR     = /home/danny/workspace/deploy/PPM_Core
# METIS Library
METISDIR   = /home/danny/workspace/deploy/Metis
#----------------------------------------------------------------------
#  Compiler settings
#----------------------------------------------------------------------
F77         = mpif90
FC          = mpif90
CC          = mpicc
INC         = -I$(PPMDIR)/include/ppm
#DEBUG       = -C -g
DEBUG       = 
#FFLAGS      = -FR -O3 $(INC) $(DEBUG) 
# change FFLAGS for debugging purpose 
FFLAGS      = -FR $(INC) $(DEBUG)  
CFLAGS      = -O3
LD          = $(F77)
LDFLAGS     = $(FFLAGS) -static
CPP         = cpp
LIBS        = -L$(PPMDIR)/lib64 -lppm -L$(METISDIR) -lmetis 
LDFLAGS     = $(LIBS)


