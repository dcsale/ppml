
#----------------------------------------------------------------------
# PATHS (need to be adjusted)
#----------------------------------------------------------------------
# PPM Library
PPMDIR     = /opt/libppm
# METIS Library
METISDIR   = /opt/metis/intel/lib
#----------------------------------------------------------------------
#  Compiler settings
#----------------------------------------------------------------------
F77         = mpif90
FC          = mpif90
CC          = mpicc
INC         = -I$(PPMDIR)/include
#DEBUG       = -C -g
DEBUG       = 
#FFLAGS      = -FR -O3 $(INC) $(DEBUG) 
# change FFLAGS for debugging purpose 
FFLAGS      = -FR $(INC) $(DEBUG)  
CFLAGS      = -O3
LD          = $(F77)
LDFLAGS     = $(FFLAGS) -static
CPP         = cpp
LIBS        = -L$(PPMDIR)/lib -lppm -L$(METISDIR) -lmetis 
LDFLAGS     = $(LIBS)


