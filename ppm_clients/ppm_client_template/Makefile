#----------------------------------------------------------------------
#  Makefile for RC client
#----------------------------------------------------------------------
SHELL       = /bin/sh
TARGET      = ppm_client
BIN         = ppm_client
#----------------------------------------------------------------------
#  Compiler settings (Linux)
#----------------------------------------------------------------------
include Make.inc
DEFINE      =
#DEFINE      = -D__MPI

#----------------------------------------------------------------------
#  Search path for RCS files                                           
#----------------------------------------------------------------------
VPATH = ./RCS

#----------------------------------------------------------------------
#  Additional suffix rules                                             
#----------------------------------------------------------------------
.SUFFIXES : .inc .inc,v .f,v .c,v
.f,v.f :
	 co $*.f

.c,v.c :
	 co $*.c

.inc,v.inc :
	 co $*.inc

#----------------------------------------------------------------------
#  Binary directory
#----------------------------------------------------------------------
bindir      = .

#----------------------------------------------------------------------
#  Default target
#----------------------------------------------------------------------
all: $(TARGET)

#----------------------------------------------------------------------
#  Object files:                                                       
#----------------------------------------------------------------------
OBJS =  client_global.o\
		pwrite.o\
		client_check_abort.o\
		client_create_particles.o\
		client_dump_parameters.o\
		client_finalize.o\
		client_module_io.o\
		client_module_realloc.o\
		client_init.o\
		client_time.o\
		client_write_diag.o\
		client_write_output.o\
		ppm_client.o

#----------------------------------------------------------------------
#  Dependencies:                                                       
#----------------------------------------------------------------------
client_global.o: client_global.f
	$(F77) $(FFLAGS)  -c client_global.f
client_module_realloc.o: client_module_realloc.f client_global.o
	$(CPP) $(DEFINE)  -P client_module_realloc.f __client_module_realloc.f
	$(F77) $(FFLAGS)  -c __client_module_realloc.f -o client_module_realloc.o
client_create_particles.o: client_create_particles.f client_global.o 
	$(CPP) $(DEFINE)  -P client_create_particles.f __client_create_particles.f
	$(F77) $(FFLAGS)  -c __client_create_particles.f -o client_create_particles.o
client_module_io.o: client_module_io.f \
	client_write_output.f \
	client_write_diag.f 
	$(CPP) $(DEFINE)  -P client_module_io.f __client_module_io.f
	$(F77) $(FFLAGS)  -c __client_module_io.f -o client_module_io.o
client_time.o: client_time.f
	$(CPP) $(DEFINE)  -P client_time.f __client_time.f
	$(F77) $(FFLAGS)  -c __client_time.f -o client_time.o
client_init.o: client_init.f client_module_io.f client_global.o
	$(CPP) $(DEFINE)  -P client_init.f __client_init.f
	$(F77) $(FFLAGS)  -c __client_init.f -o client_init.o
client_finalize.o: client_finalize.f client_global.o
	$(CPP) $(DEFINE)  -P client_finalize.f __client_finalize.f
	$(F77) $(FFLAGS)  -c __client_finalize.f -o client_finalize.o
pwrite.o: pwrite.f
	$(F77) $(FFLAGS)  -c pwrite.f
client_dump_parameters.o: client_dump_parameters.f
	$(F77) $(FFLAGS)  -c client_dump_parameters.f
client_check_abort.o: client_check_abort.f
	$(CPP) $(DEFINE)  -P client_check_abort.f __client_check_abort.f
	$(F77) $(FFLAGS)  -c __client_check_abort.f -o client_check_abort.o
ppm_client.o: ppm_client.f client_global.o
	$(CPP) $(DEFINE)  -P ppm_client.f __ppm_client.f
	$(F77) $(FFLAGS)  -c __ppm_client.f -o ppm_client.o

#----------------------------------------------------------------------
#  link                                                                
#----------------------------------------------------------------------
$(TARGET): $(OBJS)
	$(LD) $(OBJS) -o $(BIN) $(LDFLAGS) $(LIBS)

#----------------------------------------------------------------------
#  Install                                                             
#----------------------------------------------------------------------
install: $(BIN)
	(cp -f $(BIN) $(bindir))

#----------------------------------------------------------------------
#  Clean                                                               
#----------------------------------------------------------------------
new: clean ppm_client
clean:
	 rm -f __*.f
	 rm -f $(OBJS)
	 rm -f *.lst
	 rm -f *.l
	 rm -f *.L
	 rm -f *.mod
