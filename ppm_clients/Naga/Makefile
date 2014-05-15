# Directory where to move the program when running 'make install'
# DESTDIR=/data/Danny/simulation_output
DESTDIR=.

# Directories where dependencies are installaed
DIR_METIS  = /home/danny/workspace/deploy/Metis
DIR_FFTW   = /home/danny/workspace/deploy/FFTW
DIR_PPM    = /home/danny/workspace/deploy/PPM_Core
DIR_PPMNUM = /home/danny/workspace/deploy/PPM_Numerics

# Build settings (comment/uncomment)
build_Debug = true
build_Local = true
# build_Ubuntu = true

# Compiler flags and settings for my machine
FC   = mpif90
CPP  = cpp

# for some reason, Ubuntu installs to "lib" instead of "lib64" (which most distros use)
ifdef build_Ubuntu
	DIR_LIB = lib
else
	DIR_LIB = lib64
endif

# if building locally, Naga uses it's own version of PPM Numerics FFT Poisson solver -- the PPM Numerics FFT Poisson solver is not yet compatible
ifdef build_Local
	INC  = -I$(DIR_PPM)/include/ppm -I$(DIR_FFTW)/include
	LIBS = -L$(DIR_PPM)/$(DIR_LIB) -lppm -L$(DIR_FFTW)/$(DIR_LIB) -lfftw3 -L$(DIR_METIS) -lmetis
else
	INC  = -I$(DIR_PPM)/include/ppm -I$(DIR_PPMNUM)/include/ppm -I$(DIR_FFTW)/include
	LIBS = -L$(DIR_PPM)/$(DIR_LIB) -lppm -L$(DIR_PPMNUM)/$(DIR_LIB) -lppmnumerics -L$(DIR_FFTW)/$(DIR_LIB) -lfftw3 -L$(DIR_METIS) -lmetis
endif

ifdef build_Debug
	FCFLAGS = $(INC) -ffree-form -ffree-line-length-none -fimplicit-none -g -Wall -fbounds-check -fbacktrace
	# FCFLAGS     = $(INC) -ffree-form -ffree-line-length-none -fimplicit-none -g3
else
	FCFLAGS = $(INC) -ffree-form -ffree-line-length-none -fimplicit-none -O3
endif

OBJ1 = Naga.o\
       naga_time.o\
       naga_say.o\
       naga_globals.o\
       naga_defaults.o\
       naga_read_ctrl.o\
       uppercase.o\
       naga_init.o\
       naga_check_setup.o\
       naga_setup.o\
       naga_allocate_patch.o\
       naga_case_testreprojection.o\
       naga_case_testperiodic.o\
       naga_case_testring.o\
       naga_case_zero.o\
       naga_case_inittaylorgreen.o\
       naga_case_torus.o\
       naga_dump_vtk.o\
       naga_setup_topologies.o\
       naga_get_topologies.o\
       naga_copy_patch.o\
       naga_remesh_particles.o\
       naga_integrate_rk1.o\
       naga_integrate_rk2.o\
       naga_map_particles.o\
       naga_rhs.o\
       naga_catch_particles.o\
       naga_add_uinfinity.o\
       naga_diagnostics.o\
       naga_strainrate.o\
       naga_output.o\
       naga_save_fields.o\
       naga_penalize_explicit.o\
       naga_penalize_interpolation.o\
       naga_penalize_interpolation_update.o\
       naga_structure_sphere.o\
       naga_structure_none.o\
       naga_concentration_sin.o\
       naga_ghost_fields.o\
       naga_extrapolate.o\
       naga_vorticity_rhs.o\
       naga_stl.o\
       naga_stepfunction.o\
       naga_validation.o\
       naga_flowcase.o
ifdef build_Local
OBJS = $(OBJ1) \
       ppm_module_fft.o\
       ppm_module_poisson.o
else
OBJS = $(OBJ1)
endif

.SUFFIXES: .o .f90
.f90.o :
	$(CPP) $(DEFINE) -P   $< -o __$<
	$(FC) $(FCFLAGS) -c __$< -o $@
.f.o :
	$(CPP) $(DEFINE) -P   $< -o __$<
	$(FC) $(FCFLAGS) -c __$< -o $@
%.o : %.mod

Naga: $(OBJS)
	$(FC) -o $@ $(OBJS) $(LIBS)

#------------------------------------------------------------------------------
# Dependencies (additional)
#------------------------------------------------------------------------------
Naga.o :\
  naga_time.o\
  naga_say.o\
  naga_init.o\
  naga_setup.o\
  naga_integrate_rk1.o\
  naga_integrate_rk2.o\
  naga_validation.o

ifdef build_Local
ppm_module_fft.o :\
  ppm_fft_plan_3d_vec_fc2c_z.f\
  ppm_fft_plan_3d_vec_bc2c_z.f\
  ppm_fft_exec_3d_vec_c2c_z.f\
  ppm_fft_plan_3d_vec_fr2c_xy.f\
  ppm_fft_plan_3d_vec_bc2r_xy.f\
  ppm_fft_exec_3d_vec_fr2c_xy.f\
  ppm_fft_exec_3d_vec_bc2r_xy.f\
  ppm_fft_normalize_r.f\
  ppm_fft_normalize_c.f

  #fft/ppm_fft_plan_3d_vec_bc2c_xy.f\
  #fft/ppm_fft_exec_3d_vec_bc2c_xy.f\

ppm_module_poisson.o :\
  ppm_module_fft.o\
  ppm_poisson_init.f\
  ppm_poisson_solve.f\
  ppm_poisson_finalize.f\
  ppm_poisson_fd.f\
  ppm_poisson_extrapolateghost.f
endif

ifdef build_Local
naga_globals.o :\
  ppm_module_fft.o\
  ppm_module_poisson.o
else
naga_globals.o :
endif

naga_init.o :\
  naga_say.o\
  naga_globals.o\
  naga_constants.o\
  naga_defaults.o\
  naga_read_ctrl.o\
  naga_check_setup.o
naga_constants.o :\
  naga_globals.o
naga_defaults.o :\
  naga_globals.o
naga_read_ctrl.o :\
  naga_globals.o\
  naga_say.o\
  uppercase.o

ifdef build_Local
naga_setup.o :\
  naga_globals.o\
  naga_say.o\
  naga_case_test.o\
  naga_case_testreprojection.o\
  naga_case_testperiodic.o\
  naga_case_testring.o\
  naga_case_zero.o\
  naga_case_inittaylorgreen.o\
  naga_case_torus.o\
  naga_allocate_patch.o\
  naga_setup_topologies.o\
  naga_get_topologies.o\
  naga_stl.o\
  naga_stepfunction.o\
  naga_copy_patch.o\
  naga_structure_sphere.o\
  naga_structure_none.o\
  naga_concentration_sin.o\
  ppm_module_fft.o\
  ppm_module_poisson.o
else
naga_setup.o :\
  naga_globals.o\
  naga_say.o\
  naga_case_test.o\
  naga_case_testreprojection.o\
  naga_case_testperiodic.o\
  naga_case_testring.o\
  naga_case_zero.o\
  naga_case_inittaylorgreen.o\
  naga_case_torus.o\
  naga_allocate_patch.o\
  naga_setup_topologies.o\
  naga_get_topologies.o\
  naga_stl.o\
  naga_stepfunction.o\
  naga_copy_patch.o\
  naga_structure_sphere.o\
  naga_structure_none.o\
  naga_concentration_sin.o
endif

naga_check_setup.o :\
  naga_globals.o\
  naga_say.o
naga_allocate_patch.o :\
  naga_globals.o\
  naga_say.o
naga_case_test.o :\
  naga_globals.o\
  naga_say.o
naga_case_testreprojection.o :\
  naga_globals.o\
  naga_say.o
naga_case_testperiodic.o :\
  naga_globals.o\
  naga_say.o
naga_case_testring.o :\
  naga_globals.o\
  naga_say.o
naga_case_zero.o :\
  naga_globals.o\
  naga_say.o
naga_case_inittaylorgreen.o :\
  naga_globals.o\
  naga_say.o
naga_case_torus.o :\
  naga_globals.o\
  naga_say.o
naga_dump_vtk.o :\
  naga_globals.o\
  naga_say.o
naga_setup_topologies.o :\
  naga_globals.o\
  naga_say.o\
  naga_allocate_patch.o\
  naga_remesh_particles.o
naga_get_topologies.o :\
  naga_globals.o\
  naga_say.o
naga_copy_patch.o :\
  naga_globals.o\
  naga_say.o
naga_remesh_particles.o :\
  naga_globals.o\
  naga_say.o
naga_integrate_rk1.o :\
  naga_globals.o\
  naga_say.o\
  naga_rhs.o\
  naga_map_particles.o\
  naga_add_uinfinity.o
naga_integrate_rk2.o :\
  naga_globals.o\
  naga_say.o\
  naga_rhs.o\
  naga_map_particles.o\
  naga_add_uinfinity.o
naga_map_particles.o :\
  naga_globals.o\
  naga_catch_particles.o\
  naga_say.o
naga_rhs.o :\
  naga_globals.o\
  naga_say.o\
  naga_add_uinfinity.o\
  naga_diagnostics.o\
  naga_strainrate.o\
  naga_penalize_explicit.o\
  naga_penalize_interpolation.o\
  naga_penalize_interpolation_update.o\
  naga_save_fields.o\
  naga_ghost_fields.o\
  naga_vorticity_rhs.o\
  naga_flowcase.o
naga_catch_particles.o :\
  naga_globals.o\
  naga_say.o
naga_flowcase.o :\
  naga_globals.o\
  naga_say.o
naga_add_uinfinity.o :\
  naga_globals.o\
  naga_say.o
naga_diagnostics.o :\
  naga_output.o\
  naga_globals.o\
  naga_say.o
naga_strainrate.o :\
  naga_globals.o\
  naga_say.o
naga_save_fields.o :\
  naga_globals.o\
  naga_say.o\
  naga_dump_vtk.o
naga_penalize_interpolation.o :\
  naga_globals.o\
  naga_say.o
naga_penalize_interpolation_update.o :\
  naga_globals.o\
  naga_say.o
naga_penalize_explicit.o :\
  naga_globals.o\
  naga_say.o
naga_structure_sphere.o :\
  naga_globals.o\
  naga_say.o
naga_structure_none.o :\
  naga_globals.o\
  naga_say.o
naga_ghost_fields.o :\
  naga_globals.o\
  naga_extrapolate.o\
  naga_say.o
naga_vorticity_rhs.o :\
  naga_globals.o\
  naga_say.o
naga_extrapolate.o :\
  naga_globals.o\
  naga_say.o
naga_concentration_sin.o :\
  naga_globals.o
naga_stl.o :\
  naga_globals.o\
  naga_stepfunction.o\
  naga_stl_read.f90\
  naga_stl_init.f90\
  naga_stl_inout.f90
naga_stepfunction.o :\
  naga_globals.o\
  naga_say.o
naga_validation.o :\
  naga_globals.o\
  naga_say.o\
  naga_dump_vtk.o
naga_output.o :\
  naga_globals.o\
  naga_say.o
naga_time.o :\
  naga_globals.o\
  naga_say.o

new: veryclean Naga

veryclean: clean
	rm Naga

clean:
	rm -f *.o
	rm -f __*.f90
	rm -f __*.f
	rm -f *.mod
	rm -f ppm_log_0.out
                               
#                        @@@@+   
#                       '@# @@                          
#                       #@   @@                                    @@##             
#                       @@    @#                                 @@@@@@#            
#                       ##    #@                               ;@@'   #@ `          
#                       @+     @@                             '@+      @'           
#                       @:    `#@+                           @@`       @@           
#                       @       @@                          @@`        #@           
#                       @        @#                        @@,         `@           
#                       @        @@'                      @@:           @#          
#                      `@         @@                     @@             #@          
#                      '@         #@                   `@@              '@          
#                      +@          @@                  @@               :@          
#                     `+@          #@@@@@@@@@@@@@@@@@@@@;               `@          
#                      @#                           `##.                 @          
#                      @#                                                @          
#                     @@                                                 @          
#                   #@@                                                  @          
#                 @@@@;                                                  @          
#               @@@@                                                    `@          
#             @@@#                                                      ,@          
#           `@@                                                         #@          
#           @@                                                          @#          
#          #@`                                                          @:          
#          @#                                                           @`          
#         @@                                                            @+          
#         @@`       `                                                   #@'         
#         @        ####`                                                 @@         
#        +@       '##@@#                                                 @@         
#        @@       @.@@@@ `                                               '@         
#        @        `@#@@#.                                                `@+        
#       @@       @ @@@@@.                                                 @@`       
#      #@        +.#@@@.#               .@@@@+`                           #@        
#      @@         ##@@+'               +@,,@@@@                           ,@        
#     `@`          ++;#               :@``#@@@@@                          `@        
#     +@             +                @: ;@@@@@@+                          @+       
#     @@            '.                @. '@@@@@@:                          @#       
#     @@           ,:                 @@ ;@@@@@'                           @@       
#    `@@          :`                 ``@@,@@@@#                            @@       
#    ;@          '`                    `,':;,`                             #@       
#    :@`       +``                                                         #@       
#    :@      +#@@@@@#:                                                     +@       
#    @+     +'@@@@@@@@@;                                                   :@       
#    @.    , @@@@@@@@@@@@                                                  ,@       
#  ``@     # @,@@@@@@@@@@;                                                 :@+      
#   #@     ` ` ;@@@# @@@@                                                   @@      
#   @#    ,    @@@@: '@@@                                                   ;@      
#   @,    ,  #@@@@@@  '+,                                                    @@     
#   @`    '  @@@@@@@+,                                                       #@     
#   @     '  :@@@@@@@@                                                        @     
#  ;@     +   #@@@@@                                                          +@    
#  :@     :     `                                                              @    
#  #@      '                                                                   #'   
#  @@      #                                                                   ;@   
#  @@       `                                                                  `@`  
#  @@       #@##@@@@,                                                           @#  
#  #@        :    ``,;@                                                         @@  
#  `@+       @         .@:                                                      ;@' 
#   @@                    :+  ,#                                                 @@ 
#    @'       #               ``                                                 #@ 
#    @@        #'                                                                `@+
#    #@           :#';'##`                                                        @@
#     @                                                                           ,@
#     @#                                                                           `
#     @@                                                                            
#     #@`                                                                           
#      @@                                                                           
#      @@                                                                           
#      @#                                                                           
#     @@                                                                            
#    #@                                                                             
#    @@                                                                             
#    @+                                                                             
#  `+@                                                                              
#   @@                                                                              
#  @@:                                                                              
# '@#                                                                               
# @@                                                                                
