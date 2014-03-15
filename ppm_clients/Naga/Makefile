DIR_BLD=/home/danny/ppm/ppm_clients/Naga
DIR_LINK=/home/danny/deploy

# some user parameters (comment/uncomment)
build_Parallel = true
# build_Debug    = true
build_Local    = true

ifdef build_Parallel
	FC = mpif90
else
	FC = gfortran
endif
CPP = cpp

ifdef build_Debug
	FCFLAGS     = $(INC) -ffree-form -ffree-line-length-none -fimplicit-none -g -Wall -fbounds-check -fbacktrace
	# FCFLAGS     = $(INC) -ffree-form -ffree-line-length-none -fimplicit-none -g3
else
	FCFLAGS     = $(INC) -ffree-form -ffree-line-length-none -fimplicit-none -O3
endif

INC_PPM_CORE = -I$(DIR_LINK)/PPM_Core/include/ppm
INC_PPM_NUM  = -I$(DIR_LINK)/PPM_Numerics/include/ppm
INC_FFTW     = -I$(DIR_LINK)/FFTW/include
# INC_PPM      = -I$(DIR_BLD)/ppm
# INC_FFT      = -I$(DIR_BLD)/ppm/fft
# INC_POISSON  = -I$(DIR_BLD)/ppm/poisson
# INC          = $(INC_PPM_CORE) $(INC_PPM_NUM) $(INC_POISSON) $(INC_FFT) $(INC_FFTW)
# INC          = $(INC_PPM_CORE) $(INC_FFTW) $(INC_POISSON) $(INC_FFT)
INC          = $(INC_PPM_CORE) $(INC_FFTW)
# INC          = $(INC_PPM_CORE) $(INC_FFTW) $(INC_PPM)
# INC          = $(INC_PPM_CORE) $(INC_FFTW) $(INC_PPM) $(INC_FFT) $(INC_POISSON)
# INC          = $(INC_PPM_CORE) $(INC_PPM_NUM) $(INC_FFTW)
# INC          = $(INC_PPM_CORE) $(INC_FFTW)
# INC          = $(INC_PPM_CORE) $(INC_PPM_NUM)

LIBS_METIS    = -L$(DIR_LINK)/Metis -lmetis
LIBS_PPM_CORE = -L$(DIR_LINK)/PPM_Core/lib64 -lppm
LIBS_PPM_NUM  = -L$(DIR_LINK)/PPM_Numerics/lib64 -lppmnumerics
LIBS_FFTW     = -L$(DIR_LINK)/FFTW/lib64 -lfftw3
# LIBS          = $(LIBS_METIS) $(LIBS_PPM_CORE) $(LIBS_PPM_NUM) $(LIBS_FFTW)
# LIBS          = $(LIBS_PPM_CORE) $(LIBS_FFTW)
LIBS          = $(LIBS_METIS) $(LIBS_PPM_CORE) $(LIBS_FFTW)
# LIBS          = $(LIBS_METIS) $(LIBS_PPM_CORE) $(LIBS_PPM_NUM)
# LIBS          = $(LIBS_PPM_CORE) $(LIBS_PPM_NUM)
# LIBS          = $(LIBS_PPM_CORE)

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
