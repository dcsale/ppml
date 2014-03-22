./objects/ppm_module_poisson.o objects/ppm_module_poisson.d: src/ppm_module_poisson.f \
 src/poisson/ppm_poisson_init.f src/poisson/ppm_poisson_solve.f \
 src/poisson/ppm_poisson_fd.f src/poisson/ppm_poisson_extrapolateghost.f\
# end of source dependencies for .o and .d files
./objects/ppm_module_poisson.o: \
./objects/ppm_module_fft.o \
./objects/ppm_module_substart.o \
./objects/ppm_module_substop.o \
./objects/ppm_module_write.o \
./objects/ppm_module_data.o \
./objects/ppm_module_mktopo.o \
./objects/ppm_module_topo_get.o \
./objects/ppm_module_mesh_define.o \
./objects/ppm_module_map_field.o \
./objects/ppm_module_map_field_global.o \
./objects/ppm_module_map.o \
./objects/ppm_module_map_field.o \
./objects/ppm_module_map_field_global.o \
./objects/ppm_module_map.o \
./objects/ppm_module_topo_get.o \
./objects/ppm_module_topo_get.o \
# end of module dependencies for .o file
