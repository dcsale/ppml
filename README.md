## Parallel Particle-Mesh Library
This repo contains an install script that compiles all the PPM
dependencies and then compiles a few examples of PPM clients.

The Parallel Particle Mesh (PPM) library provides a transparent parallelization middleware for 
particle mesh simulations.
PPM is a Fortran 90 software layer between the Message Passing Interface (MPI) and client 
applications for simulations of physical systems using Particle-Mesh methods. The PPM library 
runs on single and multi-processor architectures, and handles 2D and 3D problems.

The newest PPML is a domain-specific programming language for parallel particle-mesh simulations 
on distributed-memory computers. 
It implements the PPM abstractions and allows rapid prototyping and automatic code generation.


More information about the PPML domain-specific programming language for parallel particle 
simulations, including code examples, and information about the new object-oriented PPM core 
design and its multi-core and GPU support can be found on the [PPM 
website](http://mosaic.mpi-cbg.de/?q=downloads/ppm_lib)

## Ongoing work - May 2014
### Naga PPM client
* adding ADMESH and Shapes codes for manipulation of STL files, this will be used to add moving solid boundaries and interactive analysis
* adding control files for generating 3D wings with NACA profiles
* adding VTK writers from PPM Core v1.2.2 - hopefully this will improve the VisIt / Paraview output capability
* adding python scripting to control VisIt and Paraview sessions - automation of creating figures and analysis
