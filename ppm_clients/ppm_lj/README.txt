This example client simulates particles interacting through the Lennard-Jones potential

Compiling the client:
1) Edit the Makefile and put the absolute paths to the ppmcore, ppmnumerics and Metis libraries.
2) make

Running the client:
1) mpirun -n 3 ./lennardjones 
will run the programm on 3 processors using default parameters.
Other parameter values can be supplied through command-line arguments 
Ex: ./lennardjones -N 300 -n 2000 -f 15 
will run the programm with 300 particles during 2000 steps and output data every 15 steps.
(see ./lennardjones --help for the full list of available parameters)

