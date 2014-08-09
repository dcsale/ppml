#!/bin/bash
# Description: 
# * A bash script to compile the Parallel-Particle Mesh (PPM) library and dependencies from scratch
# 
# REQUIREMENTS: (these things are usually best installed through your OS package manager)
# * compilers: C, C++, Fortran
# * packages: make, BLAS, LAPACK, curl, OpenMPI
# * other: Ruby, Python

# test: how can you check if a variable has been defined in bash <4.2 and >=4.2?
# in bash >= 4.2 use
# if [[ -v foobar ]]
#
# in bash < 4.2
# declare -p foobar $>/dev/null
#
# or can try (this seems to work just fine in bash 4.1 - on Hyak computer):
#
#foobar=false
#echo $foobar
#if [ -n "${foobar+1}" ] && $foobar
#  then
#  echo "foobar is defined. foobar is true."
#else
#  echo "foobar is not defined as true"
#fi

#install_foobar=true
#echo $install_foobar
#if [ -n "${install_foobar+1}" ] && $install_foobar
#  then
#  echo "foobar will be installed"
#else
#  echo "foobar will NOT be installed"
#fi

# =========================================================================== #
# User Parameters - build & install directories
# =========================================================================== #
DIR_DEPLOY=/home/danny/workspace/deploy_ppml
DIR_BLD=$PWD

# =========================================================================== #
# User Parameters - Specify which software to install (set comment/uncomment)
# =========================================================================== #
# dependencies
# INSTALL_MPI=true
# INSTALL_RUBY=true
# INSTALL_FFTW=true
# INSTALL_METIS=true
# INSTALL_PPMCORE=true
# INSTALL_PPMNUMERICS=true

# PPM clients
# INSTALL_exClient=true
# INSTALL_LJ=true
INSTALL_NAGA=true
# INSTALL_GRAY=true

# =========================================================================== #
# User Parameters - build settings (set true / false)
# =========================================================================== #
buildParallel=true
buildDebug=false
runTests=true

# =========================================================================== #
# User Parameters - source directories for PPM clients
# =========================================================================== #
# easy to include client sources here:
DIR_CLIENT=$DIR_BLD/ppm_clients

SRC_exClient=$DIR_CLIENT/ppm_client_template
SRC_LJ=$DIR_CLIENT/ppm_lj
SRC_NAGA=$DIR_CLIENT/Naga

# =========================================================================== #
# User Parameters - source directories (comment / uncomment)
# =========================================================================== #
# easy to include all external libraries and source here:
DIR_EXT=$DIR_BLD/external

# MPI
SRC_MPI=$DIR_EXT/openmpi-1.6.1
# SRC_MPI=$DIR_EXT/openmpi-1.6.5

# FFTW
# SRC_FFTW=$DIR_EXT/fftw-3.3.2
SRC_FFTW=$DIR_EXT/fftw-3.3.3

# Metis
SRC_METIS=$DIR_EXT/metis-4.0.3
# SRC_METIS=$DIR_EXT/metis4/metis-4.0

# PPM Core
# SRC_PPMCORE=$DIR_EXT/libppm-1.2.1
# SRC_PPMCORE=$DIR_EXT/libppm-1.2.1-naga-debug
SRC_PPMCORE=$DIR_EXT/ppmcore

# PPM Numerics
# SRC_PPMNUMERICS=$DIR_EXT/libppmnumerics-r1036
# SRC_PPMNUMERICS=$DIR_EXT/libppmnumerics-r1036+
SRC_PPMNUMERICS=$DIR_EXT/ppmnumerics

# =========================================================================== #
# Specify the deployment directories
# =========================================================================== #
DIR_DEPLOY_MPI=$DIR_DEPLOY/OpenMPI
DIR_DEPLOY_FFTW=$DIR_DEPLOY/FFTW
DIR_DEPLOY_METIS=$DIR_DEPLOY/Metis
DIR_DEPLOY_PPMCORE=$DIR_DEPLOY/PPM_Core
DIR_DEPLOY_PPMNUMERICS=$DIR_DEPLOY/PPM_Numerics

# =========================================================================== #
# Specify the correct compiler if doing parallel builds (only OpenMPI implemented for now)
# =========================================================================== #
if $buildParallel; then 	
	MPI=openmpi
	c_CC=mpicc
	c_CXX=mpicxx
	c_FC=mpif90
else
	MPI=no 					# leave this as 'no' if not using any MPI
	c_CC=gcc
	c_CXX=g++
	c_FC=gfortran
fi

# initialize all git submodules (like PPM Core and PPM Numerics)
git submodule update --init --recursive

# =========================================================================== #
# echo '*******************************************************************************'
# echo '             _____________________________________________________             '
# echo '            /                                                    \             '
# echo '           |    _____________________________________________     |            '
# echo '           |   |                                             |    |            '
# echo '           |   | $ setup complete...                         |    |            '
# echo '           |   |                                             |    |            '
# echo '           |   | $ let's compile this beast...               |    |            '
# echo '           |   |                                             |    |            '
# echo '           |   |                                             |    |            '
# echo '           |   |                                             |    |            '
# echo '           |   |                                             |    |            '
# echo '           |   |                                             |    |            '
# echo '           |   |                                             |    |            '
# echo '           |   |                                             |    |            '
# echo '           |   |                                             |    |            '
# echo '           |   |                                             |    |            '
# echo '           |   |_____________________________________________|    |            '
# echo '           |                                                      |            '
# echo '            \_____________________________________________________/            '
# echo '                   \_______________________________________/                   '
# echo '                _______________________________________________                '
# echo '             _-'    .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.  --- `-_             '
# echo '          _-'.-.-. .---.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.--.  .-.-.`-_          '
# echo '       _-'.-.-.-.a.-s-.d.-f-.-.-.-.-.-.-.-.-.h.-j-.k.-l-`__`. .-.-.-.`-_       '
# echo '    _-'.-.-.-.-. .-----.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-----. .-.-.-.-.`-_    '
# echo ' _-'.-.-.-.-.-. .---.-. .-----------------------------. .-.---. .---.-.-.-.`-_ '
# echo ':-----------------------------------------------------------------------------:'
# echo '`---._.-----------------------------------------------------------------._.---''
# echo '*******************************************************************************'
# =========================================================================== #

# =========================================================================== #
# Ruby is only needed for the new PPML Fortran 2003+ versions
if [ -n "${INSTALL_RUBY+1}" ] && $INSTALL_RUBY; then
#if [[ -v INSTALL_RUBY ]]; then
	echo '*********************************************'
	echo '* ______      _           '
	echo '* | ___ \    | |          '
	echo '* | |_/ /   _| |__  _   _ '
	echo '* |    / | | | |_ \| | | |'
	echo '* | |\ \ |_| | |_) | |_| |'
	echo '* \_| \_\__,_|_.__/ \__, |'
	echo '*                    __/ |'
	echo '*                   |___/ '
	echo '*'
	echo '*prepare for rubies and gems and yeah...'
	echo '*********************************************'
	curl -L https://get.rvm.io | bash -s stable --ruby
	source ~/.rvm/scripts/rvm
	rvm pkg install libyaml
	rvm reinstall 1.9.3 --with-libyaml-dir=~/.rvm/usr
	rvm reinstall all --force
	rvm use 1.9.3
	gem install cucumber rspec rake thor configatron antlr3 fortran funit
fi
# =========================================================================== #

# =========================================================================== #
# NOTE: I think most clusters will already have MPI installed and configured optimally, probably better to use system MPI instead of building?
if [ -n "${INSTALL_MPI+1}" ] && $INSTALL_MPI; then
#if [[ -v INSTALL_MPI ]]; then
	echo '*********************************************'
	echo '*  _____                 ___  _________ _____ '
	echo '* |  _  |                |  \/  || ___ \_   _|'
	echo '* | | | |_ __   ___ _ __ | .  . || |_/ / | |  '
	echo '* | | | | '_ \ / _ \ '_ \| |\/| ||  __/  | |  '
	echo '* \ \_/ / |_) |  __/ | | | |  | || |    _| |_ '
	echo '*  \___/| .__/ \___|_| |_\_|  |_/\_|    \___/ '
	echo '*       | |                                   '
	echo '*       |_|                                   '
	echo '*'
	echo '*making MPI from: '$SRC_MPI
	echo '*********************************************'
	rm -rf $DIR_DEPLOY_MPI
	mkdir -p $DIR_DEPLOY_MPI
	
	make clean
	cd $SRC_MPI
	./configure --prefix=$DIR_DEPLOY_MPI CC=gcc CPP=cpp CXX=g++ FC=gfortran F90=gfortran
	make
	make install

	# Add environmental variables (by making changes to .bashrc)
	echo 'export PATH='"$DIR_DEPLOY_MPI"'/bin:$PATH' >> ~/.bashrc
	echo 'export LD_LIBRARY_PATH='"$DIR_DEPLOY_MPI"'/lib64:$LD_LIBRARY_PATH' >> ~/.bashrc

	# you might also need to create the symbolic link 
	# ln -s $DIR_DEPLOY_MPI/etc/openmpi-default-hostfile /usr/lib64/mpi/gcc/openmpi/etc/openmpi-default-hostfile

	# if you encounter an error message like
	# ./conftest: error while loading shared libraries: libmpi_cxx.so.1: cannot open shared object file: No such file or directory
	# try the following solutions from:
	# http://stackoverflow.com/questions/11368215/loading-shared-library-in-open-mpi-mpi-run
	# in summary the problem is caused by
	# LD_LIBRARY_PATH is not exported automatically to MPI processes, spawned by mpirun. You should use
	# mpirun -x LD_LIBRARY_PATH ...
	# to push the value of LD_LIBRARY_PATH. Also make sure that the specified path exists on all nodes in the cluster and that libarmadillo.so.3 is available everywhere.
	# or if you get the Open RTE cannot open hostfile error, try:
	# sudo ln -s /etc /usr/lib64/mpi/gcc/openmpi/etc
fi
# =========================================================================== #

# =========================================================================== #
# most clusters will probably have FFTW installed, but we sometimes need to compile and include sources
if [ -n "${INSTALL_FFTW+1}" ] && $INSTALL_FFTW; then
#if [[ -v INSTALL_FFTW ]]; then
	echo '*********************************************'
	echo '*____________ _____ _    _ '
	echo '*|  ___|  ___|_   _| |  | |'
	echo '*| |_  | |_    | | | |  | |'
	echo '*|  _| |  _|   | | | |/\| |'
	echo '*| |   | |     | | \  /\  /'
	echo '*\_|   \_|     \_/  \/  \/ '
	echo '*'
	echo '*making FFTW3 from: '$SRC_FFTW
	echo '*********************************************'
	rm -rf $DIR_DEPLOY_FFTW
	mkdir -p $DIR_DEPLOY_FFTW
	cd $SRC_FFTW
	
	if $buildParallel; then
		# ./configure --prefix=$DIR_DEPLOY_FFTW --exec-prefix=$DIR_DEPLOY_FFTW CC=$c_CC F77=$c_FC --enable-mpi
		./configure --prefix=$DIR_DEPLOY_FFTW CC=$c_CC F77=$c_FC --enable-mpi
	else
		# ./configure --prefix=$DIR_DEPLOY_FFTW --exec-prefix=$DIR_DEPLOY_FFTW CC=$c_CC F77=$c_FC
		./configure --prefix=$DIR_DEPLOY_FFTW CC=$c_CC F77=$c_FC
	fi
	make clean
	make
	# FFTW tests can take a very long time to run
	# if $runTests; then
	# 	make check
	# fi
	make install
	make clean
fi
# =========================================================================== #

# =========================================================================== #
# Metis is used to partition graphs, PPM lib is only compatible with v4.0.3 (or the patched Metis available from PPM devs)
if [ -n "${INSTALL_METIS+1}" ] && $INSTALL_METIS; then
#if [[ -v INSTALL_METIS ]]; then
	echo '*********************************************'
	echo '*  ___  ___     _   _     '
	echo '*  |  \/  |    | | (_)    '
	echo '*  | .  . | ___| |_ _ ___ '
	echo '*  | |\/| |/ _ \ __| / __|'
	echo '*  | |  | |  __/ |_| \__ \'
	echo '*  \_|  |_/\___|\__|_|___/'
	echo '*'
	echo '*making METIS from: '$SRC_METIS
	echo '*********************************************'
	rm -rf $DIR_DEPLOY_METIS
	mkdir -p $DIR_DEPLOY_METIS
	cd $SRC_METIS
	
	make realclean
	make
	cp -r $SRC_METIS/libmetis.a $DIR_DEPLOY_METIS
	# cp -r $SRC_METIS/* $DIR_DEPLOY_METIS
	if $runTests; then
		cd Graphs
		echo $SRC_METIS/Graphs
		../kmetis 4elt.graph 40file
		../onmetis 4elt.graph
		../pmetis test.mgraph 2
		../kmetis test.mgraph 2
		../kmetis test.mgraph 5
		../partnmesh metis.mesh 10
		../partdmesh metis.mesh 10
		../mesh2dual metis.mesh
		../kmetis metis.mesh.dgraph 10
	fi
	make realclean
fi
# =========================================================================== #

# =========================================================================== #
if [ -n "${INSTALL_PPMCORE+1}" ] && $INSTALL_PPMCORE; then
#if [[ -v INSTALL_PPMCORE ]]; then
	echo '*********************************************'
	echo '*_______________  ___      _____ ___________ _____ '
	echo '*| ___ \ ___ \  \/  |     /  __ \  _  | ___ \  ___|'
	echo '*| |_/ / |_/ / .  . |     | /  \/ | | | |_/ / |__  '
	echo '*|  __/|  __/| |\/| |     | |   | | | |    /|  __| '
	echo '*| |   | |   | |  | |     | \__/\ \_/ / |\ \| |___ '
	echo '*\_|   \_|   \_|  |_/      \____/\___/\_| \_\____/ '
	echo '*'
	echo '*making PPM Core Library from: '$SRC_PPMCORE
	echo '*********************************************'
	rm -rf $DIR_DEPLOY_PPMCORE
	mkdir -p $DIR_DEPLOY_PPMCORE
	cd $SRC_PPMCORE

	# for the newest PPML, checkout the correct branch for the client
	# if $branchDevelop; then
		# git checkout develop
	# else
		git checkout master
	# fi

	make clean
	if $buildParallel; then
		if $buildDebug; then
			./configure --prefix=$DIR_DEPLOY_PPMCORE --enable-linux CC=$c_CC CXX=$c_CXX FC=$c_FC LDFLAGS=-L$DIR_DEPLOY_METIS --enable-mpi=$MPI --enable-debug
		else
			./configure --prefix=$DIR_DEPLOY_PPMCORE --enable-linux CC=$c_CC CXX=$c_CXX FC=$c_FC LDFLAGS=-L$DIR_DEPLOY_METIS --enable-mpi=$MPI
		fi
	else
		if $buildDebug; then
			./configure --prefix=$DIR_DEPLOY_PPMCORE --enable-linux CC=$c_CC CXX=$c_CXX FC=$c_FC LDFLAGS=-L$DIR_DEPLOY_METIS --enable-debug
		else
			./configure --prefix=$DIR_DEPLOY_PPMCORE --enable-linux CC=$c_CC CXX=$c_CXX FC=$c_FC LDFLAGS=-L$DIR_DEPLOY_METIS
		fi	
	fi
	make
	if $runTests; then
		# NOTE: the tests will crash unless compiled with MPI
		# make test
		make ftest
		# FUNIT_FLAGS="--procs=1,3" make ftest
	fi
	make install
	# I find that unit testing requires to copy some more files
	cp -r $SRC_PPMCORE/utils $DIR_DEPLOY_PPMCORE/utils
	cp -r $SRC_PPMCORE/src $DIR_DEPLOY_PPMCORE/src

	make clean
fi
# =========================================================================== #

# =========================================================================== #
if [ -n "${INSTALL_PPMNUMERICS+1}" ] && $INSTALL_PPMNUMERICS; then
#if [[ -v INSTALL_PPMNUMERICS ]]; then
	echo '*********************************************'
	echo '*_______________  ___   __   ___   ____  ___ ___________ _____ _____  _____ '
	echo '*| ___ \ ___ \  \/  |   | \ | | | | |  \/  ||  ___| ___ \_   _/  __ \/  ___|'
	echo '*| |_/ / |_/ / .  . |   |  \| | | | | .  . || |__ | |_/ / | | | /  \/\  \__ '
	echo '*|  __/|  __/| |\/| |   | . \`| | | | |\/| ||  __||    /  | | | |     \___ \'
	echo '*| |   | |   | |  | |   | |\  | |_| | |  | || |___| |\ \ _| |_| \__/\/\__/ /'
	echo '*\_|   \_|   \_|  |_/   \_| \_/\___/\_|  |_/\____/\_| \_|\___/ \____/\____/ '
	echo '*'
	echo '*making PPM Numerics Library from: '$SRC_PPMNUMERICS
	echo '*********************************************'
	rm -rf $DIR_DEPLOY_PPMNUMERICS
	mkdir -p $DIR_DEPLOY_PPMNUMERICS
	cd $SRC_PPMNUMERICS

	# for the newest PPML, checkout the correct branch for the client
	# if $branchDevelop; then
		# git checkout develop
	# else
		git checkout master
	# fi

	make clean
	if $buildParallel; then
		if $buildDebug; then
			./configure --prefix=$DIR_DEPLOY_PPMNUMERICS --enable-linux CC=$c_CC CXX=$c_CXX FC=$c_FC --with-ppm=$DIR_DEPLOY_PPMCORE FCFLAGS=-I$DIR_DEPLOY_FFTW/include LDFLAGS=-L$DIR_DEPLOY_FFTW/lib64 --enable-mpi=$MPI --enable-debug
		else
			./configure --prefix=$DIR_DEPLOY_PPMNUMERICS --enable-linux CC=$c_CC CXX=$c_CXX FC=$c_FC --with-ppm=$DIR_DEPLOY_PPMCORE FCFLAGS=-I$DIR_DEPLOY_FFTW/include LDFLAGS=-L$DIR_DEPLOY_FFTW/lib64 --enable-mpi=$MPI
		fi
	else
		if $buildDebug; then
			./configure --prefix=$DIR_DEPLOY_PPMNUMERICS --enable-linux CC=$c_CC CXX=$c_CXX FC=$c_FC --with-ppm=$DIR_DEPLOY_PPMCORE FCFLAGS=-I$DIR_DEPLOY_FFTW/include LDFLAGS=-L$DIR_DEPLOY_FFTW/lib64 --enable-debug
		else
			./configure --prefix=$DIR_DEPLOY_PPMNUMERICS --enable-linux CC=$c_CC CXX=$c_CXX FC=$c_FC --with-ppm=$DIR_DEPLOY_PPMCORE FCFLAGS=-I$DIR_DEPLOY_FFTW/include LDFLAGS=-L$DIR_DEPLOY_FFTW/lib64
		fi	
	fi
	
	make
	if $runTests; then
		make ftest
		# FUNIT_FLAGS="--procs=1,3" make ftest
	fi
	make install
	# I find that unit testing requires to copy some more files
	# cp -r $SRC_PPMNUMERICS/utils $DIR_DEPLOY_PPMNUMERICS/utils
	cp -r $SRC_PPMNUMERICS/src $DIR_DEPLOY_PPMNUMERICS/src

	make clean
fi
# =========================================================================== #


### test some example PPM clients just to make sure everything compiled okay
### Example PPM clients - v1.2.1
# * client template
# * molecular gas dynamics - Lennard-Jones Potential
if [ -n "${INSTALL_LJ+1}" ] && $INSTALL_LJ; then
#if [[ -v INSTALL_LJ ]]; then 
	# make the client
	cd $SRC_LJ
	make clean
	# make
	make 

	# run the client
	if $buildParallel; then 	
		# mpirun -np 1 ./lennardjones --debug 1
		mpirun -np 1 ./lennardjones
	else
		# ./lennardjones
		./lennardjones -N 300 -n 2000 -f 15
	fi
fi
# * diffusion in complex spaces - Particle Strength Exchange
### Example PPM clients - v1.2.2
# PPML Client Generator examples
# * web client for automatic code generation
# * gray-scott client
# if [[ -v INSTALL_GRAY ]]; then
# 
# fi


# =========================================================================== #
if [ -n "${INSTALL_NAGA+1}" ] && $INSTALL_NAGA; then
#if [[ -v INSTALL_NAGA ]]; then
	echo '*********************************************'
	echo '* _   _   ___  _____   ___  '
	echo '*| \ | | / _ \|  __ \ / _ \ '
	echo '*|  \| |/ /_\ \ |  \// /_\ \'
	echo '*| . \ ||  _  | | __ |  _  |'
	echo '*| |\  || | | | |_\ \| | | |'
	echo '*\_| \_/\_| |_/\____/\_| |_/'
	echo '*'
	echo '*making Naga from: '$SRC_NAGA
	echo '*********************************************'

	# make the client
	cd $SRC_NAGA
	make clean
	make

    # copy the compiled Naga to the directory for output
    cp Naga $DIR_DEPLOY

	# run the client
	# if $buildParallel; then 	
	# 	mpirun -np 1 ./Naga CTRL
	# else
	# 	./Naga CTRL
	# fi

	# if you need, try the -x LD_LIBRARY_PATH command to export the LD_LIBRARY_PATH to all connected nodes, like:
	# mpirun -np 1 -x LD_LIBRARY_PATH


	make clean
fi
# =========================================================================== #











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
