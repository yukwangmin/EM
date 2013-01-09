TARGET=PEM


CC = g++


CXX=mpicxx
#CXX=mpicxx -pedantic -fopenmp
#CXX=mpicxx -pedantic
BGCXX=mpixlcxx_r
#BGCXX=mpixlcxx_r -O3 -qarch=450d -qtune=450 -qnosmp
#         -qarch=<suboption>
#                Specifies the general processor architecture for which the code
#                (instructions) should be generated.  For any given -qarch
#                setting, the compiler defaults to a specific, matching -qtune
#                setting, which can provide additional performance improvements.
#                The suboptions are:
#
#                440
#                     Generates code for a single floating-point unit (FPU) only,
#                     on a PowerPC 440 hardware platform.
#                440d
#                     Generates parallel instructions for the 440d Double Hummer
#                     dual FPU. Note that if you encounter problems with code
#                     generation, try resetting this option to -qarch=440.
#                450
#                     Generates code for a single floating-point unit (FPU) only,
#                     on a PowerPC 450 hardware platform.
#                450d
#                     Generates parallel instructions for the 450d Double Hummer
#                     dual FPU. Note that if you encounter problems with code
#                     generation, try resetting this option to -qarch=450.
#
#                Default: -qarch=450d

#         -O[<level>]
#                Optimizes code at a choice of levels during compilation. This is
#                equivalent to -qoptimize[=<level>]. <level> can be:
#
#                0
#                     Performs only quick local optimizations such as constant
#                     folding and elimination of local common subexpressions.
#                2
#                     Performs optimizations that the compiler developers
#                     considered the best combination for compilation speed and
#                     runtime performance. The optimizations may change from
#                     product release to release.
#                3
#                     Performs some memory and compile-time intensive
#                     optimizations in addition to those executed with -O2. The
#                     -O3 specific optimizations have the potential to alter the
#                     semantics of a program. The compiler guards against these
#                     optimizations at -O2 and the option -qstrict is provided at
#                     -O3 to turn off these aggressive optimizations.
#                     Specifying -O3 implies -qhot=level=0.
#                4
#                     This option is the same as -O3, but also:
#                       o sets the -qarch and -qtune options to the architecture
#                       of the compiling machine unless -qnoautoconfig is
#                       specified. -qnoautoconfig is set by default on Blue
#                       Gene/P.
#                       o sets the -qcache option most appropriate to the
#                       characteristics of the compiling unless -qnoautoconfig is
#                       specified. -qnoautoconfig is set by default on Blue
#                       Gene/P.
#                       o sets the -qipa option.
#                       o sets the -qhot option to level=1.
#                5
#                     Equivalent to -O4 -qipa=level=2.

#         -qtune={440|450}
#                Specifies the architecture system for which the executable
#                program is optimized.
#
#                440
#                     Optimizes object code for the 440 family of processors.
#                     This is the default for -qarch=440 and -qarch=440d.
#                450
#                     Optimizes object code for the 450 family of processors.
#                     This is the default for -qarch=450 and -qarch=450d.







ifeq ($(TARGET),BGPEM)

XERCES_ROOT = ~/lib/xerces-c-3.1.1
#XERCES_ROOT = ~/lib/xerces-c-3.1.1-powerpc-aix-xlc-7.0
# export LIBPATH=$LIBPATH:~/lib/xerces-c-3.1.1-powerpc-aix-xlc-7.0/lib

MPARSER_ROOT = ~/projects/mparser


CFLAGS =  -O3 -qarch=450d -qtune=450 -qsmp=omp
LFLAGS = 

INCLUDE = -I ./ -I $(XERCES_ROOT)/src -I $(MPARSER_ROOT)
LIBS = -L ./ -L $(XERCES_ROOT)/src/.libs -L $(MPARSER_ROOT) -lxerces-c -lmparser

XERCES_DYNAMIC_LINK = 


else


XERCES_ROOT = ~/lib/xerces-c-3.1.1-x86-linux-gcc-3.4
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/lib/xerces-c-3.1.1-x86-linux-gcc-3.4/lib

MPARSER_ROOT = ~/projects/mparser

#Note: Change -O to -g when using debugger.
#CFLAGS = -g 
#CFLAGS = -g -Wall
CFLAGS = -O -fopenmp
#CFLAGS = -O -fopenmp -D RELEASE 

LFLAGS = 
#LFLAGS = -static

INCLUDE = -I ./ -I $(XERCES_ROOT)/include -I $(MPARSER_ROOT)
LIBS = -L ./ -L $(XERCES_ROOT)/lib -L $(MPARSER_ROOT) -lxerces-c -lmparser


#XERCES_DYNAMIC_LINK = -Wl,-rpath,$(XERCES_ROOT)/lib
XERCES_DYNAMIC_LINK = 


endif





#.PHONY: all clean distclean 
all: ${TARGET}

SEM: serial_em.cpp serial_solver.cpp time_controller.cpp resultview.cpp serial_mover.cpp particle.cpp em_error.cpp initial_setting.cpp user_specific_initial_particle_setting.cpp solver.h serial_solver.h time_controller.h resultview.h particle_mover.h serial_mover.h particle.h em_error.h initial_setting.h user_specific_initial_particle_setting.h
	$(CC) $(CFLAGS) $(LFLAGS) serial_em.cpp serial_solver.cpp time_controller.cpp resultview.cpp serial_mover.cpp particle.cpp em_error.cpp initial_setting.cpp user_specific_initial_particle_setting.cpp $(INCLUDE) $(LIBS) -o $@ $(XERCES_DYNAMIC_LINK)

PEM: parallel_em.cpp parallel_solver.cpp time_controller.cpp resultview.cpp parallel_mover.cpp particle.cpp em_error.cpp initial_setting.cpp user_specific_initial_particle_setting.cpp solver.h parallel_solver.h time_controller.h resultview.h particle_mover.h parallel_mover.h particle.h em_error.h initial_setting.h user_specific_initial_particle_setting.h
	${CXX} $(CFLAGS) $(LFLAGS) parallel_em.cpp parallel_solver.cpp time_controller.cpp resultview.cpp parallel_mover.cpp particle.cpp em_error.cpp initial_setting.cpp user_specific_initial_particle_setting.cpp $(INCLUDE) $(LIBS) -o $@ $(XERCES_DYNAMIC_LINK)


BGPEM: parallel_em.cpp parallel_solver.cpp time_controller.cpp resultview.cpp parallel_mover.cpp particle.cpp em_error.cpp initial_setting.cpp user_specific_initial_particle_setting.cpp solver.h parallel_solver.h time_controller.h resultview.h particle_mover.h parallel_mover.h particle.h em_error.h initial_setting.h user_specific_initial_particle_setting.h
	${BGCXX} $(CFLAGS) $(LFLAGS) parallel_em.cpp parallel_solver.cpp time_controller.cpp resultview.cpp parallel_mover.cpp particle.cpp em_error.cpp initial_setting.cpp user_specific_initial_particle_setting.cpp $(INCLUDE) $(LIBS) -o $@



# export OMP_NUM_THREADS=4
# --dt=value : value must be between 0 and 100. This means the percentage of dt of CFL condition limit. The default value is 75%
# example : --dt=70
# Have to execute "export LD_LIBRARY_PATH=$(XERCES_ROOT)/lib"
# Execute "export OMP_NUM_THREADS=8" in order to speed up. The number "8" should be modified according to your CPU numbers.
srun:	
	./${TARGET} setting_test_cooling.xml

#srun:	
#	./${TARGET} --start_time=0 --end_time=0.000001 --left_x=0 --right_x=0.256 --dx=0.008 --left_y=0 --right_y=16.384 --dy=0.008 --left_z=0 --right_z=0.256 --dz=0.008 --vis_interval=2000 --vis_left_x=0 --vis_right_x=0.256 --vis_left_y=0 --vis_right_y=16.384 --vis_left_z=0 --vis_right_z=0.256 --errormode=0

#srun:	
#	./${TARGET} --start_time=0 --end_time=0.000001 --left_x=0 --right_x=16.384 --dx=0.008 --left_y=0 --right_y=0.256 --dy=0.008 --left_z=0 --right_z=0.256 --dz=0.008 --vis_interval=2000 --vis_left_x=0 --vis_right_x=16.384 --vis_left_y=0 --vis_right_y=0.256 --vis_left_z=0 --vis_right_z=0.256 --errormode=0

#srun:
#	./${TARGET} --start_time=0 --end_time=0.000001 --left_x=0 --right_x=0.256 --dx=0.008 --left_y=0 --right_y=0.256 --dy=0.008 --left_z=0 --right_z=0.512 --dz=0.008 --vis_interval=10 --vis_left_x=0 --vis_right_x=0.256 --vis_left_y=0 --vis_right_y=0.256 --vis_left_z=0 --vis_right_z=0.512 --errormode=0

#srun:	
#	./${TARGET} --start_time=0 --end_time=0.000001 --left_x=0 --right_x=0.256 --dx=0.008 --left_y=0 --right_y=0.256 --dy=0.008 --left_z=0 --right_z=16.384 --dz=0.008 --vis_interval=2000 --vis_left_x=0 --vis_right_x=0.256 --vis_left_y=0 --vis_right_y=0.256 --vis_left_z=0 --vis_right_z=16.384 --errormode=0


# export OMP_NUM_THREADS=4
prun:
	mpirun -np 8 ./${TARGET} --procNumX=2 --procNumY=2 --procNumZ=2 --setting=setting_MPI_test_cooling.xml

#prun:
#	mpirun -np 1 -x OMP_NUM_THREADS ./${TARGET} --procNumX=1 --procNumY=1 --procNumZ=1 --start_time=0 --end_time=0.000000001 --left_x=0 --right_x=0.256 --dx=0.001 --left_y=0 --right_y=0.256 --dy=0.001 --left_z=0 --right_z=0.256 --dz=0.001 --vis_interval=-1 --vis_left_x=-0.128 --vis_right_x=0.128 --vis_left_y=-0.128 --vis_right_y=0.128 --vis_left_z=-0.128 --vis_right_z=0.128 --errormode=-2

#prun:
#	mpirun -np 8 -x OMP_NUM_THREADS ./${TARGET} --procNumX=2 --procNumY=2 --procNumZ=2 --start_time=0 --end_time=0.1 --left_x=-1 --right_x=1 --dx=0.01 --left_y=-1 --right_y=1 --dy=0.01 --left_z=-1 --right_z=1 --dz=0.01




clean:
	rm $(TARGET) *.o *.txt

tar:
	tar -cvf EM.tar *.cpp *.h Makefile *.xml *.ll *.out 

