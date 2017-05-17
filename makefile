#======================================================================#
#                                                                      #
#      Makefile -- build the qubit function library                    #
#                                                                      #
#======================================================================#

#
# --- COMMON CONFIG
#

# COMPILER options: GNU, INTEL
COMPILER = INTEL
EXE = demo
MY_FILE_NAME = basicTemplate
USE_MPI=0
# Do not use MPI and GPU at the same time
USE_GPU=1
QUEST_DIR = QUEST

#
# --- compiler
#

ifneq ($(USE_MPI), 0)
	ifeq ($(COMPILER), GNU)
		# GCC compilers
		CC         = mpicc
		CFLAGS	   = -c
		CLFLAGS    = -O2 -std=c99 -mavx -Wall
		CFLAGS_OMP = -fopenmp
	else ifeq ($(COMPILER), INTEL)
		# Mvapich2
		CC         = mpicc
		CFLAGS	   = -c 
		CLFLAGS     = -O2 -std=c99
		CFLAGS_OMP = -openmp
	else 
		$(error " *** error: invalid compiler")
	endif
else
	ifneq ($(USE_GPU), 0)
		CC	   = nvcc
		CFLAGS	   = -dc
		CLFLAGS	   = -O2 -arch=compute_30 -code=sm_30 -lineinfo
		CFLAGS_OMP = -Xcompiler -fopenmp 
	else ifeq ($(COMPILER), GNU)
		# GCC compilers
		CC         = gcc
		CFLAGS	   = -c
		CLFLAGS     = -O2 -std=c99 -mavx -Wall
		CFLAGS_OMP = -fopenmp
	else ifeq ($(COMPILER), INTEL)
		# Intel compilers
		CC         = icc
		CFLAGS	   = -c
		CLFLAGS     = -O2 -std=c99 -Wall -xAVX -axCORE-AVX2 -restrict
		CFLAGS_OMP = -openmp
	else 
		$(error " *** error: invalid compiler")
	endif

endif

#
# --- libraries
#
LIBS = -lm -lgomp


#
# --- targets
#
OBJ = $(MY_FILE_NAME).o qubits.o
ifneq ($(USE_MPI), 0)
	OBJ += qubits_env_mpi.o
else
	ifneq ($(USE_GPU), 0)
		OBJ += qubits_env_localGPU.o
	else
		OBJ += qubits_env_local.o
	endif
endif

#
# --- rules
#
%.o: %.c
	$(CC) $(CFLAGS) $(CLFLAGS) $(CFLAGS_OMP) $<
%.o: %.cu
	$(CC) $(CFLAGS) $(CLFLAGS) $(CFLAGS_OMP) $<
%.o: %.cpp
	$(CC) $(CFLAGS) $(CLFLAGS) $(CFLAGS_OMP) $<

%.o: $(QUEST_DIR)/%.c
	$(CC) $(CFLAGS) $(CLFLAGS) $(CFLAGS_OMP) $<

%.o: $(QUEST_DIR)/%.cu
	$(CC) $(CFLAGS) $(CLFLAGS) $(CFLAGS_OMP) $<
%.o: $(QUEST_DIR)/%.cpp
	$(CC) $(CFLAGS) $(CLFLAGS) $(CFLAGS_OMP) $<


#
# --- build
#
default:	demo

demo:		$(OBJ)
		$(CC) $(CLFLAGS) $(CFLAGS_OMP) -o $(EXE) $(OBJ) $(LIBS)

.PHONY:		clean veryclean
clean:
		/bin/rm -f *.o demo
veryclean:	clean
		/bin/rm -f *.h~ *.c~ makefile~
