#======================================================================#
#                                                                      #
#      Makefile														   # 
#      -- build the QuEST function library and link user sources       #
#                                                                      #
#======================================================================#

#
# --- USER CONFIG
#

# COMPILER options: {GNU, INTEL}
COMPILER = GNU

# EXECUTABLE TO GENERATE
EXE = demo

# USER SOURCE FILES
# This makefile expects all user sources to be in the root directory.
# If using more than one source, separate by spaces eg 
# MY_C_SOURCES = tutorialExample helperFunctions
MY_C_SOURCES = tutorialExample

# GPU architecture dependent. To find the compute capability for your
# device see https://developer.nvidia.com/cuda-gpus
# Remove . before entering number ie for compute capability 3.0 write 30
GPU_COMPUTE_CAPABILITY = 30

# PATH TO QUEST LIBRARY SOURCES FROM ROOT DIRECTORY
QUEST_DIR = QuEST





#======================================================================#
#                                                                      #
#      Makefile execution 											   #
#                                                                      #
#======================================================================#

INCLUDE = -I$(QUEST_DIR)

## ------------------- COMPILE USER SOURCES ---------------------------

#
# --- C COMPILER
#

ifeq ($(COMPILER), GNU)
	# COMPILER = GNU
	QuEST_CC = gcc
else ifeq ($(COMPILER), INTEL)
	# COMPILER = INTEL
	QuEST_CC = icc
else 
	$(error " *** build error: invalid compiler")
endif

#
# --- C COMPILER FLAGS
#


ifeq ($(COMPILER), GNU)
	# GCC compilers
	QuEST_CFLAGS     = -x c -O2 -std=c99 -mavx -Wall
else ifeq ($(COMPILER), INTEL)
	# Intel compilers
	QuEST_CFLAGS     = -x c -O2 -std=c99 -fprotect-parens -Wall -xAVX -axCORE-AVX2 -diag-disable cpu-dispatch
else 
	$(error " *** error: invalid compiler")
endif

#
## --- GPU compiler
#

QuEST_CPP 		   = nvcc
QuEST_CPPFLAGS       = -dc -O2 -arch=compute_$(GPU_COMPUTE_CAPABILITY) -code=sm_$(GPU_COMPUTE_CAPABILITY) -ccbin $(QuEST_CC)
QuEST_CPPLFLAGS      = -O2 -arch=compute_$(GPU_COMPUTE_CAPABILITY) -code=sm_$(GPU_COMPUTE_CAPABILITY) -lineinfo -ccbin $(QuEST_CC)

#
# --- targets
#
USER_OBJ = $(addsuffix .o, $(MY_C_SOURCES))
OBJ = qubits.o qubits_env_localGPU.o mt19937ar.o
OBJ += $(USER_OBJ)

all: $(OBJ)
	$(QuEST_CPP) $(INCLUDE) $(QuEST_CPPLFLAGS) -o $(EXE) $(OBJ)


%.o: %.c
	$(QuEST_CC) $(INCLUDE) $(QuEST_CFLAGS) -c $<

%.o: %.cu
	$(QuEST_CPP) $(QuEST_CPPFLAGS) $<
%.o: %.cpp
	$(QuEST_CPP) $(QuEST_CPPFLAGS) $<

%.o: $(QUEST_DIR)/%.cu
	$(QuEST_CPP) $(QuEST_CPPFLAGS) $<
%.o: $(QUEST_DIR)/%.cpp
	$(QuEST_CPP) $(QuEST_CPPFLAGS) $<


.PHONY:		clean veryclean
clean:
	/bin/rm -f *.o $(EXE)
veryclean:	clean
	/bin/rm -f *.h~ *.c~ makefile~


