#======================================================================#
#                                                                      #
#      Makefile -- build the qubit function library                    #
#                                                                      #
#======================================================================#

#
# --- COMMON CONFIG
#

CC = nvcc
EXE = demo
MY_FILE_NAME = tutorialExample
QUEST_DIR = QuEST

#
# --- compiler
#

CFLAGS	   = -dc
CLFLAGS	   = -O2 -arch=compute_30 -code=sm_30 -lineinfo
	
#
# --- targets
#
OBJ = $(MY_FILE_NAME).o qubits.o qubits_env_localGPU.o mt19937ar.o

#
# --- rules
#
%.o: %.cu
	$(CC) $(CFLAGS) $(CLFLAGS) $(CFLAGS_OMP) $<
%.o: %.cpp
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
		$(CC) $(CLFLAGS) $(CFLAGS_OMP) -o $(EXE) $(OBJ) 

.PHONY:		clean veryclean
clean:
		/bin/rm -f *.o demo
veryclean:	clean
		/bin/rm -f *.h~ *.c~ makefile~
