TARGET=""
ifeq "${TARGET}" "MPI"
EXE				    = basic_lib_MPI.exe
COMPILATOR 			= mpic++
ifneq "$(strip $(findstring MINGW32, $(shell uname)))" ""
COMPILATOR = g++
endif
COMPILATOR_OPTIONS 	= -O3 -std=c++0x -DUSE_MPI
L1 					= $(NOMAD_HOME)/lib/nomad.MPI.a
LIBS 				= $(L1) -lm -lmpi
OBJS 				= basic_lib_MPI.o
else
EXE                 = basic_lib.exe
COMPILATOR          = g++
COMPILATOR_OPTIONS  = -O3 -std=c++0x
L1                  = $(NOMAD_HOME)/lib/nomad.a
LIBS                = $(L1) -lm
OBJS                = basic_lib.o
endif

INCLUDE             = -I$(NOMAD_HOME)/src -I.
COMPILE             = $(COMPILATOR) $(COMPILATOR_OPTIONS) $(INCLUDE) -c

.PHONY: all clean check_env
.SILENT:
.SUFFIXES:

scalar: $(EXE)
mpi: $(EXE)
	make TARGET=MPI
all: scalar mpi

$(EXE): check_env $(L1) $(OBJS)
	$(COMPILATOR) $(OBJS) $(LIBS) $(COMPILATOR_OPTIONS) -o $(EXE)
	strip $(EXE)

%.o: basic_lib.cpp
	$(COMPILE) -o $@ $<

check_env:
ifndef NOMAD_HOME
	$(error NOMAD_HOME is undefined)
endif

clean:
	rm -f *.o *.exe
