SWARMSRC=swarm/swarm.cpp swarm/particle.cpp swarm/neldermead.cpp 
BFGSSRC=bfgs/bfgs.h bfgs/linesearch.h # bfgs_test.cpp
UTILITYSRC=utilities/RTree.h utilities/rrect.h utilities/function.h utilities/tsqueue.h utilities/matrix.h utilities/vector_ops.h
PCH=utilities/vector_ops.h.gch utilities/matrix.h.gch utilities/tsqueue.h.gch
SRC=${SWARMSRC} main.cpp
F=utilities/dgemm.o utilities/dgemv.o utilities/dger.o utilities/dgetf2.o utilities/dgetrf.o utilities/dgetri.o utilities/dlamch.o utilities/dlaswp.o utilities/dswap.o utilities/dtrmm.o utilities/dtrmv.o utilities/dtrsm.o utilities/dtrti2.o utilities/dtrtri.o utilities/idamax.o utilities/ieeeck.o utilities/ilaenv.o utilities/iparmq.o utilities/lsame.o utilities/xerbla.o
OBJS=${SRC:.cpp=.o} ${F}
EXE=particle
CXXFLAGS=-std=c++0x -I. -I${HOME}/openmpi/env/include -Wall -Wextra -g -march=native -fopenmp -D_GLIBCXX_PARALLEL #-fsanitize=thread -fPIE
LDFLAGS=-g -fopenmp -Llbfgsb -llbfgsb -lgfortran #-ltsan -fsanitize=thread -pie

all: ${OBJS} ${PCH} fortran
	/usr/bin/g++-4.8 ${OBJS} -o ${EXE} ${LDFLAGS}

fortran:
	$(MAKE) -C lbfgsb library 

%.o: %.cpp $(wildcard %.h)
	/usr/bin/g++-4.8 ${CXXFLAGS} -c -o $@ $<

%.o: %.f
	gfortran -O3 -Wall -fbounds-check -Wno-uninitialized -c $< -o $@	

%.gch: %
	/usr/bin/g++-4.8 ${CXXFLAGS} $<

clean:
	rm -f ${OBJS} ${EXE}

fullclean:
	rm -f ${OBJS} ${EXE} ${PCH}

archive:
	tar caf source.tar.gz Makefile ${SRC} ${SRC:.cpp=.h} ${UTILITYSRC}
