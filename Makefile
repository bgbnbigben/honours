TLD=${CURDIR}
CC=$(shell which g++)
MPICC=$(shell which mpicc)
SWARMSRC=swarm/swarm.cpp swarm/particle.cpp swarm/neldermead.cpp 
BFGSSRC=bfgs/bfgs.h bfgs/linesearch.h # bfgs_test.cpp
UTILITYSRC=utilities/RTree.h utilities/rrect.h utilities/function.h utilities/tsqueue.h utilities/matrix.h utilities/vector_ops.h
PCH=utilities/vector_ops.h.gch utilities/matrix.h.gch utilities/tsqueue.h.gch
SRC=${SWARMSRC} main.cpp
F=utilities/dgemm.o utilities/dgemv.o utilities/dger.o utilities/dgetf2.o utilities/dgetrf.o utilities/dgetri.o utilities/dlamch.o utilities/dlaswp.o utilities/dswap.o utilities/dtrmm.o utilities/dtrmv.o utilities/dtrsm.o utilities/dtrti2.o utilities/dtrtri.o utilities/idamax.o utilities/ieeeck.o utilities/ilaenv.o utilities/iparmq.o utilities/lsame.o utilities/xerbla.o
OBJS=${SRC:.cpp=.o} ${F}
EXE=particle
CXXFLAGS=-std=c++0x -I. -I${HOME}/openmpi/env/include -Ispatialindex/include -Wall -Wextra -g -march=native -fopenmp -D_GLIBCXX_PARALLEL #-fsanitize=thread -fPIE
LDFLAGS=-g -fopenmp -Llbfgsb -Lspatialindex/lib -llbfgsb -lgfortran -lspatialindex #-ltsan -fsanitize=thread -pie

.PHONY: all fortran libspatialindex clean fullclean archive

all: fortran libspatialindex ${OBJS} ${PCH}
	$(CC) ${OBJS} -o ${EXE} ${LDFLAGS}

libspatialindex: fortran ${TLD}/spatialindex/lib/libspatialindex.a

${TLD}/spatialindex/lib/libspatialindex.a:
	mkdir -p spatialindex
	if [ ! -e ${TLD}/libspatialindex/Makefile ]; then \
		if [ ! -e ${TLD}/libspatialindex/configure ]; then \
			cd ${TLD}/libspatialindex && ./autogen.sh; \
		fi; \
		cd ${TLD}/libspatialindex && ./configure --prefix=${TLD}/spatialindex; \
	fi;
	$(MAKE) -C libspatialindex all install

fortran:
	$(MAKE) -C lbfgsb library 

main.o: main.cpp
	$(MPICC) ${CXXFLAGS} -c -o $@ $<

%.o: %.cpp $(wildcard %.h)
	$(CC) ${CXXFLAGS} -c -o $@ $<

%.o: %.f
	gfortran -O3 -Wall -fbounds-check -Wno-uninitialized -c $< -o $@	

%.gch: %
	$(CC) ${CXXFLAGS} $<

clean:
	$(MAKE) -C lbfgsb fullclean
	rm -f ${OBJS} ${EXE}
	if [ -e ${TLD}/libspatialindex/Makefile ]; then \
		$(MAKE) -C libspatialindex distclean; \
	fi
	rm -rf spatialindex/

fullclean:
	rm -f ${OBJS} ${EXE} ${PCH}

archive:
	tar caf source.tar.gz Makefile ${SRC} ${SRC:.cpp=.h} ${UTILITYSRC}
