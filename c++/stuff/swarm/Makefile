TLD=${CURDIR}
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

all: libspatialindex ${OBJS} ${PCH} fortran
	/usr/bin/g++-4.8 ${OBJS} -o ${EXE} ${LDFLAGS}

libspatialindex: ${TLD}/spatialindex/lib/libspatialindex.a

${TLD}/spatialindex/lib/libspatialindex.a:
	mkdir -p spatialindex
	if [ ! -e ${TLD}/libspatialindex/Makefile ]; then \
		cd libspatialindex && ./configure --prefix=${TLD}/spatialindex; \
	fi
	$(MAKE) -C libspatialindex all install

fortran:
	$(MAKE) -C lbfgsb library 

%.o: %.cpp $(wildcard %.h)
	/usr/bin/g++-4.8 ${CXXFLAGS} -c -o $@ $<

%.o: %.f
	gfortran -O3 -Wall -fbounds-check -Wno-uninitialized -c $< -o $@	

%.gch: %
	/usr/bin/g++-4.8 ${CXXFLAGS} $<

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
