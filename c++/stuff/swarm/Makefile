SWARMSRC=swarm/swarm.cpp swarm/particle.cpp swarm/neldermead.cpp 
BFGSSRC=bfgs/bfgs.h bfgs/linesearch.h # bfgs_test.cpp
UTILITYSRC=utilities/RTree.h utilities/rrect.h utilities/function.h utilities/tsqueue.h utilities/matrix.h utilities/vector_ops.h
PCH=utilities/vector_ops.h.gch utilities/matrix.h.gch utilities/tsqueue.h.gch
SRC=${SWARMSRC} main.cpp
OBJS=${SRC:.cpp=.o}
EXE=particle
CXXFLAGS=-std=c++0x -I. -Wall -g -march=native -fopenmp #-fsanitize=thread -fPIE -D_GLIBCXX_PARALLEL 
LDFLAGS=-g -fopenmp #-ltsan -fsanitize=thread -pie

all: ${OBJS} ${PCH}
	/usr/bin/g++-4.8 ${OBJS} -o ${EXE} ${LDFLAGS}

%.o: %.cpp $(wildcard %.h)
	/usr/bin/g++-4.8 ${CXXFLAGS} -c -o $@ $<

%.gch: %
	/usr/bin/g++-4.8 ${CXXFLAGS} $<

clean:
	rm -f ${OBJS} ${EXE} ${PCH}

archive:
	tar caf source.tar.gz Makefile ${SRC} ${SRC:.cpp=.h} ${UTILITYSRC}
