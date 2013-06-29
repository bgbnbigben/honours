SRC=swarm.cpp particle.cpp neldermead.cpp
OBJS=${SRC:.cpp=.o}
EXE=particle
CXXFLAGS=-std=c++0x -Wall -g -march=native -fopenmp #-fsanitize=thread -fPIE -D_GLIBCXX_PARALLEL 
LDFLAGS=-g -fopenmp #-ltsan -fsanitize=thread -pie

all: ${OBJS}
	/usr/bin/g++-4.8 ${OBJS} -o ${EXE} ${LDFLAGS}

%.o: %.cpp %.h
	/usr/bin/g++-4.8 ${CXXFLAGS} -c -o $@ $<

clean:
	rm -f ${OBJS} ${EXE}

archive:
	tar caf source.tar.gz Makefile ${SRC} ${SRC:.cpp=.h} RTree.h rrect.h
