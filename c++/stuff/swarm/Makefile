SRC=swarm.cpp particle.cpp
OBJS=${SRC:.cpp=.o}
EXE=particle

all: ${OBJS}
	g++ ${OBJS} -o ${EXE}

%.o: %.cpp %.h
	g++ -std=c++0x -Wall -c -o $@ $<

clean:
	rm -f ${OBJS} ${EXE}
