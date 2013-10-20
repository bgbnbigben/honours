FC = gfortran
AR = ar
RM = rm -f

FFLAGS = -g -Wall -fbounds-check -Wno-uninitialized 

LBFGSB  = lbfgsb.f
LINPACK = linpack.f
BLAS    = blas.f
TIMER   = timer.f

OBJS    = lbfgsb.o linpack.o blas.o timer.o

LIBNAME = liblbfgsb.a
.PHONY: all library clean fullclean

all: library

clean:
	$(RM) x.* *.o

fullclean: clean
	$(RM) $(LIBNAME)

library: $(LIBNAME)

%.a: $(OBJS)
	$(AR) rvs $(LIBNAME) $(OBJS)

%.o: %.f
	$(FC) $(FFLAGS) -c $<
