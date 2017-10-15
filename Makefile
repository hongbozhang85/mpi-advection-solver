.SUFFIXES:
.PRECIOUS: %.o

HDRS=serAdvect.h parAdvect.h
OBJS=serAdvect.o parAdvect.o
PROG=testAdvect
CCFLAGS=-O3

all: $(PROG) 

%: %.o $(OBJS)
	mpicc -o $* $*.o $(OBJS) -lm
#	mpicc -DBLOCK -o $* $*.o $(OBJS) -lm
%.o: %.c $(HDRS)
	mpicc -Wall $(CCFLAGS) -c $*.c
#	mpicc -DBLOCK -Wall $(CCFLAGS) -c $*.c
clean:
	rm -f *.o $(PROG)
