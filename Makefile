
# Makefile for the sequential version

CC = gcc

CFLAGS = -Wall -g -fno-builtin-cexp

LIBRARIES = -lm

EXECNAME = comus3

OBJS = comus.o comus_streec.o  comus_phylo.o rand1.o comus_evolve.o gamma.o model.o twister.o global.o nucmodels.o  aamodels.o eigen.o treefile.o

HFILE = comus.h comus_tree.h

all: $(EXECNAME)

$(EXECNAME) : $(OBJS) $(HFILE)
	$(CC) $(CFLAGS) -o $(EXECNAME) $(OBJS) $(LIBRARIES)

.c.o:
	$(CC) $(CFLAGS) -c $*.c

comus.o: comus.c
	$(CC) $(CFLAGS) -c comus.c

comus_streec.o: comus_streec.c
	$(CC) $(CFLAGS) -c comus_streec.c

comus_phylo.o: comus_phylo.c
	$(CC) $(CFLAGS) -c comus_phylo.c

comus_evolve.o: comus_evolve.c
	$(CC) $(CFLAGS) -c comus_evolve.c

gamma.o : gamma.c gamma.h
	$(CC) $(CFLAGS) -c gamma.c

model.o : model.c model.h
	$(CC) $(CFLAGS) -c model.c

twister.o : twister.c twister.h
	  $(CC) $(CFLAGS) -c twister.c

global.o : global.c global.h
	 $(CC) $(CFLAGS) -c global.c

nucmodels.o : nucmodels.c nucmodels.h
	    $(CC) $(CFLAGS) -c nucmodels.c

rand1.o: rand1.c
	$(CC) $(CFLAGS) -c rand1.c

clean:
	rm -f $(EXECNAME)
	rm -f $(OBJS)
