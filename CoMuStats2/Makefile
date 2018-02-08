
# Makefile for the sequential version

CC = g++ 

CFLAGS = -std=c++11 -Wall  -g

LIBRARIES = -lm -lsequence -lz -static

INCLUDE = -I/usr/local/include/Sequence

EXECNAME = CoMuStats

OBJS = comustats.o

HFILE = 

all: $(EXECNAME)

CoMuStats : $(OBJS) $(HFILE)
	$(CC) $(CFLAGS) -o $(EXECNAME) $(OBJS)  $(INCLUDE)  $(LIBRARIES)

comustats.o: comustats.cc
	$(CC) $(CFLAGS) -c comustats.cc
clean:
	rm -f $(EXECNAME)
	rm -f $(OBJS)
