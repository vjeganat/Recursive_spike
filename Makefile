# Adapted from https://sites.google.com/lbl.gov/cs267-spr2018/hw-2?authuser=0

CC = g++
MPCC = mpicxx
CFLAGS = -O3 -std=c++11
LIBS = 


TARGETS = serial mpi

all:	$(TARGETS)

serial: serial.o common.o
mpi: mpi.o common.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi.o common.o -lm

mpi.o: mpi.cpp common.h
	$(MPCC) -c $(CFLAGS) mpi.cpp

serial.o: serial.cpp common.h

common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt
