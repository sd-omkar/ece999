MKLDIR = /usr/local/sbel/intel/2013sp1/mkl/include/
MKLLIB = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
#MKLLIB = -lmkl_solver_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
#MKLLIB = -mkl -lm
MICLIB = /tmp/

CC = icc
CFLAGS = -O3 -m64 -Wall -I. -I$(MKLDIR) -vec-report1 
#LDFLAGS = -L$(MICLIB) $(MKLLIB)
LDFLAGS = $(MKLLIB)

all : pardiso 

pardiso : pardiso.o mmio.o
	$(CC) -o pardiso pardiso.o mmio.o $(LDFLAGS) 

pardiso.o : pardiso.c
	$(CC) -c $(CFLAGS) pardiso.c

mmio.o : mmio.h mmio.c
	$(CC) -c $(CFLAGS) mmio.c

clean:
	rm -vf *.o pardiso


