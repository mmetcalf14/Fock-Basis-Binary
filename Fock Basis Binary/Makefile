ifeq "$(ARCH)" ""
ARCH = x86_64
endif

ifeq "$(OS)" ""
OS = $(shell uname -s)
endif

ifeq ("$(OS)", "Darwin")
	CC = clang++ -O3 -m64 -std=c++11 -stdlib=libc++ -I./
	LAPACK = -lblas -llapack -lm
else ifeq ("$(OS)", "Linux")
	MKLROOT= /data2/intel/composer_xe_2015.2.164/mkl
	LAPACK = $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
	$(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group \
	$(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
	$(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a \
	-Wl,--end-group -lpthread -lm
	CC = icpc -O3 -Wall -std=c++11 -I./ -I$(HDF5ROOT)/include -DMKL
endif

LIB := Hamiltonian_Template.o Lanczos_Algorithm.o

all: ed.app

ed.app: main.o $(LIB)
	$(CC) -o $@ $< $(LIB)

main.o: main.cpp
	$(CC) -c -o $@ $<

Hamiltonian_Template.o: Hamiltonian_Template.cpp
	$(CC) -c -o $@ $<

Lanczos_Algorithm.o: Lanczos_Algorithm.cpp
	$(CC) -c -o $@ $<

.PHONY: clean

clean:
	rm -f *.o *.app
