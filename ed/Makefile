
ifeq "$(ARCH)" ""
ARCH = x86_64
endif

ifeq "$(OS)" ""
OS = $(shell uname -s)
endif

ifeq ("$(OS)", "Darwin")
	CC = clang++ -O3 -m64 -std=c++11 -stdlib=libc++ -I./
	LAPACK = -lblas -llapack -lm
	HDF5LIB = -lhdf5 -lhdf5_cpp
else ifeq ("$(OS)", "Linux")
	HDF5DIR = /home/share/apps/HDF_Group/HDF5/1.8.14
	CC =icpc -O3 -Wall -std=c++11 -I./ -I$(HDF5DIR)/include -DMKL
	MKLROOT=/home/share/apps/intel/composer_xe_2013_sp1.3.174/mkl
	LAPACK = $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a  $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a  -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
	$(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
	HDF5LIB = -L$(HDF5DIR)/lib -lhdf5 -lhdf5_cpp
endif

SRC = fbasis.o fhamiltonian.o lanczos.o

fed.app:fmain.o fbasis.o fhamiltonian.o lanczos.o
	$(CC) -o $@ $< ${SRC}  $(LAPACK)

fmain.o: fermion/main.cpp
	$(CC) -c -o $@ $<

node.o: LinkedList/node.cpp
	$(CC) -c -o $@ $<

fbasis.o: fermion/basis.cpp
	$(CC) -c -o $@ $<

fhamiltonian.o: fermion/hamiltonian.cpp
	$(CC) -c -o $@ $<

eigen.app: EigenMatrix.o
	$(CC) -o $@ $<

EigenMatrix.o: matrix/EigenMatrix.cpp
	$(CC) -c -o $@ $<

search.app: fermion/search.cpp
	$(CC) -o $@ $<

search.o: fermion/search.cpp
	$(CC) -c -o $@ $<

lanczos.o: lanczos/lanczos.cpp
	$(CC) -c -o $@ $<

.PHONY:
	clean

clean:
	rm -f *.o *.app
