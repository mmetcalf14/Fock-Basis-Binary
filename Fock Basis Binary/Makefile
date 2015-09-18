ifeq "$(ARCH)" ""
ARCH = x86_64
endif

ifeq "$(OS)" ""
OS = $(shell uname -s)
endif

ifeq "$(NODE)" ""
NODE = $(shell uname -n)
endif

ifeq ("$(NODE)", "kagome.rcc.ucmerced.edu")
	EIGENINC = /usr/local/include
	HDF5ROOT = /data2/hdf5/HDF_Group/HDF5/1.8.15-patch1
	MKLROOT = /data2/intel/composer_xe_2015.2.164/mkl
else ifeq ("$(NODE)", "edgestate.ucmerced.edu")
	EIGENINC = /usr/local/include
	HDF5ROOT = /opt/HDF_Group/HDF5/1.8.15-patch1
	MKLROOT = /opt/intel/compilers_and_libraries_2016.0.109/linux/mkl
else ifneq (, $(filter "$(NODE)", "comet-ln1.sdsc.edu" "comet-ln2.sdsc.edu" "comet-ln3.sdsc.edu" "comet-ln4.sdsc.edu"))
	EIGENINC = /opt/eigen/include
	HDF5ROOT = /home/chenyen/hdf5/HDF_Group/HDF5/1.8.15-patch1
	MKLROOT = /opt/intel/composer_xe_2015.2.164/mkl
endif

ifeq ("$(OS)", "Darwin")
	OPENMP = #-fopenmp
	CC = clang++ $(OPENMP) -O3 -m64 -std=c++11 -stdlib=libc++ -I./
	HDF5LIB = -lhdf5 -lhdf5_cpp
	LAPACK = -lblas -llapack -lm
	LAPACK_OMP = -lblas -llapack -lm
else ifeq ("$(OS)", "Linux")
	OPENMP = -openmp
	#NOTE: the order of linker matters!
	HDF5LIB =  $(HDF5ROOT)/lib/libhdf5_cpp.a $(HDF5ROOT)/lib/libhdf5.a $(HDF5ROOT)/lib/libz.a
	LAPACK = $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
	$(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group \
	$(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
	$(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a \
	-Wl,--end-group -lpthread -lm
	LAPACK_OMP = $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
	$(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group \
	${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
	${MKLROOT}/lib/intel64/libmkl_core.a \
	${MKLROOT}/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm
	CC = icpc $(OPENMP) -O3 -Wall -std=c++11 -I./ -I$(HDF5ROOT)/include -I$(EIGENINC) -DMKL
endif

LIB := Hamiltonian_Template.o Lanczos_Algorithm.o

all: ed.app

ed.app: main.o $(LIB)
	$(CC) -o $@ $< $(LIB)

main.o: main_static.cpp
	$(CC) -c -o $@ $<

Hamiltonian_Template.o: Hamiltonian_Template.cpp
	$(CC) -c -o $@ $<

Lanczos_Algorithm.o: Lanczos_Algorithm.cpp
	$(CC) -c -o $@ $<

.PHONY: clean

clean:
	rm -f *.o *.app
