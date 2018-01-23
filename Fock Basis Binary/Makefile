#Make File for Static Fock Basis Binary
#Mekena Metcalf 081715

#Macros
OBJS =  main_dynamic.o Hamiltonian.o Lanczos.o Basis.o Observables.o main_functions.o
TARGET = ActAro_ED_Static.app
# Compile Macros
## This is for Mac OS
CC = clang++ -O3 -m64 -std=c++11 -stdlib=libc++ -I./ -I/usr/local/include
## This is GCC on most Linux
#CC = g++ -std=c++11 -lstdc++ -lm -I/usr/local/include
##For Merced Cluster
#CC = g++ -std=c++11 -lstdc++ -lm -I/home/mmcgrew/Library/Eigen
## This is Intel C/C++ compiler, which is under test.
#CC = icpc -O3 -Wall -std=c++11 -I/usr/local/include
##For Merced Cluster
#CC = icpc -O3 -Wall -std=c++11 -I/home/mmcgrew/Library/Eigen

#Rules
all: $(TARGET)

%.o: %.cpp
	$(CC) -c -o $@ $<

ED_Static.app: $(OBJS)
	$(CC) -o $@ $(OBJS)

ED_Dynamic.app: $(OBJS)
	$(CC) -o $@ $(OBJS)
ED_Dynamic_SSCorr.app: $(OBJS)
	$(CC) -o $@ $(OBJS)
ED_Dynamic_Corr.app: $(OBJS)
	$(CC) -o $@ $(OBJS)
ActAro_ED_Static.app: $(OBJS)
	$(CC) -o $@ $(OBJS)
ED_ProbePhase.app: $(OBJS)
	$(CC) -o $@ $(OBJS)
.PHONY: clean

clean:
	rm *.o
