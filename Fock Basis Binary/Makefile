#Make File for Static Fock Basis Binary
#Mekena Metcalf 081715

#Macros
OBJS =  main_dynamic.o Hamiltonian.o Lanczos.o Basis.o
TARGET = ED_Dynamic.app
# Compile Macros
CC = clang++ -O3 -m64 -std=c++11 -stdlib=libc++ -I./ -I/usr/local/include
# CC = g++ -lstdc++ -lm -I/usr/local/include
#CFLAGS = -o $@ -c
LIB = -lstdc++

#Rules
all: $(TARGET)

%.o: %.cpp
	$(CC) -c -o $@ $<

ED_Static.app: $(OBJS)
	$(CC) -o $@ $(OBJS)

ED_Dynamic.app: $(OBJS)
	$(CC) -o $@ $(OBJS)

clean:
	rm *.o
