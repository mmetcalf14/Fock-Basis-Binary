#Make File for Static Fock Basis Binary
#Mekena Metcalf 081715

#Macros
OBJS =  main_static.o Hamiltonian_Template.o Lanczos_Algorithm.o
TARGET = ED_Static.app
# Compile Macros
CC = g++ -g -I./
#CFLAGS = -o $@ -c
LIB = -lstdc++

#Rules 
all: ${TARGET}

*.o: *.cpp
	${CC} -c -o $@ $<

ED_Static.app: $(OBJS)
	${CC} -o $@ $(OBJS)

ED_Dynamic.app: $(OBJS)
	${CC} -o $@ $(OBJS)

clean:
	rm *.o
