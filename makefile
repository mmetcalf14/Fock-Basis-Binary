#Make File for Static Fock Basis Binary
#Mekena Metcalf 081715

#Macros
OBJS =  main.o Hamiltonian.o Lanczos.o
TARGET = ED.app
# Compile Macros
CC = g++ -g
CFLAGS = -o $@ -c
LIB = -lstdc++

#Rules 
all: ${TARGET}

Hamiltonian.o : Hamiltonian_Template.cpp Hamiltonian_Template.h
		${CC} ${CFLAGS}  $<

Lanczos.o : Lanczos_Algorithm.cpp Hamiltonian_Template.h
	    ${CC} ${CFLAGS}  $<

main.o : main.cpp Hamiltonian_Template.h
	${CC} ${CFLAGS}  $<

${TARGET} : ${OBJS}
	 ${CC} -o $@ $< ${LIB}
