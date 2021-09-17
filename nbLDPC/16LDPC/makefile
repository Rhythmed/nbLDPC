# -Wextra 
LIBS = -lm
# -static -static-libgcc -static-libstdc++
CCC = g++ -O3 -Wall -DLINUX_C

OBJ_FILES_TOOLS = nbLDPC.o
OBJ_FILES = ${OBJ_FILES_TOOLS} main.o

H_FILES_TOOLS = nbLDPC.h rand.h
H_FILES = ${H_FILES_TOOLS}

all: 16LDPC

nbLDPC.o: nbLDPC.cpp ${H_FILES}
	${CCC} -c nbLDPC.cpp -o nbLDPC.o
main.o: main.cpp ${H_FILES}
	${CCC} -c main.cpp -o main.o
16LDPC: ${OBJ_FILES}
	${CCC} ${OBJ_FILES} ${LIBS} -o 16ldpc


doxygen:
	doxygen conf.doxyfile
latex: doxygen
	cd doxygen/latex;make
clean:
	rm -rf *o 16LDPC doxygen