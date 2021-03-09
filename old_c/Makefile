COMPILER=gcc
FLAGS = -O3 -std=c++17 -fopenmp -Wall 
FLAGS2 = -Wall -std=c99
BIN_DIR = bin/

all: directories marzone

directories:
	mkdir -p $(BIN_DIR)

marzone.o: marzone.c marzone.h
	$(COMPILER) $(FLAGS2) -c marzone.c input.c

marzone: marzone.o
	$(COMPILER) marzone.o -o $(BIN_DIR)$@ 

