COMPILER=g++
FLAGS = -O3 -std=c++17 -fopenmp -Wall -fpermissive
BIN_DIR = bin/

all: directories marzone

directories:
	mkdir -p $(BIN_DIR)

marzone: marzone.c marzone.h
	$(COMPILER) $(FLAGS) marzone.c input.c output.c -o $(BIN_DIR)$@

