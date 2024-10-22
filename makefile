# Project name
PROJNAME=vagmd

# Definition of the compiler to be used
CC=gcc

# Binary path
BINPATH=./bin/$(PROJNAME)

# Definition of C files
CFILES=$(wildcard ./src/*.c) $(wildcard ./src/*/*.c)

# Definition of the object files
OBJ=$(subst .c,.o,$(CFILES))

# Run the binary
run:
	@ mkdir -p results
	@ time ./bin/$(PROJNAME) \
	-membrane_area 12.96 \
	-vacuum_pressure -81325.0 \
	-number_channels 6 \
	-feed_mass_flow_rate 0.083333333333333333333333 \
	-cool_mass_flow_rate 0.083333333333333333333333
	@ cat ./results/report.csv

# Build the binary
build: binfolder $(BINPATH)

$(BINPATH): $(OBJ)
	@ $(LINK.C) -s -O3 -o $@ $^ $(LDLIBS)
	@ rm -rf ./src/*.o
	@ rm -rf ./src/*/*.o

# Create bin folder
binfolder:
	@ mkdir -p bin

# Clean all files in the bin and results folders
cleanall:
	@ rm -rf ./bin ./graphs
	@ rm -rf ./*.o ./src/*.o ./src/*/*.o

# Inclusion of PETSc config information
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test