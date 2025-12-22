# Makefile for Course Exercise 1 !Lecture 6

#--- Compiler and Flags ---#
# Define the gfortran as Fortran Compiler (FC)
FC = gfortran

# Define *general* flags (debug, warnings)
FFLAGS = -g -Wall -fcheck=all -mmacosx-version-min=11.1

# Define compile-only flags, based on lecture example 
FCFLAGS = $(FFLAGS) -c

# Define linker flags
LDFLAGS = $(FFLAGS)


#--- Source and Object Files ---#
# Define the list of object files (.o) explicitly
OBJS = geometry.o particle.o ex1.o

# Define the name of the final program
TARGET = ex1


#--- Build Rules ---#
# All depends on the target executable
all: $(TARGET)

# Rule to link the final executable
$(TARGET): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

# Pattern rule to compile source files into object files
%.o: %.f90
	$(FC) $(FCFLAGS) $< -o $@


#--- Module Dependencies ---#
# Dependencies between files
particle.o: geometry.mod
ex1.o: geometry.mod particle.mod

# Rules to generate .mod files
geometry.mod: geometry.o
particle.mod: particle.o


#--- Cleaning Rule ---#
# Rule to remove generated files
clean:
	rm -f $(TARGET) $(OBJS) *.mod

# Make Targets that are not actual files 
.PHONY: all clean