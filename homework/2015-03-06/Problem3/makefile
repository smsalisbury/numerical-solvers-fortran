# This makefile handles the compiling of the source code
# Start of the makefile

# ======================
# DECLARATIONS
# ======================

# The compiler
FC = gfortran

# Directories
MOD_DIR = /home/spencer/numerical-solvers-fortran/modules
OBJ_DIR = obj
INCLUDES = -I. -I$(MOD_DIR)

# Compiling flags
FCFLAGS = -g -O2 -fbounds-check

# Program name
PROGRAM = a.out
MOD_OBJECTS = config.o functions.o file_io.o derivate1D.o roots1D.o initial_value_1D.o main.o

# ======================
# MAKE COMMANDS
# ======================

# Default make command ($ make)
all: $(PROGRAM)

# Make the main program
$(PROGRAM): $(MOD_OBJECTS)
	$(FC) -o $(PROGRAM) $(MOD_OBJECTS)

# Make the objects
main.o:
	$(FC) $(FCFLAGS) -c $(FILENAME) main.f90

config.o:
	$(FC) $(FCFLAGS) -c $(FILENAME) config.f90
	
functions.o:
	$(FC) $(FCFLAGS) -c $(FILENAME) functions.f90

initial_value_1D.o:
	$(FC) $(FCFLAGS) -c $(FILENAME) $(MOD_DIR)/initial_value_1D.f90

file_io.o:
	$(FC) $(FCFLAGS) -c $(FILENAME) $(MOD_DIR)/file_io.f90

roots1D.o:
	$(FC) $(FCFLAGS) -c $(FILENAME) $(MOD_DIR)/roots1D.f90

derivate1D.o:
	$(FC) $(FCFLAGS) -c $(FILENAME) $(MOD_DIR)/derivate1D.f90

# Special - make clean
clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAM)