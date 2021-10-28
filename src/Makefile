# **********************************
# Makefile for SWE HLLC code
#
# Use  make both for windows and linux
#  tested for gcc 5.2.0 and 11.1.0        
#
# **********************************

# Compiler
FC = gfortran

# Flags and debug option
DEBUG ?= 0
ifeq ($(DEBUG), 1)
    FLAGS = -ggdb -ffpe-trap='invalid','zero','overflow' -fcheck=all -Wall -Wno-tabs
else
    FLAGS = -ffpe-trap='invalid','zero','overflow'
endif


# REMOVE = del
REMOVE = rm

# Name of executable:
TARGET = SWE1D_HLLC

# Main program located in:
PRIN = main.f95

# Variables and parameters:
VARS = \
variables.f95 \
parameters.f95

# Modules with subroutines:
MODS = \
tests.f95 \
bc.f95 \
print_results.f95 \
initialise.f95 \
finalise.f95 \
equations.f95 \
time_integrate.f95 \
solver.f95 \
tests_bugs.f95

$(TARGET): $(PRIN) $(MODS) $(VARS)

	$(FC) $(FLAGS) $(VARS) $(MODS) $(PRIN) -o $(TARGET) && $(REMOVE) *.mod

clean:

	$(REMOVE) $(TARGET)
