# Compiler
FC = gfortran

# Compiler flags
# FFLAGS = -O2
FFLAGS = -g

# BLAS library
LIBS = -L/usr/lib/x86_64-linux-gnu -lopenblas

# Target executable
TARGET = test_blas

# Source files
SRC = test_blas.f95

# Build the program
$(TARGET): $(SRC)
	$(FC) $(FFLAGS) $(SRC) $(LIBS) -o $(TARGET)

# Clean up generated files
.PHONY: clean
clean:
	rm -f $(TARGET)
